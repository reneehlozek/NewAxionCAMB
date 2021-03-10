    module DarkEnergyFluid
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyFluid
        !comoving sound speed is always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyFluid_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyFluid_SelfPointer
    procedure :: Init =>TDarkEnergyFluid_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    end type TDarkEnergyFluid

    !Example implementation of fluid model using specific analytic form
    !(approximate effective axion fluid model from arXiv:1806.10608, with c_s^2=1 if n=infinity (w_n=1))
    !This is an example, it's not supposed to be a rigorous model!  (not very well tested)
    type, extends(TDarkEnergyModel) :: TAxionEffectiveFluid
        real(dl) :: w_n = 1._dl !Effective equation of state when oscillating
        real(dl) :: Om = 0._dl !Omega of the early DE component today (assumed to be negligible compared to omega_lambda)
        real(dl) :: a_c  !transition scale factor
        real(dl) :: theta_i = const_pi/2 !Initial value
        real(dl), private :: pow, omL, acpow, freq, n !cached internally
    contains
    procedure :: ReadParams =>  TAxionEffectiveFluid_ReadParams
    procedure, nopass :: PythonClass => TAxionEffectiveFluid_PythonClass
    procedure, nopass :: SelfPointer => TAxionEffectiveFluid_SelfPointer
    procedure :: Init => TAxionEffectiveFluid_Init
    procedure :: w_de => TAxionEffectiveFluid_w_de
    procedure :: grho_de => TAxionEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TAxionEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionEffectiveFluid_PerturbationEvolve
    end type TAxionEffectiveFluid

    ! RH added new modified DE fluid
    type, extends(TDarkEnergyModel) :: TModDeltaEffectiveFluid
        integer ::  max_abin=2 ! number of bins allowed
        real(dl), dimension(2) :: deltavals  !Values of the delta bin
        real(dl), dimension(2) :: lnavals ! values of ln a
        real(dl) :: smoothscale = 0.3_dl !Smoothing scale
        
        real(dl) :: grhornomass 
        real(dl) :: grhocrit
        real(dl) :: grhov
        real(dl) :: grho
        real(dl) :: grhoc
        real(dl) :: grhom
        real(dl) :: grhog
        real(dl) :: grhob
        real(dl) :: grhok
        real(dl) :: a_c  !transition scale factor
        real(dl), private :: pow, omL, acpow, freq, n !cached internally
    contains
    procedure :: ReadParams =>  TModDeltaEffectiveFluid_ReadParams
    procedure, nopass :: PythonClass => TModDeltaEffectiveFluid_PythonClass
    procedure, nopass :: SelfPointer => TModDeltaEffectiveFluid_SelfPointer
    procedure :: Init => TModDeltaEffectiveFluid_Init
    procedure :: w_de => TModDeltaEffectiveFluid_w_de
    procedure :: grho_de => TModDeltaEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TModDeltaEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TModDeltaEffectiveFluid_PerturbationEvolve
    procedure :: TotalDensity => TModDeltaEffectiveFluid_TotalDensity
    end type TModDeltaEffectiveFluid

    contains


    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergyFluid_ReadParams


    function TDarkEnergyFluid_PythonClass()
    character(LEN=:), allocatable :: TDarkEnergyFluid_PythonClass

    TDarkEnergyFluid_PythonClass = 'DarkEnergyFluid'

    end function TDarkEnergyFluid_PythonClass

    subroutine TDarkEnergyFluid_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyFluid), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyFluid_SelfPointer

    subroutine TDarkEnergyFluid_Init(this, State)
    use classes
    class(TDarkEnergyFluid), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State

    call this%TDarkEnergyEqnOfState%Init(State)

    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        if (this%use_tabulated_w) then
            if (any(this%equation_of_state%F<-1)) &
                error stop 'Fluid dark energy model does not allow w crossing -1'
        elseif (this%wa/=0 .and. &
            ((1+this%w_lam < -1.e-6_dl) .or. 1+this%w_lam + this%wa < -1.e-6_dl)) then
            error stop 'Fluid dark energy model does not allow w crossing -1'
        end if
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyFluid_Init


    subroutine TDarkEnergyFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    if (this%no_perturbations) then
        dgrhoe=0
        dgqe=0
    else
        dgrhoe = ay(w_ix) * grhov_t
        dgqe = ay(w_ix + 1) * grhov_t * (1 + w)
    end if
    end subroutine TDarkEnergyFluid_PerturbedStressEnergy


    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, loga

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - w) *  (y(w_ix) + (1 + w) * Hv3_over_k) &
        -  (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
    if (this%use_tabulated_w) then
        !account for derivatives of w
        loga = log(a)
        if (loga > this%equation_of_state%Xmin_interp .and. loga < this%equation_of_state%Xmax_interp) then
            ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)* Hv3_over_k
        end if
    elseif (this%wa/=0) then
        ayprime(w_ix) = ayprime(w_ix) + Hv3_over_k*this%wa*adotoa*a
    end if
    !velocity
    if (abs(w+1) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
            k * this%cs2_lam * y(w_ix) / (1 + w)
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TDarkEnergyFluid_PerturbationEvolve



    subroutine TAxionEffectiveFluid_ReadParams(this, Ini)
    use IniObjects
    class(TAxionEffectiveFluid) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyModel%ReadParams(Ini)
    this%w_n  = Ini%Read_Double('AxionEffectiveFluid_w_n')
    this%om  = Ini%Read_Double('AxionEffectiveFluid_om')
    this%a_c  = Ini%Read_Double('AxionEffectiveFluid_a_c')
    call Ini%Read('AxionEffectiveFluid_theta_i', this%theta_i)

    end subroutine TAxionEffectiveFluid_ReadParams


    function TAxionEffectiveFluid_PythonClass()
    character(LEN=:), allocatable :: TAxionEffectiveFluid_PythonClass

    TAxionEffectiveFluid_PythonClass = 'AxionEffectiveFluid'
    end function TAxionEffectiveFluid_PythonClass

    subroutine TAxionEffectiveFluid_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TAxionEffectiveFluid), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TAxionEffectiveFluid_SelfPointer

    subroutine TAxionEffectiveFluid_Init(this, State)
    use classes
    class(TAxionEffectiveFluid), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State
    real(dl) :: grho_rad, F, p, mu, xc, n

    select type(State)
    class is (CAMBdata)
        this%is_cosmological_constant = this%om==0
        this%pow = 3*(1+this%w_n)
        this%omL = State%Omega_de - this%om !Omega_de is total dark energy density today
        this%acpow = this%a_c**this%pow
        this%num_perturb_equations = 2
        if (this%w_n < 0.9999) then
            ! n <> infinity
            !get (very) approximate result for sound speed parameter; arXiv:1806.10608  Eq 30 (but mu may not exactly agree with what they used)
            n = nint((1+this%w_n)/(1-this%w_n))
            !Assume radiation domination, standard neutrino model; H0 factors cancel
            grho_rad = (kappa/c**2*4*sigma_boltz/c**3*State%CP%tcmb**4*Mpc**2*(1+3.046*7._dl/8*(4._dl/11)**(4._dl/3)))
            xc = this%a_c**2/2/sqrt(grho_rad/3)
            F=7./8
            p=1./2
            mu = 1/xc*(1-cos(this%theta_i))**((1-n)/2.)*sqrt((1-F)*(6*p+2)*this%theta_i/n/sin(this%theta_i))
            this%freq =  mu*(1-cos(this%theta_i))**((n-1)/2.)* &
                sqrt(const_pi)*Gamma((n+1)/(2.*n))/Gamma(1+0.5/n)*2.**(-(n**2+1)/(2.*n))*3.**((1./n-1)/2)*this%a_c**(-6./(n+1)+3) &
                *( this%a_c**(6*n/(n+1.))+1)**(0.5*(1./n-1))
            this%n = n
        end if
    end select

    end subroutine TAxionEffectiveFluid_Init


    function TAxionEffectiveFluid_w_de(this, a)
    class(TAxionEffectiveFluid) :: this
    real(dl) :: TAxionEffectiveFluid_w_de
    real(dl), intent(IN) :: a
    real(dl) :: rho, apow, acpow

    apow = a**this%pow
    acpow = this%acpow
    rho = this%omL+ this%om*(1+acpow)/(apow+acpow)
    TAxionEffectiveFluid_w_de = this%om*(1+acpow)/(apow+acpow)**2*(1+this%w_n)*apow/rho - 1

    end function TAxionEffectiveFluid_w_de

    function TAxionEffectiveFluid_grho_de(this,grhok, grhoc, grhob, grhog, grhornomass, grhov, a)
!this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TAxionEffectiveFluid),intent(inout) :: this
    real(dl) :: TAxionEffectiveFluid_grho_de, apow
    real(dl), intent(IN) :: grhok, grhoc, grhob, grhog, grhornomass, grhov, a

    if(a == 0.d0)then
        TAxionEffectiveFluid_grho_de = 0.d0
    else
        apow = a**this%pow
        TAxionEffectiveFluid_grho_de = (this%omL*(apow+this%acpow)+this%om*(1+this%acpow))*a**4 &
            /((apow+this%acpow)*(this%omL+this%om))

!        write(86,*) a, TAxionEffectiveFluid_grho_de
    endif

    end function TAxionEffectiveFluid_grho_de

    subroutine TAxionEffectiveFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TAxionEffectiveFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv, apow, acpow, cs2, fac

    if (this%w_n < 0.9999) then
        fac = 2*a**(2-6*this%w_n)*this%freq**2
        cs2 = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
    else
        cs2 = 1
    end if
    apow = a**this%pow
    acpow = this%acpow
    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k

    ! dw/dlog a/(1+w)

    ! fix the axion fluid too

    deriv  = (acpow**2*(this%om+this%omL)+this%om*acpow-apow**2*this%omL)*this%pow &
        /((apow+acpow)*(this%omL*(apow+acpow)+this%om*(1+acpow)))


    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + &
        k * cs2 * y(w_ix)

    end subroutine TAxionEffectiveFluid_PerturbationEvolve


    subroutine TAxionEffectiveFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TAxionEffectiveFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TAxionEffectiveFluid_PerturbedStressEnergy

     !!! RH  modified below to include the delta

    subroutine TModDeltaEffectiveFluid_ReadParams(this, Ini)
      use IniObjects
      integer i
      class(TModDeltaEffectiveFluid) :: this
      class(TIniFile), intent(in) :: Ini
      
      call this%TDarkEnergyModel%ReadParams(Ini)
      
      do i=1,this%max_abin
         this%lnavals(i) = Ini%Read_Double_Array('ModDeltaEffectiveFluid_lnavals', i)
         this%deltavals(i) = Ini%Read_Double_Array('ModDeltaEffectiveFluid_deltavals', i)
      enddo
      this%smoothscale = Ini%Read_Double('ModDeltaEffectiveFluid_tau')


    end subroutine TModDeltaEffectiveFluid_ReadParams


    function TModDeltaEffectiveFluid_PythonClass()
      character(LEN=:), allocatable :: TModDeltaEffectiveFluid_PythonClass
      TModDeltaEffectiveFluid_PythonClass = 'ModDeltaEffectiveFluid'
    end function TModDeltaEffectiveFluid_PythonClass


    subroutine TModDeltaEffectiveFluid_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TModDeltaEffectiveFluid), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TModDeltaEffectiveFluid_SelfPointer

    subroutine TModDeltaEffectiveFluid_Init(this, State)
    use classes

    class(TModDeltaEffectiveFluid), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State
    real(dl) :: grho_rad, F, p, mu, xc, n

    select type(State)
    class is (CAMBdata)
        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 2
        this%grhornomass = State%grhornomass
        this%grhoc = State%grhoc
        this%grhob = State%grhob
        this%grhov = State%grhov
        this%grhog = State%grhog
        this%grhocrit = State%grhocrit
        this%grhom = State%grhoc + State%grhob
    end select

    end subroutine TModDeltaEffectiveFluid_Init


    function deltaStep(this,a)
      class(TModDeltaEffectiveFluid),intent(IN) :: this
      real(dl) ::  deltaStep
      real(dl), intent(IN) :: a
      real(dl) :: rho, apow, acpow, lnaval, dlndeltadlna
      real(dl) :: dx, splin_lndelta, lna, Q, ddelta_sum, d2delta_sum, wbg
      real(dl) :: bigtermA, bigtermB
      integer i, j,k
      
      Q = TModDeltaQ(this,a) 
      if(a .gt.  (this%grhom/this%grhov)**(1/3.)) then
         wbg=-1
      else if( (a .lt.  (this%grhom/this%grhov)**(1/3.)) .and. (a .gt. (this%grhornomass+this%grhog)/this%grhom)) then
         wbg=0
      else
         wbg=1/3
      end if

      ! if u is ouside the x() interval take a boundary value (left or right)
      deltaStep= 0
      bigtermA = 0
      bigtermB = 0
      do i =1,this%max_abin-1
         bigtermA = 1 + (a/exp(this%lnavals(i+1)))**(1/this%smoothscale)
         bigtermB = 1 + (a/exp(this%lnavals(i)))**(1/this%smoothscale)
         deltaStep = deltaStep+ this%deltavals(i)*(1./bigtermA -1./bigtermB)
      enddo
      

    end function deltaStep


    function deltaStepFirstDeriv(this,a)
      class(TModDeltaEffectiveFluid),intent(IN) :: this
      real(dl) :: deltaStepFirstDeriv, deltaStep
      real(dl), intent(IN) :: a
      real(dl) :: rho, apow, acpow, lnaval, dlndeltadlna
      real(dl) :: dx, splin_lndelta, lna, Q, ddelta_sum, d2delta_sum, wbg
      real(dl) :: bracketA, bracketB
      integer i, j,k
      
      Q = TModDeltaQ(this,a) !(this%grhom+this%grhornomass)/this%grhov
      if(a .gt.  (this%grhom/this%grhov)**(1/3.)) then
         wbg=-1
      else if( (a .lt.  (this%grhom/this%grhov)**(1/3.)) .and. (a .gt. (this%grhornomass+this%grhog)/this%grhom) ) then
         wbg=0
      else
         wbg=1/3
      end if

      ! if u is ouside the x() interval take a boundary value (left or right)
      deltaStepFirstDeriv= 0
      bracketA = 0
      bracketB = 0
      do i =1,this%max_abin-1
         bracketA = a/exp(this%lnavals(i))
         bracketB = a/exp(this%lnavals(i+1))
         
         deltaStepFirstDeriv= deltaStepFirstDeriv +  a*( bracketA**(-1+1/this%smoothscale)/( (1+bracketA**(1/this%smoothscale))**2*exp(this%lnavals(i))*this%smoothscale) - bracketB**(-1+1/this%smoothscale)/((1+bracketB**(1/this%smoothscale) )**2*exp(this%lnavals(i+1))*this%smoothscale)   )
         
         
      enddo

    end function deltaStepFirstDeriv
    

    function deltaStepSecondDeriv(this,a)

      class(TModDeltaEffectiveFluid),intent(IN) :: this
      real(dl) :: deltaStepSecondDeriv
      real(dl), intent(IN) :: a
      real(dl) :: rho, apow, acpow, lnaval, dlndeltadlna
      real(dl) :: dx, splin_lndelta, lna, Q, ddelta_sum, d2delta_sum, wbg
      real(dl) :: bracketA, bracketB
      integer i, j,k
      
      Q = (this%grhom+this%grhornomass)/this%grhov
      if(a .gt.  (this%grhom/this%grhov)**(1/3.)) then
         wbg=-1
      else if( (a .lt.  (this%grhom/this%grhov)**(1/3.)) .and. (a .gt. (this%grhornomass+this%grhog)/this%grhom) ) then
         wbg=0
      else
         wbg=1/3
      end if
      
      ! if u is ouside the x() interval take a boundary value (left or right)
      deltaStepSecondDeriv= 0
      bracketA = 0
      bracketB = 0
      do i =1,this%max_abin-1
         bracketA = a/exp(this%lnavals(i))
         bracketB = a/exp(this%lnavals(i+1))
         
         deltaStepSecondDeriv= deltaStepSecondDeriv + a*(a*(-2*bracketA**(-2+2/this%smoothscale)/((1+bracketA**(1/this%smoothscale))**3*exp(this%lnavals(i))**2*this%smoothscale**2 ) + 2*bracketB**(-2+2/this%smoothscale)/((1+bracketB**(1/this%smoothscale))**3*exp(this%lnavals(i+1))**2*this%smoothscale**2) + bracketA**(-2+1/this%smoothscale)*(-1+1/this%smoothscale)/((1+bracketA**(1/this%smoothscale))**2*exp(this%lnavals(i))*2*this%smoothscale) - bracketB**(-2+(1/this%smoothscale))*(-1+1/this%smoothscale)/((1+bracketB**(1/this%smoothscale))**2*exp(this%lnavals(i+1))**2*this%smoothscale)  ) + bracketA**(-1+(1/this%smoothscale))/((1+bracketA**(1/this%smoothscale))**2*exp(this%lnavals(i))*this%smoothscale )  -bracketB**(-1+(1/this%smoothscale))/((1+bracketB**(1/this%smoothscale) )**2*exp(this%lnavals(i+1))*this%smoothscale)   )
                  
      enddo

    end function deltaStepSecondDeriv
    
    function TModDeltaEffectiveFluid_w_de(this, a)
    class(TModDeltaEffectiveFluid) :: this
    real(dl) :: TModDeltaEffectiveFluid_w_de
    real(dl), intent(IN) :: a
    real(dl) :: rho, apow, acpow, lnaval, dlndeltadlna, dx, splin_lndelta, delta_sum, lna, Q, ddelta_sum, d2delta_sum, wbg
    real(dl) :: bigtermA, bigtermB
    integer i, j,k
    lnaval = dlog(a)
    
    ! need to fix this - not radiation included
    Q = (this%grhom+this%grhornomass)/this%grhov
    ! need to set the bg equation of state
    !   write(*,*) ' a, (this%grhom/this%grhov)**(1/3), this%grhornomass/this%grhom
    
    if(a .gt.  (this%grhom/this%grhov)**(1/3.)) then 
       wbg=-1
    else if( (a .lt.  (this%grhom/this%grhov)**(1/3.)) .and. (a .gt. (this%grhornomass+this%grhog)/this%grhom)) then 
       wbg=0
    else 
       wbg=1/3
    end if
    

    delta_sum=deltaStep(this,a)    
    ddelta_sum=deltaStepFirstDeriv(this,a)    
    d2delta_sum=deltaStepSecondDeriv(this,a)    
    TModDeltaEffectiveFluid_w_de=  (Q*delta_sum/(1+delta_sum*(1+Q)))*(1+wbg) - (1/3.)*((1+Q)/(1+delta_sum*(1+Q)))*ddelta_sum - 1
    
    write(99,*) a, delta_sum, ddelta_sum, d2delta_sum, TModDeltaEffectiveFluid_w_de
  end function TModDeltaEffectiveFluid_w_de

    function TModDeltaEffectiveFluid_grho_de(this,grhok, grhoc, grhob, grhog, grhornomass, grhov, a)  !relative density (8 pi G a^4 rho_de /grhov)
      class(TModDeltaEffectiveFluid),intent(inout) :: this
!    class(TModDeltaEffectiveFluid) :: this
    real(dl) :: TModDeltaEffectiveFluid_grho_de, delta_sum, omega_de_bg, grho_fid, modH2, TModDeltaEffectiveFluid_grho_de1, H2
    real(dl), intent(IN) :: grhok, grhoc, grhob, grhog, grhornomass, grhov, a 
    integer i

    if(a < 1e-10) then
       TModDeltaEffectiveFluid_grho_de = 0.d0
    else

       ! grho_v = grhocrit*omega_de

       delta_sum=deltaStep(this,a) ! need to fix the step

       grho_fid = grhok/(a**2)+ (grhoc + grhob)/(a**3) + (grhog + grhornomass)/(a**4) +  grhov

       ! equation (8) in 1208.4845
       TModDeltaEffectiveFluid_grho_de = (grhov + delta_sum*grho_fid)/grhov !normalizing as fraction 

       TModDeltaEffectiveFluid_grho_de = a**4*TModDeltaEffectiveFluid_grho_de !sending as a^4grho

! a*a*a*a*(omega_de_bg+delta_sum)/(1+delta_sum)
!       write(88,*) a, delta_sum, TModDeltaEffectiveFluid_grho_de, grho_fid, grhov*a**4

    endif

    end function TModDeltaEffectiveFluid_grho_de


    function TModDeltaEffectiveFluid_TotalDensity(this,grhok, grhoc, grhob, grhog, grhornomass, grhov, a)
      class(TModDeltaEffectiveFluid),intent(inout) :: this
      real(dl) :: TModDeltaEffectiveFluid_TotalDensity, delta_sum
      real(dl), intent(IN) :: grhok, grhoc, grhob, grhog, grhornomass, grhov, a
      real(dl) :: grhov_de
      integer i

         delta_sum=deltaStep(this,a)
         ! equation (1) in 1304.3724
         grhov_de = grhov * this%grho_de(grhok, grhoc, grhob, grhog, grhornomass, grhov, a) / (a **4)
         ! 8*pi*G*rho*a**4
         TModDeltaEffectiveFluid_TotalDensity =  grhok/(a**2)+ (grhoc + grhob)/(a**3) + (grhog + grhornomass)/(a**4) +  grhov_de

         TModDeltaEffectiveFluid_TotalDensity = TModDeltaEffectiveFluid_TotalDensity*(a**4)


    end function TModDeltaEffectiveFluid_TotalDensity


    function TModDeltaQ(this, a)  
    class(TModDeltaEffectiveFluid) :: this
    real(dl) :: TModDeltaQ
    real(dl), intent(IN) :: a
    integer i


    TModDeltaQ = (this%grhom/a**3+this%grhornomass/a**4+this%grhog/a**4)/(this%grhov)

    end function TModDeltaQ


    function TModDeltaQderiv(this, a)  
    class(TModDeltaEffectiveFluid) :: this
    real(dl) :: TModDeltaQderiv, TModDeltaQm, TModDeltaQp, am, ap, step
    real(dl), intent(IN) :: a
    integer i

    step=1e-5
    am=a*(1-step)
    ap=a*(1+step)

    TModDeltaQm = TModDeltaQ(this,am)
    TModDeltaQp = TModDeltaQ(this,ap)

    TModDeltaQderiv=a*(TModDeltaQp-TModDeltaQm)/(2*step) ! logarithmic derivative
    
    end function TModDeltaQderiv

    subroutine TModDeltaEffectiveFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TModDeltaEffectiveFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv, cs2, Q,Qderiv, fac, delta_sum, ddelta_sum, d2delta_sum, grhov, wbg

    if(a .gt.  (this%grhom/this%grhov)**(1/3.)) then
       wbg=-1
    else if( (a .lt.  (this%grhom/this%grhov)**(1/3.)) .and. (a .gt. (this%grhornomass+this%grhog)/this%grhom) ) then
       wbg=0
    else
       wbg=1/3
    end if

    cs2=1


    Q = TModDeltaQ(this,a)
    Qderiv = TModDeltaQDeriv(this,a)

    delta_sum=deltaStep(this,a)
    ddelta_sum=deltaStepFirstDeriv(this,a)
    d2delta_sum=deltaStepSecondDeriv(this,a)

    ! dw/dlog a
    deriv = (1+wbg)*(Qderiv*delta_sum/(1+delta_sum*(1+Q)) + Q*ddelta_sum*(1/(1+delta_sum*(1+Q))) -Q*delta_sum/(1+delta_sum*(1+Q))**2*(ddelta_sum*(1+Q)+ delta_sum*Qderiv))  -(1/3.)*(Qderiv*(1/(1+delta_sum*(1+Q)))*ddelta_sum  -(1+Q)/(1+delta_sum*(1+Q))**2*(ddelta_sum*(1+Q) + delta_sum*Qderiv )*ddelta_sum + ((1+Q)/(1+delta_sum*(1+Q)))*d2delta_sum ) 

    ! RH hardcoded this to test
!    deriv=0


    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + &
        k * cs2 * y(w_ix)

 !   write(*,*) Q, Qderiv, a, w, delta_sum, ddelta_sum, d2delta_sum, deriv, ayprime(w_ix), ayprime(w_ix + 1)
 !   write(*,*) 'Q, Qderiv, a, w, delta_sum, ddelta_sum, d2delta_sum, deriv, ayprime(w_ix), ayprime(w_ix + 1)'


    ! RH hardcoded 
!    ayprime(w_ix + 1) = 0
!    ayprime(w_ix) = 0
    end subroutine TModDeltaEffectiveFluid_PerturbationEvolve


    subroutine TModDeltaEffectiveFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TModDeltaEffectiveFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TModDeltaEffectiveFluid_PerturbedStressEnergy


    end module DarkEnergyFluid
