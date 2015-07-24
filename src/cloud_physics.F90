module cloud_physics

    use ui
    use advection_solver
    use netcdf

    implicit none

    private

    public cloud_physics_init
    public cloud_physics_run
    public cloud_physics_final
    public cloud_physics_output

    ! --------------------------------------------------------------------------
    !                               Constants
    ! --------------------------------------------------------------------------
    real(8), parameter :: kb = 1.3806488d-23 ! Boltzmann's constant [J K-1].
    real(8), parameter :: NA = 6.0221413d23 ! Avogadro's number [mol-1].
    real(8), parameter :: Ma = 28.9644d0 ! Dry air molecular weight [g mol-1].
    real(8), parameter :: Mw = 18.0160d0 ! Water molecular weight [g mol-1].
    real(8), parameter :: Rg = kb*NA ! Universal gas constant [J K-1 mol-1]
    real(8), parameter :: Rd = Rg/Ma ! Dry air constant [J g-1 K-1].
    real(8), parameter :: Rv = Rg/Mw ! Water vapor constant [J g-1 K-1].
    real(8), parameter :: T_freeze = 273.15d0 ! Freezing point temperature of water [K].
    real(8), parameter :: T_triple = 273.16d0 ! Triple point temperature of water [K].
    real(8), parameter :: Rho_water = 1.0d0 ! Water density [g cm-3].

    ! "Constants" that is approximatively constant within an ordinary situation.
    real(8), parameter :: cpv = 1.870d0 ! Specific heat capacity of water vapor at constant pressure [J g-1 K-1].
    real(8), parameter :: cvv = 1.410d0 ! Specific heat capacity of water vapor at constant volume [J g-1 K-1].

    real(8), parameter :: PI = atan(1.0d0)*4.0d0 ! Just pie.

    ! Unit conversion constants.
    real(8), parameter :: cal_to_J = 4.1868d0

    ! --------------------------------------------------------------------------
    !                          Enumerate constants
    ! --------------------------------------------------------------------------
    integer, parameter :: saturation_vapor_pressure_formula_first_approximation = 1
    integer, parameter :: saturation_vapor_pressure_formula_more_accurate = 2
    integer, parameter :: saturation_vapor_pressure_formula_unknown = 3

    integer, parameter :: vaporization_latent_heat_formula_unknown = 1

    integer, parameter :: vapor_molecular_diffusivity_formula_hall_and_pruppacher_1976 = 1

    integer, parameter :: dry_air_thermal_conductivity_formula_beard_pruppacher_1971 = 1

    integer, parameter :: vapor_thermal_conductivity_formula_beard_pruppacher_1971 = 1

    integer, parameter :: moist_air_thermal_conductivity_formula_mason_sexena = 1

    ! --------------------------------------------------------------------------
    !                              Variables
    ! --------------------------------------------------------------------------
    real(8), allocatable :: r(:)
    real(8), allocatable :: dr(:)
    real(8), allocatable :: f(:)

    ! --------------------------------------------------------------------------
    !                           Control Parameters
    ! --------------------------------------------------------------------------
    character(30) advection_scheme
    integer saturation_vapor_pressure_formula
    integer vaporization_latent_heat_formula
    integer vapor_molecular_diffusivity_formula
    integer dry_air_thermal_conductivity_formula
    integer vapor_thermal_conductivity_formula
    integer moist_air_thermal_conductivity_formula

    integer num_bin
    real(8) time_step_size

    namelist /cloud_physics_params/ &
        advection_scheme, num_bin, time_step_size, &
        saturation_vapor_pressure_formula, &
        vaporization_latent_heat_formula, &
        vapor_molecular_diffusivity_formula, &
        dry_air_thermal_conductivity_formula, &
        vapor_thermal_conductivity_formula, &
        moist_air_thermal_conductivity_formula

contains

    ! --------------------------------------------------------------------------
    !                            Public Interfaces
    ! --------------------------------------------------------------------------

    subroutine cloud_physics_init(namelist_file_path, T, p, qv, qc)

        character(*), intent(in) :: namelist_file_path
        real(8), intent(in) :: T  ! Ambient temperature [K].
        real(8), intent(in) :: p  ! Ambient pressure [hPa].
        real(8), intent(in) :: qv ! Ambient vapor mixing ratio [g g-1].
        real(8), intent(in) :: qc ! Droplet mixing ratio [g g-1].

        integer i

        open(10, file=namelist_file_path)
        read(10, nml=cloud_physics_params)
        close(10)

        allocate(r(num_bin))
        allocate(dr(num_bin-1))
        allocate(f(num_bin))

        call set_initial_bins(num_bin, r)
        call set_initial_droplet_spectrum(num_bin, qc, Rho_water, r, f)

        do i = 1, num_bin-1
            dr(i) = r(i+1)-r(i)
        end do

        call advection_solver_init(num_bin, r, dr, f)

    end subroutine cloud_physics_init

    subroutine cloud_physics_run(T, p, qv, qc)

        real(8), intent(in) :: T  ! Ambient temperature [K].
        real(8), intent(in) :: p  ! Ambient pressure [hPa].
        real(8), intent(in) :: qv ! Ambient vapor mixing ratio [g g-1].
        real(8), intent(in) :: qc ! Droplet mixing ratio [g g-1].

        real(8) es ! Ambient saturation vapor pressure [hPa].
        real(8) S  ! Ambient saturation ratio [1].

        es = calc_saturation_vapor_pressure(T)
        S = qv/calc_saturation_vapor_mixing_ratio(p, es)

        call advection_solver_run(T, p, qv, es, S, time_step_size, calc_condensation_growth_rate)

    end subroutine cloud_physics_run

    subroutine cloud_physics_final()

        deallocate(r)
        deallocate(dr)
        deallocate(f)

    end subroutine cloud_physics_final

    subroutine cloud_physics_output(file_path)

        character(*), intent(in) :: file_path

        integer ncid, ierr

        ierr = nf90_create(file_path, nf90_clobber, ncid)

        ierr = nf90_put_att(ncid, nf90_global, "time_step_size", time_step_size)

        ierr = nf90_enddef(ncid)

        call advection_solver_output(ncid)

        ierr = nf90_close(ncid)

    end subroutine cloud_physics_output

    ! --------------------------------------------------------------------------
    !                           Private Interfaces
    ! --------------------------------------------------------------------------

    subroutine diagnose(T, p, qv, r)

        real(8), intent(in) :: T  ! Temperature [K].
        real(8), intent(in) :: p  ! Pressure [hPa].
        real(8), intent(in) :: qv ! Vapor mixing ratio [1].
        real(8), intent(in) :: r  ! Droplet radius [cm].

        real(8) es
        real(8) qvs
        real(8) S

        es = calc_saturation_vapor_pressure(T)
        qvs = calc_saturation_vapor_mixing_ratio(p, es)
        S = qv/qvs

        print *, "es =", es
        print *, "qv =", qv
        print *, "S  =", S
        print *, "Lv =", calc_vaporization_latent_heat(T)
        print *, "Ka =", calc_dry_air_thermal_conductivity(T)
        print *, "Kv =", calc_vapor_thermal_conductivity(T)
        print *, "Dv =", calc_vapor_molecular_diffusivity(T, p, r)

    end subroutine diagnose

    function calc_saturation_vapor_pressure(T) result(es)

        real(8), intent(in) :: T ! Temperature [K].
        real(8) es ! Saturation vapor pressure [hPa].

        ! The stauration vapor pressure is calculated from Clausius-Clapeyron
        ! equation.
        select case (saturation_vapor_pressure_formula)
        case (saturation_vapor_pressure_formula_first_approximation)
            es = 2.53d9*exp(-5.42d3/T)
        case (saturation_vapor_pressure_formula_more_accurate)
            ! TODO: Figure out the formula for es in this branch.
            call report_error("This saturation vapor formula is not available!", &
                __FILE__, __LINE__)
        case (saturation_vapor_pressure_formula_unknown)
            es = 6.107d0*exp(17.15d0*(T-T_freeze)/(T-38.33d0))
        case default
            call report_error("Invalid saturation_vapor_pressure_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_saturation_vapor_pressure

    function calc_saturation_vapor_mixing_ratio(p, es) result(qvs)

        real(8), intent(in) :: p  ! Ambient pressure [hPa].
        real(8), intent(in) :: es ! Saturation vapor pressure [hPa].
        real(8) qvs ! Saturation vapor mixing ratio [g g-1].

        qvs = 0.622d0*es/(p-0.378*es)

    end function calc_saturation_vapor_mixing_ratio

    function calc_vaporization_latent_heat(T) result(Lv)

        real(8), intent(in) :: T ! Temperature [K].
        real(8) Lv ! Vaporization latent heat [J g-1].

        select case (vaporization_latent_heat_formula)
        case (vaporization_latent_heat_formula_unknown)
            Lv = 2.5d3-(T-T_freeze)*2.371d0
        case default
            call report_error("Invalid vaporization_latent_heat_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_vaporization_latent_heat

    subroutine calc_condensation_growth_rate(T, p, es, S, r, drdt, div_drdt)

        real(8), intent(in) :: T  ! Ambient temperature [K].
        real(8), intent(in) :: p  ! Ambient pressure [hPa].
        real(8), intent(in) :: es ! Ambient saturation vapor pressure [hPa].
        real(8), intent(in) :: S  ! Ambient saturation ratio [1].
        real(8), intent(in) :: r  ! Droplet radius [cm].
        real(8), intent(out) :: drdt ! Condensation growth rate [cm s-1].
        real(8), intent(out) :: div_drdt ! Condensation growth rate divergence [s-1].

        real(8) Lv ! Vaporization latent heat [J g-1].
        real(8) K  ! Air thermal conductivity [J s-1 cm-1 K-1].
        real(8) Dv ! Vapor molecular diffusivity [cm2 s-1].
        real(8) Fk ! Heat conduction term.
        real(8) Fd ! Vapor diffusivity term.

        Lv = calc_vaporization_latent_heat(T)
        K  = calc_dry_air_thermal_conductivity(T)
        Dv = calc_vapor_molecular_diffusivity(T, p, r)
        Fk = (Lv/(Rv*T)-1)*(Lv*Rho_water)/(K*T)
        Fd = (Rho_water*Rv*T)/(Dv*es)*1.0d4

        drdt = (S-1.0d0)/(Fk+Fd)/r
        div_drdt = -drdt/r

    end subroutine calc_condensation_growth_rate

    function calc_vapor_molecular_diffusivity(T, p, r) result(Dv)

        real(8), intent(in) :: T  ! Ambient temperature [K].
        real(8), intent(in) :: p  ! Ambient pressure [Pa].
        real(8), intent(in) :: r  ! Droplet radius [cm].
        real(8) Dv ! Vapor molecular diffusivity [cm2 s-1]

        real(8) p0 ! Reference pressure [hPa].
        real(8) alpha_c ! TODO: What is this?
        real(8) delta_v ! Water vapor jump length [cm].
        real(8) Tr ! Temperature at drop surface [K].

        select case (vapor_molecular_diffusivity_formula)
        case (vapor_molecular_diffusivity_formula_hall_and_pruppacher_1976)
#ifdef DEBUG
            if (T > T_freeze+40.0d0 .or. T < -T_freeze-40.0d0) then
                call report_error("Temperature is out of range [-40,40] degree!", &
                    __FILE__, __LINE__)
            end if
#endif
            p0 = 1013.25d0
            Dv = 0.211d0*(T/T_freeze)**1.94d0*p0/p
            alpha_c = 0.04d0 ! TODO: What is this?
            delta_v = 1.3*6.6d-6*1013.25d0/293.15*T/p ! TODO: Where does this come from?
            Tr = T ! TODO: Is this OK?
            Dv = Dv/(r/(r+delta_v)+Dv/(r*alpha_c)*(2.0d0*PI/Rv/Tr)**0.5d0)
        case default
            call report_error("Invalid vapor_molecular_diffusivity_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_vapor_molecular_diffusivity

    function calc_dry_air_thermal_conductivity(T) result(Ka)

        real(8), intent(in) :: T ! Temperature [K].
        real(8) Ka ! Dry air thermal conductivity [W cm-1 K-1].

        select case (dry_air_thermal_conductivity_formula)
        case (dry_air_thermal_conductivity_formula_beard_pruppacher_1971)
            Ka = (5.69d0+0.017d0*(T-T_freeze))*1.0d-5*cal_to_J
        case default
            call report_error("Invalid dry_air_thermal_conductivity_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_dry_air_thermal_conductivity

    function calc_vapor_thermal_conductivity(T) result(Kv)

        real(8), intent(in) :: T ! Temperature [K].
        real(8) Kv ! Water vapor thermal conductivity [W cm-1 K-1].

        select case (vapor_thermal_conductivity_formula)
        case (vapor_thermal_conductivity_formula_beard_pruppacher_1971)
            Kv = (3.78d0+0.020d0*T)*1.0d-5*cal_to_J
        case default
            call report_error("Invalid vapor_thermal_conductivity_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_vapor_thermal_conductivity

    function calc_moist_air_thermal_conductivity(T, Ka, Kv) result(Km)

        real(8), intent(in) :: T ! Temperature [K].
        real(8), intent(in), optional :: Ka ! Dry air thermal conductivity [W m-1 K-1].
        real(8), intent(in), optional :: Kv ! Water vapor thermal conductivity [W m-1 K-1].
        real(8) Km ! Moist air thermal conductivity [W cm-1 K-1].

        real(8), parameter :: gamma1 = 1.17d0
        real(8), parameter :: gamma2 = 1.02d0
        real(8) Ka_, Kv_
        real(8) xv ! Mole fraction of water vapor [1]. TODO: How to set this quantity?

        if (present(Ka) .and. present(Kv)) then
            Ka_ = Ka
            Kv_ = Kv
        else
            Ka_ = calc_dry_air_thermal_conductivity(T)
            Kv_ = calc_vapor_thermal_conductivity(T)
        end if

        select case (moist_air_thermal_conductivity_formula)
        case (moist_air_thermal_conductivity_formula_mason_sexena)
            Km = Ka_*(1.0d0-(gamma1-gamma2*Kv_/Ka_)*xv)
        case default
            call report_error("Invalid moist_air_thermal_conductivity_formula!", &
                __FILE__, __LINE__)
        end select

    end function calc_moist_air_thermal_conductivity

    subroutine set_initial_bins(num_bin, r)

        integer, intent(in) :: num_bin
        real(8), intent(out) :: r(num_bin) ! [cm]

        real(8), parameter :: a = 2.0d0**0.5d0
        real(8), parameter :: b = 4.0d0/3.0d0*PI
        real(8) m(num_bin) ! Droplet mass per bin [g].
        integer i

        m(1) = 4.0d-15
        r(1) = (m(1)/b)**(1.0d0/3.0d0)
        do i = 2, num_bin
            m(i) = m(i-1)*a
            r(i) = (m(i)/b)**(1.0d0/3.0d0)
        end do

    end subroutine set_initial_bins

    subroutine set_initial_droplet_spectrum(num_bin, q, rho, r, f)

        integer, intent(in) :: num_bin
        real(8), intent(in) :: q   ! Mixing ratio [g g-1]
        real(8), intent(in) :: rho ! Droplet density [g cm-3]
        real(8), intent(in) :: r(num_bin)
        real(8), intent(out) :: f(num_bin)

        real(8), parameter :: Nt = 200.0d0 ! Droplet number concentration [cm-3].
        real(8), parameter :: alpha = 2.0d0 ! Shape parameter.
        real(8), parameter :: gamma1 = 1.999964d0
        real(8), parameter :: gamma2 = 119.9964d0
        real(8) lambda
        real(8) N0
        integer i

        lambda = (gamma2*Nt*PI/6.0d0/(gamma1*1.0d-3*rho*q))**(1.0d0/3.0d0)
        N0 = Nt*lambda**(1+alpha)/gamma1
        do i = 1, num_bin
            f(i) = N0*r(i)**3*exp(-lambda*r(i))*dlog(2.0d0)/2.0d0/3.0d0
            if (f(i) < 1.0d-30) f(i) = 0.0d0
        end do

    end subroutine set_initial_droplet_spectrum

end module cloud_physics
