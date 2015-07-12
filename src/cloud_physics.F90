module cloud_physics

    use advection_solver

    implicit none

    private

    public cloud_physics_init
    public cloud_physics_run
    public cloud_physics_final

    ! --------------------------------------------------------------------------
    !                               Constants
    ! --------------------------------------------------------------------------
    real(8), parameter :: Rd = ...
    real(8), parameter :: Cp = ...

    ! --------------------------------------------------------------------------
    !                           Control Parameters
    ! --------------------------------------------------------------------------
    character(30) advection_scheme
    character(30) saturated_vapor_pressure_formula

    namelist /cloud_physics/ &
        advection_scheme, saturated_vapor_pressure_formula

contains

    ! --------------------------------------------------------------------------
    !                            Public Interfaces
    ! --------------------------------------------------------------------------

    subroutine cloud_physics_init(...)

    end subroutine cloud_physics

    subroutine cloud_physics_run(...)

    end subroutine cloud_physics_run

    subroutine cloud_physics_final()

    end subroutine cloud_physics_final

    ! --------------------------------------------------------------------------
    !                           Private Interfaces
    ! --------------------------------------------------------------------------

    subroutine calc_saturated_vapor_pressure(...)

    end subroutine calc_saturated_vapor_pressure

    subroutine calc_condensation_rate

    end subroutine calc_condensation_rate

end module cloud_physics
