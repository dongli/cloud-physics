module advection_solver

#ifdef CLOUD_MICROPHYSICS_USE_LAGRANGIAN_SOLVER
    use lagrangian_solver
#endif

    implicit none

contains

    ! --------------------------------------------------------------------------
    !                            Public Interfaces
    ! --------------------------------------------------------------------------

    subroutine advection_solver_init(...)

    end subroutine advection_solver_init

    subroutine advection_solver_run(...)

    end subroutine advection_solver_run

    subroutine advection_solver_final()

    end subroutine advection_solver_final

end module advection_solver
