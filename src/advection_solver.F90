module advection_solver

    use lagrangian_solver

    implicit none

    private

    public advection_solver_init
    public advection_solver_run
    public advection_solver_final
    public advection_solver_output

contains

    ! --------------------------------------------------------------------------
    !                            Public Interfaces
    ! --------------------------------------------------------------------------

    subroutine advection_solver_init(num_bin, r, dr, f)

        integer, intent(in) :: num_bin ! Initial discretized bin number.
        real(8), intent(in) :: r(num_bin) ! Initial droplet radius left side [cm].
        real(8), intent(in) :: dr(num_bin) ! Initial bin interval [cm].
        real(8), intent(in) :: f(num_bin) ! Initial drop-size distribution [per unit air mass?].

        call lagrangian_solver_init(num_bin, r, dr, f)

    end subroutine advection_solver_init

    subroutine advection_solver_run(T, p, qv, es, S, dt, growth_rate)

        real(8), intent(in) :: T  ! Ambient temperature [K].
        real(8), intent(in) :: p  ! Ambient pressure [hPa].
        real(8), intent(in) :: qv ! Ambient saturation vapor mixing ratio [g g-1].
        real(8), intent(in) :: es ! Ambient saturation vapor pressure [hPa].
        real(8), intent(in) :: S  ! Ambient saturation ratio [1].
        real(8), intent(in) :: dt ! Time step size [s].

        interface
            subroutine growth_rate(T, p, es, S, r, drdt, div_drdt)
                real(8), intent(in) :: T  ! Ambient temperature [K].
                real(8), intent(in) :: p  ! Ambient pressure [hPa].
                real(8), intent(in) :: es ! Ambient saturation vapor pressure [hPa].
                real(8), intent(in) :: S  ! Ambient saturation ratio [1].
                real(8), intent(in) :: r  ! Droplet radius [cm].
                real(8), intent(out) :: drdt ! Droplet growth rate [cm s-1].
                real(8), intent(out) :: div_drdt ! Droplet growth rate divergence [s-1].
            end subroutine growth_rate
        end interface

        call lagrangian_solver_run(T, p, qv, es, S, dt, growth_rate)

    end subroutine advection_solver_run

    subroutine advection_solver_final()

    end subroutine advection_solver_final

    subroutine advection_solver_output(ncid)

        integer, intent(in) :: ncid

        call lagrangian_solver_output(ncid)

    end subroutine advection_solver_output

end module advection_solver
