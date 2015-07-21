module lagrangian_solver

    use netcdf

    implicit none

    private

    public lagrangian_solver_init
    public lagrangian_solver_run
    public lagrangian_solver_final
    public lagrangian_solver_output

    integer, parameter :: max_num_bin = 1000

    type bin
        real(8) r
        real(8) dr
        real(8) f
    end type bin

    integer num_bin
    type(bin), allocatable :: bins(:)

contains

    ! --------------------------------------------------------------------------
    !                            Public Interfaces
    ! --------------------------------------------------------------------------

    subroutine lagrangian_solver_init(num_bin_, r, dr, f)

        integer, intent(in) :: num_bin_ ! Initial discretized bin number.
        real(8), intent(in) :: r(num_bin) ! Initial droplet radius left side [cm].
        real(8), intent(in) :: dr(num_bin) ! Initial bin interval [cm].
        real(8), intent(in) :: f(num_bin) ! Initial drop-size distribution [per unit air mass?].

        integer i

        allocate(bins(max_num_bin))

        num_bin = 0
        do i = 1, num_bin_-1
            if (f(i) > 0) then
                bins(i)%r = r(i)+dr(i)*0.5
                bins(i)%dr = dr(i)
                bins(i)%f = f(i)
                num_bin = num_bin+1
            end if
        end do

    end subroutine lagrangian_solver_init

    subroutine lagrangian_solver_run(T, p, qv, es, S, dt, growth_rate)

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

        real(8) drdt1, k1
        real(8) drdt2, k2
        real(8) drdt3, k3
        real(8) drdt4, k4
        real(8) r0, dr0
        real(8) r1, dr1
        real(8) div_drdt
        real(8) total_droplet_num
        integer i

        total_droplet_num = 0
        do i = 1, num_bin
            r0 = bins(i)%r
            dr0 = bins(i)%dr
            ! Update central droplet radius r.
            ! Update bin interval dr.
            call growth_rate(T, p, es, S, r0, drdt1, div_drdt)
            r1 = r0+drdt1*dt*0.5d0
            k1 = dr0*div_drdt
            dr1 = dr0+k1*dt*0.5d0
            call growth_rate(T, p, es, S, r0, drdt2, div_drdt)
            r1 = r0+drdt2*dt*0.05d0
            k2 = dr1*div_drdt
            dr1 = dr0+k2*dt*0.5d0
            call growth_rate(T, p, es, S, r0, drdt3, div_drdt)
            r1 = r0+drdt3*dt
            k3 = dr1*div_drdt
            dr1 = dr0+k3*dt
            call growth_rate(T, p, es, S, r0, drdt4, div_drdt)
            bins(i)%r = r0+(drdt1+2.0d0*drdt2+2.0d0*drdt3+drdt4)/6.0d0*dt
            k4 = dr1*div_drdt
            bins(i)%dr = dr0+(k1+2.0d0*k2+2.0d0*k3+k4)/6.0d0*dt
            bins(i)%f = bins(i)%f*dr0/bins(i)%dr
            total_droplet_num = total_droplet_num+bins(i)%dr*bins(i)%f
        end do
        print "(F30.15)", total_droplet_num
        print *, bins(1)%r, bins(2)%r

    end subroutine lagrangian_solver_run

    subroutine lagrangian_solver_final()

        deallocate(bins)

    end subroutine lagrangian_solver_final

    subroutine lagrangian_solver_output(ncid)

        integer, intent(in) :: ncid

        integer num_bin_dimid, r_varid, dr_varid, f_varid
        integer i, ierr

        real(8), allocatable :: var(:)

        allocate(var(num_bin))

        ierr = nf90_redef(ncid)

        ierr = nf90_def_dim(ncid, "num_bin", num_bin, num_bin_dimid)

        ierr = nf90_def_var(ncid, "r", nf90_double, [num_bin_dimid], r_varid)

        ierr = nf90_def_var(ncid, "dr", nf90_double, [num_bin_dimid], dr_varid)

        ierr = nf90_def_var(ncid, "f", nf90_double, [num_bin_dimid], f_varid)

        ierr = nf90_enddef(ncid)

        do i = 1, num_bin
            var(i) = bins(i)%r
        end do
        ierr = nf90_put_var(ncid, r_varid, var)

        do i = 1, num_bin
            var(i) = bins(i)%dr
        end do
        ierr = nf90_put_var(ncid, dr_varid, var)

        do i = 1, num_bin
            var(i) = bins(i)%f
        end do
        ierr = nf90_put_var(ncid, f_varid, var)

        deallocate(var)

    end subroutine lagrangian_solver_output

end module lagrangian_solver
