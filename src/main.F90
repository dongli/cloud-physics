program main

    use cloud_physics

    implicit none

    ! Test variables.
    real(8) :: T = 283.0d0                ! [K]
    real(8) :: p = 1000.0d0               ! [hPa]
    real(8) :: qv = 1.104969761952449e-2  ! [g g-1]
    real(8) :: qc = 1.6d-6                ! [g g-1]

    character(30) file_path
    integer i

    call cloud_physics_init("./namelist", T, p, qv, qc)

    write(file_path, "('output.', I5.5, '.nc')") 0
    call cloud_physics_output(file_path)
    do i = 1, 6000
        call cloud_physics_run(T, p, qv)
        write(file_path, "('output.', I5.5, '.nc')") i
        call cloud_physics_output(file_path)
    end do

    call cloud_physics_final()

end program main
