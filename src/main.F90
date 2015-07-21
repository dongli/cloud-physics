program main

    use cloud_physics

    implicit none

    ! Test variables.
    real(8) :: T = 278.15d0
    real(8) :: p = 800.0d0
    real(8) :: qv = 1.10497d-02

    character(30) file_path
    integer i

    call cloud_physics_init("./namelist")

    write(file_path, "('output.', I5.5, '.nc')") 0
    call cloud_physics_output(file_path)
    do i = 1, 1000
        call cloud_physics_run(T, p, qv)
        write(file_path, "('output.', I5.5, '.nc')") i
        call cloud_physics_output(file_path)
    end do

    call cloud_physics_final()

end program main