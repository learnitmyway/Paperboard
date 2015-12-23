    ! rm a.out & gfortran -Wall -Wextra -Wconversion -pedantic input.f90 umat.f90 && ./a.out

program input

    implicit none  

    ! timer step
    integer :: n

    ! number of increments per time interval
    integer, parameter :: i = 5

    real (kind=8):: t_data(2), eps33_data(2), eps13_data(2), eps23_data(2)

    ! total number of intervals
    integer, parameter :: ii = size(t_data)+(i-1)*(size(t_data)-1)

    ! number of state variables
    integer, parameter :: nstatv = 13
    
    ! other variables
    integer :: remainder
    real (kind=8):: &
        time(2,ii), eps33(ii), eps13(ii), eps23(ii), &
        stran(6,ii), dstran(6,ii), &
        ddsdde(6,6), stress(6,ii), statev(nstatv,ii)

    ! print *, 'ii: ', ii

    t_data = (/0.0, 5.0/) ! (/0.0, 5.0, 15.0, 25.0, 35.0/)

    ! strain input
    eps33_data = (/0.0, 0.005/) ! (/0.0, 0.005, -0.005, 0.005, -0.005/)
    eps13_data = (/0.0, 0.0/) ! (/0.0, 0.0, 0.0, 0.0, 0.0/) 
    eps23_data = (/0.0, 0.0/)

    ! linear interpolation
    time = 0.0

    call linear_interpolation(t_data, eps33_data, i, time(1,:), eps33, size(t_data))
    call linear_interpolation(t_data, eps13_data, i, time(1,:), eps13, size(t_data))
    call linear_interpolation(t_data, eps23_data, i, time(1,:), eps23, size(t_data))

    ! Initialisation (t = 0)
    stran = 0.0
    dstran = 0.0
    stress = 0.0
    statev = 0.0

    ! print *, 'stran(initial): ', stran
    ! print *, 
    ! print *, 'dstran(initial): ', dstran
    ! print *, 
    ! print *, 'statev(initial): ', statev
    ! print *, 

    do n = 1, ii - 1

        ! strain 
        stran(3,n+1) = eps33(n+1)
        stran(5,n+1) = 2 * eps13(n+1)
        stran(6,n+1) = 2 * eps23(n+1)

        ! print *, 'stran(input): ', stran(:,n+1)

        ! strain increment
        dstran(:,n+1) = stran(:,n+1) - stran(:,n)

        ! print *, 'dstran(input): ', dstran(:,n+1)

        ! fill arrays
        stress(:,n+1) = stress(:,n)
        statev(1:13,n+1) = statev(1:13,n)

        ! print *, 'n: ', n
        ! print *, 
        ! print *, 'stran(before): ', stran(:,n+1)
        ! print *, 
        ! print *, 'dstran(before): ', dstran(:,n+1)
        ! print *, 
        ! print *, 'statev(before): ', statev(:,n+1)
        ! print *, 
        ! print *, 'stress(before): ', stress(:,n+1)
        ! print *, 
           
        ! stress vector
        call umat(stress(:,n+1),statev(:,n+1),ddsdde,&
        stran(:,n+1),dstran(:,n+1),time(:,n+1), nstatv)

        remainder = modulo(n,i)
        ! print *, 'remainder(n+1): ', remainder
        ! print *, 
        ! if (remainder == 0) then
        !     print *, 'remainder(n+1): ', remainder
        !     print *, 
            ! print *, 'time(n+1): ', time(:,n+1)
            ! print *, 
            ! print *, 'stran(n+1): ', stran(:,n+1)
            ! print *, 
            ! print *, 'dstran(n+1): ', dstran(:,n+1)
            ! print *, 
            ! print *, 'statev(n+1): ', statev(:,n+1)
            ! print *, 
            ! print *, 'stress(n+1): ', stress(:,n+1)
        ! end if

    end do

end program input

subroutine linear_interpolation(x_data, y_data, n, x, y, dim)
    implicit none

    integer, intent(in) :: n, dim
    real (kind=8), intent(in) :: x_data(dim), y_data(dim)
    integer :: i, j, i_max, j_max
    real (kind=8):: dx, dy
    real (kind=8), intent(out) :: &
        x(size(x_data)+(n-1)*(size(x_data)-1)), & 
        y(size(y_data)+(n-1)*(size(y_data)-1)) 

    ! Initialisation
    x(1) = x_data(1)
    y(1) = y_data(1)

    ! Linear interpolation 
    i_max = size(x_data) - 1
    j_max = n + 1

    do i=1, i_max

        dx = (x_data(i+1) - x_data(i)) / n
        dy = (y_data(i+1) - y_data(i)) / n
        
        do j=2, j_max
            x((i - 1) * n + j) = x((i - 1) * n + (j - 1)) + dx;
            y((i - 1) * n + j) = y((i - 1) * n + (j - 1)) + dy;
        end do

    end do

end subroutine linear_interpolation
