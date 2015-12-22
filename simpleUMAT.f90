subroutine umat(stress,statev,ddsdde,sse,spd,scd,  &
        rpl,ddsddt,drplde,drpldt,stran,dstran,  &
        time,dtime,temp,dtemp,predef,dpred,materl,ndi,nshr,ntens,  &
        nstatv,props,nprops,coords,drot,pnewdt,celent,  &
        dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)

    implicit none

    integer :: k1, k2, k, Newton, nshr
    integer, intent(in) :: ndi, ntens, nprops, nstatv
    real (kind=8), intent(in) :: &
        time(2), stran(ntens), props(nprops), dstran(ntens)
    real (kind=8):: & 
        diffstran(ntens), &
        eps_e(ntens), eps_p(ntens), kap, tol, &
        E_const, E_const33, nu, Sy0, H, &
        ebulk3, eg2, eg, eg3, elam, C_mat(3,3), &
        Dstress(6), S_trial(3), Sy, &
        yieldFunc(1), S(3), &
        S1, S2, S3, v1, &
        f_S(3), r(3), &
        r_S(3,3), &
        Sy_kap, &
        f_Sy(1), &
        f_kap(1), &
        r_Sy(3), &
        r_kap(3), &
        invC(3,3), Dlam, dellam(1), a_vec(3), b_vec, invA(4,4), A_mat(4,4), &
        DSkap(4,1), DS(3), Dkap, Deps_p(3), &
        C_mod(3,3), C_alg(3,3), C_algd(1)
    real (kind=8), intent(out) :: ddsdde(ntens,ntens)
    real (kind=8), intent(inout) :: stress(ntens), statev(nstatv)
    real (kind=8), parameter :: one = 1.0d0, two = 2.0d0, three = 3.0d0

    ! unused variables
    real (kind=8)                    :: sse
    real (kind=8)                    :: spd
    real (kind=8)                    :: scd
    real (kind=8)                    :: rpl
    real (kind=8)                    :: ddsddt(ntens)
    real (kind=8)                    :: drplde(ntens)
    real (kind=8)                    :: drpldt
    real (kind=8)                    :: dtime
    real (kind=8)                    :: temp
    real (kind=8)                    :: dtemp
    real (kind=8)                    :: predef(1)
    real (kind=8)                    :: dpred(1)
    CHARACTER (LEN=80)       :: materl

    real (kind=8)                    :: coords(3)
    real (kind=8)                    :: drot(3,3)
    real (kind=8)                    :: pnewdt
    real (kind=8)                    :: celent
    real (kind=8)                    :: dfgrd0(3,3)
    real (kind=8)                    :: dfgrd1(3,3)
    INTEGER                  :: noel
    INTEGER                  :: npt
    INTEGER                  :: kslay
    INTEGER                  :: kspt
    INTEGER                  :: kstep
    INTEGER                  :: kinc

    print *, 'time: ', time(1)
    print *, 'statev(k): ', statev
    print *, 'dstran(k): ', dstran
    print *, 'stran(k): ', stran
    print *, 'stress(k): ', stress

    !    RECOVER ELASTIC AND PLASTIC STRAINS AND ROTATE FORWARD
    !     ALSO RECOVER EQUIVALENT PLASTIC STRAIN
    ! CALL rotsig(statev( 1), drot, eps_e, 2, ndi, nshr)
    ! CALL rotsig(statev(ntens+1), drot, eps_p, 2, ndi, nshr)
    ! eps_e = eps_e + dstran

    ! initialise elastic and plastic strains
    eps_e = statev(1:ntens) + dstran
    eps_p = statev(1+ntens:6+ntens)

    ! print *, 'eps_e: ', eps_e
    ! print *, 'eps_p: ', eps_p
    
    ! initialise accumulated plastic strain
    kap = statev(13)

    ! tolerance
    tol = 1.d-6

    ! max number of iterations
    Newton = 50

    ! material parameters 
    E_const = 1000
    nu = 0.0
    Sy0 = 1.0
    H = 100.0

    ! elastic properties
    IF(nu > 0.4999.AND.nu < 0.5001) nu=0.499
    ebulk3=E_const/(one-two*nu)
    eg2=E_const/(one+nu)
    eg=eg2/two
    eg3=three*eg
    elam=(ebulk3-eg2)/three

    ! algorithmic modulus
    DO  k1=1,ntens
      DO  k2=1,ntens
        ddsdde(k2,k1)=0.0
      END DO
    END DO

    DO  k1=1,ndi
      DO  k2=1,ndi
        ddsdde(k2,k1)=elam
      END DO
      ddsdde(k1,k1)=eg2+elam
    END DO
    DO  k1=ndi+1,ntens
      ddsdde(k1,k1)=eg
    END DO

    ! Elasticity matrix
    C_mat(1,1) = ddsdde(3,3)
    C_mat(2,2) = ddsdde(5,5)
    C_mat(3,3) = ddsdde(6,6)

    print *,'C_mat: ', C_mat(1,:)
    print *, C_mat(2,:)
    print *, C_mat(3,:)

    ! stress
    Dstress = matmul(ddsdde, dstran) 

    stress = stress + Dstress
    S_trial(1) = stress(3)
    S_trial(2) = stress(5)
    S_trial(3) = stress(6)
     
    ! Yield stress
    Sy = Sy0 + H * kap

    ! Yield function
    yieldFunc = norm2(S_trial) - Sy

    print *, 'stress = ', stress
    print *, 'yieldFunc = ', yieldfunc(1)
    
    ! Plastic corrector step (Return-Mapping Algorithm)

if (yieldFunc(1) > tol) then
    
    ! Initialisation (k = 0)
    S = S_trial
    ! f_S
    call initialise(f_S,3,1)
    ! r_S
    call initialise(r_S,3,3)
    ! r_kap
    call initialise(r_kap,3,1)
    ! invA
    call initialise(invA,4,4)
    ! A_mat
    call initialise(A_mat,4,4)
    Dlam = 0.0

    ! Newton-Raphson
    
    do k = 1, Newton

        ! Check for convergence
        if (abs(yieldFunc(1)) < tol .and. norm2(a_vec) < tol) then
            exit
        end if

        ! Symbolic variables
        S1 = S(1)
        S2 = S(2)
        S3 = S(3)

        print *, 'S1: ', S1
        print *, 'S2: ', S2
        print *, 'S3: ', S3

        ! Substitute variables
        v1 = S1**2.0 + S2**2.0 + S3**2.0

        print *, 'v1: ', v1

        ! gradient of yield function with respect to S
        f_S(1) = S1/(v1)**(1.0/2.0)
        f_S(2) = S2/(v1)**(1.0/2.0)
        f_S(3) = S3/(v1)**(1.0/2.0)

        print *, 'f_S: ', f_S(1)
        print *, f_S(2)
        print *, f_S(3)


        ! flow vector (normalised f_S)
        r = f_S/norm2(f_S)

        print *, 'r: ', r(1)
        print *, r(2)
        print *, r(3)

        ! gradient of flow vector with respect to S

        r_S(1,1) = 1.0D0/sqrt(S1**2/(v1)+S2**2/(v1)+S3**2/(v1))*1.0D0/sqrt(v1)-S1**2*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0)+S1*1.0D0/(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S1*(-2.0D0))/(v1)+S1**3*1.0D0/(v1)**2*2.0D0+&
            S1*S2**2*1.0D0/(v1)**2*2.0D0+S1*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0) 

        r_S(1,2) = S1*1.0D0/(S1**2/(v1)+S2**2/(v1)+S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S2*(-2.0D0))/(v1)+&
            S2**3*1.0D0/(v1)**2*2.0D0+S1**2*S2*1.0D0/(v1)**2*2.0D0+&
            S2*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S1*S2*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(1,3) = S1*1.0D0/(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S3*(-2.0D0))/(v1)+S3**3*1.0D0/(v1)**2*2.0D0+&
            S1**2*S3*1.0D0/(v1)**2*2.0D0+S2**2*S3*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S1*S3*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(2,1) = S2*1.0D0/(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S1*(-2.0D0))/(v1)+S1**3*1.0D0/(v1)**2*2.0D0+&
            S1*S2**2*1.0D0/(v1)**2*2.0D0+S1*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S1*S2*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(2,2) = 1.0D0/sqrt(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))*1.0D0/sqrt(v1)-S2**2*1.0D0/sqrt(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0)+S2*1.0D0/(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S2*(-2.0D0))/(v1)+&
            S2**3*1.0D0/(v1)**2*2.0D0+S1**2*S2*1.0D0/(v1)**2*2.0D0+&
            S2*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0) 

        r_S(2,3) = S2*1.0D0/(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S3*(-2.0D0))/(v1)+S3**3*1.0D0/(v1)**2*2.0D0+&
            S1**2*S3*1.0D0/(v1)**2*2.0D0+S2**2*S3*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S2*S3*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(3,1) = S3*1.0D0/(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S1*(-2.0D0))/(v1)+S1**3*1.0D0/(v1)**2*2.0D0+&
            S1*S2**2*1.0D0/(v1)**2*2.0D0+S1*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S1*S3*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(3,2) = S3*1.0D0/(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S2*(-2.0D0))/(v1)+S2**3*1.0D0/(v1)**2*2.0D0+&
            S1**2*S2*1.0D0/(v1)**2*2.0D0+S2*S3**2*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)-S2*S3*1.0D0/sqrt(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0) 

        r_S(3,3) = 1.0D0/sqrt(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))*1.0D0/sqrt(v1)-S3**2*1.0D0/sqrt(S1**2/(v1)+S2**2/(v1)+&
            S3**2/(v1))*1.0D0/(v1)**(3.0D0/2.0D0)+S3*1.0D0/(S1**2/(v1)+&
            S2**2/(v1)+S3**2/(v1))**(3.0D0/2.0D0)*1.0D0/sqrt(v1)*((S3*(-2.0D0))/(v1)+&
            S3**3*1.0D0/(v1)**2*2.0D0+S1**2*S3*1.0D0/(v1)**2*2.0D0+&
            S2**2*S3*1.0D0/(v1)**2*2.0D0)*(1.0D0/2.0D0)

        print *,'r_S: ', r_S(1,:)
        print *, r_S(2,:)
        print *, r_S(3,:)

        ! derivative of Sy with respect to kap
        Sy_kap = H
        
        ! derivative of yield function with respect to kap
        
        f_kap = f_Sy * Sy_kap

        print *,'f_kap: ', f_kap

        ! gradient of normalised flow vector with respect to Sy
        r_Sy(1) = 0.0
        r_Sy(2) = 0.0
        r_Sy(3) = 0.0

        ! gradient of normalised flow vector with respect to kap

        r_kap = r_Sy * Sy_kap

        print *,'r_kap: ', r_kap(1)
        print *, r_kap(2)
        print *, r_kap(3)
        
        ! inverse of C 

        call inverse(C_mat,invC)

        print *, 'invC: ', invC(1,:)
        print *, invC(2,:)
        print *, invC(3,:)

        ! increment in plasticity parameter
        
        a_vec = matmul(invC, (S - S_trial)) + Dlam * r
        
        print *, 'a_vec: ', a_vec

        b_vec = -kap + statev(13) + Dlam
        
        print *, 'b_vec: ', b_vec

        invA(1:3,1:3) = invC + Dlam*r_S
        invA(1:3,4) = Dlam*r_kap
        invA(4,4) = -1.0d0
        
        print *,'invA: ', invA(1,:)
        print *, invA(2,:)
        print *, invA(3,:)
        print *, invA(4,:)

        ! A_mat = inverse of invA 

        call M44INV (invA, A_mat)
        
        print *, 'A_mat: ', A_mat(1,:)
        print *, A_mat(2,:)
        print *, A_mat(3,:)
        print *, A_mat(4,:)

        dellam = (yieldFunc - matmul(matmul((/reshape((/f_S/),(/1,3/)), f_kap/), A_mat), reshape((/a_vec, b_vec/), (/4,1/)))) / &
           (matmul(matmul((/reshape((/f_S/),(/1,3/)), f_kap/), A_mat), reshape((/r, 1.0d0/), (/4,1/))))

        print *, 'dellam: ', dellam

        ! Obtain increments in stress and internal variables
        DSkap = - matmul(A_mat, reshape((/a_vec, b_vec/), (/4,1/))) - dellam(1) * matmul(A_mat, reshape((/r, 1.0d0/), (/4,1/)))
        DS = DSkap(1:3,1)
        Dkap = DSkap(4,1)

        print *, 'DSkap: ', DSkap

        ! Update plastic strain and internal variables
        
        Deps_p = matmul(invC, DS)
        
        print *, 'Deps_p', Deps_p

        eps_p(3) = eps_p(3) - Deps_p(1)
        eps_p(5) = eps_p(5) - Deps_p(2)
        eps_p(6) = eps_p(6) - Deps_p(3)
        
        eps_e(3) = eps_e(3) + Deps_p(1)
        eps_e(5) = eps_e(5) + Deps_p(2)
        eps_e(6) = eps_e(6) + Deps_p(3)

        kap = kap + Dkap
        
        Dlam = Dlam + dellam(1)
        
        S = S + DS
        
        Sy = Sy0 + H * kap

        yieldFunc = norm2(S) - Sy

        print *, 'eps_p(k+1): ', eps_p
        print *, 'eps_e(k+1)', eps_e
        print *, 'kap: ', kap
        print *, 'Dlam: ', Dlam
        print *, 'S: ', S
        print *, 'yieldFunc: ', yieldFunc

    end do! for k=1 : Newton
    
    ! If it does not converge
    
    if (k > Newton - 1 .AND. yieldFunc(1) > tol) then
        print*, 'Warning: did not converge'
        pause
    end if
        ! Update stress state variables if plastic
        stress(3) = S(1)
        stress(5) = S(2)
        stress(6) = S(3)
        statev(1:ntens) = eps_e
        statev(1+ntens:6+ntens) = eps_p
        statev(1+2*ntens) = kap   

        ! Update algorithmic modulus if plastic
        C_mod = A_mat(1:3,1:3)

        C_algd = (matmul(reshape((/f_S/),(/1,3/)), matmul(C_mod,r)) - f_kap)
        C_alg = C_mod - ((matmul(reshape(matmul(C_mod,r), (/3,1/)), matmul(reshape((/f_S/),(/1,3/)),C_mod)))) / C_algd(1)

        ! print *, 'C_alg: ', C_alg(1,:)
        ! print *, C_alg(2,:)
        ! print *, C_alg(3,:)

        ddsdde(3,3) = C_alg(1,1)
        ddsdde(3,5) = C_alg(1,2)
        ddsdde(3,6) = C_alg(1,3)
        ddsdde(5,3) = C_alg(2,1)
        ddsdde(5,5) = C_alg(2,2)
        ddsdde(5,6) = C_alg(2,3)
        ddsdde(6,3) = C_alg(3,1)
        ddsdde(6,5) = C_alg(3,2)
        ddsdde(6,6) = C_alg(3,3)

else

    ! Update stress state variables if elastic
    S = S_trial
    stress(3) = S(1)
    stress(5) = S(2)
    stress(6) = S(3)
    statev(1:ntens) = eps_e
    statev(1+ntens:6+ntens) = eps_p
    statev(1+2*ntens) = kap
    
end if! if yieldFunc > tol (plasticity)
    
!     print *, 'ddsdde(k+1): '
!     print *, ddsdde(1,:)
!     print *, ddsdde(2,:)
!     print *, ddsdde(3,:)
!     print *, ddsdde(4,:)
!     print *, ddsdde(5,:)
!     print *, ddsdde(6,:)
    print *, 'statev(k+1): ', statev
    print *, ''
    print *, 'stress(k+1): ', stress
    print *, ''
    do k1=1,13
      if (isnan(statev(k1))) then
        pause
      end if
    end do
    
end subroutine umat

subroutine initialise(output,m,n)

    implicit none
    integer, intent(in) :: m, n
    real (kind=8), intent(out) :: output(m,n)
    integer :: k1, k2

    ! initialise an m x n array to 0
    DO  k1=1,m
      DO  k2=1,n
        output(k2,k1)=0.0
      END DO
    END DO

end subroutine initialise

subroutine inverse(input,output)

  implicit none 
  real (kind=8), intent(in) :: input(3,3)
  real (kind=8), intent(out) :: output(3,3)
  real (kind=8) :: determinant
  integer :: i, j
   
  ! calculate the inverse of a 3x3 matrix
   
  determinant = -input(1,3)*input(2,2)*input(3,1)+input(1,2)*input(2,3)*input(3,1)+input(1,3)*input(2,1)*input(3,2)-&
  input(1,1)*input(2,3)*input(3,2)-input(1,2)*input(2,1)*input(3,3)+input(1,1)*input(2,2)*input(3,3)

  output(1,1) = (-input(2,3)*input(3,2)+input(2,2)*input(3,3)) / determinant
  output(1,2) = (+input(1,3)*input(3,2)-input(1,2)*input(3,3)) / determinant
  output(1,3) = (-input(1,3)*input(2,2)+input(1,2)*input(2,3)) / determinant

  output(2,1) = (+input(2,3)*input(3,1)-input(2,1)*input(3,3)) / determinant
  output(2,2) = (-input(1,3)*input(3,1)+input(1,1)*input(3,3)) / determinant
  output(2,3) = (+input(1,3)*input(2,1)-input(1,1)*input(2,3)) / determinant

  output(3,1) = (-input(2,2)*input(3,1)+input(2,1)*input(3,2)) / determinant
  output(3,2) = (+input(1,2)*input(3,1)-input(1,1)*input(3,2)) / determinant
  output(3,3) = (-input(1,2)*input(2,1)+input(1,1)*input(2,2)) / determinant

end subroutine inverse

SUBROUTINE M44INV (A, AINV)
! http://web.hku.hk/~gdli/UsefulFiles/matrix/m44inv_f90.txt
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
  DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: AINV

  DOUBLE PRECISION :: DET
  DOUBLE PRECISION, DIMENSION(4,4) :: COFACTOR


  DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
         A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
         A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

  COFACTOR(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
  COFACTOR(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
  COFACTOR(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
  COFACTOR(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
  COFACTOR(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
  COFACTOR(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
  COFACTOR(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
  COFACTOR(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
  COFACTOR(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
  COFACTOR(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
  COFACTOR(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
  COFACTOR(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
  COFACTOR(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
  COFACTOR(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  COFACTOR(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
  COFACTOR(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

  AINV = TRANSPOSE(COFACTOR) / DET

  RETURN

  END SUBROUTINE M44INV

