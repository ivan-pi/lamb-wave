module lamb

  implicit none

  integer, parameter :: wp = kind(1.0d0)

  real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

  complex(wp), parameter :: j_ = complex(0.0_wp,1.0_wp)

contains

  subroutine sym(n,h,kt,kl,gamma,a,b,c,d)
    integer, intent(in) :: n
    real(wp), intent(in) :: h, kt, kl, gamma
    complex(wp), intent(out), dimension(n,n) :: a, b, c, d

    real(wp) :: alpha(n), beta(n)
    integer :: i, j

    do i = 1, n
      alpha(i) = (real(i,wp) - 1.0_wp)*pi/h
      beta(i) = (real(i,wp) - 0.5_wp)*pi/h
    end do

    a = 0
    b = 0
    c = 0
    d = 0

    do i = 1, n
      d(i,i) = alpha(i)**2/gamma - kl**2
      b(i,i) = gamma*beta(i)**2 - kt**2
    end do

    do i = 1, n
      a(i,1) = sqrt(2.0_wp)*j_*(2.0_wp - gamma)*(-1)**i
    end do
    do j = 2, n
      do i = 1, n
        a(i,j) = 2*j_*(-1)**(i+j)*(alpha(j)**2 + (gamma - 2)*beta(i)**2) / &
                                  (beta(i)**2 - alpha(j)**2)
      end do
    end do

    a = a/(h*gamma)

    c = -gamma * transpose(a)

  end subroutine

  subroutine full_eig(n,h,kt,kl,gamma,k)
    integer, intent(in) :: n
    real(wp), intent(in) :: h, kt, kl, gamma
    complex(wp), intent(out) :: k(4*n)

    ! Local computations
    complex(wp), dimension(n,n) :: a, b, c, d
    integer :: i

    ! Eigenvalue computation
    complex(wp), allocatable :: M(:,:)
    complex(wp), allocatable :: zvl(:,:), zvr(:,:)
    complex(wp), allocatable :: work(:), rwork(:)
    integer :: s, lwork, info

    s = 4*n
    allocate(m(s,s))

    lwork = 2*s
    allocate(work(lwork),rwork(lwork))

    call sym(n,h,kt,kl,gamma,a,b,c,d)

    M = 0
    do i = 1, 2*n
      M(i,2*n+i) = 1.0_wp
    end do

    M(2*n+1:3*n,1:n) = -b
    M(3*n+1:4*n,n+1:2*n) = -d

    M(2*n+1:3*n,3*n+1:4*n) = -a
    M(3*n+1:4*n,2*n+1:3*n) = -c

    call zgeev('N','N',4*n,M,4*n,k,zvl,1,zvr,1,work,lwork,rwork,info)

    if (info /= 0) then
      print *, "{full_eig} failed with info = ", info
      stop
    end if

  end subroutine

  subroutine reduced_eig(n,h,kt,kl,gamma,k)
    integer, intent(in) :: n
    real(wp), intent(in) :: h, kt, kl, gamma
    complex(wp), intent(out) :: k(2*n)

    ! Local computations
    complex(wp), dimension(n,n) :: a, b, c, d, f2, g2
    integer :: i

    ! Eigenvalue computation
    complex(wp), allocatable, target :: M(:,:)
    complex(wp), allocatable :: zldum(:,:), zrdum(:,:)
    complex(wp), allocatable :: work(:), rwork(:)
    integer :: s, lwork

    ! Linear system solve
    complex(wp), allocatable :: aba(:,:)
    integer, allocatable :: ipiv(:)
    integer :: info

    s = 2*n
    allocate(m(s,s),aba(n,n),ipiv(n))

    lwork = 2*s
    allocate(work(lwork),rwork(lwork))

    call sym(n,h,kt,kl,gamma,a,b,c,d)

    ! Pre-compute any expressions using matrix A, before it is factorized
    aba = matmul(b,a)
    f2 = d - matmul(c,a)

    ! Solve Ax = BA, solution is A^{-1}BA
    call zgesv(n,n,a,n,ipiv,aba,n,info)

    if (info /= 0) then
      print *, "{full_eig} zgesv failed with info = ", info
      stop
    end if

    ! Remaining operations
    f2 = f2 + aba
    g2 = matmul(d,aba)

    ! Zero upper half
    M(1:n,:) = 0.0_wp

    ! Bottom left block
    M(n+1:2*n,1:n) = -g2

    ! Upper right block
    do i = 1, n
      M(i,n+i) = 1.0_wp
    end do 

    ! Bottom right block
    M(n+1:2*n,n+1:2*n) = -f2

    !
    ! Eigenvalue computation
    !
    call zgeev('N','N',s,M,s,k,zldum,1,zrdum,1,work,lwork,rwork,info)

    if (info /= 0) then
      print *, "{full_eig} zgeev failed with info = ", info
      stop
    end if

    k = sqrt(k)

  end subroutine

  subroutine disp(a,text)
    complex(wp), intent(in) :: a(:,:)
    character(len=*), intent(in), optional :: text

    integer :: i

    if (present(text)) then
      write(*,'(A)') text
    end if
    do i = 1, size(a,1)
      print *, a(i,:)
    end do
  end subroutine


  function positive_eigvals(k) result(kpos)
    complex(wp), intent(in) :: k(:)
    complex(wp), allocatable :: kpos(:)

    kpos = pack(k, aimag(k) >= 0.0_wp)
  end function

end module

program main

  use lamb
  implicit none

  integer, parameter :: n = 32
  real(wp) :: omega, vl, vs, kt, kl, gamma, h

  integer, parameter :: s = 2*n
  complex(wp), allocatable :: w(:), ym(:)
  integer :: i

  allocate(w(s))

  h = 0.01_wp

  omega = 215.0e3_wp

  vl = 4170.0_wp
  vs = 2150.0_wp

  kt = omega/vs
  kl = omega/vl
  gamma = (vl/vs)**2

  !call full_eig(n,h,kt,kl,gamma,w)
  call reduced_eig(n,h,kt,kl,gamma,w)

  w = positive_eigvals(w)

  block
    integer :: unit
    open(newunit=unit,file='eigvals.txt')
    
    do i = 1, size(w)
      write(unit,*) real(w(i))*h, aimag(w(i))*h
    end do

    write(unit,*)
    write(unit,*)

    allocate(ym(s))
    call merkulov(s,ym)

    do i = 1, s
      write(unit,*) real(ym(i)), aimag(ym(i))
    end do
    close(unit)
  end block

  print *, "kt = ", kt

contains

  subroutine merkulov(n,y)
    integer, intent(in) :: n
    complex(wp), intent(out) :: y(n)

    real(wp) :: nr, re, im
    integer :: i

    do i = 1, n
      nr = real(i,wp) + 0.5_wp
      re = 0.5_wp*log(2*pi*nr)
      im = 0.5_wp*(pi*nr - log(2*pi*nr)/(pi*nr))
      y(i) = complex(re,im)
    end do
  end subroutine

end program