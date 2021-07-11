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
      b(i,i) = beta(i)**2/gamma - kl**2
      d(i,i) = gamma*alpha(i)**2 - kt**2
    end do

    do i = 1, n
      a(1,i) = sqrt(2.0_wp)*j_*(-1)**i
    end do
    do i = 2, n
      do j = 1, n
        a(i,j) = 2*j_*(-1)**(i+j)*((gamma - 2)*alpha(j)**2 + beta(i)**2) / &
                                  (alpha(i)**2 - beta(j)**2)
      end do
    end do

    a = a/(h*gamma)

    c = -gamma * transpose(a)

  end subroutine

end module

program main

  use lamb
  implicit none

  integer, parameter :: n = 8
  real(wp) :: omega, vl, vs, kt, kl, gamma, h

  complex(wp), allocatable, dimension(:,:) :: a, b, c, d

  integer, parameter :: s = 4*n
  complex(wp), allocatable :: m(:,:), w(:), ym(:)
  complex(wp), allocatable :: zvl(:,:), zvr(:,:)
  integer, parameter :: lwork = 2*s
  complex(wp), allocatable :: work(:), rwork(:)
  integer :: info, i

  external :: zgeev

  allocate(a(n,n),b(n,n),c(n,n),d(n,n))

  allocate(m(s,s),w(s))

  allocate(work(lwork),rwork(lwork))


  h = 0.01_wp

  omega = 215.0e3_wp

  vl = 4170.0_wp
  vs = 2150.0_wp

  kt = omega/vs
  kl = omega/vl
  gamma = (vl/vs)**2


  call sym(n,h,kt,kl,gamma,a,b,c,d)

  print *, a

  print *, "b = "
  print *, b

  print *, "d = "
  print *, d

  M = 0
  do i = 1, 2*n
    M(i,2*n+i) = 1.0_wp
  end do

  M(2*n+1:3*n,1:n) = -b
  M(3*n+1:4*n,n+1:2*n) = -d

  M(2*n+1:3*n,3*n+1:4*n) = -a
  M(3*n+1:4*n,2*n+1:3*n) = -c

  call zgeev('N','N',4*n,M,4*n,w,zvl,1,zvr,1,work,lwork,rwork,info)
  print *, "info = ", info

  print *, "Eigenvalues: "
  do i = 1, s
    print *, w(i)
  end do

  block
    integer :: unit
    open(newunit=unit,file='eigvals.txt')
    do i = 1, s
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