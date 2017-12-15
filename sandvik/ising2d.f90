!----------------------!
 module systemvariables
 implicit none

 integer :: l                        ! system length
 integer :: n                        ! number of spins (n=l*l)
 real(8) :: pflip(-1:1,-4:4)         ! flip probabilities
 integer, allocatable :: spin(:)     ! spin array

 end module systemvariables
!--------------------------!

!-------------------------------------------------------------------!
!Metropolis-algorithm simulation of the two-dimensional Ising model.!
!Reads the following from a file 'ising.in':                        !
! l,t,h = system length, temperature (T/J), field (h.J)             !
! bins, binsteps = number of bins, MC steps per bin                 !
!Writes these bin averages to the file 'bindata.dat':               ! 
! <E>/N, <E**2>/N**2, <|M|>/N, <M**2>/N**2                          !
!These results should be processed by the program average.f90.      !
!A single random numbers seeds (3) are read from a file 'seed.in'.
!-------------------------------------------------------------------!
 program ising2d
!----------------------------------!
 use systemvariables; implicit none

 integer :: i,j,bins,binsteps,seed
 real(8) :: t,h

 open(10,file='ising.in',status='old')
 read(10,*)l,t,h
 read(10,*)bins,binsteps
 close(10)

 call initialize(t,h)
 do i=1,binsteps
    call mcstep
 end do
 do j=1,bins
    call resetdatasums
    do i=1,binsteps
       call mcstep
       call measure(h)
    end do
    call writebindata(binsteps)
 end do
 deallocate(spin)
     
!-------------------!
 end program ising2d
!-------------------!

!--------------------------------------------!
!Carries out one Monte Carlo step, defined as! 
!n flip attempts of randomly selected spins. !
!--------------------------------------------!
 subroutine mcstep
!----------------------------------!
 use systemvariables; implicit none

 integer :: i,s,x,y,s1,s2,s3,s4
 real(8), external :: ran

 do i=1,n
    s=int(ran()*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l)+y*l)
    s2=spin(x+mod(y+1,l)*l)
    s3=spin(mod(x-1+l,l)+y*l)
    s4=spin(x+mod(y-1+l,l)*l)
    if (ran()<pflip(spin(s),s1+s2+s3+s4)) spin(s)=-spin(s)
 end do

 end subroutine mcstep
!---------------------!

!--------------------------------------------------------!
!Measures the energy (enrg1), the squared energy (enrg2),! 
!the absolute value of the magnetization (magn1), and the!
!squared magnetization (magn2).                          !
!--------------------------------------------------------!
 subroutine measure(h)
!----------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 integer :: s,x,y,e,m
 real(8) :: h

 e=0
 do s=0,n-1
    x=mod(s,l); y=s/l
    e=e-spin(s)*(spin(mod(x+1,l)+y*l)+spin(x+mod(y+1,l)*l))
 enddo
 m=sum(spin)
 enrg1=enrg1+dble(e-m*h)
 enrg2=enrg2+dble(e-m*h)**2

 magn1=magn1+dble(abs(m))
 magn2=magn2+dble(m)**2

 end subroutine measure
!----------------------!

!-------------------------------------------------!
!Writes accumulated data to the file 'bindata.dat'!
!-------------------------------------------------!
 subroutine writebindata(steps)
!------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 integer :: steps
 
 open(1,file='bindata.dat',status='unknown',position='append')
 enrg1=enrg1/(dble(steps)*dble(n))
 enrg2=enrg2/(dble(steps)*dble(n)**2)
 magn1=magn1/(dble(steps)*dble(n))
 magn2=magn2/(dble(steps)*dble(n)**2)
 write(1,1)enrg1,enrg2,magn1,magn2
 1 format(4f18.12)
 close(1)

 end subroutine writebindata
!---------------------------!

!--------------------------------------------!
!Sets the data accumulation variables to zero!
!--------------------------------------------!
 subroutine resetdatasums
!------------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2
 common/measurments/enrg1,enrg2,magn1,magn2

 enrg1=0.d0
 enrg2=0.d0
 magn1=0.d0
 magn2=0.d0

 end subroutine resetdatasums
!----------------------------!

!------------------------------------------------------------------!
!Various initializations: system size n=l*l, flipping probabilities!
!pflip, random number generator, initial random spin state.        !
!------------------------------------------------------------------!
 subroutine initialize(t,h)
!--------------------------!
 use systemvariables
 implicit none

 integer :: i,j,ns
 real(8) :: t,h

 integer, allocatable :: seed(:)
 real(8), external :: ran

 n=l*l
 do i=-4,4
    pflip(-1,i)=exp(+2.*(i+h)/t)
    pflip(+1,i)=exp(-2.*(i+h)/t)
 enddo

 call initran(1)

 allocate (spin(0:n-1))
 do i=0,n-1
    spin(i)=2*int(2.d0*ran())-1
 enddo

 end subroutine initialize
!-------------------------!

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!
