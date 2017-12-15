!----------------!
 program averages
!----------------!
 implicit none

 integer :: i,l,bins
 real(8) :: temp,av,er
 real(8),allocatable :: data(:,:)

 open(1,file='ising.in',status='old')
 read(1,*)l,temp
 close(1)

 open (1,file='bindata.dat',status='old')
 bins=0
 do 
   read(1,*,end=1)av
   bins=bins+1
 end do
 1 close(1)
 allocate(data(bins,4))

 open (1,file='bindata.dat',status='old')
 do i=1,bins
    read(1,*)data(i,1),data(i,2),data(i,3),data(i,4)
 end do
 close(1)

 print*,'---------------------------------------'
 print*,'L = ',l
 print*,'T = ',temp
 print*,'Number of bins: ',bins
 print*,'---------------------------------------'
 call averageanderror1(data(:,1),bins,av,er)
 write(*,2)' <E>/N  : ',av,er
 call averageanderror2(data(:,1),data(:,2),bins,l**2,temp**2,av,er)
 write(*,2)' <C>/N  : ',av,er
 call averageanderror1(data(:,3),bins,av,er)
 write(*,2)' <|m|>  : ',av,er
 call averageanderror1(data(:,4),bins,av,er)
 write(*,2)' <m**2> : ',av,er
 call averageanderror2(data(:,3),data(:,4),bins,l**2,temp,av,er)
 write(*,2)' <X>    : ',av,er
 2 format(a,2f15.7) 

 end program averages
!--------------------!

!--------------------------------------------!
 subroutine averageanderror1(data,bins,av,er)
!--------------------------------------------!
 implicit none

 integer :: i,bins
 real(8) :: data(bins),av,er

 av=sum(data)/dble(bins)
 er=sum(data**2)/dble(bins)
 er=sqrt((er-av**2)/dble(bins-1))

 end subroutine averageanderror1
!-------------------------------!

!----------------------------------------------------------!
 subroutine averageanderror2(data1,data2,bins,n,temp,av,er)
!----------------------------------------------------------!
 implicit none

 integer :: i,n,bins
 real(8) :: data1(bins),data2(bins),temp,diff,av,er

 av=sum(data2-data1**2)/dble(bins)
 er=sum((data2-data1**2)**2)/dble(bins)
 er=sqrt((er-av**2)/dble(bins-1))
 av=av*dble(n)/temp
 er=er*dble(n)/temp

 end subroutine averageanderror2
!-------------------------------!
