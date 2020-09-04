module random
implicit none
save
integer i1,i2,ir(670)
integer sem(670),see(33)
real(8) cmax
real(8) r 
real susu
integer hour
end module random


module integral
implicit none
save
real(8)::sig=1d0
real(8) x,y
real(8)::Pi=3.1419
end module integral

!!!!!!!!!!!!!!!!!!!!
program gaussina
use random
use integral

implicit none
integer i
 
call initialize_ramdom
open(1,file='GaussianDistr.dat')
do i=1,1000000
   call gauss
   write(1,*) y
enddo
stop
end program gaussina

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!C RUTINA  Hecha en Mar del plata para distribucion gaussiana   C


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


recursive SUBROUTINE GAUSS
  use  integral
  use random
  implicit none
  REAL(8) PHI,RO
  real(8)::x0=0,y0=0 !x0=4,y0=4
  CALL rand_int
  X=R
  CALL RANd_int
  Y=R
  PHI=2.D0*PI*X
  RO=DSQRT(-DLOG(Y)*2.D0/SIG)
  X=RO*DCOS(PHI)+x0
  Y=RO*DSIN(PHI)+y0
  RETURN
END SUBROUTINE GAUSS



!--------------------------------------------------------------
subroutine initialize_ramdom
use random
implicit none
integer i

i1=1
i2=104
cmax=2*(2**30-1)+1
call system_clock(hour)    !llama al reloj
see=hour       
call random_seed(put=see)  !semilla de la realizacion 
call random_seed(get=see)
do I=1,670 
   call random_number(susu)
   sem(i)=susu*cmax 
   ir(i)=i+1 
end do
ir(670)=1           
return
end subroutine initialize_ramdom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine rand_int
use random
implicit none

i1=ir(i1)
i2=ir(i2)
sem(i1)=xor(sem(i1),sem(i2))
r=dfloat(sem(i1))/cmax
if(r ==0.d0) call rand_int
return  
end subroutine rand_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
