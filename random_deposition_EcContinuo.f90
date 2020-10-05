module globals
implicit none
save
real(8),allocatable,dimension(:)::h
end module globals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module random
implicit none
save
INTEGER I1,I2,IR(670)
INTEGER SEM(670),SEE(40)
real(8) cmax
REAL(8) R 
real susu
integer hour
end module random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module temporal
implicit none
save
integer tptos
integer xmax
real(8),allocatable,dimension(:)::time
end module temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module integral
implicit none
save
real(8)::sig=1d0
real(8) x,y
real(8)::Pi=3.141592653589793d0
end module integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use globals
use random
use temporal
use integral

implicit none

REAL(8),allocatable,dimension(:):: rugcuad,media  

INTEGER:: i,j,k,nrea,rea
real(8) t,dt
real(8) h1,h2
integer Lx
character(len=80) output

WRITE(*,*) 'Nro de experimentos.'
READ(*,*) nrea
write(*,*) 'xmax potencia'
read(*,*) xmax
write(*,*) 'tamanio Lx'
read(*,*) Lx
write(*,*) 'dt'
read(*,*) dt
write(*,*) 'output file?'
read(*,*) output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call initialize_ramdom
call tiempos

allocate(rugcuad(0:tptos),h(Lx),media(0:tptos))

rugcuad =0d0
media=0
avg: do rea=1,nrea
   !print*,rea
   t=0.d0
   h=0d0
   timeinterval:do j=0,tptos
      evol:do while (t < time(j))
         do i=1,Lx
            call gauss
            h(i)=h(i)+dsqrt(2*dt)*x
         enddo
         t=t+dt
      enddo evol
      h1=sum(1.d0*h)/Lx
      h2=sum(1.d0*h*h)/Lx

      rugcuad(j) = rugcuad(j) + sqrt(h2 - h1**2)

      media(j) = media(j) + h1
   enddo timeinterval
   if(mod(rea,10)==0) then
      !print*,rea
      open(2,file=output)
      write(2,*) '#L=',Lx,'rea=',rea
      do j=0,tptos
         write(2,100) time(j),rugcuad(j)/rea,media(j)/rea
         write(*,100) time(j),rugcuad(j)/rea,media(j)/rea
      enddo
      close(2)
   endif
enddo avg
100 format(1x,f17.7,10(1x,f25.8))
1000 format(80A)
stop
end program




!**********************
subroutine initialize_ramdom
use random
implicit none
integer i

i1=1
i2=104
CMAX=2*(2**30-1)+1
call SYSTEM_clock(hour)       		!llama al reloj
see=hour       
call RANDOM_seed(put=see)!semilla de la realizacion 
call RANDOM_seed(get=see)

do I=1,670 
   call random_number(SUSU)
   SEM(I)=SUSU*CMAX 
   IR(I)=I+1 
enddo
IR(670)=1           
return
end subroutine initialize_ramdom

!***********************************
SUBROUTINE RAND
use random
IMPLICIT NONE

I1=IR(I1)
I2=IR(I2)
SEM(I1)=XOR(SEM(I1),SEM(I2))
R=DFLOAT(SEM(I1))/CMAX
RETURN
END SUBROUTINE RAND


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********
subroutine tiempos
!use pasos
use temporal
implicit none
integer i
real(8) xmin
real(8)::dt1=0.1

xmin=0
tptos=(xmax-xmin)/dt1+1
allocate(time(0:tptos))
do i=0,tptos
   time(i)=10**xmin
   xmin=xmin+dt1
enddo
time(tptos)=10**xmax
return
end subroutine tiempos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gauss
  use  integral
  use random

  implicit none
!! En general tomaremos sig=1
  real(8) phi,ro
  !real(8)::x0=0,y0=0 !x0=4,y0=4
  CALL rand
  x=r
  CALL rand
  do while(r==0.d0)
     call rand
  enddo
  y=r
  phi=2.d0*pi*x
  ro=dsqrt(-dlog(y)*2d0/sig)
  x=ro*dcos(phi)!+x0
  y=ro*dsin(phi)!+x0
  return
end subroutine gauss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

