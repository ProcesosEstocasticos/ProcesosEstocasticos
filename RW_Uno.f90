module temporal
implicit none
save
integer tptos
!integer tmax
real(8),allocatable,dimension(:)::time
end module temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
implicit none
integer,allocatable,dimension(:)::der,izq
end module globals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module random
implicit none
save
INTEGER I1,I2,IR(670)
INTEGER SEM(670),SEE(33)
real(8) cmax
REAL(8) R 
real susu
integer hour
end module random


program rw
use globals
use temporal
use random

implicit none
real(8),allocatable,dimension(:)::R2
integer j
integer x,paso
character(len=80) output
integer jt, N, L
integer t
integer pos,pos0,sel
integer:: rea,nrea




write(*,*) 'entre el numero de pasos'
read(*,*) N
write(*,*) 'entre tamanio red'
read(*,*) L
write(*,*) 'entre el numero de realizaciones'
read(*,*) nrea
write(*,*) 'entre nombre de archivo salida'
read(*,*) output

call initialize_random !para el generador
call bc(L) !condiciones de contorno periodicas
call tiempo(N) !tiempo equiespaciado pot 10

allocate(R2(0:tptos))

do rea=1,nrea
   print*,rea
   pos=L/2
   pos0=pos
   t=0
   do jt=0,tptos
      do while(t<=time(jt))
         call rand !!!llamo generador
         sel=r*2+1
         select case (sel)
         case(1)
            pos=der(pos)
         case(2)
            pos=izq(pos)
         end select
         t=t+1
         paso=abs(pos-pos0)
      enddo
      R2(jt)=R2(jt)+paso
   enddo
enddo
open(1,file=output)
write(1,*) '#Npasos=',N,'# caminantes=',nrea
do j=0,tptos
    write(1,*) time(j),R2(j)/nrea
enddo
stop
end program rw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bc(L)
use Globals
implicit none
integer i,L
allocate(der(L),izq(L))
do i=1,L
    der(i)=i+1
    izq(i)=i-1
enddo
der(L)=1
izq(1)=L
return
end subroutine bc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tiempo(N)
use temporal
integer N
!equiespaciodo en potencias de 10

tptos=log10(1.d0*N)
allocate(time(0:tptos))
do i=0,tptos
   time(i)=10**i
enddo
return
end subroutine tiempo

!**********************
subroutine initialize_random
use random
implicit none
integer i

i1=1
i2=104
CMAX=2*(2**30-1)+1
CALL SYSTEM_CLOCK(hour)       		!llama al reloj
see=hour       
CALL RANDOM_SEED(put=see)!semilla de la realizacion 
CALL RANDOM_SEED(get=see)

DO I=1,670 
   CALL random_number(SUSU)
   SEM(I)=SUSU*CMAX 
   IR(I)=I+1 
ENDDO
IR(670)=1           
return
end subroutine initialize_random

!***********************************
RECURSIVE SUBROUTINE  RAND
use random
IMPLICIT NONE

I1=IR(I1)
I2=IR(I2)
SEM(I1)=XOR(SEM(I1),SEM(I2))
R=DFLOAT(SEM(I1))/CMAX
if(R==1.d0) call Rand
RETURN
END SUBROUTINE RAND

