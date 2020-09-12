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
use random

implicit none
real(8),allocatable,dimension(:)::rho_t
real(8),allocatable,dimension(:)::position
integer i,j
integer x
character(len=80) output
integer jt, N, L
real(8) t
integer pos,sel,nw
integer:: rea,nrea
real(8) rho




write(*,*) 'entre el numero de pasos'
read(*,*) N
write(*,*) 'entre tamanio red'
read(*,*) L
write(*,*) 'entre densidad'
read(*,*) rho
write(*,*) 'entre el numero de realizaciones'
read(*,*) nrea
write(*,*) 'entre nombre de archivo salida'
read(*,*) output

call initialize_random !para el generador
call bc(L) !condiciones de contorno periodicas

allocate(position(L),rho_t(0:N))
rho_t=0d0

do rea=1,nrea
   
   nw=0
   position=0
   do i=1,L
       call rand
       if(r < rho) then
           nw=nw+1
           position(nw)=i
       endif
   enddo
   print*,rea,nw

   t=0
   jt=0
   
   evol:do while(t<=N)
      call rand !!!vamos a seleccionar al azar un caminante      
      sel=r*nw+1
      x=position(sel)
      call rand !!!vamos a seleccionar al azar una direccion
      sel=r*2+1
      select case (sel)
      case(1)
         pos=der(x)
      case(2)
         pos=izq(x)
      end select
      
      position(sel) =pos !actualizamos la posicion del caminante
      t=t+1d0/nw
      if(t>=jt)then
         rho_t(jt)=rho_t(jt)+nw/dble(L)
         jt=jt+1
      endif
   enddo evol
   
enddo
open(1,file=output)
write(1,*) '#Npasos=',N,'# caminantes=',nrea
do j=0,N
    write(1,*) j,rho_t(j)/nrea
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
print*,'semilla',hour    
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
Recursive SUBROUTINE RAND
use random
IMPLICIT NONE

1 I1=IR(I1)
I2=IR(I2)
SEM(I1)=XOR(SEM(I1),SEM(I2))
R=DFLOAT(SEM(I1))/CMAX
if(R==1d0) call rand
RETURN
END SUBROUTINE RAND

