module temporal
  implicit none
  save
  integer tptos
  integer tmax
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program rw
  use globals
  use random
  use temporal
  
  implicit none
  real(8),allocatable,dimension(:)::rho_t
  real(8),allocatable,dimension(:)::position
  integer,allocatable,dimension(:)::Lat
  integer i,j
  integer x
  character(len=80) output
  integer jt, N, L
  real(8) t
  integer pos,sel,dx,nw
  integer:: rea,nrea
  real(8) rho

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  write(*,*) 'entre el tiempo maximo'
  read(*,*) tmax
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
  call tiempo(tmax) !tiempo tomado lineal en este programa

  allocate(position(L),Lat(L),rho_t(0:tptos))
  rho_t=0d0

  do rea=1,nrea
     nw=0
     position=0
     lat=0
     do i=1,L
        call rand
        if(r < rho) then
           nw=nw+1
           position(nw)=i
           Lat(i)=nw !en la posicion i-esima se encuentra el caminante nw-esimo
        endif
     enddo

     t=0
     do jt=0,tptos
        evol:do while(t<=time(jt)) !aca está lineal
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
           if(Lat(pos)==0)then
              position(sel) =pos !actualizamos la posicion del caminante
              Lat(pos)=sel !el lugar "pos" es ocupado por el caminante "sel"
              Lat(x)=0 !el lugar "x" ya no es ocupado por el caminante "sel"
           endif
           t=t+1d0/nw  ! el tiempo avanza como 1/nw, cuando elegí nw el tiempo avanza en uno
        enddo evol
        rho_t(jt)=rho_t(jt)+nw/dble(L)
     enddo
  enddo

  open(1,file=output)
  write(1,*) '#tiempos=',tptos+1,'# caminantes=',nrea
  do j=0,tptos
     write(1,*) time(j),rho_t(j)/nrea
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
  !print*,'semilla',hour    
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
recursive SUBROUTINE RAND
  use random
  IMPLICIT NONE

1 I1=IR(I1)
  I2=IR(I2)
  SEM(I1)=XOR(SEM(I1),SEM(I2))
  R=DFLOAT(SEM(I1))/CMAX
  if(R==1d0) call rand
  RETURN
END SUBROUTINE RAND
!*****************************************

!*************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tiempo(N)
  use temporal
  integer N,i
  !equiespaciado a eleccion

  tptos=N!! log10(1.d0*N)
  allocate(time(0:tptos))
  do i=0,tptos
     time(i)=i!! 10**i
  enddo
  return
end subroutine tiempo

!**********************
