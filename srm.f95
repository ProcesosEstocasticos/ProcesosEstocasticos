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
  real(8) xmax
  real(8),allocatable,dimension(:)::time
end module temporal

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$module integral
!!$implicit none
!!$save
!!$real(8)::sig=1d0
!!$real(8) x,y
!!$real(8)::Pi=3.1419
!!$end module integral
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module pasos
  implicit none
  save
  integer,allocatable,dimension(:)::der,izq
  integer Lx
end module pasos



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program srm
  use globals
  use random
  use temporal
  use pasos

  implicit none

  REAL(8),allocatable,dimension(:):: rugcuad,media  

  INTEGER:: i,j,k,nrea,rea
  real(8) t,dt
  real(8) h1
  character(len=80) output

  WRITE(*,*) 'Nro de experimentos.' !! 100
  READ(*,*) nrea           
  write(*,*) 'xmax potencia'  ! t > t_s~ Lˆ2 
  read(*,*) xmax
  write(*,*) 'tamanio Lx'
  read(*,*) Lx                  !128 y 256
  write(*,*) 'output file?'
  read(*,*) output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call initialize_ramdom
  call tiempos
  call bc
  allocate(rugcuad(0:tptos),h(Lx)) !alloco

  rugcuad =0d0  !inicializo acumulador
  avg: do rea=1,nrea
     !print*,rea
     t=0.d0
     h=0
     timeinterval:do j=0,tptos
        evol:do while (t < time(j))
           i=Lx*r+1
           call rand
           RULES:if(h(i)<= h(izq(i)).and.h(i)<=h(der(i))) then
              h(i)=h(i)+1
           else if(h(i)< h(izq(i)).and.h(i)>h(der(i))) then
              h(der(i))=h(der(i))+1
           else if(h(i)<h(der(I)).and.h(i)>h(izq(i))) then
              h(izq(i))=h(izq(i))+1
           else if(h(der(i))==h(izq(i))) then
              call rand
              if(r < 0.5) then
                 h(der(i))=h(der(i))+1
              else
                 h(izq(i))=h(izq(i))+1
              endif
           endif RULES
           t=t+1.d0/Lx
        enddo evol
        h1=sum(1.d0*h)/Lx
        do i=1,Lx
           rugcuad(j) = rugcuad(j) +(1.d0*h(i) - h1)**2 
        enddo
     enddo timeinterval
     if(mod(rea,10)==0) then
        print*,rea
        open(2,file=output)
        write(2,*) '#L=',Lx,'rea=',rea
        do j=0,tptos
           write(2,*) time(j),sqrt(rugcuad(j)/(rea*lx))
           !write(*,100) time(j),sqrt(rugcuad(j)/rea)
        enddo
        close(2)
     endif
  enddo avg
100 format(1x,f17.7,10(1x,f25.8))
1000 format(80A)
  stop
end program srm




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
recursive SUBROUTINE RAND
  use random
  IMPLICIT NONE

  I1=IR(I1)
  I2=IR(I2)
  SEM(I1)=XOR(SEM(I1),SEM(I2))
  R=DFLOAT(SEM(I1))/CMAX
  if(r==1.d0) call rand
  RETURN
END SUBROUTINE RAND


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********
subroutine tiempos
  use pasos
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

!*******************************

subroutine bc
  use pasos
  implicit none
  integer i
  allocate(der(Lx),izq(Lx))

  do i=1,Lx
     der(i)=i+1
     izq(i)=i-1
  enddo
  der(Lx)=1
  izq(1)=Lx
  return
end subroutine bc
