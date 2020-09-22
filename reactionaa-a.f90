module random !modulo con las variables para generar numeros al azar 
implicit none
save
INTEGER I1,I2,IR(670)
INTEGER SEM(670),SEE(40)
real(8) cmax
REAL(8) R 
real susu
integer hour
end module random

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
module globals !van ls condiciones de contorno
implicit none
save
integer,allocatable,dimension(:)::der,izq
integer L
end module globals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module temporal ! variable temporal
implicit none
save
integer tptos !importante
integer tmax  
real(8),allocatable,dimension(:)::time
end module temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program reaction
use random !
use globals! 
use temporal ! NO OLVIDAR DE PONER

implicit none !importante
integer,allocatable,dimension(:)::Lat,position  !el vector position almacenar'a una lista de los caminantes aun "vivos", siempre enumerados de manera creciente
real(8),allocatable,dimension(:)::rho_t
integer rea,nrea
integer i,sel,x,pos,dx
integer nw
integer j
real(8) rho,t !rho es la densidad inicial
character(80)::output

write(*,*) 'entre L'
read(*,*) L
write(*,*) 'entre densidad'
read(*,*) rho
!write(*,*) 'tmax'
!read(*,*) tmax
write(*,*) 'entre # de realizaciones'
read(*,*) nrea
write(*,*) 'output file'
read(*,*) output

call initialize_ramdom
call bc
call tiempos !chequear como da con ese esquiespaciado

allocate(Lat(L),position(L),rho_t(0:tptos))
rho_t=0 ! inicializo el acumulador
do rea=1,nrea
   nw=0
   t=0d0
   Lat=0
   position=0
   do i=1,L
      call rand
      if(r < rho) then
         nw=nw+1
         position(nw)=i
         lat(i)=nw
      endif
   enddo
   !print*,'tptos=',tptos
   do j=0,tptos
      !print*,'j=',j
      do while(t <  time(j))
         !print*,'t=',t
         call rand
         sel=r*nw+1 !elijo al rw sel-esimo
         x=position(sel) ! que esta en la posicion x
         call rand
         dx=2*r+1 !para que lado difunde?
         select case(dx)
         case(1)
            pos=der(x)
         case(2)
            pos=izq(x)
         end select
         if(Lat(pos) ==0 ) then !si la posicion a la que va esta desocupada entonces
            Lat(x)=0   ! posicion x queda vacia
            Lat(pos)=sel !en esa posicion pongo al rw seleccionado
            position(sel)=pos !la posicion del que rw # sel se actualiza 
         else
            Lat(x)=0 !aca desocupo al que elegi, se evapora            
            if(sel/=nw)then
               x=position(nw) !miramos  la coordenada del ultimo rw  
               position(sel)=x !intercambio la lista (por ejemplo) si la posicion del ultimo es x=8, hay 4 rw, y se selecciono al 2        
               Lat(x)=sel  !Latt(8)=2               
            endif
            nw=nw-1  !el new disminuye en uno (ya que uno se evaporo)
         endif
         t=t+1d0/nw !t avanza en 1/nw. Cuando elegí a todos los nw, el tiempo avanza en 1, simultaniedad
         !do i=1,L
         !   write(*,'(I4)',advance='no')Lat(i)
         !enddo
         !write(*,*)
         !print*,"............................."
         !read(*,*)
      enddo
      rho_t(j)=rho_t(j)+1.d0*nw !mido la densidad al paso j-esimo
   enddo
   if(mod(rea,10)==0) then !guardo cada 10 realizaciones
      print*,rea
      open(1,file=output)
      write(1,*) '#L=',L,'rho=',rho,'rea=',rea !el # es para graficar con el grace 
      do i=0,tptos
         write(1,*) time(i),rho_t(i)/L/rea !promedio
      enddo
      close(1) !DEBO CERRARLO PORQUE LO VUELO A ABRIR CADA 10 EXPERIMENTOS
   endif
enddo

stop
end program reaction


!************************
!*******************************
subroutine bc
use globals
implicit none
integer i
allocate(der(L),izq(L))

do i=1,L
	der(i)=i+1
	izq(i)=i-1
enddo
der(L)=1
izq(1)=L
return
end

!****************************
subroutine tiempos
use temporal
implicit none
integer i
real(8) xmin,xmax

xmin=-2 !potencia minima
xmax=6! !potencia maxima 10**xmax deber ser > que 1/(rho**2)
tptos=(xmax-xmin)/0.1d0 !numero de puntos con equiespaciado log con equiespaciado 0.1
allocate(time(0:tptos))
do i=0,tptos
   time(i)=10**xmin
   xmin=xmin+0.1d0
enddo
return
end subroutine tiempos


!**********************
subroutine initialize_ramdom
use random
implicit none
integer i

i1=1
i2=104
CMAX=2*(2**30-1)+1
CALL SYSTEM_CLOCK(hour)       		!llama al reloj
hour=17065204
print*,hour
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
end subroutine initialize_ramdom

!***********************************
SUBROUTINE RAND
use random
IMPLICIT NONE

221 I1=IR(I1)
I2=IR(I2)
SEM(I1)=XOR(SEM(I1),SEM(I2))
R=DFLOAT(SEM(I1))/CMAX
if(R==1d0) goto 221
RETURN
END SUBROUTINE RAND



