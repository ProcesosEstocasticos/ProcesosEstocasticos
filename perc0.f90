! Este programa escibe las coordenadas del cluster de percolacion
program perc0
INTEGER see(1),hour
integer,allocatable,dimension(:,:)::latt
integer rea,nrea,L
integer i,j,k
real(8) r,p

write(*,*) 'L?'!longitud de la red
read(*,*) L
write(*,*) 'p?'!percolation probability
read(*,*) p

allocate(latt(L,L))

!call seed(hour,see)
open(1,file='percolacion.dat')
latt=0
do i=1,L
   do j=1,L
      call random_number(r)
      if(r <=p) then
         latt(i,j)=1
         write(1,*) i,j
      endif
   enddo
enddo
close(1)
stop
end



subroutine seed(hour,see)
INTEGER see(1),hour
CALL SYSTEM_CLOCK(hour)       		!llama al reloj
see=hour       
CALL RANDOM_SEED(put=see)!semilla de la realizacion 
CALL RANDOM_SEED(get=see)
return
end


