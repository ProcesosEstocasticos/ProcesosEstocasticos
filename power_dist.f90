module globals

implicit none
save
INTEGER I1,I2,IR(670)
INTEGER SEM(670),SEE(33)
real(8) cmax
REAL(8) R 
real susu
integer hour

end module globals

program power_dist
use globals
implicit none
integer i,k,N
real(8)::lambda = 3.5, z
real(8),allocatable, dimension(:)::P

call initialize_ramdom


write(*,*) 'entre N'
read(*,*) N
allocate(P(N))
do i =1, N
   call random_number(r) !!usar otro de abajo
   z=1-lambda
   k=2*r**(1/z)
   if( k > N) k= N 
   P(k)=P(k)+1
enddo
open(1,file ='PowerDist')
write(1,*)'# N=', N, 'lambda=', lambda
do i=1, N
	if(P(i) /= 0) write(1,*) i, P(i)
enddo
stop
end

!**********************
subroutine initialize_ramdom
use globals
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
end

!***********************************
SUBROUTINE RAND
use globals
IMPLICIT NONE

I1=IR(I1)
I2=IR(I2)
SEM(I1)=XOR(SEM(I1),SEM(I2))
R=DFLOAT(SEM(I1))/CMAX
RETURN
END






