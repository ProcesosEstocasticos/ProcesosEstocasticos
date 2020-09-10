program test_derived_type

module globals

implicit none
save
real(8) r
end module globals

!!!!!!!!!!!!!!!!!!!!!!!!!
use globals
implicit none

type::random !! 1 valor entero, 1 real
   integer::value
   real(8)::hist
end type random
type(random),allocatable,dimension(:)::dist
integer::i,k,L
character(80) output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fin declaracion de variables

write(*,*) 'entre L' !numero de veces que llamo al generador
read(*,*) L
write(*,*) ' output file'
read(*,*) output

allocate(dist(L))
do i=1,L
    call random_mumber(r) !generador nativo de fortran, no es muy bueno
    k=r*L+1 !selecciono 1 valor entero entre L (sin el cero)
    dist(k)%value=k
    dist(k)%hist=dist(k)%value+1d0
enddo
open(1,file=output)
do i=1,L
    k=dist(i)%hist
    if(k /= 0) write(1,*)  dist(i)%value,dist(i)%hist
enddo
deallocate(dist)
100 format(80A)
stop
end program test_derived_type
