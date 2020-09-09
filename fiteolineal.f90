program fiteo_lineal

implicit none
real(8) x(100),y(100)
real(8) x_med,y_med,xy_med,x2_med
real(8) m,b
integer i,n_datos
character(3) n_c

y_med=0.d0
x_med=0.d0
xy_med=0.d0
x2_med=0.d0

n_datos=1
do 
   write(*,*) 'entre x e y'
   read(*,*) x(n_datos),y(n_datos)
   write(*,*) 'ingrese otro dato? si o no'
   read(*,*) n_c
   if(n_c.ne.'si') exit
   n_datos=n_datos+1
enddo
do i=1,n_datos
   y_med=y_med+y(i)
   x_med=x_med+x(i)
enddo
y_med=y_med/n_datos
x_med=x_med/n_datos
do i=1,n_datos
   xy_med=xy_med+x(i)*y(i)
   x2_med=x2_med+x(i)*x(i)
enddo
m=x2_med-n_datos*x_med*x_med
m=(xy_med-y_med*x_med*n_datos)/m
b=y_med-m*x_med
do i=1,n_datos
   write(*,*) x(i),b+m*x(i)
enddo
write(*,*) b,m
stop
end


