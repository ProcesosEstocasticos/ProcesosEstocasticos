module globals
  implicit none
  save

  integer nptos
  integer kmin,kmax
  integer N_node  
  integer cluster_number,max_mass,ss2,gc

  real(8) lambda
  real(8) r

  integer,allocatable,dimension(:)::edge
  integer,allocatable,dimension(:)::node
  integer,allocatable,dimension(:)::kk
  integer,allocatable,dimension(:)::cluster,mass
  
  real(8),allocatable,dimension(:)::Pk
  real(8),allocatable,dimension(:)::frac
  real(8),allocatable,dimension(:)::Pinf,S2
    
end module globals

module random
	save
	integer::q1=1,q2=104
	integer ir(670)
	integer::sem(670)

	real(8)::nmax=2*(2**30-1)+1
	real(8) sig
end module random





Program Perc
  use globals
  use random

  implicit none

  integer i,j,rea,nrea
  real(8) product
  real(8) pmin,deltp

  N_node   =100000

  kmin     =0
  kmax     =20
  lambda   =4d0


  nrea     =100
  nptos    =100
  pmin     =0d0
  deltp    =0.01d0


  allocate(kk(N_node),node(N_node))
  allocate(Pk(kmin:kmax))
  allocate(cluster(N_node),mass(N_node))
  allocate(Pinf(nptos),S2(nptos))
  allocate(frac(nptos))  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Pk           =0d0
  Pk(0)        =exp(-lambda)
  product      =1d0
  do i=1,kmax
     product=product*i         
     Pk(i)=1.d0*exp(-lambda)*lambda**(i)/dble(product)
  enddo

  Pk=Pk/sum(Pk)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  frac=0d0
  do i=1,nptos
     frac(i)=pmin+i*deltp
  enddo

    
  call initialize_random

  Pinf  =0d0
  S2    =0d0

  do rea=1,nrea
     print*,rea
     call redmr
     call bombing_nodes
     deallocate(edge)
     
     open(1,file='Pinf.dat')
     open(2,file='S2.dat')
     write(1,*)'#',rea
     write(2,*)'#',rea
     do i=1,nptos
        write(1,*) frac(i),Pinf(i)/dble(rea)/dble(N_node)
        write(2,*) frac(i),S2(i)/dble(rea)
     enddo
     close(1)
     close(2)     
  enddo

end Program Perc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Percolacion de nodos
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bombing_nodes
  use globals
  use random
  implicit none
  integer i,j,k
  integer pto,sel,nodo1,vec
  integer N_neto
  
  real(8) pp
  
  integer, dimension(N_node)::ListNod

  do i=1,N_node
     ListNod(i)=i
  enddo

  N_neto =N_node
  do pto=nptos,1,-1
     pp   =frac(pto)
     elim: do while(dble(N_neto)/dble(N_node)>pp)
        call rand
        sel   =r*N_neto+1
        nodo1 =ListNod(sel)
        ListNod(sel)=ListNod(N_neto)
        N_neto=N_neto-1
        
        !actualizamos el vector edge y kk
        do k=1,kk(nodo1)
           vec  =edge(node(nodo1)+k-1)
           do j=1,kk(vec)
              if(edge(node(vec)+j-1)==nodo1)then
                 edge(node(vec)+j-1)=edge(node(vec)+kk(vec)-1)
                 kk(vec) =kk(vec)-1
                 exit
              endif
           enddo
        enddo
        kk(nodo1)=-1
     enddo elim
     
     call cluster_id
     Pinf(pto)=Pinf(pto)+max_mass
     S2(pto)  =S2(pto)+ss2  
  enddo

end subroutine bombing_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Creaci'on
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine redmr
  use globals
  use random
  implicit none
  integer i,j,k
  integer suma   
  integer m,n,mayor,menor
  integer cuento, mlinks
  integer link,AuxN1,AuxN2

  real(8) ws

  integer, allocatable, dimension(:)::listaAguj,kkAux

  real(8),dimension(kmin:kmax)::IntegralPk


  allocate(kkAux(n_node))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Realizamos la integral del Pk
  
  IntegralPk   =0
  ws           =0d0
  do i=kmin,kmax
     ws          = ws+Pk(i)
     IntegralPk(i)   = ws
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!asignamos las conexiones a los nodos
2 kk=0 
  do i=1,N_node
     call rand
     do j=kmin,kmax
        if(IntegralPk(j-1)<=r.and.r<=IntegralPk(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo

  suma=sum(kk)

  if(mod(suma,2).ne.0) go to 2
  
  allocate(listaAguj(suma))  
  allocate(edge(suma))

  mlinks        =suma

  listaAguj     =0
  link          =0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!armamos la lista de agujitas

  do i=1,n_Node
     if(kk(i)==0) cycle
     do k=1,kk(i)
        link                =link+1
        listaAguj(link)    =i
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!armamos el vector node

  node(1)       =1
  do i=1,n_node-1
     node(i+1)    = node(i)+kk(i)
  enddo

  edge        =0
  kkAux     =0
  cuento      =0
  do while(mlinks>0)
3    if(cuento==100) then
	deallocate(listaAguj,edge)
	go to 2
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Elegimos dos agujitas
     
     call rand
     AuxN1    =r*mlinks+1
     call rand
     AuxN2    =r*mlinks+1

     m        =listaAguj(AuxN1)
     n        =listaAguj(AuxN2)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Vamos a comprobar si podemos conectar las agujitas o no

     if(m.ne.n) then
	do k=1,kkAux(m)
           if(edge(node(m)+k-1)==n)then
              cuento     =cuento+1
              go to 3
           endif
	enddo
	edge(node(m)+kkAux(m))   =n
	edge(node(n)+kkAux(n))   =m
	call rand
	kkAux(m)                 =kkAux(m)+1
	kkAux(n)                 =kkAux(n)+1
	mayor                      =max0(AuxN1,AuxN2)
	menor                      =min0(AuxN1,AuxN2)

	if(mayor==mlinks) then
           listaAguj(mayor)       =listaAguj(mlinks)
           listaAguj(menor)       =listaAguj(mlinks-1)
	else 
           if (mayor==mlinks-1) then
              listaAguj(menor)    =listaAguj(mlinks)
           else
              listaAguj(mayor)    =listaAguj(mlinks)
              listaAguj(menor)    =listaAguj(mlinks-1)
           endif
	endif
	mlinks           =mlinks-2
	cuento           =0

     else
        cuento           =cuento+1
        go to 3
     endif
  enddo
  deallocate(kkAux,listaAguj)
  
  
end subroutine redmr

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Cluster identification
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_id
use Globals
implicit none

integer M,N,i,j,z,k
integer nburn,cs

integer,dimension(N_node)::ocupp
integer,dimension(0:N_node+1)::w


ss2=0
ocupp=1
cluster=0
cluster_number=0
mass=0
max_mass=0
do i=1,N_node
   if(kk(i)==-1) cycle
   if(cluster(i) .eq. 0.and.ocupp(i)==1) then 
      cluster_number=cluster_number+1
      w(0)=i
      nburn=1
      cs=0
      do while(nburn.ne.0)
         nburn=nburn-1
         j=w(nburn)
         cluster(j)=cluster_number
         cs=cs+1
         do k=0,kk(j)-1
            m=edge(node(j)+k)
            if(cluster(m) == 0.and.ocupp(m)==1) then
               ocupp(m)=0  
               w(nburn)=m
               nburn=nburn+1
            endif
         enddo
      enddo
      mass(cluster_number)=cs
      if(cs > max_mass) then
         max_mass=cs
         gc=cluster_number
      endif
      if(cs > ss2.and.cs < max_mass) ss2=cs
    endif
enddo
return
end subroutine cluster_id


!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)			!llama al reloj
  !print*,hour
  !hour=1142413639
  see=hour
  CALL RANDOM_SEED(put=see)		!semilla de la realizaci\F3n
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)   =r*nmax
     ir(i)    =i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random

!****************************************
subroutine rand
  use globals
  use random
  implicit none

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


