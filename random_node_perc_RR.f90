module globals
implicit none
save
type link_info
   integer :: nn
   real(8) :: energy
   integer :: value
end type link_info
type(link_info),allocatable,dimension(:)::edge
integer, allocatable,dimension(:)::node
integer,allocatable,dimension(:)::kk
integer::gc
integer,allocatable,dimension(:)::ocupp
integer,allocatable,dimension(:)::cluster
integer,allocatable,dimension(:)::Histo_Link,Histo_original,Smax
real(8),allocatable,dimension(:)::fraction,s2,Pinf,P_v,s2_v
integer ss2
integer nptos,rea
integer :: factor=1
character(len=3):: option,option2

end module globals

module random

implicit none
save
real(8) :: CMAX=2*(2**30-1)+1
INTEGER::I1=1,I2=104,IR(670)
INTEGER SEM(670),SEE(20)
integer hour,azar 
real susu
real(8) r
end module random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use globals
use random
implicit none 

integer,allocatable,dimension(:)::Histo_mass,Histo_max_mass

integer,allocatable,dimension(:)::Histo_all_length,Histo_opt_length
integer,allocatable,dimension(:)::mass

integer Ntarget
integer N_node
integer M,N,i,j,jj,k,kpart,k1,k2,l,z,origin,endpoint

real(8) lamda
integer n_rea
integer max_mass,cluster_number,cs
integer label
real(8) norm2
character(len=80)::output,output2,outputmass
real(8) p,q
real(8) kkk
integer M_number
integer lambda
integer::k_min=2
real(8) pmin,pmax,dp

character(7) var1 
character(6) var2
character(6) var3
character(6) var4
character(6) var5
character(6) var7
character(6) var8




write(*,*) 'Number of nodes'
read(*,*)  N_node
write(*,*) '<z>'
read(*,*)  lambda
write(*,*) 'pmin'
read(*,*)  pmin
write(*,*) 'pmax'
read(*,*)  pmax
write(*,*) 'dp'
read(*,*)  dp
write(*,*) 'Number of realizations'
read(*,*)  n_rea

 
write(var1,'(I7)') N_node            
call zeros(var1)
write(var4,'(I6)') lambda
call zeros(var4)

allocate (mass(N_node))
allocate (cluster(N_node))
allocate (kk(N_node))
allocate(Node(N_node))
allocate(ocupp(N_node))

nptos=(pmax-pmin)/dp

allocate(fraction(0:nptos),Pinf(0:nptos),s2(0:nptos))!,P_v(0:,nptos))!,P_all(0:nptos),norm(0:nptos),Smax(0:nptos))

do i=0,nptos
   fraction(i)=pmin
   pmin=pmin+dp
enddo
if(pmax==1d0) fraction(nptos)=1.d0

call initialize_random
s2=0d0
Pinf=0d0

!M_number=N_node*kkk/2
!Main loop in realizations
do rea=1,n_rea
   !call net(N_node,M_number) !draw the net
   call redmr(n_node,lambda,k_min)
   call bombing_nodes(N_node,max_mass)
   !!temporary output
   if(mod(rea,10)==0) then
      print*,rea
      open(2,file='Pinf_node_RR_s2_N_'//var1//'_z_'//var4//'.dat')
      write(2,*) '# N=',N_node, 'z=', lambda, 'exp=',rea	 
      do i=0,nptos
         write(2,*) fraction(i),Pinf(i)/rea/N_node,s2(i)/rea
         write(*,*) fraction(i),Pinf(i)/rea/N_node,s2(i)/rea
      enddo
      close(2)
   endif
   deallocate(edge)
enddo !realizations),Pinf(i)/rea

100 format(80A)
200 format('#Nodes=',I8,1x,'<k>=',f10.2,1x,'realizations=',I8)
1000  format(I10,5(1x,I15))
2000  format('+',I10,I10,4(1x,f20.8),1x,I10,I10)
3000  format(I10,6(1x,f25.12))
4000  format(5(1x,f15.8))
stop
end
!!!!!!!!!!!!!!!!!!
subroutine bombing_links(N_node,max_mass)
use Globals
use random

implicit none
integer,allocatable,dimension(:)::nodo1,nodo2,target
real(8),allocatable,dimension(:)::w
integer i,ii,j,jj,k,m,ik,max_mass,Ntarget,N_node
real(8) p
integer Nsel1,Nsel2,sumk
integer total
integer,dimension(N_node)::mass
integer cluster_number,ptos,Node_ongc



total=sum(kk)/2

allocate(nodo1(total),nodo2(total))

Ntarget=0
edge%value=0
n1:do i=1,N_node
   n2:do k=0,kk(i)-1
      if(edge(node(i)+k)%value==1) cycle
      edge(node(i)+k)%value=1
      m=edge(node(i)+k)%nn
      do j=0,kk(m)-1
         if(edge(node(m)+j)%nn == i) then
            edge(node(m)+j)%value=1
            Ntarget=Ntarget+1
            nodo1(Ntarget)=i
            nodo2(Ntarget)=m
         endif
      enddo
   enddo n2
enddo n1

sumk=Ntarget
ptos=nptos
do ptos=nptos,0,-1
   p=fraction(ptos)
   kill:do while(p <1.d0*Ntarget/total)
      call rand
      i=r*Ntarget+1
      do k=1,kk(nodo1(i))
         m=edge(node(nodo1(i))+k-1)%nn
         if(m == nodo2(i)) then
            edge(node(nodo1(i))+k-1)%nn= edge(node(nodo1(i))+kk(nodo1(i))-1)%nn
            kk(nodo1(i))=kk(nodo1(i))-1
            do j=1,kk(nodo2(i))
               if(edge(node(nodo2(i))+j-1)%nn==nodo1(i)) then
                  edge(node(nodo2(i))+j-1)%nn= edge(node(nodo2(i))+kk(nodo2(i))-1)%nn
                  kk(nodo2(i))=kk(nodo2(i))-1
                  exit
               endif
            enddo
            exit
         endif
      enddo
      nodo1(i)=nodo1(Ntarget)
      nodo2(i)=nodo2(Ntarget)
      Ntarget=Ntarget-1
   enddo kill
   call cluster_id(N_node,cluster_number,mass,max_mass)
   Pinf(ptos)=Pinf(ptos)+1.d0*max_mass
   s2(ptos)=s2(ptos)+ss2
end do
return
end subroutine bombing_links

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive SUBROUTINE RAND
  use Globals
  use random

  IMPLICIT NONE
  
  I1=IR(I1)
  I2=IR(I2)
  SEM(I1)=XOR(SEM(I1),SEM(I2))
  R=DFLOAT(SEM(I1))/CMAX
  if(r == 1d0) call rand
  RETURN
END SUBROUTINE RAND

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Cluster identification
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_id(N_node,cluster_number,mass,max_mass)
use Globals
implicit none


integer N_node
integer M,N,i,j,z,k,l,flag,origin
integer,dimension(0:N_node+1)::w
integer,dimension(N_node)::mass
integer nburn,cluster_number,cs,max_mass

ss2=0
ocupp=1
cluster=0
cluster_number=0
mass=0
max_mass=0
do i=1,N_node
   if(kk(i)==0) cycle
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
            m=edge(node(j)+k)%nn
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


!!!!!!!!!!!



!!!!!!!!!!!!!!!!
!Draw the net
!!!!!!!!!!!!!!!!

subroutine net(N_node,M_number)
use Globals
use random

implicit none
 
integer N_node,M_number
integer M,N,i,j,k,flag

kk=0
allocate(edge(0:50*N_node))
edge%nn=0
Node(1)=0
do i=2,N_node
   Node(i)=Node(i-1)+50 
enddo
M=0
do while(M < M_number)
   flag=0
   CALL rand
   i=r*N_node+1
   CALL rand
   j=r*N_node+1
   if(i.eq. j) cycle
   do k=0,kk(i)-1
      if(edge(node(i)+k)%nn == j) then
            flag=1
            exit
      endif
   enddo
   if(flag.eq. 0) then
         edge(node(i)+kk(i))%nn=j
         edge(node(j)+kk(j))%nn=i
         kk(i)=kk(i)+1
         kk(j)=kk(j)+1
         M=M+1
   endif
enddo

return
end subroutine net


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Cluster identification
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_id_perc(target,Node_ongc,cluster_number,mass,max_mass,N_node)
use Globals
implicit none


integer N_node
integer M,N,i,j,z,k,l,ii,flag,origin,Node_ongc
integer,dimension(0:Node_ongc)::w
integer,dimension(Node_ongc)::target
integer,dimension(N_node)::mass
integer nburn,cluster_number,cs,max_mass

ocupp=1
cluster=0
cluster_number=0
mass=0
max_mass=0
do ii=1,Node_ongc
   i=Target(ii)		
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
            m=edge(node(j)+k)%nn
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
    endif
enddo

return
end


!!!!!!!!!!!

subroutine bombing_nodes(N_node,max_mass)
use Globals
use random

implicit none
integer,allocatable,dimension(:)::nodo1,nodo2,pos1,pos2
real(8),allocatable,dimension(:)::w
integer i,ii,j,k,m,jj,max_mass,Ntarget,N_node
real(8) p
integer ptos
integer Nsel,sumk
integer total
integer cluster_number,S_big
integer,dimension(0:N_node)::mass
!allocate(Target(N_node),nodo(N_node))


allocate(pos1(N_node),nodo1(N_node))

Ntarget=0
n1:do i=1,N_node
   Ntarget=Ntarget+1
   nodo1(Ntarget)=Ntarget
enddo n1
S_big=Ntarget
do ptos=nptos,0,-1
   p=fraction(ptos)
   kill:do while(p <1.d0*Ntarget/S_big)
      call rand
      Nsel=r*Ntarget+1
      i=nodo1(Nsel)
      do k=1,kk(i)
         ii=edge(node(i)+k-1)%nn
         links:do j=1,kk(ii)
            if(edge(node(ii)+j-1)%nn==i) then
               edge(node(ii)+j-1)%nn=edge(node(ii)+kk(ii)-1)%nn
               kk(ii)=kk(ii)-1
               exit
            endif
         enddo links
      enddo
      kk(i)=0
      nodo1(Nsel)=nodo1(ntarget)
      Ntarget=Ntarget-1
   enddo kill
   call cluster_id(N_node,cluster_number,mass,max_mass)
   !print*,'after bombing',p,max_mass,Ntarget,sum(mass)
   !pause
   Pinf(ptos)=Pinf(ptos)+1.d0*max_mass
   s2(ptos)=s2(ptos)+ss2
enddo
return
end subroutine bombing_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine redmr(n_node,lambda,k_min)
  use globals
  use random
  implicit none

  integer i,j,k,suma
  integer m,n,mayor,menor
  integer cuento,mlinks,k_max,k_min
  integer link,nodo1,nodo2,n_node
  integer lambda
  integer,allocatable,dimension(:)::listalinks,conecto

  allocate(conecto(n_node))

  kk=lambda
  suma=sum(kk)
!  if(mod(suma,2).ne.0)  go to 2                !2: arma nuevamente la distribucion de conectividad
  allocate(listalinks(suma))
  allocate(edge(0:suma))
  listalinks=0
  mlinks=suma
  
  link=0
  do i=1,n_Node
     do k=1,kk(i)
        link=link+1
        listalinks(link)=i 
     enddo
  enddo
  
  node(1)=1
  do i=2,n_node
     node(i)=node(i-1)+kk(i-1)
  enddo
  
  !conecto los nodos
  edge%nn=0
  conecto=0
  cuento=0
  do while(mlinks > 0)   
3    if(cuento==100) then
        deallocate(listalinks,edge)
        call redmr(n_node,lambda,k_min)                                !2: arma nuevamente la distribucion de conectividad
     endif
     
     call rand           !elijo los nodos a conectar
     nodo1=r*mlinks+1
     call rand
     nodo2=r*mlinks+1
     
     m=listalinks(nodo1)
     n=listalinks(nodo2)
     
     if (m.ne.n) then     !evito multiples conexiones
        do k=1,conecto(m)      
           if (edge(node(m)+k-1)%nn==n)then
              cuento=cuento+1                 
              go to 3                          !3: elijo nuevamente que nodos voy a conectar
           endif
        enddo
        edge(node(m)+conecto(m))%nn=n
        edge(node(n)+conecto(n))%nn=m
        conecto(m)=conecto(m)+1
        conecto(n)=conecto(n)+1
      
  
       ! if (conecto(m)-1 > conectividad(m)) write(*,*)m, conectividad(m),conecto(m)-1


        !para reordenar los links
        
        mayor=max0(nodo1,nodo2)
        menor=min0(nodo1,nodo2)

        if (mayor==mlinks ) then                    !si uno de los elegidos es el ultimo
           listalinks(mayor)=listalinks(mlinks)
           listalinks(menor)=listalinks(mlinks-1)
        else if (mayor==mlinks-1)then               !si uno de los elegidos es el anteultimo
           listalinks(mayor)=listalinks(mayor)
           listalinks(menor)=listalinks(mlinks)         
        else
           listalinks(mayor)=listalinks(mlinks)     !si ninguno es ultimo o anteultimo
           listalinks(menor)=listalinks(mlinks-1)
        endif
        mlinks=mlinks-2
        cuento=0
     else
        cuento=cuento+1
        go to 3
     endif
  enddo
 ! kk=conecto
  return
end subroutine redmr

!!!!!!!!!!!!!!!
!-------------------------------------------------------------
!-----------------------------------------------------------------
recursive subroutine rand_int

use globals
use random
implicit none

i1=ir(i1)
i2=ir(i2)
sem(i1)=xor(sem(i1),sem(i2))
r=dfloat(sem(i1))/cmax
if(r ==0.d0) call rand_int
return	
end subroutine rand_int

!!!!!!!!!!!
subroutine initialize_random
use random
implicit none
integer i
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zeros(b)
  implicit none

  integer i
  character(len=1):: b(6)
  do i=1,6
     if (b(i).eq.' ') b(i)='0'
  end do
  return

end subroutine zeros



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











