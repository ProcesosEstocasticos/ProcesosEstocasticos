module globals
  implicit none
  save

  integer nptos
  integer kmin,kmax
  integer N_node  

  real(8) lambda
  real(8) r

  integer,allocatable,dimension(:)::edge
  integer,allocatable,dimension(:)::node
  integer,allocatable,dimension(:)::kk

  real(8),allocatable,dimension(:)::Pk
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

  integer i,rea,nrea
  real(8) product


  N_node   =300

  kmin     =2
  kmax     =100
  lambda   =2.5d0


  nrea     =1



  allocate(kk(N_node),node(N_node))
  allocate(Pk(kmin:kmax))


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Pk           =0d0
  !Pk(0)        =exp(-lambda)
  !product      =1d0
  do i=kmin,kmax
     !product=product*i         
     !Pk(i)=1.d0*exp(-lambda)*lambda**(i)/dble(product)
     Pk(i)=dble(i)**(-lambda)
  enddo

  Pk=Pk/sum(Pk)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  call initialize_random

  do rea=1,nrea

     call redmr
     call gephy
     deallocate(edge)

  enddo

end Program Perc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!       Gephy
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gephy
  use globals
  implicit none

  integer i,j,k
  integer L_tot
  
  character(len=23) :: str1,str2,strn
  
  integer, allocatable, dimension(:)::Listhead,Listtail


  L_tot=sum(kk)/2
  allocate(Listhead(L_tot),Listtail(L_tot))  
  call MakeList(Listhead,Listtail,L_tot)

  open(1,file='RedGephi.gexf')
  write(1,'(a)')'<?xml version="1.0" encoding="UTF-8"?>'
  write(1,'(a)')'<gexf xmlns="http://www.gexf.net/1.3" xmlns:viz="&
             http://www.gexf.net/1.3/viz" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
  write(1,'(a)')'xsi:schemaLocation="http://www.gexf.net/1.3 http://www.gexf.net/1.3/gexf.xsd" version="1.3">'
  
  
 
  write(1,'(a)') '   <graph>'
  write(1,'(a)') '       <nodes>'


  do i=1,N_node
     write(strn,'(I23)') i
     write(1,'(a)')'           <node id="'//trim(adjustl(strn))//'">'
     write(1,'(a)')'           </node>'
  enddo
  write(1,'(a)')'       </nodes>'
  write(1,'(a)')'       <edges>'

  do i=1,L_tot
     write(strn,'(I23)') i
     write(str1,'(I23)') Listhead(i)
     write(str2,'(I23)') Listtail(i)

   write(1,'(a)')'           <edge id="'//trim(adjustl(strn))//'" source="'//trim(adjustl(str1))//&
                  '" target="'//trim(adjustl(str2))//'"/>'

  enddo

  write(1,'(a)')'       </edges>'
  write(1,'(a)')'   </graph>'
  write(1,'(a)')'</gexf>'
  close(1)

112 format(a38)


end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Hacer lista de links
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeList(Listhead,Listtail,L_tot)
  use globals
  implicit none
  integer i,j,k
  integer vec
  integer posEdge
  integer L_tot,L_aux
  
  integer, dimension(L_tot)::Listhead,Listtail
  integer, allocatable, dimension(:)::edgeUsado

  allocate(edgeUsado(sum(kk)))
  
  edgeUsado=0
  Listhead =0
  Listtail =0
  L_Aux    =0
  do i=1,N_node
     if(kk(i)==0)cycle
     do k=1,kk(i)
        posEdge =node(i)+k-1
        if(edgeUsado(posEdge)==1)cycle
        edgeUsado(posEdge)=1
        vec =edge(posEdge)
        do j=1,kk(vec)
           if(edge(node(vec)+j-1)==i)then
              edgeUsado(node(vec)+j-1)=1
              exit
           endif
        enddo
        L_aux=L_aux+1
        Listhead(L_aux)=i
        Listtail(L_aux)=vec   
     enddo
  enddo
end subroutine MakeList


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

  real(8),dimension(kmin-1:kmax)::IntegralPk


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


