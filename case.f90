program case_file
  implicit none

  integer nota

  write(*,*) 'entre nota de 0 a 100'
  read(*,*) nota

  select case(nota)
  case(96:100)
     write(*,*)'A'
  case(87:95) 
     write(*,*)'B'
  case(77:86) 
     write(*,*)'C'
  case(67:76) 
     write(*,*)'D'
  case(1:66) 
     write(*,*)'B'
  end select
  stop
end program case_file
