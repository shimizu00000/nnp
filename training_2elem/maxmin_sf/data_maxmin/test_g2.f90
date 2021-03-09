program test
implicit none
integer i
integer,parameter :: nparam = 8
double precision dummy(nparam,2)
character work*8,work1*50

open(unit=1,file="g2_a_a_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_a_a_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_a_a_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_a_b_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_a_b_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_a_b_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_a_c_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_a_c_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_a_c_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo





open(unit=1,file="g2_b_a_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_b_a_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_b_a_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_b_b_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_b_b_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_b_b_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_b_c_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_b_c_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_b_c_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo





open(unit=1,file="g2_c_a_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_c_a_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_c_a_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_c_b_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_c_b_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_c_b_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g2_c_c_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)

write(*,*)" "
write(*,'(A10)')"g2_c_c_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo

write(*,*)" "
write(*,'(A10)')"g2_c_c_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo




end
