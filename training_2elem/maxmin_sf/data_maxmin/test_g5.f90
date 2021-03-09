program test
implicit none
integer i
integer,parameter :: nparam = 22
double precision dummy(nparam,2)
character work*8,work1*100

open(unit=1,file="g5_a_aa_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_aa_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_aa_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_a_ab_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_ab_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_ab_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_a_ac_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_ac_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_ac_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_a_bb_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_bb_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_bb_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_a_bc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_bc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_bc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_a_cc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_a_cc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_a_cc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo






open(unit=1,file="g5_b_aa_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_aa_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_aa_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_b_ab_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_ab_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_ab_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_b_ac_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_ac_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_ac_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_b_bb_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_bb_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_bb_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_b_bc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_bc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_bc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_b_cc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_b_cc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_b_cc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



open(unit=1,file="g5_c_aa_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_aa_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_aa_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_c_ab_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_ab_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_ab_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_c_ac_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_ac_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_ac_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_c_bb_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_bb_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_bb_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_c_bc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_bc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_bc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo


open(unit=1,file="g5_c_cc_maxmin.dat",action="read",form="unformatted")
read(1)dummy
close(1)
write(*,*)" "
write(*,'(A11)')"g5_c_cc_max"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,1)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo
write(*,*)" "
write(*,'(A11)')"g5_c_cc_min"
do i = 1, nparam
  write(work,'(i8)')i-1
  write(work1,'(f60.45)')dummy(i,2)
  write(*,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(*,'(a)')trim(adjustl(work1))
enddo



end
