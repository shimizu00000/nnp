program main_training
use allocarray
implicit none
integer i, j, k, m, dummyi
character(len=20) dummyc
character work*8,work1*100,work2*100
! Variables for NNP
!-------------------
integer nparam
double precision,allocatable,dimension(:) :: weight
character(len=4) elem_a, elem_b
integer io_nnp_a,  io_nnp_b
integer nlayer_a,  nlayer_b
integer nweight_a, nweight_b
integer nhidden_a, nhidden_b
integer,allocatable,dimension(:) :: &
  network_a, network_b
double precision,allocatable,dimension(:) :: &
  weight_a, weight_b
double precision,allocatable,dimension(:) :: &
  atomic_ene_a, atomic_ene_b
double precision,allocatable,dimension(:,:) :: &
  hidden_a, hidden_b
double precision  :: factr, pgtol
integer io_weight, percent_tst
! rc : cutoff radius
double precision rc, rn
integer neighbor
character(len=10) phi

! Variables for symmetry function
!---------------------------------
!A
integer &
  num_g2_a_a,  num_g2_a_b, &
  num_g5_a_aa, num_g5_a_ab, &
  num_g5_a_bb, &
  num_g_a(5)
!B
integer &
  num_g2_b_a,  num_g2_b_b, &
  num_g5_b_aa, num_g5_b_ab, &
  num_g5_b_bb, &
  num_g_b(5)
!A
double precision,allocatable,dimension(:,:) :: &
  g2_a_a,  g2_a_b, &
  g5_a_aa, g5_a_ab, &
  g5_a_bb
!B
double precision,allocatable,dimension(:,:) :: &
  g2_b_a,  g2_b_b, &
  g5_b_aa, g5_b_ab, &
  g5_b_bb
!A
double precision,allocatable,dimension(:,:) :: &
  g2_a_a_maxmin,  g2_a_b_maxmin, &
  g5_a_aa_maxmin, g5_a_ab_maxmin, &
  g5_a_bb_maxmin
!B
double precision,allocatable,dimension(:,:) :: &
  g2_b_a_maxmin,  g2_b_b_maxmin, &
  g5_b_aa_maxmin, g5_b_ab_maxmin, &
  g5_b_bb_maxmin
!A
double precision,allocatable,dimension(:) :: &
  eta2_a_a,  rs2_a_a, &
  eta2_a_b,  rs2_a_b, &
  eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
  eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
  eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb
!B
double precision,allocatable,dimension(:) :: &
  eta2_b_a,  rs2_b_a, &
  eta2_b_b,  rs2_b_b, &
  eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
  eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
  eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb


 !*************************
 ! READ input_nnp.dat FILE
 !*************************
   ! 2 elements
   call read_input_nnp2( &
        elem_a, &
        num_g2_a_a,  num_g2_a_b, &
        num_g5_a_aa, num_g5_a_ab, &
        num_g5_a_bb, &
        elem_b, &
        num_g2_b_a,  num_g2_b_b, &
        num_g5_b_aa, num_g5_b_ab, &
        num_g5_b_bb, &
        rc, rn, phi, neighbor, &
        io_weight, percent_tst, &
        factr, pgtol )

   num_g_a(1) = num_g2_a_a  ; num_g_a(2) = num_g2_a_b
   num_g_a(3) = num_g5_a_aa ; num_g_a(4) = num_g5_a_ab
   num_g_a(5) = num_g5_a_bb

   num_g_b(1) = num_g2_b_a  ; num_g_b(2) = num_g2_b_b
   num_g_b(3) = num_g5_b_aa ; num_g_b(4) = num_g5_b_ab
   num_g_b(5) = num_g5_b_bb


   !------------------------------
   ! Read number of hidden layers
   !------------------------------
   call read_layer_node2( nlayer_a, nlayer_b )

   allocate( network_a( nlayer_a ) )
   allocate( network_b( nlayer_b ) )

   !--------------------------------
   ! Read number of nodes per layer
   !--------------------------------
   call read_layer_node_frex2( network_a, network_b, &
                               nlayer_a,  nlayer_b, &
                               io_nnp_a,  io_nnp_b )

 !if( io_nnp_a == 0 ) deallocate( network_a )
 !if( io_nnp_b == 0 ) deallocate( network_b )

! Network structure
!-------------------------------------------------------------------------
! network(1)    network(2)      network(3)      network(4)      network(5)
!-------------------------------------------------------------------------
! Bias          Bias            Bias            Bias
! G1            H1-1            H2-1            H3-1            Ea
! G2            H1-2            H2-2            H3-2    
! ...           ...             ...             ...     
!
! #node = network(layer) + 1(bias)

!------------------------------------------------------------
! Total number of weight parameters
! Total number of nodes (excluding input nodes and output(atomic_ene) node )
!------------------------------------------------------------
!A
 if( io_nnp_a == 1 )then
   nweight_a = 0
   do i = 1, nlayer_a - 1
     nweight_a = nweight_a + ( network_a(i) + 1 ) * network_a(i+1)
   enddo
   allocate( weight_a( nweight_a ) )
 endif
!B
 if( io_nnp_b == 1 )then
   nweight_b = 0
   do i = 1, nlayer_b - 1
     nweight_b = nweight_b + ( network_b(i) + 1 ) * network_b(i+1)
   enddo
   allocate( weight_b( nweight_b ) )
 endif


 !----------------------
 ! Number of parameters
 !----------------------
   nparam = 0
   if( io_nnp_a == 1 ) nparam = nparam + nweight_a
   if( io_nnp_b == 1 ) nparam = nparam + nweight_b
   write(*,*) " Number of parameters : ", nparam
   write(*,*) " ===================================="
   if( io_nnp_a == 1 ) write(*,*) " nweight_a: ", nweight_a
   if( io_nnp_b == 1 ) write(*,*) " nweight_b: ", nweight_b
   write(*,*) " ===================================="
   write(*,*) " "

 allocate( weight( nparam ) )

 !************************************************************
 ! DEFINE INITIAL WEIGHT PARAMETERS
 !************************************************************
   ! From file
   open(unit=7,file="./data_weight/weight.dat",action="read",form="unformatted")
   read(7) weight
   close(7)

   !A,B
   if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
     do i = 1, nweight_a
       weight_a(i) = weight( i )
     enddo
     dummyi = nweight_a
     do i = 1, nweight_b
       weight_b(i) = weight( i + dummyi )
     enddo
   endif


 !***********
 ! READ DATA
 !***********
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
       call allocate_2( "g2_a_a",  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  eta2_a_b,  rs2_a_b, &
                        "g5_a_aa", eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ab", eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                        "g5_a_bb", eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                        num_g2_a_a,  num_g2_a_b, &
                        num_g5_a_aa, num_g5_a_ab, &
                        num_g5_a_bb, &
                        g2_a_a_maxmin,  g2_a_b_maxmin, &
                        g5_a_aa_maxmin, g5_a_ab_maxmin, &
                        g5_a_bb_maxmin )
       call allocate_2( "g2_b_a",  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  eta2_b_b,  rs2_b_b, &
                        "g5_b_aa", eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                        "g5_b_ab", eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                        "g5_b_bb", eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        num_g2_b_a,  num_g2_b_b, &
                        num_g5_b_aa, num_g5_b_ab, &
                        num_g5_b_bb, &
                        g2_b_a_maxmin,  g2_b_b_maxmin, &
                        g5_b_aa_maxmin, g5_b_ab_maxmin, &
                        g5_b_bb_maxmin )
 endif


 !*************
 ! LAMMPS file
 !*************
     open(unit=99,file="./lammps.nnp",action="write")
     
     call date_and_time(dummyc)
     write(99,'(2A8)')"# DATE: ",trim(adjustl(dummyc))
     write(99,*)" "
     write(99,'(A15)')"#num_elements 2"
     write(99,'(A16,F10.5)')"cutoff_distance ",rc
     write(99,*)" "
     write(99,'(A8,I5)')"num_g_a ",5
     write(99,*)" "
     write(99,'(A11,I5)')"num_g2_a_a ",  num_g_a(1)
     write(99,'(A11,I5)')"num_g2_a_b ",  num_g_a(2)
     write(99,'(A12,I5)')"num_g5_a_aa ", num_g_a(3)
     write(99,'(A12,I5)')"num_g5_a_ab ", num_g_a(4)
     write(99,'(A12,I5)')"num_g5_a_bb ", num_g_a(5)
     write(99,*)" "
     write(99,'(A8,I5)')"num_g_b ",5
     write(99,*)" "
     write(99,'(A11,I5)')"num_g2_b_a ",  num_g_b(1)
     write(99,'(A11,I5)')"num_g2_b_b ",  num_g_b(2)
     write(99,'(A12,I5)')"num_g5_b_aa ", num_g_b(3)
     write(99,'(A12,I5)')"num_g5_b_ab ", num_g_b(4)
     write(99,'(A12,I5)')"num_g5_b_bb ", num_g_b(5)
     write(99,*)" "
     write(99,'(A9,I5)')"nlayer_a ",   nlayer_a
     write(99,'(A10,I5)')"node_h1_a ",  network_a(2)
     write(99,'(A10,I5)')"node_h2_a ",  network_a(3)
     write(99,'(A11,I5)')"node_out_a ", network_a(4)
     write(99,*)" "
     write(99,'(A9,I5)')"nlayer_b ",   nlayer_b
     write(99,'(A10,I5)')"node_h1_b ",  network_b(2)
     write(99,'(A10,I5)')"node_h2_b ",  network_b(3)
     write(99,'(A11,I5)')"node_out_b ", network_b(4)
     write(99,*)" "


write(*,*)" Writing input parameters..."

write(99,*)" "
write(99,'(A12)')"param_g2_a_a"
do i = 1, num_g2_a_a
  write(99,'(2F15.10)') eta2_a_a(i),  rs2_a_a(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_a_b"
do i = 1, num_g2_a_b
  write(99,'(2F15.10)') eta2_a_b(i),  rs2_a_b(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_aa"
do i = 1, num_g5_a_aa
  write(99,'(4F15.10)') eta5_a_aa(i), theta5_a_aa(i), zeta5_a_aa(i), lambda5_a_aa(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_ab"
do i = 1, num_g5_a_ab
  write(99,'(4F15.10)') eta5_a_ab(i), theta5_a_ab(i), zeta5_a_ab(i), lambda5_a_ab(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_bb"
do i = 1, num_g5_a_bb
  write(99,'(4F15.10)') eta5_a_bb(i), theta5_a_bb(i), zeta5_a_bb(i), lambda5_a_bb(i)
enddo
write(99,*)" "


write(99,*)" "
write(99,'(A12)')"param_g2_b_a"
do i = 1, num_g2_b_a
  write(99,'(2F15.10)') eta2_b_a(i),  rs2_b_a(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_b_b"
do i = 1, num_g2_b_b
  write(99,'(2F15.10)') eta2_b_b(i),  rs2_b_b(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_aa"
do i = 1, num_g5_b_aa
  write(99,'(4F15.10)') eta5_b_aa(i), theta5_b_aa(i), zeta5_b_aa(i), lambda5_b_aa(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_ab"
do i = 1, num_g5_b_ab
  write(99,'(4F15.10)') eta5_b_ab(i), theta5_b_ab(i), zeta5_b_ab(i), lambda5_b_ab(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_bb"
do i = 1, num_g5_b_bb
  write(99,'(4F15.10)') eta5_b_bb(i), theta5_b_bb(i), zeta5_b_bb(i), lambda5_b_bb(i)
enddo
write(99,*)" "



write(*,*)" Writing max_min values..."

write(99,*)" "
write(99,'(A10)')"g2_a_a_max"
do i = 1, num_g2_a_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_a_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_a_min"
do i = 1, num_g2_a_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_a_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_b_max"
do i = 1, num_g2_a_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_b_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_b_min"
do i = 1, num_g2_a_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_b_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo


write(99,*)" "
write(99,'(A11)')"g5_a_aa_max"
do i = 1, num_g5_a_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_aa_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_aa_min"
do i = 1, num_g5_a_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_aa_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_a_ab_max"
do i = 1, num_g5_a_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ab_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_ab_min"
do i = 1, num_g5_a_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ab_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_a_bb_max"
do i = 1, num_g5_a_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bb_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_bb_min"
do i = 1, num_g5_a_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bb_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A10)')"g2_b_a_max"
do i = 1, num_g2_b_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_a_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_b_a_min"
do i = 1, num_g2_b_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_a_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_b_b_max"
do i = 1, num_g2_b_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_b_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_b_b_min"
do i = 1, num_g2_b_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_b_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo


write(99,*)" "
write(99,'(A11)')"g5_b_aa_max"
do i = 1, num_g5_b_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_aa_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_aa_min"
do i = 1, num_g5_b_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_aa_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_b_ab_max"
do i = 1, num_g5_b_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ab_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_ab_min"
do i = 1, num_g5_b_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ab_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_b_bb_max"
do i = 1, num_g5_b_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bb_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_bb_min"
do i = 1, num_g5_b_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bb_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo




write(*,*)" Writing weight parameters..."

write(99,*)" "
write(99,'(A15,I8)')"weight_a_params ", nweight_a
write(99,*)" "
write(99,'(A8)')"weight_a"
do i = 1, nweight_a
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')weight_a(i)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A15,I8)')"weight_b_params ", nweight_b
write(99,*)" "
write(99,'(A8)')"weight_b"
do i = 1, nweight_b
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')weight_b(i)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo



write(*,*)" "
write(*,*)" Finish"



end
