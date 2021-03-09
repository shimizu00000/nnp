program main_training
use allocarray
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank, dummyi
integer,allocatable,dimension(:) :: displs, recvcount
integer i, j, k, m
integer stat, counter, nfile, infile
integer imd, md_steps
double precision t0, t1
character(len=20) dummyc
character(len=120) tag, filename, dir_data, dir_sf
character(len=120),allocatable,dimension(:) :: tag_total
! num_data_trn   : number of training data
! num_data_tst   : number of test data
! num_data_total : num_data_trn + num_data_tst
! num_data_x     : number of data for each element
! list           : list of trn and tst data
integer num_data_trn, num_data_tst, num_data_total
integer num_data_trn_mpi, num_data_tst_mpi, num_data_total_mpi
integer num_data_a, num_data_b
integer num_data_a_trn, num_data_b_trn
integer num_data_a_trn_mpi, num_data_b_trn_mpi
integer,allocatable,dimension(:) :: &
  num_atom, num_atom_a, num_atom_b
integer natom_total, natom_total_mpi
integer natom_a_total, natom_ab_total
integer natom_b_total, natom_ba_total
integer natom_a_sum, natom_ab_sum
integer natom_b_sum, natom_ba_sum
integer,allocatable,dimension(:)   :: list, list_mpi
integer,allocatable,dimension(:) :: info_io, info_io_mpi
double precision seed(2)
double precision,allocatable,dimension(:) :: ransu
! nelem   : total number of elements
! natom   : total number of atoms
! natom_x : distribute natom to each element
integer nelem, natom
integer natom_a, natom_b
integer io_a, io_b
integer counter_ene, counter_io
integer counter_a, counter_b
integer counter_natom
integer counter_natom_a,   counter_natom_ab
integer counter_natom_b,   counter_natom_ba
integer counter_deriv_a,   counter_deriv_ab
integer counter_deriv_b,   counter_deriv_ba
character(len=4) elem_a, elem_b

! Variables for L-BFGS-B
!------------------------
! nparam : number of parameters
! factr, pgtol : convergence criteria
integer nparam
!integer,parameter :: memory = 5
!integer,parameter :: iprint = 1
integer,parameter :: dp = kind( 1.0d0 )
!double precision,parameter :: factr = 1.0d+7
!double precision,parameter :: pgtol = 1.0d-5
character(len=60) :: task, csave
logical           :: lsave(4)
integer           :: isave(44), memory, iprint
double precision  :: dsave(29), factr, pgtol
integer,allocatable,dimension(:) :: nbd, iwa
double precision,allocatable,dimension(:) :: weight, l, u, wa, g
! f : error function
double precision  :: f, f_mpi
! g  : diffential of error function with respect to weight parameters
! g_x: distribute g to each element
double precision,allocatable,dimension(:) :: &
  g_a, g_e_a, g_f_a, g_a_mpi, &
  g_b, g_e_b, g_f_b, g_b_mpi

! Variables for NNP
!-------------------
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
double precision energy, energy0
double precision,allocatable,dimension(:) :: &
  energy0_all, energy_all
double precision,allocatable,dimension(:,:) :: force, force0
double precision,allocatable,dimension(:) :: &
  force_all, force0_all
double precision alpha, beta
! mse_energy_trn : mean square error of energy in trn
! mse_energy_tst : mean square error of energy in tst
double precision mse_energy_trn, mse_energy_tst
double precision mse_energy_trn_mpi, mse_energy_tst_mpi
double precision mse_force_trn, mse_force_tst
double precision mse_force_trn_mpi, mse_force_tst_mpi
! io_weight   : weight from scratch or file
! percent_tst : fraction of tst data
integer io_weight, percent_tst
! rc : cutoff radius
! rn : cutoff radius + rn
! neighbor : number of neighboring atoms
! phi : activation function
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
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_a_a,  g2_deriv_a_b, &
  g5_deriv_a_aa, g5_deriv_a_ab, &
  g5_deriv_a_bb
!B
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_b_a,  g2_deriv_b_b, &
  g5_deriv_b_aa, g5_deriv_b_ab, &
  g5_deriv_b_bb
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


 !MPI INIT
 call mpi_init( ierr )
 call mpi_comm_size( mpi_comm_world, num_mpi, ierr )
 call mpi_comm_rank( mpi_comm_world, myrank, ierr )
 write(*,*) " NUM_MPI: ", myrank, num_mpi
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 ) call cpu_time( t0 )

 !*************************************************************
 ! READ input_nnp.dat FILE
 !*************************************************************
 if( myrank == 0 )then

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
        dir_data, dir_sf, &
        io_weight, percent_tst, alpha, beta, &
        factr, pgtol, memory, iprint )

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

 endif ! myrank == 0

 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_bcast( elem_a, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_b, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( rc, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_data, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_sf, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( percent_tst, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( alpha, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( beta, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g_a, 5, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g_b, 5, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( nlayer_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( nlayer_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( network_a( nlayer_a ) )
 allocate( network_b( nlayer_b ) )

 !--------------------------------
 ! Read number of nodes per layer
 !--------------------------------
 if( myrank == 0 )then
   call read_layer_node_frex2( network_a, network_b, &
                               nlayer_a,  nlayer_b,&
                               io_nnp_a,  io_nnp_b )
 endif ! myrank == 0

 call mpi_bcast( network_a, nlayer_a, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( network_b, nlayer_b, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( io_nnp_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( io_nnp_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

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
   allocate(      g_a( nweight_a ), g_e_a( nweight_a ), g_f_a( nweight_a ), &
              g_a_mpi( nweight_a ) )
   nhidden_a = 0
   do i = 2, nlayer_a - 1
     nhidden_a = nhidden_a + ( network_a(i) + 1 )
   enddo
 endif
!B
 if( io_nnp_b == 1 )then
   nweight_b = 0
   do i = 1, nlayer_b - 1
     nweight_b = nweight_b + ( network_b(i) + 1 ) * network_b(i+1)
   enddo
   allocate( weight_b( nweight_b ) )
   allocate(      g_b( nweight_b ), g_e_b( nweight_b ), g_f_b( nweight_b ), &
              g_b_mpi( nweight_b ) )
   nhidden_b = 0
   do i = 2, nlayer_b - 1
     nhidden_b = nhidden_b + ( network_b(i) + 1 )
   enddo
 endif
 call mpi_barrier( mpi_comm_world, ierr )

 !********************
 ! PARAMETER FOR BFGS
 !********************
 !----------------------
 ! Number of parameters
 !----------------------
 if( myrank == 0 )then
   nparam = 0
   if( io_nnp_a == 1 ) nparam = nparam + nweight_a
   if( io_nnp_b == 1 ) nparam = nparam + nweight_b
   write(*,*) " Number of parameters : ", nparam
   write(*,*) " ===================================="
   if( io_nnp_a == 1 ) write(*,*) " nweight_a: ", nweight_a
   if( io_nnp_b == 1 ) write(*,*) " nweight_b: ", nweight_b
   write(*,*) " ===================================="
   write(*,*) " "
 endif ! myrank == 0
 call mpi_bcast( nparam, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 )then
   allocate( nbd( nparam ), l( nparam ), u( nparam ) )
   allocate( iwa( 3*nparam ) )
   allocate( wa( 2*memory*nparam + 5*nparam + 11*memory*memory + 8*memory ) )
 endif ! myrank == 0
 allocate( g( nparam ), weight( nparam ) )
 call mpi_barrier( mpi_comm_world, ierr )

 !------------------------
 ! Upper and lower bounds
 !------------------------
 if( myrank == 0 )then
   do i = 1, nparam
     nbd(i) =  0 ! 0 nbounded, 2 both bounded
     l(i)   = -1.0d2
     u(i)   =  1.0d2 
   enddo
 endif ! myrank == 0
!************************************************************
! DEFINE INITIAL WEIGHT PARAMETERS
!************************************************************
 if( myrank == 0 )then
   ! From scrach
   if( io_weight == 0 )then
     allocate ( ransu( nparam + 2 ) )
     call pre_random
     call random_number( ransu )
     do  i = 1, nparam
       weight(i) = 2.0d0 * ransu(i) - 1.0d0
     enddo
     deallocate( ransu )
   ! From file
   elseif( io_weight == 1 )then
     open(unit=7,file="./data_weight/weight.dat",action="read",form="unformatted")
     read(7) weight
     close(7)
   ! Bug
   else
     write(*,*) "Error(main): weight parameter setting."
     stop
   endif
 endif ! myrank == 0

 call mpi_bcast( weight, nparam, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


 !*************************************************************
 ! COUNTING NUMBER OF DATA from input_tag_training.dat FILE
 !*************************************************************
 if( myrank == 0 )then
   open(unit=99,file="input_tag_training.dat",action="read")
   counter = 0
   rewind(99)
   do
     read( 99, '(A5)', iostat=stat ) dummyc
     if( stat /= 0 ) exit
     if( trim( dummyc ) /= '' ) counter = counter + 1
   enddo
 endif ! myrank == 0
 call mpi_bcast( counter, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


 !--------------------------------------
 ! DO LOOP : nfile (input_tag_training)
 !--------------------------------------
 allocate( tag_total( counter ) )
 if( myrank == 0 )then
   rewind(99)
   do nfile = 1, counter
     read(99,*) tag
     tag_total( nfile ) = trim(adjustl(tag))
   enddo
 endif ! myrank == 0
 close(99) ! input_tag_training.dat

 call mpi_bcast( tag_total, counter*120, mpi_character, 0, mpi_comm_world, ierr )
 !bcast character: #data*len
 call mpi_barrier( mpi_comm_world, ierr )

 num_data_total = 0
 if( io_nnp_a == 1 ) num_data_a = 0
 if( io_nnp_b == 1 ) num_data_b = 0

 natom_total = 0
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   natom_a_total   = 0 ; natom_a_sum   = 0
   natom_ab_total  = 0 ; natom_ab_sum  = 0
   natom_b_total   = 0 ; natom_b_sum   = 0
   natom_ba_total  = 0 ; natom_ba_sum  = 0
   !total: natom_x (g2, g5)
   !sum: natom_x * natom (g2_deriv, g5_deriv)
 endif

 do nfile = 1, ceiling( dble(counter)/dble(num_mpi) )
   !infile = nfile + myrank*ceiling( dble(counter)/dble(num_mpi) )
   infile = ( myrank + 1 ) + ( ( nfile - 1 ) * num_mpi )
   if( infile <= counter )then

     tag = tag_total( infile )
     ! info (binary), READ
     !filename = '../step2_data_binary/data/info_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(tag))//'.dat'
     open(unit=1,file=filename,action="read",form="unformatted")
     !---------------------------------
     ! READ info_tag.dat FILE (BINARY)
     !---------------------------------
     ! 2 elements
     call read_io2( natom,   nelem, &
                    elem_a,  elem_b, &
                    natom_a, natom_b, &
                    io_a,    io_b )
     ! MD steps
     read(1) md_steps

     num_data_total = num_data_total + md_steps
     natom_total = natom_total + natom * md_steps
     !A
     if(     io_a == 1 .and. io_b == 0 )then
       num_data_a    = num_data_a    + md_steps
       natom_a_total = natom_a_total + natom_a * md_steps
       natom_a_sum   = natom_a_sum   + natom_a * natom * md_steps
     !B
     elseif( io_a == 0 .and. io_b == 1 )then
       num_data_b    = num_data_b    + md_steps
       natom_b_total = natom_b_total + natom_b * md_steps
       natom_b_sum   = natom_b_sum   + natom_b * natom * md_steps
     !A,B
     elseif( io_a == 1 .and. io_b == 1 )then
       num_data_a     = num_data_a     + md_steps
       natom_a_total  = natom_a_total  + natom_a * md_steps
       natom_a_sum    = natom_a_sum    + natom_a * natom * md_steps
       natom_ab_total = natom_ab_total + natom_a * md_steps
       natom_ab_sum   = natom_ab_sum   + natom_a * natom * md_steps
       num_data_b     = num_data_b     + md_steps
       natom_b_total  = natom_b_total  + natom_b * md_steps
       natom_b_sum    = natom_b_sum    + natom_b * natom * md_steps
       natom_ba_total = natom_ba_total + natom_b * md_steps
       natom_ba_sum   = natom_ba_sum   + natom_b * natom * md_steps
     else
       write(*,*)"Error(main): counting data io."
       stop
     endif
  
     close(1) ! info(binary)

   endif ! infile <= counter
 enddo ! nfile
 call mpi_barrier( mpi_comm_world, ierr )
 !--------------
 ! END COUNTING 
 !--------------

 !***********
 ! READ DATA
 !***********
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   call alloc_2( num_atom, num_data_total, &
                 num_atom_a, num_atom_b, &
                 num_data_a, num_data_b, &
                 natom_a_total, natom_ab_total, &
                 natom_b_total, natom_ba_total, &
                 natom_a_sum,   natom_ab_sum, &
                 natom_b_sum,   natom_ba_sum, &
                 "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                 "g2_a_b",  g2_a_b,  g2_deriv_a_b,  g2_a_b_maxmin,  num_g2_a_b, &
                 "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                 "g5_a_ab", g5_a_ab, g5_deriv_a_ab, g5_a_ab_maxmin, num_g5_a_ab, &
                 "g5_a_bb", g5_a_bb, g5_deriv_a_bb, g5_a_bb_maxmin, num_g5_a_bb, &
                 "g2_b_a",  g2_b_a,  g2_deriv_b_a,  g2_b_a_maxmin,  num_g2_b_a, &
                 "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                 "g5_b_aa", g5_b_aa, g5_deriv_b_aa, g5_b_aa_maxmin, num_g5_b_aa, &
                 "g5_b_ab", g5_b_ab, g5_deriv_b_ab, g5_b_ab_maxmin, num_g5_b_ab, &
                 "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb )
 endif

 ! 2 elements system
 allocate( info_io( num_data_total*2 ) )
 allocate( energy0_all( num_data_total ) )
 allocate( force0_all( natom_total*3 ) )

 !--------------------------------------
 ! DO LOOP : nfile (input_tag_training)
 !--------------------------------------
 if( myrank == 0 )then
   write(*,*)" Reading data"
   write(*,*)" "
 endif ! myrank == 0

 counter_io = 0
 counter_ene = 0
 if( io_nnp_a == 1 ) counter_a = 0
 if( io_nnp_b == 1 ) counter_b = 0

 counter_natom = 0
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   counter_natom_a   = 0 ; counter_deriv_a   = 0
   counter_natom_ab  = 0 ; counter_deriv_ab  = 0
   counter_natom_b   = 0 ; counter_deriv_b   = 0
   counter_natom_ba  = 0 ; counter_deriv_ba  = 0
   ! counter_deriv_x: natom_x * natom (g2_deriv & g5_deriv)
 endif

 do nfile = 1, ceiling( dble(counter)/dble(num_mpi) )
   !infile = nfile + myrank*ceiling( dble(counter)/dble(num_mpi) )
   infile = ( myrank + 1 ) + ( ( nfile - 1 ) * num_mpi )
   if( infile <= counter )then

     tag = tag_total( infile )
     ! info (binary), READ
     !filename = '../step2_data_binary/data/info_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(tag))//'.dat'
     open(unit=1,file=filename,action="read",form="unformatted")
     ! energies (binary), READ
     !filename = '../step2_data_binary/data/energies_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/energies_'//trim(adjustl(tag))//'.dat'
     open(unit=3,file=filename,action="read",form="unformatted")
     ! forces (binary), READ
     !filename = '../step2_data_binary/data/forces_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/forces_'//trim(adjustl(tag))//'.dat'
     open(unit=4,file=filename,action="read",form="unformatted")
     ! sf (binary), READ
     !write( dummyc, * ) int(rc)
     !filename = '../step3_sf_binary/data_sf/rc_'//trim(adjustl(dummyc))//&
     !           '/sf_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_sf))//'/sf_'//trim(adjustl(tag))//'.dat'
     open(unit=11,file=filename,action="read",form="unformatted")
     ! sf_deriv (binary), READ
     !filename = '../step3_sf_binary/data_sf/rc_'//trim(adjustl(dummyc))//&
     !           '/sf_deriv_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_sf))//'/sf_deriv_'//trim(adjustl(tag))//'.dat'
     open(unit=12,file=filename,action="read",form="unformatted")

     !---------------------------------
     ! READ info_tag.dat FILE (BINARY)
     !---------------------------------
     ! 2 elements system
     call read_io2( natom,   nelem, &
                    elem_a,  elem_b, &
                    natom_a, natom_b, &
                    io_a,    io_b )
     ! Force
     allocate( force0( natom, 3 ) )

     ! MD steps
     read(1) md_steps
     do imd = 1, md_steps

       read(3) energy0
       counter_ene = counter_ene + 1
       energy0_all( counter_ene ) = energy0
       num_atom( counter_ene ) = natom

       read(4) force0
       do i = 1, natom
         force0_all( counter_natom + 3*(i-1) + 1 ) = force0(i,1)
         force0_all( counter_natom + 3*(i-1) + 2 ) = force0(i,2)
         force0_all( counter_natom + 3*(i-1) + 3 ) = force0(i,3)
       enddo
       counter_natom = counter_natom + 3*natom
  
       info_io( counter_io + 1 ) = io_a
       info_io( counter_io + 2 ) = io_b
       counter_io = counter_io + 2

       ! 1 element
       !A
       if(     io_a == 1 .and. io_b == 0 )then
         call read_data_1( num_atom_a, num_data_a, natom_a_total, natom_a_sum, &
                           counter_a, counter_natom_a, counter_deriv_a, &
                           imd, natom_a, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa )
       !B
       elseif( io_a == 0 .and. io_b == 1 )then
         call read_data_1( num_atom_b, num_data_b, natom_b_total, natom_b_sum, &
                           counter_b, counter_natom_b, counter_deriv_b, &
                           imd, natom_b, natom, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb )
       !A,B
       elseif( io_a == 1 .and. io_b == 1 )then
         call read_data_2( num_atom_a, num_atom_b, &
                           num_data_a, num_data_b, &
                           natom_a_total, natom_ab_total, &
                           natom_b_total, natom_ba_total, &
                           natom_a_sum, natom_ab_sum, &
                           natom_b_sum, natom_ba_sum, &
                           counter_a, counter_b, &
                           counter_natom_a, counter_natom_ab, &
                           counter_natom_b, counter_natom_ba, &
                           counter_deriv_a, counter_deriv_ab, &
                           counter_deriv_b, counter_deriv_ba, &
                           imd, natom_a, natom_b, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g2_a_b",  g2_a_b,  g2_deriv_a_b,  g2_a_b_maxmin,  num_g2_a_b, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                           "g5_a_ab", g5_a_ab, g5_deriv_a_ab, g5_a_ab_maxmin, num_g5_a_ab, &
                           "g5_a_bb", g5_a_bb, g5_deriv_a_bb, g5_a_bb_maxmin, num_g5_a_bb, &
                           "g2_b_a",  g2_b_a,  g2_deriv_b_a,  g2_b_a_maxmin,  num_g2_b_a, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g5_b_aa", g5_b_aa, g5_deriv_b_aa, g5_b_aa_maxmin, num_g5_b_aa, &
                           "g5_b_ab", g5_b_ab, g5_deriv_b_ab, g5_b_ab_maxmin, num_g5_b_ab, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb )
       else
         write(*,*)"Error(main): read data."
         stop
       endif

     enddo ! imd

     deallocate( force0 )

     if( io_a == 1 ) counter_a = counter_a + md_steps
     if( io_b == 1 ) counter_b = counter_b + md_steps

     close(1)  ! info
     close(3)  ! energies
     close(4)  ! forces
     close(11) ! sf
     close(12) ! sf_deriv

   endif ! infile <= counter
 enddo ! nfile
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 )then
   write(*,*)" Finish reading data"
   write(*,*)" "
 endif ! myrank == 0
 !---------------
 ! END READ DATA
 !---------------

 !---------------------
 ! List of trn and tst
 !---------------------
 allocate( list( num_data_total ) )
 list = 0
 num_data_trn = 0
 if( io_nnp_a == 1 ) num_data_a_trn = 0
 if( io_nnp_b == 1 ) num_data_b_trn = 0
 call pre_random

 counter_io = 0
 do i = 1, num_data_total
   call random_number( seed )
   seed(2) = seed(2) * 100.0d0
   ! list = 1: training, list = 0: test
   if( seed(2) > percent_tst )then
     list(i) = 1
     num_data_trn = num_data_trn + 1
     if( info_io( counter_io + 1 ) == 1 ) num_data_a_trn = num_data_a_trn + 1
     if( info_io( counter_io + 2 ) == 1 ) num_data_b_trn = num_data_b_trn + 1
   endif
   counter_io = counter_io + 2
 enddo ! i

 call mpi_reduce( natom_total, natom_total_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_total, num_data_total_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_trn, num_data_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_a_trn, num_data_a_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_b_trn, num_data_b_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_data_total_mpi, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


 ! DISPLS for MPI_GATHERV
 allocate( displs( num_mpi ), recvcount( num_mpi ) )
 displs = 0 ; recvcount = 0
 if( myrank == 0 )then
   m = 0
   do j = 1, num_mpi
     do i = 1, ceiling( dble(counter)/dble(num_mpi) )
       m = j + ( ( i - 1 ) * num_mpi )
       !m = m + 1
       if( m <= counter ) recvcount(j) = recvcount(j) + 1
     enddo
   enddo
   write(*,*) " Data devided into (MPI): "
   write(*,*) recvcount
   write(*,*) " "
   do i = 1, num_mpi
     j = num_mpi + 1 - i
     do k = 1, j - 1
       displs(j) = displs(j) + recvcount(k)
     enddo
   enddo
   displs(1) = 0
 endif ! myrank == 0
 call mpi_bcast( displs, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( recvcount, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( list_mpi( num_data_total_mpi ) )
 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_gatherv( list, num_data_total, mpi_integer, &
                   list_mpi, recvcount, displs, &
                   mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 deallocate( displs, recvcount )


 ! DISPLS for MPI_GATHERV
 allocate( displs( num_mpi ), recvcount( num_mpi ) )
 displs = 0 ; recvcount = 0
 if( myrank == 0 )then
   m = 0
   do j = 1, num_mpi
     do i = 1, ceiling( dble(counter)/dble(num_mpi) )
       m = j + ( ( i - 1 ) * num_mpi )
       !m = m + 1
       if( m <= counter ) recvcount(j) = recvcount(j) + 2
     enddo
   enddo
   do i = 1, num_mpi
     j = num_mpi + 1 - i
     do k = 1, j - 1
       displs(j) = displs(j) + recvcount(k)
     enddo
   enddo
   displs(1) = 0
 endif ! myrank == 0
 call mpi_bcast( displs, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( recvcount, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( info_io_mpi( num_data_total_mpi*2 ) )
 call mpi_gatherv( info_io, 2*num_data_total, mpi_integer, &
                   info_io_mpi, recvcount, displs, &
                   mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 deallocate( displs, recvcount )


 if( myrank == 0 )then
   num_data_tst_mpi = num_data_total_mpi - num_data_trn_mpi
   write(*,*) " Total number of data: ",    num_data_total_mpi
   write(*,*) " ======================================"
   write(*,*) " Number of training data: ", num_data_trn_mpi
   write(*,*) " Number of test data: ",     num_data_tst_mpi
   write(*,*) " ======================================"
   write(*,*) percent_tst, "% is allocated for test set."
   write(*,*) " "
1112 format(I10,I3,A3,A4,I2,A4,I2)
   open(unit=98,file="list.dat",action="write")
   write(98,*) "# total number of data"
   write(98,*) num_data_total_mpi
   write(98,*) "# 1: training, 0: test, A, B 1: on, 0:off"
   counter_io = 0
   do i = 1, num_data_total_mpi
     write(98,1112) i, list_mpi(i), " : ", &
                    "  A ", info_io_mpi(counter_io+1), "  B ", info_io_mpi(counter_io+2)
     counter_io = counter_io + 2
   enddo
   close(98)

   open(unit=98,file="tag_out.dat",action="write")
   do i = 1, num_mpi
     do j = 1, ceiling(dble(counter)/dble(num_mpi))
       infile = i + ( ( j - 1 ) * num_mpi ) 
       if( infile <= counter )then
         write(98,*) trim(adjustl( tag_total(infile)) )
       endif
     enddo
   enddo
   close(98)

 endif ! myrank == 0


 !***************
 ! L-BFGS-B LOOP
 !***************
 if( myrank == 0 )then
   open(unit=13,file="./mse_training.dat",action="write")
   open(unit=14,file="./mse_test.dat",action="write")
   write(13,*) "# MSE energy (eV^2/atom^2) & force in training (eV^2/ang^2)"
   write(14,*) "# MSE energy (eV^2/atom^2) & force in test (eV^2/ang^2)"
 endif ! myrank == 0

 task = 'START'
 do while( task(1:2) .eq. 'FG' .or. task .eq. 'NEW_X' .or. task .eq. 'START' )
   !---------------
   ! L-BFGS-B code
   !---------------
   if( myrank == 0 )then
     call setulb( nparam, memory, weight, l, u, nbd, f_mpi, g, factr, pgtol, &
                  wa, iwa, task, iprint, csave, lsave, isave, dsave )
   endif ! myrank == 0
   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( weight, nparam, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_bcast( task, 60, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   if( task(1:2) .eq. 'FG' )then

     ! Evaluate f and g
     f = 0.0d0 ; g = 0.0d0
     if( io_nnp_a == 1 ) g_a = 0.0d0 
     if( io_nnp_b == 1 ) g_b = 0.0d0 

     ! Evaluate RMSE of energy
     mse_energy_trn = 0.0d0 ; mse_energy_tst = 0.0d0
     mse_force_trn  = 0.0d0 ; mse_force_tst  = 0.0d0

     !------------------------------
     ! Distribute weight parameters
     !------------------------------
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

     !*************************************************************
     ! DO LOOP MD_STEPS
     !*************************************************************
     counter_io = 0
     if( io_nnp_a == 1 ) counter_a = 0
     if( io_nnp_b == 1 ) counter_b = 0
     counter_natom = 0
     !A,B
     if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
       counter_natom_a   = 0 ; counter_deriv_a   = 0
       counter_natom_ab  = 0 ; counter_deriv_ab  = 0
       counter_natom_b   = 0 ; counter_deriv_b   = 0
       counter_natom_ba  = 0 ; counter_deriv_ba  = 0
       ! g2_deriv & g5_deriv
     endif

     do imd = 1, num_data_total

       energy0 = energy0_all( imd )
       io_a = info_io( counter_io + 1 )
       io_b = info_io( counter_io + 2 )
       counter_io = counter_io + 2

       !---------------
       ! TRAINING DATA
       !---------------
       if( list( imd ) == 1 )then
         !------------------------------------------------------------
         ! Energy(f) &
         ! Diffential of energy with respsct to weight parameters (g)
         !------------------------------------------------------------
         !*** A,B NNP only ***
         !A,B
         !elseif( io_a == 1 .and. io_b == 1 )then
         natom_a = 0; natom_b = 0
         if( io_a == 1 )then
           counter_a = counter_a + 1
           natom_a   = num_atom_a( counter_a )
           allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
         endif
         if( io_b == 1 )then
           counter_b = counter_b + 1
           natom_b   = num_atom_b( counter_b )
           allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
         endif

         natom = num_atom( imd )

         allocate( force( natom, 3 ), force0( natom, 3 ) )
         do i = 1, natom
           force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
           force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
           force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
         enddo

         call train2ab( io_a, io_b, &
                        natom, natom_a, natom_b, &
                        natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                        natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                        natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                        natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                        g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                        g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                        g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                        g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                        g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                        network_a, nlayer_a, num_g_a, &
                        weight_a, nweight_a, hidden_a, nhidden_a, &
                        g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                        g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                        g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                        g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                        g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                        network_b, nlayer_b, num_g_b, &
                        weight_b, nweight_b, hidden_b, nhidden_b, &
                        beta, g_e_a, g_f_a, g_e_b, g_f_b, &
                        energy, energy0, force, force0 )

         if( io_a == 1 )then
           deallocate( hidden_a, atomic_ene_a )
           g_a =   g_a &
                 + alpha*( g_e_a/dble(natom_a) ) + beta*( g_f_a/dble(3*natom_a) )
         endif
         if( io_b == 1 )then
           deallocate( hidden_b, atomic_ene_b )
           g_b =   g_b &
                 + alpha*( g_e_b/dble(natom_b) ) + beta*( g_f_b/dble(3*natom_b) )
         endif
         !-----------
         ! END f & g
         !-----------

         f = f + alpha*( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           f = f + beta*( force(i,1)-force0(i,1) )**2/dble(3*natom)
           f = f + beta*( force(i,2)-force0(i,2) )**2/dble(3*natom)
           f = f + beta*( force(i,3)-force0(i,3) )**2/dble(3*natom)
         enddo

         mse_energy_trn =   mse_energy_trn &
                          + ( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           mse_force_trn = mse_force_trn + ( ( force(i,1)-force0(i,1) )**2/dble(3*natom) )
           mse_force_trn = mse_force_trn + ( ( force(i,2)-force0(i,2) )**2/dble(3*natom) )
           mse_force_trn = mse_force_trn + ( ( force(i,3)-force0(i,3) )**2/dble(3*natom) )
         enddo

       !-----------
       ! TEST DATA
       !-----------
       elseif( list( imd ) == 0 )then

         !--------
         ! Energy
         !--------
         !A,B
         natom_a = 0; natom_b = 0
         if( io_a == 1 )then
           counter_a = counter_a + 1
           natom_a   = num_atom_a( counter_a )
           allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
         endif
         if( io_b == 1 )then
           counter_b = counter_b + 1
           natom_b   = num_atom_b( counter_b )
           allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
         endif

         natom = num_atom( imd )

         allocate( force( natom, 3 ), force0( natom, 3 ) )
         do i = 1, natom
           force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
           force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
           force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
         enddo

         call valid2ab( io_a, io_b, &
                        natom, natom_a, natom_b, &
                        natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                        natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                        natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                        natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                        g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                        g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                        g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                        g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                        g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                        network_a, nlayer_a, num_g_a, &
                        weight_a, nweight_a, hidden_a, nhidden_a, &
                        g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                        g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                        g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                        g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                        g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                        network_b, nlayer_b, num_g_b, &
                        weight_b, nweight_b, hidden_b, nhidden_b, &
                        energy, force )

         if( io_a == 1 ) deallocate( hidden_a, atomic_ene_a )
         if( io_b == 1 ) deallocate( hidden_b, atomic_ene_b )

         mse_energy_tst =   mse_energy_tst &
                          + ( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           mse_force_tst = mse_force_tst + ( ( force(i,1)-force0(i,1) )**2/dble(3*natom) )
           mse_force_tst = mse_force_tst + ( ( force(i,2)-force0(i,2) )**2/dble(3*natom) )
           mse_force_tst = mse_force_tst + ( ( force(i,3)-force0(i,3) )**2/dble(3*natom) )
         enddo

       endif ! list, training or test

       counter_natom = counter_natom + 3*natom
       deallocate( force, force0 )

     enddo ! do imd = 1, num_data_total ! md_steps
     call mpi_barrier( mpi_comm_world, ierr )
     !------------
     ! END do imd
     !------------

     call mpi_reduce( f, f_mpi, 1, mpi_double_precision, mpi_sum, 0, &
                      mpi_comm_world, ierr )
     call mpi_reduce( mse_energy_trn, mse_energy_trn_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_energy_tst, mse_energy_tst_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_force_trn, mse_force_trn_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_force_tst, mse_force_tst_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_barrier( mpi_comm_world, ierr )


     if( myrank == 0 )then

       f_mpi = f_mpi / dble( num_data_trn_mpi )

       mse_energy_trn_mpi = mse_energy_trn_mpi / dble( num_data_trn_mpi )
       mse_force_trn_mpi  = mse_force_trn_mpi  / dble( num_data_trn_mpi )
       write(13,fmt='(2F25.12)') mse_energy_trn_mpi, mse_force_trn_mpi ! mse_training

       mse_energy_tst_mpi = mse_energy_tst_mpi / dble( num_data_tst_mpi )
       mse_force_tst_mpi  = mse_force_tst_mpi  / dble( num_data_tst_mpi )
       write(14,fmt='(2F25.12)') mse_energy_tst_mpi, mse_force_tst_mpi ! mse_test

     endif ! myrank == 0


     !----------------------------
     ! Concatenate g_x parameters
     !----------------------------
     call mpi_reduce( g_a, g_a_mpi, nweight_a, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( g_b, g_b_mpi, nweight_b, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_barrier( mpi_comm_world, ierr )

     if( myrank == 0 )then
       !A,B
       if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
         if( num_data_a_trn_mpi /= 0 )then
           g_a_mpi = g_a_mpi / dble( num_data_a_trn_mpi )
           do i = 1, nweight_a
             g(i) = g_a_mpi(i)
           enddo
         endif
         if( num_data_b_trn_mpi /= 0 )then
           g_b_mpi = g_b_mpi / dble( num_data_b_trn_mpi )
           do i = 1, nweight_b
             g( nweight_a + i ) = g_b_mpi(i)
           enddo
         endif
       endif
     endif ! myrank == 0

   endif ! if( task(1:2) .eq. 'FG' )then

   ! Temporary write weight parameters
   if( myrank == 0 )then
     open(unit=7,file="./data_weight/weight.dat",action="write",form="unformatted")
     write(7) weight
     close(7)
   endif ! myrank == 0

 enddo !end of loop do while
 !---------------
 ! END BFGS LOOP
 !---------------


 if( myrank == 0 )then

   mse_energy_trn_mpi = dsqrt( mse_energy_trn_mpi )
   write(*,*)  "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
   write(13,*) "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
   mse_force_trn_mpi = dsqrt( mse_force_trn_mpi )
   write(*,*)  "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"
   write(13,*) "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"

   mse_energy_tst_mpi = dsqrt( mse_energy_tst_mpi )
   write(*,*)  "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
   write(14,*) "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
   mse_force_tst_mpi = dsqrt( mse_force_tst_mpi )
   write(*,*)  "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"
   write(14,*) "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"

   close(13) ! mse_trn.dat
   close(14) ! mse_tst.dat

   !-------------------------
   ! Write weight parameters
   !-------------------------
   open(unit=7,file="./data_weight/weight.dat",action="write",form="unformatted")
   write(7) weight
   close(7)

 endif ! myrank == 0


 !***************
 ! FINALIZE DATA
 !***************
 !------------------------------
 ! Distribute weight parameters
 !------------------------------
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   do i = 1, nweight_a
     weight_a(i) = weight(i)
   enddo
   dummyi = nweight_a
   do i = 1, nweight_b
     weight_b(i) = weight( dummyi + i )
   enddo
 endif

 !******************
 ! DO LOOP MD_STEPS
 !******************
 counter_io = 0
 if( io_nnp_a == 1 ) counter_a = 0
 if( io_nnp_b == 1 ) counter_b = 0
 counter_natom = 0
 !A,B
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   counter_natom_a   = 0 ; counter_deriv_a   = 0
   counter_natom_ab  = 0 ; counter_deriv_ab  = 0
   counter_natom_b   = 0 ; counter_deriv_b   = 0
   counter_natom_ba  = 0 ; counter_deriv_ba  = 0
   ! g2_deriv & g5_deriv
 endif

 allocate( energy_all( num_data_total ) )
 allocate( force_all( natom_total*3 ) )

 do imd = 1, num_data_total

   energy0 = energy0_all( imd )

   io_a = info_io( counter_io + 1 )
   io_b = info_io( counter_io + 2 )
   counter_io = counter_io + 2

   !A,B
   natom_a = 0; natom_b = 0
   if( io_a == 1 )then
     counter_a = counter_a + 1
     natom_a   = num_atom_a( counter_a )
     allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
   endif
   if( io_b == 1 )then
     counter_b = counter_b + 1
     natom_b   = num_atom_b( counter_b )
     allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
   endif

   natom = num_atom( imd )
   if( natom /= natom_a + natom_b ) write(*,*)" natom does not match"
   if( natom /= natom_a + natom_b ) stop

   allocate( force( natom, 3 ), force0( natom, 3 ) )
   do i = 1, natom
     force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
     force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
     force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
   enddo

   call valid2ab( io_a, io_b, &
                  natom, natom_a, natom_b, &
                  natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                  natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                  natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                  natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                  g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                  g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                  g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                  g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                  g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                  network_a, nlayer_a, num_g_a, &
                  weight_a, nweight_a, hidden_a, nhidden_a, &
                  g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                  g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                  g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                  g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                  g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                  network_b, nlayer_b, num_g_b, &
                  weight_b, nweight_b, hidden_b, nhidden_b, &
                  energy, force )

   if( io_a == 1 ) deallocate( hidden_a, atomic_ene_a )
   if( io_b == 1 ) deallocate( hidden_b, atomic_ene_b )

   energy_all(  imd ) = energy  / dble( natom )
   energy0_all( imd ) = energy0 / dble( natom )

   do i = 1, natom
     force_all( counter_natom + 3*(i-1) + 1 ) = force(i,1)
     force_all( counter_natom + 3*(i-1) + 2 ) = force(i,2)
     force_all( counter_natom + 3*(i-1) + 3 ) = force(i,3)
   enddo

   counter_natom = counter_natom + 3*natom

   deallocate( force, force0 )

 enddo ! do imd = 1, num_data_total
 !-----------------
 ! END DO LOOP imd
 !-----------------
 call mpi_barrier( mpi_comm_world, ierr )


 write(filename,'("energy_training",i4.4,".dat")') myrank
 open(unit=13,file=filename,action="write")
 write(filename,'("energy_test",i4.4,".dat")') myrank
 open(unit=14,file=filename,action="write")
 write(13,*) "# DFT(eV/atom) NNP(eV/atom) Natom"
 write(14,*) "# DFT(eV/atom) NNP(eV/atom) Natom"

 write(filename,'("force_training",i4.4,".dat")') myrank
 open(unit=15,file=filename,action="write")
 write(filename,'("force_test",i4.4,".dat")') myrank
 open(unit=16,file=filename,action="write")
 write(15,*) "# DFT(eV/ang) NNP(eV/ang) Natom"
 write(16,*) "# DFT(eV/ang) NNP(eV/ang) Natom"


 counter_natom = 0

 do imd = 1, num_data_total

   natom = num_atom( imd )

   !---------------
   ! TRAINING DATA
   !---------------
   if( list( imd ) == 1 )then

     write(13,*) energy0_all( imd ), energy_all( imd )
     do i = 1, natom
       write(15,fmt='(6F16.9,I5)') &
       force0_all( counter_natom + 3*(i-1) + 1 ), &
       force0_all( counter_natom + 3*(i-1) + 2 ), &
       force0_all( counter_natom + 3*(i-1) + 3 ), &
       force_all(  counter_natom + 3*(i-1) + 1 ), &
       force_all(  counter_natom + 3*(i-1) + 2 ), &
       force_all(  counter_natom + 3*(i-1) + 3 ), &
       imd
     enddo

   !-----------
   ! TEST DATA
   !-----------
   elseif( list( imd ) == 0 )then

     write(14,*) energy0_all( imd ), energy_all( imd )
     do i = 1, natom
       write(16,fmt='(6F16.9,I5)') &
       force0_all( counter_natom + 3*(i-1) + 1 ), &
       force0_all( counter_natom + 3*(i-1) + 2 ), &
       force0_all( counter_natom + 3*(i-1) + 3 ), &
       force_all(  counter_natom + 3*(i-1) + 1 ), &
       force_all(  counter_natom + 3*(i-1) + 2 ), &
       force_all(  counter_natom + 3*(i-1) + 3 ), &
       imd
     enddo

   else
     write(*,*)" Error(main): Finalize"
     stop
   endif ! list

   counter_natom = counter_natom + 3*natom

 enddo ! do imd = 1, num_data_total
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 ) write(13,*) "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
 if( myrank == 0 ) write(14,*) "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
 if( myrank == 0 ) write(15,*) "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"
 if( myrank == 0 ) write(16,*) "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"
 close(13) ! energy_training.dat
 close(14) ! energy_test.dat
 close(15) ! force_training.dat
 close(16) ! force_test.dat


 !***************************
 ! END FINALIZE TRANING DATA
 !***************************

 if( myrank == 0 )then
   write(*,*) " Finish calculations."
   call cpu_time( t1 )
   write(*,*) " Total time: ", t1-t0, "seconds"
 endif ! myrank == 0

 call mpi_finalize( ierr )





end
