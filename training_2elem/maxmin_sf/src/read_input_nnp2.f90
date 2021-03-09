subroutine read_input_nnp2( &
           elem_a, &
           num_g2_a_a,  num_g2_a_b, &
           num_g5_a_aa, num_g5_a_ab, &
           num_g5_a_bb, &
           elem_b, &
           num_g2_b_a,  num_g2_b_b, &
           num_g5_b_aa, num_g5_b_ab, &
           num_g5_b_bb, &
           rc, rn, phi, neighbor, &
           dir_data, dir_sf, io_mix )

implicit none

character(len=100) tag
character(len=100) cread
integer            iread
double precision   dread

integer i, stat, counter

integer node_in_a, node_in_b

!A-*
integer num_g2_a_a,  num_g2_a_b, &
        num_g5_a_aa, num_g5_a_ab, &
        num_g5_a_bb
!B-*
integer num_g2_b_a,  num_g2_b_b, &
        num_g5_b_aa, num_g5_b_ab, &
        num_g5_b_bb

integer   nelem, natom

integer   natom_a, natom_b
integer   natom_1, natom_2

character(len=4) elem_a, elem_b
character(len=4) elem_1, elem_2

integer md_steps

integer neighbor

double precision rc, rn

character(len=10) phi

character(len=120) dir_data, dir_sf

integer io_mix

character dummyc


!***********************
!Read input_nnp.dat file
!***********************
!input_nnp.dat includes
!----------------------
!rc = 7.0d0
!rn = 1.0d0
!neighbor = 300
!phi = tanh
!
!elemen_a = Au
!node_in_a = 44
!num_g2_a_a = 8
!num_g2_a_b = 8
!num_g2_a_c = 8
!num_g5_a_aa = 22,
!num_g5_a_ab = 22
!num_g5_a_ac = 22
!num_g5_a_bb = 22
!num_g5_a_bc = 22,
!num_g5_a_cc = 22

!also b,c,d
!---
 open(unit=98,file="input_nnp.dat",action="read")
!Number of lines in input_nnp.dat
 counter = 0
 rewind(98)
 do
   read( 98, '(A1)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
   if( trim( dummyc ) /= '' ) counter = counter + 1
 enddo

 rc = 0.0d0
 rn = 0.0d0
 neighbor = 0

!A
 elem_a = "none"
 node_in_a   = 0 
 num_g2_a_a  = 0 ; num_g2_a_b  = 0
 num_g5_a_aa = 0 ; num_g5_a_ab = 0
 num_g5_a_bb = 0
!B
 elem_b = "none"
 node_in_b   = 0 
 num_g2_b_a  = 0 ; num_g2_b_b  = 0
 num_g5_b_aa = 0 ; num_g5_b_ab = 0
 num_g5_b_bb = 0

 io_mix = 0
 

print*,"COUNTER",counter

 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   !parameters for A
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "node_in_a" )then
     read(cread,*) iread
     node_in_a   = iread
   endif
   if( tag == "num_g2_a_a" )then
     read(cread,*) iread
     num_g2_a_a  = iread
   endif
   if( tag == "num_g2_a_b" )then
     read(cread,*) iread
     num_g2_a_b  = iread
   endif
   if( tag == "num_g5_a_aa" )then
     read(cread,*) iread
     num_g5_a_aa = iread
   endif
   if( tag == "num_g5_a_ab" )then
     read(cread,*) iread
     num_g5_a_ab = iread
   endif
   if( tag == "num_g5_a_bb" )then
     read(cread,*) iread
     num_g5_a_bb = iread
   endif

   !parameters for B
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "node_in_b" )then
     read(cread,*) iread
     node_in_b   = iread
   endif
   if( tag == "num_g2_b_a" )then
     read(cread,*) iread
     num_g2_b_a  = iread
   endif
   if( tag == "num_g2_b_b" )then
     read(cread,*) iread
     num_g2_b_b  = iread
   endif
   if( tag == "num_g5_b_aa" )then
     read(cread,*) iread
     num_g5_b_aa = iread
   endif
   if( tag == "num_g5_b_ab" )then
     read(cread,*) iread
     num_g5_b_ab = iread
   endif
   if( tag == "num_g5_b_bb" )then
     read(cread,*) iread
     num_g5_b_bb = iread
   endif


   !other parameters
   if( tag == "rc" )then
     read(cread,*) dread
     rc = dread
   endif
   if( tag == "rn" )then
     read(cread,*) dread
     rn = dread
   endif
   if( tag == "neighbor" )then
     read(cread,*) iread
     neighbor = iread
   endif
   if( tag == "phi" )then
     phi = cread
   endif
   if( tag == "dir_data" )then
     dir_data = cread
   endif
   if( tag == "dir_sf" )then
     dir_sf = cread
   endif

   if( tag == "io_mix" )then
     read(cread,*) iread
     io_mix = iread
   endif


 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


 write(*,*) " "
 write(*,*) " Parameters for symmetry function"
 write(*,*) " ============================================================"

 !1 element
 !A
 if( elem_a /= "none" .and. elem_b == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa
   write(*,*) " "

 !B
 elseif( elem_a == "none" .and. elem_b /= "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb
   write(*,*) " "

 endif

 !2 elements
 !A,B
 if( elem_a /= "none" .and. elem_b /= "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-b: ", num_g2_a_b
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ab: ", num_g5_a_ab, "a-bb: ", num_g5_a_bb
   write(*,*) " "
   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-a: ", num_g2_b_a, "b-b: ", num_g2_b_b
   write(*,*) " num_g5"
   write(*,*) " b-aa: ", num_g5_b_aa, "b-ab: ", num_g5_b_ab, "b-bb: ", num_g5_b_bb
   write(*,*) " "

 endif


 write(*,*) " "
 write(*,*) " Activation function: ", phi
 write(*,*) " Cutoff distance:", rc, " [Angst]"
 write(*,*) " Margin distance:", rn, " [Angst]"
 write(*,*) " Neighbor list size:", neighbor, " atoms"
 write(*,*) " "
 write(*,*) " Data: ", trim(adjustl(dir_data))
 write(*,*) " SF: ", trim(adjustl(dir_sf))
 write(*,*) " "
 write(*,*) " ============================================================"


 112 format( A9, A9 )
 113 format( A11, I5 )
 114 format( A9, A7, I5, A7, I5, A7, I5 )
 115 format( A9, A7, I5, A7, I5, A7, I5, A7, I5, A7, I5, A7, I5 )
 116 format( A17, F5.2, A9, A17, F5.2, A8 )
 117 format( A20, I6, A6 )


end subroutine read_input_nnp2
