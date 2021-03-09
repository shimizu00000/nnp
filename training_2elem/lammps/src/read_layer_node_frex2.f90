subroutine read_layer_node_frex2( &
           network_a, network_b, &  
           nlayer_a, nlayer_b, &
           io_nnp_a, io_nnp_b )

implicit none

character(len=120) tag
character(len=120) cread
integer            iread
!double precision   dread

integer i, stat, counter
character(len=4) elem_a, elem_b

integer nlayer_a, nlayer_b
integer network_a( nlayer_a ), &
        network_b( nlayer_b )
!A
integer node_in_a, node_out_a
integer node_h1_a, node_h2_a, node_h3_a, node_h4_a, node_h5_a, &
        node_h6_a, node_h7_a, node_h8_a, node_h9_a, node_h10_a
!B
integer node_in_b, node_out_b, &
        node_h1_b, node_h2_b, node_h3_b, node_h4_b, node_h5_b, &
        node_h6_b, node_h7_b, node_h8_b, node_h9_b, node_h10_b

integer io_nnp_a, io_nnp_b

character dummyc


!***********************
!Read input_nnp.dat file
!***********************

 open(unit=98,file="input_nnp.dat",action="read")
!Number of lines in input_nnp.dat
 counter = 0
 rewind(98)
 do
   read( 98, '(A1)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
   if( trim( dummyc ) /= '' ) counter = counter + 1
 enddo

!A
 elem_a = "none"
 io_nnp_a = 0
 network_a = 0
 node_in_a = 0 ; node_out_a = 0
 node_h1_a = 0 ; node_h2_a = 0 ; node_h3_a = 0 ; node_h4_a = 0 ; node_h5_a = 0
 node_h6_a = 0 ; node_h7_a = 0 ; node_h8_a = 0 ; node_h9_a = 0 ; node_h10_a = 0
!B
 elem_b = "none"
 io_nnp_b = 0
 network_b = 0
 node_in_b = 0 ; node_out_b = 0 
 node_h1_b = 0 ; node_h2_b = 0 ; node_h3_b = 0 ; node_h4_b = 0 ; node_h5_b = 0 
 node_h6_b = 0 ; node_h7_b = 0 ; node_h8_b = 0 ; node_h9_b = 0 ; node_h10_b = 0


 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   ! Parameters for A
   !------------------
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "nlayer_a" )then
     read(cread,*) iread
     nlayer_a    = iread
   endif
   if( tag == "nnp_io_a" )then
     read(cread,*) iread
     io_nnp_a    = iread
   endif

   if( tag == "node_in_a" )then
     read(cread,*)  iread
     network_a(1) = iread
   endif
   if( tag == "node_out_a" )then
     read(cread,*) iread
     network_a(nlayer_a) = iread
   endif

   if( tag == "node_h1_a" )then
     read(cread,*)  iread
     network_a(2) = iread
   endif
   if( tag == "node_h2_a" )then
     read(cread,*)  iread
     network_a(3) = iread
   endif
   if( tag == "node_h3_a" )then
     read(cread,*)  iread
     network_a(4) = iread
   endif
   if( tag == "node_h4_a" )then
     read(cread,*)  iread
     network_a(5) = iread
   endif
   if( tag == "node_h5_a" )then
     read(cread,*)  iread
     network_a(6) = iread
   endif
   if( tag == "node_h6_a" )then
     read(cread,*)  iread
     network_a(7) = iread
   endif
   if( tag == "node_h7_a" )then
     read(cread,*)  iread
     network_a(8) = iread
   endif
   if( tag == "node_h8_a" )then
     read(cread,*)  iread
     network_a(9) = iread
   endif
   if( tag == "node_h9_a" )then
     read(cread,*)  iread
     network_a(10) = iread
   endif
   if( tag == "node_h10_a" )then
     read(cread,*)  iread
     network_a(11) = iread
   endif

   ! Parameters for B
   !------------------
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "nlayer_b" )then
     read(cread,*) iread
     nlayer_b    = iread
   endif
   if( tag == "nnp_io_b" )then
     read(cread,*) iread
     io_nnp_b    = iread
   endif

   if( tag == "node_in_b" )then
     read(cread,*)  iread
     network_b(1) = iread
   endif
   if( tag == "node_out_b" )then
     read(cread,*)  iread
     network_b(nlayer_b) = iread
   endif

   if( tag == "node_h1_b" )then
     read(cread,*)  iread
     network_b(2) = iread
   endif
   if( tag == "node_h2_b" )then
     read(cread,*)  iread
     network_b(3) = iread
   endif
   if( tag == "node_h3_b" )then
     read(cread,*)  iread
     network_b(4) = iread
   endif
   if( tag == "node_h4_b" )then
     read(cread,*)  iread
     network_b(5) = iread
   endif
   if( tag == "node_h5_b" )then
     read(cread,*)  iread
     network_b(6) = iread
   endif
   if( tag == "node_h6_b" )then
     read(cread,*)  iread
     network_b(7) = iread
   endif
   if( tag == "node_h7_b" )then
     read(cread,*)  iread
     network_b(8) = iread
   endif
   if( tag == "node_h8_b" )then
     read(cread,*)  iread
     network_b(9) = iread
   endif
   if( tag == "node_h9_b" )then
     read(cread,*)  iread
     network_b(10) = iread
   endif
   if( tag == "node_h10_b" )then
     read(cread,*)  iread
     network_b(11) = iread
   endif


 enddo ! i = 1, counter

 close(98) ! input_nnp.dat


 write(*,*) " "
 write(*,*) " Parameters for neural network"
 write(*,*) " ===================================="

 !A
 if( elem_a /= "none" )then
   write(*,fmt='(A10)',advance='no') " NNP(a): "
   do i = 1, nlayer_a - 1
     write(*,fmt='(I4,A2)',advance='no') network_a(i),"-"
   enddo
   write(*,fmt='(I4)',advance='no') network_a(nlayer_a)
   if( io_nnp_a == 0 ) write(*,fmt='(A4)') "OFF"
   if( io_nnp_a == 1 ) write(*,fmt='(A4)') "ON"
 endif
 !B
 if( elem_b /= "none" )then
   write(*,fmt='(A10)',advance='no') " NNP(b): "
   do i = 1, nlayer_b - 1
     write(*,fmt='(I4,A2)',advance='no') network_b(i),"-"
   enddo
   write(*,fmt='(I4)',advance='no') network_b(nlayer_b)
   if( io_nnp_b == 0 ) write(*,fmt='(A4)') "OFF"
   if( io_nnp_b == 1 ) write(*,fmt='(A4)') "ON"
 endif

 write(*,*) " ===================================="
 write(*,*) " "


end subroutine read_layer_node_frex2
