subroutine read_layer_node2( nlayer_a, nlayer_b )

implicit none
character(len=120) tag
character(len=120) cread
integer            iread
!double precision   dread

integer i, stat, counter
character(len=4) elem_a, elem_b
integer nlayer_a, nlayer_b

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
 nlayer_a = 0
!B
 elem_b = "none"
 nlayer_b = 0
 
 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   !parameters for A
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "nlayer_a" )then
     read(cread,*) iread
     nlayer_a    = iread
   endif

   !parameters for B
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "nlayer_b" )then
     read(cread,*) iread
     nlayer_b    = iread
   endif


 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


end subroutine read_layer_node2
