module allocarray
contains


!************************************************************
! ALLOCATE_2
!************************************************************
subroutine alloc_2( num_atom, num_data_total, &
                    num_atom_x, num_atom_y, &
                    num_data_x, num_data_y, &
                    natom_x_total, natom_xy_total, &
                    natom_y_total, natom_yx_total, &
                    natom_x_sum,   natom_xy_sum, &
                    natom_y_sum,   natom_yx_sum, &
                    tag2_x_x,  g2_x_x,  g2_deriv_x_x,  g2_x_x_maxmin,  num_g2_x_x, &
                    tag2_x_y,  g2_x_y,  g2_deriv_x_y,  g2_x_y_maxmin,  num_g2_x_y, &
                    tag5_x_xx, g5_x_xx, g5_deriv_x_xx, g5_x_xx_maxmin, num_g5_x_xx, &
                    tag5_x_xy, g5_x_xy, g5_deriv_x_xy, g5_x_xy_maxmin, num_g5_x_xy, &
                    tag5_x_yy, g5_x_yy, g5_deriv_x_yy, g5_x_yy_maxmin, num_g5_x_yy, &
                    tag2_y_x,  g2_y_x,  g2_deriv_y_x,  g2_y_x_maxmin,  num_g2_y_x, &
                    tag2_y_y,  g2_y_y,  g2_deriv_y_y,  g2_y_y_maxmin,  num_g2_y_y, &
                    tag5_y_xx, g5_y_xx, g5_deriv_y_xx, g5_y_xx_maxmin, num_g5_y_xx, &
                    tag5_y_xy, g5_y_xy, g5_deriv_y_xy, g5_y_xy_maxmin, num_g5_y_xy, &
                    tag5_y_yy, g5_y_yy, g5_deriv_y_yy, g5_y_yy_maxmin, num_g5_y_yy )
                        
implicit none
character(len=120) filename
character(len=6) tag2_x_x,  tag2_x_y
character(len=7) tag5_x_xx, tag5_x_xy, &
                 tag5_x_yy
character(len=6) tag2_y_x,  tag2_y_y
character(len=7) tag5_y_xx, tag5_y_xy, &
                 tag5_y_yy
integer num_data_total, num_data_x, num_data_y
integer natom_x_total, natom_xy_total
integer natom_y_total, natom_yx_total
integer natom_x_sum,   natom_xy_sum
integer natom_y_sum,   natom_yx_sum
integer,allocatable,dimension(:) :: &
  num_atom, num_atom_x, num_atom_y
integer num_g2_x_x,  num_g2_x_y, &
        num_g5_x_xx, num_g5_x_xy, &
        num_g5_x_yy
integer num_g2_y_x,  num_g2_y_y, &
        num_g5_y_xx, num_g5_y_xy, &
        num_g5_y_yy
double precision,allocatable,dimension(:,:) :: &
  g2_x_x,  g2_x_y, &
  g5_x_xx, g5_x_xy, &
  g5_x_yy
double precision,allocatable,dimension(:,:) :: &
  g2_y_x,  g2_y_y, &
  g5_y_xx, g5_y_xy, &
  g5_y_yy
double precision,allocatable,dimension(:,:) :: &
  g2_x_x_maxmin,  g2_x_y_maxmin, &
  g5_x_xx_maxmin, g5_x_xy_maxmin, &
  g5_x_yy_maxmin
double precision,allocatable,dimension(:,:) :: &
  g2_y_x_maxmin,  g2_y_y_maxmin, &
  g5_y_xx_maxmin, g5_y_xy_maxmin, &
  g5_y_yy_maxmin
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_x_x,  g2_deriv_x_y, &
  g5_deriv_x_xx, g5_deriv_x_xy, &
  g5_deriv_x_yy
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_y_x,  g2_deriv_y_y, &
  g5_deriv_y_xx, g5_deriv_y_xy, &
  g5_deriv_y_yy


 allocate( num_atom( num_data_total ) )

!X
 allocate( num_atom_x( num_data_x ) )

 allocate( g2_x_x(  natom_x_total,   num_g2_x_x ) )
 allocate( g2_x_y(  natom_xy_total,  num_g2_x_y ) )
 allocate( g5_x_xx( natom_x_total,   num_g5_x_xx ) )
 allocate( g5_x_xy( natom_xy_total,  num_g5_x_xy ) )
 allocate( g5_x_yy( natom_xy_total,  num_g5_x_yy ) )

 allocate( g2_x_x_maxmin(  num_g2_x_x, 2 ) )
 allocate( g2_x_y_maxmin(  num_g2_x_y, 2 ) )
 allocate( g5_x_xx_maxmin( num_g5_x_xx, 2 ) )
 allocate( g5_x_xy_maxmin( num_g5_x_xy, 2 ) )
 allocate( g5_x_yy_maxmin( num_g5_x_yy, 2 ) )

 allocate( g2_deriv_x_x(  natom_x_sum,   3, num_g2_x_x ) )
 allocate( g2_deriv_x_y(  natom_xy_sum,  3, num_g2_x_y ) )
 allocate( g5_deriv_x_xx( natom_x_sum,   3, num_g5_x_xx ) )
 allocate( g5_deriv_x_xy( natom_xy_sum,  3, num_g5_x_xy ) )
 allocate( g5_deriv_x_yy( natom_xy_sum,  3, num_g5_x_yy ) )

 filename = './data_maxmin/'//trim(adjustl(tag2_x_x))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_x_x_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_x_y))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_x_y_maxmin
 close(20)

 filename = './data_maxmin/'//trim(adjustl(tag5_x_xx))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_xx_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_xy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_xy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_yy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_yy_maxmin
 close(20)

!Y
 allocate( num_atom_y( num_data_y ) )

 allocate( g2_y_x(  natom_yx_total,  num_g2_y_x ) )
 allocate( g2_y_y(  natom_y_total,   num_g2_y_y ) )
 allocate( g5_y_xx( natom_yx_total,  num_g5_y_xx ) )
 allocate( g5_y_xy( natom_yx_total,  num_g5_y_xy ) )
 allocate( g5_y_yy( natom_y_total,   num_g5_y_yy ) )

 allocate( g2_y_x_maxmin(  num_g2_y_x, 2 ) )
 allocate( g2_y_y_maxmin(  num_g2_y_y, 2 ) )
 allocate( g5_y_xx_maxmin( num_g5_y_xx, 2 ) )
 allocate( g5_y_xy_maxmin( num_g5_y_xy, 2 ) )
 allocate( g5_y_yy_maxmin( num_g5_y_yy, 2 ) )

 allocate( g2_deriv_y_x(  natom_yx_sum,  3, num_g2_y_x ) )
 allocate( g2_deriv_y_y(  natom_y_sum,   3, num_g2_y_y ) )
 allocate( g5_deriv_y_xx( natom_yx_sum,  3, num_g5_y_xx ) )
 allocate( g5_deriv_y_xy( natom_yx_sum,  3, num_g5_y_xy ) )
 allocate( g5_deriv_y_yy( natom_y_sum,   3, num_g5_y_yy ) )

 filename = './data_maxmin/'//trim(adjustl(tag2_y_x))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_y_x_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_y_y))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_y_y_maxmin
 close(20)

 filename = './data_maxmin/'//trim(adjustl(tag5_y_xx))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_xx_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_xy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_xy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_yy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_yy_maxmin
 close(20)



end subroutine alloc_2


end module
