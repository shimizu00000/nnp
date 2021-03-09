subroutine energy_force_nnp2ab( &
           g2_x_x,  g2_deriv_x_x,  num_g2_x_x, &
           g2_x_y,  g2_deriv_x_y,  num_g2_x_y, &
           g5_x_xx, g5_deriv_x_xx, num_g5_x_xx, &
           g5_x_xy, g5_deriv_x_xy, num_g5_x_xy, &
           g5_x_yy, g5_deriv_x_yy, num_g5_x_yy, &
           num_g_x, &
           network_x, nlayer_x, weight_x, nweight_x, hidden_x, nhidden_x, &
           natom, natom_x, energy_x, force )
implicit none
integer i, j, k, m, di
integer n1
integer l0, l1, l2, l3, l4, l5
integer w1, w2, w3, w4
integer h1, h2, h3, h4
integer h2_0, h3_0, h4_0
integer w2_0, w3_0, w4_0

integer natom_x, natom
integer num_g2_x_x,  num_g2_x_y, &
        num_g5_x_xx, num_g5_x_xy, &
        num_g5_x_yy
integer num_g_x(5)
double precision &
  g2_x_x(  natom_x, num_g2_x_x ), &
  g2_x_y(  natom_x, num_g2_x_y ), &
  g5_x_xx( natom_x, num_g5_x_xx ), &
  g5_x_xy( natom_x, num_g5_x_xy ), &
  g5_x_yy( natom_x, num_g5_x_yy )
double precision &
  g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x ), &
  g2_deriv_x_y(  natom_x, natom, 3, num_g2_x_y ), &
  g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ), &
  g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ), &
  g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy )

integer nlayer_x, nweight_x, nhidden_x
integer network_x( nlayer_x )
double precision input_x( natom_x, 0 : network_x(1) ) 
double precision weight_x( nweight_x ), hidden_x( natom_x, nhidden_x )
double precision dhidden_x( natom_x, nhidden_x ), const( natom_x )


double precision &
  dGdx( natom_x, natom, network_x(1) ), &
  dGdy( natom_x, natom, network_x(1) ), &
  dGdz( natom_x, natom, network_x(1) )

double precision energy_x, atomic_ene_x( natom_x )
double precision force( natom, 3 ), dummy, dummy_x, dummy_y, dummy_z


 input_x = 0.0d0
 do i = 1, natom_x
   input_x(i,0) = 1.0d0
 enddo

 ! x_x 1
 do j = 1, num_g2_x_x
   do i = 1, natom_x
     input_x(i,j) = g2_x_x(i,j)
   enddo
 enddo
 ! x_y 2
 do j = 1, num_g2_x_y
   k = j + num_g_x(1)
   do i = 1, natom_x
     input_x(i,k) = g2_x_y(i,j)
   enddo
 enddo
 ! x_xx 3
 do j = 1, num_g5_x_xx
   k = j + sum( num_g_x(1:2) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_xx(i,j)
   enddo
 enddo 
 ! x_xy 4
 do j = 1, num_g5_x_xy
   k = j + sum( num_g_x(1:3) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_xy(i,j)
   enddo
 enddo
 ! x_yy 5
 do j = 1, num_g5_x_yy
   k = j + sum( num_g_x(1:4) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_yy(i,j)
   enddo
 enddo


 hidden_x     = 0.0d0


!----------------
! 1 hidden layer
!----------------
!-----------------
! 2 hidden layers
!-----------------
 if( nlayer_x == 4 )then

   !------------------
   ! 1st hidden layer
   !------------------
   do i = 1, natom_x
     hidden_x(i,1) = 1.0d0
   enddo
   do l2 = 1, network_x(2)
     h1 = 1 + l2
     do l1 = 0, network_x(1)
       w1 = ( l2 - 1 )*( network_x(1) + 1 ) + ( 1 + l1 )
       dummy = weight_x(w1)
       do i = 1, natom_x
         hidden_x(i,h1) = hidden_x(i,h1) + ( dummy * input_x(i,l1) )
       enddo
     enddo
     do i = 1, natom_x
       hidden_x(i,h1) = dtanh( hidden_x(i,h1) )
     enddo
   enddo
   !------------------
   ! 2nd hidden layer
   !------------------
   w2_0 = ( network_x(1) + 1 ) * network_x(2)
   h2_0 = network_x(2) + 1
   do i = 1, natom_x
     hidden_x(i,h2_0+1) = 1.0d0 ! bias
   enddo
   do l2 = 1, network_x(3)
     h2 = h2_0 + 1 + l2 ! "+1" bias
     do l1 = 0, network_x(2)
       w2 = w2_0 + ( l2 - 1 )*( network_x(2) + 1 ) + ( 1 + l1 )
       h1 = 1 + l1
       dummy = weight_x(w2)
       do i = 1, natom_x
         hidden_x(i,h2) = hidden_x(i,h2) + ( dummy * hidden_x(i,h1) )
       enddo
     enddo
     do i = 1, natom_x
       hidden_x(i,h2) = dtanh( hidden_x(i,h2) )
     enddo
   enddo
   !----------------------------
   ! Atomic energy (last) layer
   !----------------------------
   atomic_ene_x = 0.0d0
   w3_0 = w2_0 + ( network_x(2) + 1 )*network_x(3)
   do l1 = 0, network_x(3)
     w3 = w3_0 + ( l1 + 1 )
     h2 = h2_0 + 1 + l1
     dummy = weight_x(w3)
     do i = 1, natom_x
       atomic_ene_x(i) = atomic_ene_x(i) + ( dummy * hidden_x(i,h2) )
     enddo
   enddo

 endif ! 2 hidden layers

!-----------------
! 3 hidden layers
!-----------------


!---------------------------
!Summation for total energy
!---------------------------
 energy_x = 0.0d0
 do i = 1, natom_x
   energy_x = energy_x + atomic_ene_x(i)
 enddo






!************************************************************
! FORCE CALCULATION
!************************************************************

 dGdx = 0.0d0
 dGdy = 0.0d0
 dGdz = 0.0d0

 ! x_x 1
 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,k) = g2_deriv_x_x(i,j,1,k)
       dGdy(i,j,k) = g2_deriv_x_x(i,j,2,k)
       dGdz(i,j,k) = g2_deriv_x_x(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_y 2
 do k = 1, num_g2_x_y
   m = k + num_g_x(1)
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g2_deriv_x_y(i,j,1,k)
       dGdy(i,j,m) = g2_deriv_x_y(i,j,2,k)
       dGdz(i,j,m) = g2_deriv_x_y(i,j,3,k)
     enddo
   enddo
 enddo
 ! x_xx 3
 do k = 1, num_g5_x_xx
   m = k + sum( num_g_x(1:2) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xx(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xx(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xx(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_xy 4
 do k = 1, num_g5_x_xy
   m = k + sum( num_g_x(1:3) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_yy 5
 do k = 1, num_g5_x_yy
   m = k + sum( num_g_x(1:4) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_yy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_yy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_yy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k




! ( 1.0d0 - hidden_x**2 )
 do j = 1, nhidden_x
   do i = 1, natom_x
     dhidden_x(i,j) = 1.0d0 - hidden_x(i,j)**2
   enddo
 enddo


!****************
! 1 hidden layer
!****************
!*****************
! 2 hidden layers
!*****************
 if( nlayer_x == 4 )then

   w2_0 = ( 1 + network_x(1) ) * network_x(2)
   w3_0 = w2_0 + ( 1 + network_x(2) ) * network_x(3)
   h2_0 = 1 + network_x(2)

   do l1 = 1, network_x(1)
   do l2 = 1, network_x(2)
     w1 = ( 1 + network_x(1) )*( l2 - 1 ) + ( 1 + l1 )
     h1 = 1 + l2
   do l3 = 1, network_x(3)
     w2 = w2_0 + ( 1 + network_x(2) )*( l3 - 1 ) + ( 1 + l2 )
     h2 = h2_0 + ( 1 + l3 )
   do l4 = 1, network_x(4)
     w3 = w3_0 + ( 1 + network_x(3) )*( l4 - 1 ) + ( 1 + l3 )

     dummy = weight_x(w1) * weight_x(w2) * weight_x(w3)
     do n1 = 1, natom_x
       const(n1) = dummy * dhidden_x(n1,h1) * dhidden_x(n1,h2)
     enddo

     do i = 1, natom !_x

       dummy_x = 0.0d0
       dummy_y = 0.0d0
       dummy_z = 0.0d0

       do n1 = 1, natom_x
         !F_x
         dummy_x = dummy_x + ( dGdx(n1,i,l1) * const(n1) )
         !force(i,1) =   force(i,1) &
         !             - dGdx(n1,i,l1) &
         !               * const(n1)
                        !*weight_x(w1) &
                        !*( 1.0d0 - hidden_x(n1,h1)**2 ) &
                        !*weight_x(w2) &
                        !*( 1.0d0 - hidden_x(n1,h2)**2 ) &
                        !*weight_x(w3)
         !F_y
         dummy_y = dummy_y + ( dGdy(n1,i,l1) * const(n1) )
         !force(i,2) =   force(i,2) &
         !             - dGdy(n1,i,l1) &
         !               * const(n1)
                        !*weight_x(w1) &
                        !*( 1.0d0 - hidden_x(n1,h1)**2 ) &
                        !*weight_x(w2) &
                        !*( 1.0d0 - hidden_x(n1,h2)**2 ) &
                        !*weight_x(w3)
         !F_z
         dummy_z = dummy_z + ( dGdz(n1,i,l1) * const(n1) )
         !force(i,3) =   force(i,3) &
         !             - dGdz(n1,i,l1) &
         !               * const(n1)
                        !*weight_x(w1) &
                        !*( 1.0d0 - hidden_x(n1,h1)**2 ) &
                        !*weight_x(w2) &
                        !*( 1.0d0 - hidden_x(n1,h2)**2 ) &
                        !*weight_x(w3)
       enddo !n1

       force(i,1) = force(i,1) - dummy_x
       force(i,2) = force(i,2) - dummy_y
       force(i,3) = force(i,3) - dummy_z

     enddo !i

   enddo !l4
   enddo !l3
   enddo !l2
   enddo !l1

 endif ! 2 hidden layers

!*****************
! 3 hidden layers
!*****************


end subroutine energy_force_nnp2ab
