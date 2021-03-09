program test
implicit none
integer i
double precision g2(6,2),g5(18,2)

!A
open(unit=1,file="g2_a_a_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_a_a"
print*,g2
open(unit=1,file="g2_a_b_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_a_b"
print*,g2

open(unit=1,file="g5_a_aa_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_a_aa"
print*,g5
open(unit=1,file="g5_a_ab_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_a_ab"
print*,g5
open(unit=1,file="g5_a_bb_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_a_bb"
print*,g5


!B
open(unit=1,file="g2_b_a_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_b_a"
print*,g2
open(unit=1,file="g2_b_b_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_b_b"
print*,g2

open(unit=1,file="g5_b_aa_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_b_aa"
print*,g5
open(unit=1,file="g5_b_ab_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_b_ab"
print*,g5
open(unit=1,file="g5_b_bb_maxmin.dat",action="read",form="unformatted")
read(1)g5
print*,"g5_b_bb"
print*,g5




end


