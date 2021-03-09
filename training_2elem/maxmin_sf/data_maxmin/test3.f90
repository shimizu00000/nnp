program test
implicit none
integer i
double precision g2(8)

open(unit=1,file="g2_a_a_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_a_a"
print*,g2
open(unit=1,file="g2_a_b_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_a_b"
print*,g2
open(unit=1,file="g2_a_c_maxmin.dat",action="read",form="unformatted")
read(1)g2
print*,"g2_a_c"
print*,g2



end

