c
c
c       this file was generated automatically
c       it contains subroutines which load
c       symmetry info for volume code tables
c
c




       subroutine loadsymsc(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(3,*), iflip(3,*)
 
      iref( 1)  =  2
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      idimp(3, 1) =  3
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iflip(3, 1) =  1
      iref( 2)  =  3
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      idimp(3, 2) =  3
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iflip(3, 2) =  1
      iref( 3)  =  2
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      idimp(3, 3) =  3
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iflip(3, 3) =  1
      iref( 4)  =  3
      idimp(1, 4) =  2
      idimp(2, 4) =  1
      idimp(3, 4) =  3
      iflip(1, 4) =  1
      iflip(2, 4) =  1
      iflip(3, 4) =  1
      iref( 5)  =  4
      idimp(1, 5) =  1
      idimp(2, 5) =  2
      idimp(3, 5) =  3
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iflip(3, 5) =  1
      iref( 6)  =  3
      idimp(1, 6) =  2
      idimp(2, 6) =  1
      idimp(3, 6) =  3
      iflip(1, 6) = -1
      iflip(2, 6) =  1
      iflip(3, 6) =  1
      iref( 7)  =  2
      idimp(1, 7) =  1
      idimp(2, 7) =  2
      idimp(3, 7) =  3
      iflip(1, 7) =  1
      iflip(2, 7) = -1
      iflip(3, 7) =  1
      iref( 8)  =  3
      idimp(1, 8) =  1
      idimp(2, 8) =  2
      idimp(3, 8) =  3
      iflip(1, 8) =  1
      iflip(2, 8) = -1
      iflip(3, 8) =  1
      iref( 9)  =  2
      idimp(1, 9) =  1
      idimp(2, 9) =  2
      idimp(3, 9) =  3
      iflip(1, 9) = -1
      iflip(2, 9) = -1
      iflip(3, 9) =  1
      iref(10)  =  3
      idimp(1,10) =  3
      idimp(2,10) =  1
      idimp(3,10) =  2
      iflip(1,10) =  1
      iflip(2,10) =  1
      iflip(3,10) =  1
      iref(11)  =  4
      idimp(1,11) =  1
      idimp(2,11) =  3
      idimp(3,11) =  2
      iflip(1,11) =  1
      iflip(2,11) =  1
      iflip(3,11) =  1
      iref(12)  =  3
      idimp(1,12) =  3
      idimp(2,12) =  1
      idimp(3,12) =  2
      iflip(1,12) = -1
      iflip(2,12) =  1
      iflip(3,12) =  1
      iref(13)  =  4
      idimp(1,13) =  2
      idimp(2,13) =  3
      idimp(3,13) =  1
      iflip(1,13) =  1
      iflip(2,13) =  1
      iflip(3,13) =  1
      iref(14)  =  1
      idimp(1,14) =  1
      idimp(2,14) =  2
      idimp(3,14) =  3
      iflip(1,14) =  1
      iflip(2,14) =  1
      iflip(3,14) =  1
      iref(15)  =  4
      idimp(1,15) =  2
      idimp(2,15) =  3
      idimp(3,15) =  1
      iflip(1,15) = -1
      iflip(2,15) =  1
      iflip(3,15) =  1
      iref(16)  =  3
      idimp(1,16) =  3
      idimp(2,16) =  1
      idimp(3,16) =  2
      iflip(1,16) =  1
      iflip(2,16) = -1
      iflip(3,16) =  1
      iref(17)  =  4
      idimp(1,17) =  1
      idimp(2,17) =  3
      idimp(3,17) =  2
      iflip(1,17) =  1
      iflip(2,17) = -1
      iflip(3,17) =  1
      iref(18)  =  3
      idimp(1,18) =  3
      idimp(2,18) =  1
      idimp(3,18) =  2
      iflip(1,18) = -1
      iflip(2,18) = -1
      iflip(3,18) =  1
      iref(19)  =  2
      idimp(1,19) =  1
      idimp(2,19) =  2
      idimp(3,19) =  3
      iflip(1,19) =  1
      iflip(2,19) =  1
      iflip(3,19) = -1
      iref(20)  =  3
      idimp(1,20) =  1
      idimp(2,20) =  2
      idimp(3,20) =  3
      iflip(1,20) =  1
      iflip(2,20) =  1
      iflip(3,20) = -1
      iref(21)  =  2
      idimp(1,21) =  1
      idimp(2,21) =  2
      idimp(3,21) =  3
      iflip(1,21) = -1
      iflip(2,21) =  1
      iflip(3,21) = -1
      iref(22)  =  3
      idimp(1,22) =  2
      idimp(2,22) =  1
      idimp(3,22) =  3
      iflip(1,22) =  1
      iflip(2,22) =  1
      iflip(3,22) = -1
      iref(23)  =  4
      idimp(1,23) =  1
      idimp(2,23) =  2
      idimp(3,23) =  3
      iflip(1,23) =  1
      iflip(2,23) =  1
      iflip(3,23) = -1
      iref(24)  =  3
      idimp(1,24) =  2
      idimp(2,24) =  1
      idimp(3,24) =  3
      iflip(1,24) = -1
      iflip(2,24) =  1
      iflip(3,24) = -1
      iref(25)  =  2
      idimp(1,25) =  1
      idimp(2,25) =  2
      idimp(3,25) =  3
      iflip(1,25) =  1
      iflip(2,25) = -1
      iflip(3,25) = -1
      iref(26)  =  3
      idimp(1,26) =  1
      idimp(2,26) =  2
      idimp(3,26) =  3
      iflip(1,26) =  1
      iflip(2,26) = -1
      iflip(3,26) = -1
      iref(27)  =  2
      idimp(1,27) =  1
      idimp(2,27) =  2
      idimp(3,27) =  3
      iflip(1,27) = -1
      iflip(2,27) = -1
      iflip(3,27) = -1
 
       return
       end
 
 
       subroutine loadsymsbtos(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(3,*), iflip(3,*)
 
      iref( 1)  =  1
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      idimp(3, 1) =  3
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iflip(3, 1) =  1
      iref( 2)  =  2
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      idimp(3, 2) =  3
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iflip(3, 2) =  1
      iref( 3)  =  2
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      idimp(3, 3) =  3
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iflip(3, 3) =  1
      iref( 4)  =  1
      idimp(1, 4) =  1
      idimp(2, 4) =  2
      idimp(3, 4) =  3
      iflip(1, 4) = -1
      iflip(2, 4) =  1
      iflip(3, 4) =  1
      iref( 5)  =  2
      idimp(1, 5) =  2
      idimp(2, 5) =  1
      idimp(3, 5) =  3
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iflip(3, 5) =  1
      iref( 6)  =  3
      idimp(1, 6) =  1
      idimp(2, 6) =  2
      idimp(3, 6) =  3
      iflip(1, 6) =  1
      iflip(2, 6) =  1
      iflip(3, 6) =  1
      iref( 7)  =  3
      idimp(1, 7) =  1
      idimp(2, 7) =  2
      idimp(3, 7) =  3
      iflip(1, 7) = -1
      iflip(2, 7) =  1
      iflip(3, 7) =  1
      iref( 8)  =  2
      idimp(1, 8) =  2
      idimp(2, 8) =  1
      idimp(3, 8) =  3
      iflip(1, 8) = -1
      iflip(2, 8) =  1
      iflip(3, 8) =  1
      iref( 9)  =  2
      idimp(1, 9) =  2
      idimp(2, 9) =  1
      idimp(3, 9) =  3
      iflip(1, 9) =  1
      iflip(2, 9) = -1
      iflip(3, 9) =  1
      iref(10)  =  3
      idimp(1,10) =  1
      idimp(2,10) =  2
      idimp(3,10) =  3
      iflip(1,10) =  1
      iflip(2,10) = -1
      iflip(3,10) =  1
      iref(11)  =  3
      idimp(1,11) =  1
      idimp(2,11) =  2
      idimp(3,11) =  3
      iflip(1,11) = -1
      iflip(2,11) = -1
      iflip(3,11) =  1
      iref(12)  =  2
      idimp(1,12) =  2
      idimp(2,12) =  1
      idimp(3,12) =  3
      iflip(1,12) = -1
      iflip(2,12) = -1
      iflip(3,12) =  1
      iref(13)  =  1
      idimp(1,13) =  1
      idimp(2,13) =  2
      idimp(3,13) =  3
      iflip(1,13) =  1
      iflip(2,13) = -1
      iflip(3,13) =  1
      iref(14)  =  2
      idimp(1,14) =  1
      idimp(2,14) =  2
      idimp(3,14) =  3
      iflip(1,14) =  1
      iflip(2,14) = -1
      iflip(3,14) =  1
      iref(15)  =  2
      idimp(1,15) =  1
      idimp(2,15) =  2
      idimp(3,15) =  3
      iflip(1,15) = -1
      iflip(2,15) = -1
      iflip(3,15) =  1
      iref(16)  =  1
      idimp(1,16) =  1
      idimp(2,16) =  2
      idimp(3,16) =  3
      iflip(1,16) = -1
      iflip(2,16) = -1
      iflip(3,16) =  1
      iref(17)  =  2
      idimp(1,17) =  3
      idimp(2,17) =  1
      idimp(3,17) =  2
      iflip(1,17) =  1
      iflip(2,17) =  1
      iflip(3,17) =  1
      iref(18)  =  3
      idimp(1,18) =  1
      idimp(2,18) =  3
      idimp(3,18) =  2
      iflip(1,18) =  1
      iflip(2,18) =  1
      iflip(3,18) =  1
      iref(19)  =  3
      idimp(1,19) =  1
      idimp(2,19) =  3
      idimp(3,19) =  2
      iflip(1,19) = -1
      iflip(2,19) =  1
      iflip(3,19) =  1
      iref(20)  =  2
      idimp(1,20) =  3
      idimp(2,20) =  1
      idimp(3,20) =  2
      iflip(1,20) = -1
      iflip(2,20) =  1
      iflip(3,20) =  1
      iref(21)  =  3
      idimp(1,21) =  2
      idimp(2,21) =  3
      idimp(3,21) =  1
      iflip(1,21) =  1
      iflip(2,21) =  1
      iflip(3,21) =  1
      iref(22)  =  3
      idimp(1,22) =  2
      idimp(2,22) =  3
      idimp(3,22) =  1
      iflip(1,22) = -1
      iflip(2,22) =  1
      iflip(3,22) =  1
      iref(23)  =  3
      idimp(1,23) =  2
      idimp(2,23) =  3
      idimp(3,23) =  1
      iflip(1,23) =  1
      iflip(2,23) = -1
      iflip(3,23) =  1
      iref(24)  =  3
      idimp(1,24) =  2
      idimp(2,24) =  3
      idimp(3,24) =  1
      iflip(1,24) = -1
      iflip(2,24) = -1
      iflip(3,24) =  1
      iref(25)  =  2
      idimp(1,25) =  3
      idimp(2,25) =  1
      idimp(3,25) =  2
      iflip(1,25) =  1
      iflip(2,25) = -1
      iflip(3,25) =  1
      iref(26)  =  3
      idimp(1,26) =  1
      idimp(2,26) =  3
      idimp(3,26) =  2
      iflip(1,26) =  1
      iflip(2,26) = -1
      iflip(3,26) =  1
      iref(27)  =  3
      idimp(1,27) =  1
      idimp(2,27) =  3
      idimp(3,27) =  2
      iflip(1,27) = -1
      iflip(2,27) = -1
      iflip(3,27) =  1
      iref(28)  =  2
      idimp(1,28) =  3
      idimp(2,28) =  1
      idimp(3,28) =  2
      iflip(1,28) = -1
      iflip(2,28) = -1
      iflip(3,28) =  1
      iref(29)  =  2
      idimp(1,29) =  3
      idimp(2,29) =  1
      idimp(3,29) =  2
      iflip(1,29) =  1
      iflip(2,29) =  1
      iflip(3,29) = -1
      iref(30)  =  3
      idimp(1,30) =  1
      idimp(2,30) =  3
      idimp(3,30) =  2
      iflip(1,30) =  1
      iflip(2,30) =  1
      iflip(3,30) = -1
      iref(31)  =  3
      idimp(1,31) =  1
      idimp(2,31) =  3
      idimp(3,31) =  2
      iflip(1,31) = -1
      iflip(2,31) =  1
      iflip(3,31) = -1
      iref(32)  =  2
      idimp(1,32) =  3
      idimp(2,32) =  1
      idimp(3,32) =  2
      iflip(1,32) = -1
      iflip(2,32) =  1
      iflip(3,32) = -1
      iref(33)  =  3
      idimp(1,33) =  2
      idimp(2,33) =  3
      idimp(3,33) =  1
      iflip(1,33) =  1
      iflip(2,33) =  1
      iflip(3,33) = -1
      iref(34)  =  3
      idimp(1,34) =  2
      idimp(2,34) =  3
      idimp(3,34) =  1
      iflip(1,34) = -1
      iflip(2,34) =  1
      iflip(3,34) = -1
      iref(35)  =  3
      idimp(1,35) =  2
      idimp(2,35) =  3
      idimp(3,35) =  1
      iflip(1,35) =  1
      iflip(2,35) = -1
      iflip(3,35) = -1
      iref(36)  =  3
      idimp(1,36) =  2
      idimp(2,36) =  3
      idimp(3,36) =  1
      iflip(1,36) = -1
      iflip(2,36) = -1
      iflip(3,36) = -1
      iref(37)  =  2
      idimp(1,37) =  3
      idimp(2,37) =  1
      idimp(3,37) =  2
      iflip(1,37) =  1
      iflip(2,37) = -1
      iflip(3,37) = -1
      iref(38)  =  3
      idimp(1,38) =  1
      idimp(2,38) =  3
      idimp(3,38) =  2
      iflip(1,38) =  1
      iflip(2,38) = -1
      iflip(3,38) = -1
      iref(39)  =  3
      idimp(1,39) =  1
      idimp(2,39) =  3
      idimp(3,39) =  2
      iflip(1,39) = -1
      iflip(2,39) = -1
      iflip(3,39) = -1
      iref(40)  =  2
      idimp(1,40) =  3
      idimp(2,40) =  1
      idimp(3,40) =  2
      iflip(1,40) = -1
      iflip(2,40) = -1
      iflip(3,40) = -1
      iref(41)  =  1
      idimp(1,41) =  1
      idimp(2,41) =  2
      idimp(3,41) =  3
      iflip(1,41) =  1
      iflip(2,41) =  1
      iflip(3,41) = -1
      iref(42)  =  2
      idimp(1,42) =  1
      idimp(2,42) =  2
      idimp(3,42) =  3
      iflip(1,42) =  1
      iflip(2,42) =  1
      iflip(3,42) = -1
      iref(43)  =  2
      idimp(1,43) =  1
      idimp(2,43) =  2
      idimp(3,43) =  3
      iflip(1,43) = -1
      iflip(2,43) =  1
      iflip(3,43) = -1
      iref(44)  =  1
      idimp(1,44) =  1
      idimp(2,44) =  2
      idimp(3,44) =  3
      iflip(1,44) = -1
      iflip(2,44) =  1
      iflip(3,44) = -1
      iref(45)  =  2
      idimp(1,45) =  2
      idimp(2,45) =  1
      idimp(3,45) =  3
      iflip(1,45) =  1
      iflip(2,45) =  1
      iflip(3,45) = -1
      iref(46)  =  3
      idimp(1,46) =  1
      idimp(2,46) =  2
      idimp(3,46) =  3
      iflip(1,46) =  1
      iflip(2,46) =  1
      iflip(3,46) = -1
      iref(47)  =  3
      idimp(1,47) =  1
      idimp(2,47) =  2
      idimp(3,47) =  3
      iflip(1,47) = -1
      iflip(2,47) =  1
      iflip(3,47) = -1
      iref(48)  =  2
      idimp(1,48) =  2
      idimp(2,48) =  1
      idimp(3,48) =  3
      iflip(1,48) = -1
      iflip(2,48) =  1
      iflip(3,48) = -1
      iref(49)  =  2
      idimp(1,49) =  2
      idimp(2,49) =  1
      idimp(3,49) =  3
      iflip(1,49) =  1
      iflip(2,49) = -1
      iflip(3,49) = -1
      iref(50)  =  3
      idimp(1,50) =  1
      idimp(2,50) =  2
      idimp(3,50) =  3
      iflip(1,50) =  1
      iflip(2,50) = -1
      iflip(3,50) = -1
      iref(51)  =  3
      idimp(1,51) =  1
      idimp(2,51) =  2
      idimp(3,51) =  3
      iflip(1,51) = -1
      iflip(2,51) = -1
      iflip(3,51) = -1
      iref(52)  =  2
      idimp(1,52) =  2
      idimp(2,52) =  1
      idimp(3,52) =  3
      iflip(1,52) = -1
      iflip(2,52) = -1
      iflip(3,52) = -1
      iref(53)  =  1
      idimp(1,53) =  1
      idimp(2,53) =  2
      idimp(3,53) =  3
      iflip(1,53) =  1
      iflip(2,53) = -1
      iflip(3,53) = -1
      iref(54)  =  2
      idimp(1,54) =  1
      idimp(2,54) =  2
      idimp(3,54) =  3
      iflip(1,54) =  1
      iflip(2,54) = -1
      iflip(3,54) = -1
      iref(55)  =  2
      idimp(1,55) =  1
      idimp(2,55) =  2
      idimp(3,55) =  3
      iflip(1,55) = -1
      iflip(2,55) = -1
      iflip(3,55) = -1
      iref(56)  =  1
      idimp(1,56) =  1
      idimp(2,56) =  2
      idimp(3,56) =  3
      iflip(1,56) = -1
      iflip(2,56) = -1
      iflip(3,56) = -1
 
       return
       end
 
 
       subroutine loadsymsstob(iref,idimp,iflip)
       implicit real *8 (a-h,o-z)
 
       dimension iref(*), idimp(3,*), iflip(3,*)
 
      iref( 1)  =  1
      idimp(1, 1) =  1
      idimp(2, 1) =  2
      idimp(3, 1) =  3
      iflip(1, 1) =  1
      iflip(2, 1) =  1
      iflip(3, 1) =  1
      iref( 2)  =  2
      idimp(1, 2) =  1
      idimp(2, 2) =  2
      idimp(3, 2) =  3
      iflip(1, 2) =  1
      iflip(2, 2) =  1
      iflip(3, 2) =  1
      iref( 3)  =  2
      idimp(1, 3) =  1
      idimp(2, 3) =  2
      idimp(3, 3) =  3
      iflip(1, 3) = -1
      iflip(2, 3) =  1
      iflip(3, 3) =  1
      iref( 4)  =  1
      idimp(1, 4) =  1
      idimp(2, 4) =  2
      idimp(3, 4) =  3
      iflip(1, 4) = -1
      iflip(2, 4) =  1
      iflip(3, 4) =  1
      iref( 5)  =  2
      idimp(1, 5) =  2
      idimp(2, 5) =  1
      idimp(3, 5) =  3
      iflip(1, 5) =  1
      iflip(2, 5) =  1
      iflip(3, 5) =  1
      iref( 6)  =  3
      idimp(1, 6) =  1
      idimp(2, 6) =  2
      idimp(3, 6) =  3
      iflip(1, 6) =  1
      iflip(2, 6) =  1
      iflip(3, 6) =  1
      iref( 7)  =  3
      idimp(1, 7) =  1
      idimp(2, 7) =  2
      idimp(3, 7) =  3
      iflip(1, 7) = -1
      iflip(2, 7) =  1
      iflip(3, 7) =  1
      iref( 8)  =  2
      idimp(1, 8) =  2
      idimp(2, 8) =  1
      idimp(3, 8) =  3
      iflip(1, 8) = -1
      iflip(2, 8) =  1
      iflip(3, 8) =  1
      iref( 9)  =  2
      idimp(1, 9) =  2
      idimp(2, 9) =  1
      idimp(3, 9) =  3
      iflip(1, 9) =  1
      iflip(2, 9) = -1
      iflip(3, 9) =  1
      iref(10)  =  3
      idimp(1,10) =  1
      idimp(2,10) =  2
      idimp(3,10) =  3
      iflip(1,10) =  1
      iflip(2,10) = -1
      iflip(3,10) =  1
      iref(11)  =  3
      idimp(1,11) =  1
      idimp(2,11) =  2
      idimp(3,11) =  3
      iflip(1,11) = -1
      iflip(2,11) = -1
      iflip(3,11) =  1
      iref(12)  =  2
      idimp(1,12) =  2
      idimp(2,12) =  1
      idimp(3,12) =  3
      iflip(1,12) = -1
      iflip(2,12) = -1
      iflip(3,12) =  1
      iref(13)  =  1
      idimp(1,13) =  1
      idimp(2,13) =  2
      idimp(3,13) =  3
      iflip(1,13) =  1
      iflip(2,13) = -1
      iflip(3,13) =  1
      iref(14)  =  2
      idimp(1,14) =  1
      idimp(2,14) =  2
      idimp(3,14) =  3
      iflip(1,14) =  1
      iflip(2,14) = -1
      iflip(3,14) =  1
      iref(15)  =  2
      idimp(1,15) =  1
      idimp(2,15) =  2
      idimp(3,15) =  3
      iflip(1,15) = -1
      iflip(2,15) = -1
      iflip(3,15) =  1
      iref(16)  =  1
      idimp(1,16) =  1
      idimp(2,16) =  2
      idimp(3,16) =  3
      iflip(1,16) = -1
      iflip(2,16) = -1
      iflip(3,16) =  1
      iref(17)  =  2
      idimp(1,17) =  3
      idimp(2,17) =  1
      idimp(3,17) =  2
      iflip(1,17) =  1
      iflip(2,17) =  1
      iflip(3,17) =  1
      iref(18)  =  3
      idimp(1,18) =  1
      idimp(2,18) =  3
      idimp(3,18) =  2
      iflip(1,18) =  1
      iflip(2,18) =  1
      iflip(3,18) =  1
      iref(19)  =  3
      idimp(1,19) =  1
      idimp(2,19) =  3
      idimp(3,19) =  2
      iflip(1,19) = -1
      iflip(2,19) =  1
      iflip(3,19) =  1
      iref(20)  =  2
      idimp(1,20) =  3
      idimp(2,20) =  1
      idimp(3,20) =  2
      iflip(1,20) = -1
      iflip(2,20) =  1
      iflip(3,20) =  1
      iref(21)  =  3
      idimp(1,21) =  2
      idimp(2,21) =  3
      idimp(3,21) =  1
      iflip(1,21) =  1
      iflip(2,21) =  1
      iflip(3,21) =  1
      iref(22)  =  3
      idimp(1,22) =  2
      idimp(2,22) =  3
      idimp(3,22) =  1
      iflip(1,22) = -1
      iflip(2,22) =  1
      iflip(3,22) =  1
      iref(23)  =  3
      idimp(1,23) =  2
      idimp(2,23) =  3
      idimp(3,23) =  1
      iflip(1,23) =  1
      iflip(2,23) = -1
      iflip(3,23) =  1
      iref(24)  =  3
      idimp(1,24) =  2
      idimp(2,24) =  3
      idimp(3,24) =  1
      iflip(1,24) = -1
      iflip(2,24) = -1
      iflip(3,24) =  1
      iref(25)  =  2
      idimp(1,25) =  3
      idimp(2,25) =  1
      idimp(3,25) =  2
      iflip(1,25) =  1
      iflip(2,25) = -1
      iflip(3,25) =  1
      iref(26)  =  3
      idimp(1,26) =  1
      idimp(2,26) =  3
      idimp(3,26) =  2
      iflip(1,26) =  1
      iflip(2,26) = -1
      iflip(3,26) =  1
      iref(27)  =  3
      idimp(1,27) =  1
      idimp(2,27) =  3
      idimp(3,27) =  2
      iflip(1,27) = -1
      iflip(2,27) = -1
      iflip(3,27) =  1
      iref(28)  =  2
      idimp(1,28) =  3
      idimp(2,28) =  1
      idimp(3,28) =  2
      iflip(1,28) = -1
      iflip(2,28) = -1
      iflip(3,28) =  1
      iref(29)  =  2
      idimp(1,29) =  3
      idimp(2,29) =  1
      idimp(3,29) =  2
      iflip(1,29) =  1
      iflip(2,29) =  1
      iflip(3,29) = -1
      iref(30)  =  3
      idimp(1,30) =  1
      idimp(2,30) =  3
      idimp(3,30) =  2
      iflip(1,30) =  1
      iflip(2,30) =  1
      iflip(3,30) = -1
      iref(31)  =  3
      idimp(1,31) =  1
      idimp(2,31) =  3
      idimp(3,31) =  2
      iflip(1,31) = -1
      iflip(2,31) =  1
      iflip(3,31) = -1
      iref(32)  =  2
      idimp(1,32) =  3
      idimp(2,32) =  1
      idimp(3,32) =  2
      iflip(1,32) = -1
      iflip(2,32) =  1
      iflip(3,32) = -1
      iref(33)  =  3
      idimp(1,33) =  2
      idimp(2,33) =  3
      idimp(3,33) =  1
      iflip(1,33) =  1
      iflip(2,33) =  1
      iflip(3,33) = -1
      iref(34)  =  3
      idimp(1,34) =  2
      idimp(2,34) =  3
      idimp(3,34) =  1
      iflip(1,34) = -1
      iflip(2,34) =  1
      iflip(3,34) = -1
      iref(35)  =  3
      idimp(1,35) =  2
      idimp(2,35) =  3
      idimp(3,35) =  1
      iflip(1,35) =  1
      iflip(2,35) = -1
      iflip(3,35) = -1
      iref(36)  =  3
      idimp(1,36) =  2
      idimp(2,36) =  3
      idimp(3,36) =  1
      iflip(1,36) = -1
      iflip(2,36) = -1
      iflip(3,36) = -1
      iref(37)  =  2
      idimp(1,37) =  3
      idimp(2,37) =  1
      idimp(3,37) =  2
      iflip(1,37) =  1
      iflip(2,37) = -1
      iflip(3,37) = -1
      iref(38)  =  3
      idimp(1,38) =  1
      idimp(2,38) =  3
      idimp(3,38) =  2
      iflip(1,38) =  1
      iflip(2,38) = -1
      iflip(3,38) = -1
      iref(39)  =  3
      idimp(1,39) =  1
      idimp(2,39) =  3
      idimp(3,39) =  2
      iflip(1,39) = -1
      iflip(2,39) = -1
      iflip(3,39) = -1
      iref(40)  =  2
      idimp(1,40) =  3
      idimp(2,40) =  1
      idimp(3,40) =  2
      iflip(1,40) = -1
      iflip(2,40) = -1
      iflip(3,40) = -1
      iref(41)  =  1
      idimp(1,41) =  1
      idimp(2,41) =  2
      idimp(3,41) =  3
      iflip(1,41) =  1
      iflip(2,41) =  1
      iflip(3,41) = -1
      iref(42)  =  2
      idimp(1,42) =  1
      idimp(2,42) =  2
      idimp(3,42) =  3
      iflip(1,42) =  1
      iflip(2,42) =  1
      iflip(3,42) = -1
      iref(43)  =  2
      idimp(1,43) =  1
      idimp(2,43) =  2
      idimp(3,43) =  3
      iflip(1,43) = -1
      iflip(2,43) =  1
      iflip(3,43) = -1
      iref(44)  =  1
      idimp(1,44) =  1
      idimp(2,44) =  2
      idimp(3,44) =  3
      iflip(1,44) = -1
      iflip(2,44) =  1
      iflip(3,44) = -1
      iref(45)  =  2
      idimp(1,45) =  2
      idimp(2,45) =  1
      idimp(3,45) =  3
      iflip(1,45) =  1
      iflip(2,45) =  1
      iflip(3,45) = -1
      iref(46)  =  3
      idimp(1,46) =  1
      idimp(2,46) =  2
      idimp(3,46) =  3
      iflip(1,46) =  1
      iflip(2,46) =  1
      iflip(3,46) = -1
      iref(47)  =  3
      idimp(1,47) =  1
      idimp(2,47) =  2
      idimp(3,47) =  3
      iflip(1,47) = -1
      iflip(2,47) =  1
      iflip(3,47) = -1
      iref(48)  =  2
      idimp(1,48) =  2
      idimp(2,48) =  1
      idimp(3,48) =  3
      iflip(1,48) = -1
      iflip(2,48) =  1
      iflip(3,48) = -1
      iref(49)  =  2
      idimp(1,49) =  2
      idimp(2,49) =  1
      idimp(3,49) =  3
      iflip(1,49) =  1
      iflip(2,49) = -1
      iflip(3,49) = -1
      iref(50)  =  3
      idimp(1,50) =  1
      idimp(2,50) =  2
      idimp(3,50) =  3
      iflip(1,50) =  1
      iflip(2,50) = -1
      iflip(3,50) = -1
      iref(51)  =  3
      idimp(1,51) =  1
      idimp(2,51) =  2
      idimp(3,51) =  3
      iflip(1,51) = -1
      iflip(2,51) = -1
      iflip(3,51) = -1
      iref(52)  =  2
      idimp(1,52) =  2
      idimp(2,52) =  1
      idimp(3,52) =  3
      iflip(1,52) = -1
      iflip(2,52) = -1
      iflip(3,52) = -1
      iref(53)  =  1
      idimp(1,53) =  1
      idimp(2,53) =  2
      idimp(3,53) =  3
      iflip(1,53) =  1
      iflip(2,53) = -1
      iflip(3,53) = -1
      iref(54)  =  2
      idimp(1,54) =  1
      idimp(2,54) =  2
      idimp(3,54) =  3
      iflip(1,54) =  1
      iflip(2,54) = -1
      iflip(3,54) = -1
      iref(55)  =  2
      idimp(1,55) =  1
      idimp(2,55) =  2
      idimp(3,55) =  3
      iflip(1,55) = -1
      iflip(2,55) = -1
      iflip(3,55) = -1
      iref(56)  =  1
      idimp(1,56) =  1
      idimp(2,56) =  2
      idimp(3,56) =  3
      iflip(1,56) = -1
      iflip(2,56) = -1
      iflip(3,56) = -1
 
       return
       end
 
 
