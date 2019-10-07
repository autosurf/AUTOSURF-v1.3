PROGRAM evaluatePES

use dynamic_parameters
!-----------------------------------------------------------------------------------
implicit none
 character (len=40) :: NAME1
 real*8 :: xi(4),V,pii,th1,th2,phi,CONVE

 pii=acos(-1d0)
 CONVE=349.75d0 ! conversion factor: switch energies from kcal/mol to wave numbers

 ! name of the PES-file
 NAME1='PES-CO2dimer-1064'


 ! ---------------------------------------------------------------------------------
 ! ***  evaluate the PES  ***
 ! ---------------------------------------------------------------------------------
 !!  xi(1) is R, the distance between centers of mass (in Angstroms).             !!
 !!  xi(2) and xi(3) are cos(theta1) and cos(theta2) and range from (-1,1).       !!
 !!  xi(4) is the dihedral angle, in radians, with range: (0,2pi).                !!

 ! coordinate R
 xi(1)=5.64d0     
 ! coordinate theta1
 th1=0.d0        
 xi(2)=dcos(th1*pii/180.d0)
 ! coordinate theta2
 th2=180.d0      
 xi(3)=dcos(th2*pii/180.d0)
 ! coordinate phi
 phi=0d0         
 xi(4)=phi*pii/180.d0

 call PES(xi,V,NAME1)
 V=V*CONVE

 write(*,501)xi
 501 format(' coordinates [R;cos(th1);cos(th2);phi] -> ',4(F6.2))
 write(6,*) '                            potential -> ',V

end program evaluatePES
