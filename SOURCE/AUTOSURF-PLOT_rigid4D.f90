!*********************************  A U T O S U R F  *******************************
!===================================================================================
!-----------------------------------------------------------------------------------
!-                                                                                 -
!-        AUTOSURF Package: A set of programs for the automated construction       -
!-              of Potential Energy Surfaces on van der Waals systems              -
!-                                                                                 -
!-----------------------------------------------------------------------------------
!===================================================================================
!***********************************************************************************
!-      "AUTOSURF-PLOT_rigid4D": PROGRAM for ...                                   -
!-----------------------------------------------------------------------------------
!-       Input files: 
!-        "input-AUTOSURF-PES.dat" (main input file used ) 
!-        & "PES-file" (output file from program AUTOSURF-PES_rigid4D)             -
!-                                                                                 -
!***********************************************************************************


PROGRAM AUTOSURF_PLOT

use dynamic_parameters
!-----------------------------------------------------------------------------------
implicit none
 integer :: i,j,flag,i2,i3,k,nline,input,ncount,gnup1
 real*8 :: tampon,V,pii,poten2,qmat(3,6),inert(3,3),bohr,range,dx,tampon1
 real*8 :: coord1(101),coord2(101),temp,temp3,temp2(5)
 real*8 :: djac(4),shift,atobohr,xass
 real*8 :: cart3(15),EGrid2(101,101),RGrid2(101,101)
 real*8 :: Rtamp
 character (len=40) :: NAME1,NAME2
 character (len=4) :: charid
 character (len=100) :: LINE
 character (len=100),dimension(6) :: LINE1D
 real*8 :: xrmin, xrmax, th1min, th1max, th2min, th2max, phimin, phimax, dd,dd1,dd2
 real*8 :: R,th1,th2,phi,Emin,Emax,dE
 real*8,dimension(4) :: ran_vec,xi,xxrmin,xxrmax
 integer :: npts,gnup,nitv
 character(len=10) :: bdate1,bdate2,bdate3
 character(len=5) :: charmyid
 character(len=8) :: xchar,xchar1
 character(len=15),dimension(4) :: charunitsE
 character(len=15),dimension(2) :: charunitsD
 integer :: summyid,ncont,ncont1,numpoints_low_count,old_numprocs
 integer :: xpass,numpoints_high_count,n_test,xunitE,npts1,npts2
 integer,dimension(1) :: xseed
 integer,dimension(8) :: date_time
 integer,parameter :: xunitD=2, xnum_intervals_R=20
 real*8,allocatable :: xgrad(:),xR(:),xMAPot(:)
 real*8,allocatable :: xrms_R1(:),xrms_R2(:),xrms_R3(:)
 real*8,allocatable :: xR_count1(:),xR_count2(:),xR_count3(:),xcor(:)
 real*8 :: xpot,CONVE,CONVD,xdR,minR,maxR,xtemp,x1,x2,x3,x4,x5,jac3(4)
 real*8 :: Min_E,Max_R,Max_E3,bot_seed,xGlob_min,dR,xmastot,xmas1,xmas2
 logical :: logica1, logica2
 real*8,parameter :: hart2wn=219474.6313702d0, hart2meV=27211.38602d0
 real*8,parameter :: ang2bohr=1.0d0/0.529177249d0, bohr2ang=0.529177249d0
 real*8,parameter :: hart2kcl=4.359744650D-18/4184*6.022140857D23
 !real*8,parameter :: hart2kcl=627.5095d0
 real*8,parameter :: kcal2wn=349.757D0
INTERFACE! 
  FUNCTION func(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi     
    REAL(SP) :: func
  END FUNCTION func
end interface
INTERFACE! 
  FUNCTION func1(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi     
    REAL(SP) :: func1
  END FUNCTION func1
end interface
INTERFACE!
  FUNCTION dfunc_actual_anal1(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_anal1
  END FUNCTION dfunc_actual_anal1
end interface
INTERFACE! 
  FUNCTION dfunc_actual_seed(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_seed
  END FUNCTION dfunc_actual_seed
end interface
INTERFACE! 
  FUNCTION dfunc_actual_anal2(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_anal2
  END FUNCTION dfunc_actual_anal2
end interface
INTERFACE! 
  FUNCTION func_actual(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual
  END FUNCTION func_actual
end interface
INTERFACE! 
  FUNCTION func_actual_min(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_min
  END FUNCTION func_actual_min
end interface
INTERFACE! 
  FUNCTION func_actual_seed(xi)  
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_seed
  END FUNCTION func_actual_seed
end interface
INTERFACE!
  FUNCTION dfunc(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)) :: dfunc
  END FUNCTION dfunc
end interface
!-----------------------------------------------------------------------------------

 pii=acos(-1d0)
 atobohr=1.0d0/0.529177249d0 
 ugrad=atobohr*hart2kcl
 bohr=0.529177249d0
 CONVE=kcal2wn
 CONVE=219474.6313702d0/(4.359744650D-18/4184*6.022140857D23)
! hart2kcl=627.5095d0
! CONVE=1.d0/627.5095d0*219474.6313702d0
 xGlob_min=-86995.895+86992.503d0
 call system('clear')

 ! find out the name of the PES-file
 write(*,*)'Enter the name of the file containing PES-data:'
 read(*,*)NAME1
 NAME1=trim(adjustl(NAME1))
 nline=scan(NAME1,' ')-1
 inquire(file=NAME1(1:nline),exist=logica1)
 if(.not.logica1)then
   write(*,*)'The file "'//NAME1(1:nline)//'" does not exist'
   write(*,*)'Enter the name of the file:'
   read(*,*)NAME1
   NAME1=trim(adjustl(NAME1))
   nline=scan(NAME1,' ')-1
   inquire(file=NAME1(1:nline),exist=logica2)
   if(.not.logica2)then
     write(*,*)'The file "'//NAME1(1:nline)//'" does not exist either!'
     write(*,*)'Check the name and try again...'
     stop 
   endif
 endif
 ! initialize the PES...
 jac3=0.d0
 call PES(jac3,V,NAME1)

 ! ALLOCATE
 allocate(xR(xnum_intervals_R+1))
 allocate(xrms_R1(xnum_intervals_R+1),xrms_R2(xnum_intervals_R+1))
 allocate(xR_count1(xnum_intervals_R+1),xR_count2(xnum_intervals_R+1))
 allocate(xrms_R3(xnum_intervals_R+1),xR_count3(xnum_intervals_R+1))
 allocate(xMAPot(xnum_intervals_R+1))


 ! *********************************************************************************
 ! ---------------------------------------------------------------------------------
 ! ---        STARTING POINT
 ! ---------------------------------------------------------------------------------
 ! *********************************************************************************

 100 continue
 write(*,*)
 write(*,*)'(1) = PES summary'
 write(*,*)'(2) = evaluate the PES'
 write(*,*)'(3) = 1D cut of the PES'
 write(*,*)'(4) = 2D cut of the PES'
 write(*,'(A81)')'(5) = R-optimized 2D cut of the PES: "U vs. extended angles" (planar geometries)'
 write(*,*)'(6) = Prepare AUTOSURF_XGW.dat'
 write(*,*)'(7) = EXIT'
 101 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=102)input
 if (input==1) goto 200
 if (input==2) goto 1000
 if (input==3) goto 2000
 if (input==4) goto 3000
 if (input==5) goto 4000
 if (input==6) goto 5000
 if (input==7) goto 666
 102 write(*,*)'enter a valid integer: 1, 2, 3, 4, 5 or 6'
 goto 101


 200 continue
 ! ---------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------
 ! ***  PES summary  ***
 ! ---------------------
 write(*,*)
 write(*,*)'(1) = PES info'
 write(*,*)'(2) = Report of errors'
 write(*,*)'(3) = BACKUP (save ab initio data used to fit)'
! write(*,*)'(3) = BACKUP (save ab initio data and expansion coefficients)'
 write(*,*)'(4) = Restart'
 write(*,*)'(5) = EXIT'
 201 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=202)input
 if (input==1) goto 210
 if (input==2) goto 211
 if (input==3) goto 212
 if (input==4) then
   call system('clear')
   goto 100
 endif
 if (input==5) goto 666
 202 write(*,*)'enter a valid integer: 1, 2, 3, 4, or 5'
 goto 201

 210 continue
 write(*,*)
 write(*,*)
 write(*,*)'PES info...'

 write(*,'(A30,I6)')' High-degree basis size:      ',basis_1
 write(*,'(A30,I6)')' Lower-degree basis size:     ',basis_2
 write(*,'(A30,I6)')' Size of minimal basis:       ',basis_3
 write(*,*)
 write(*,'(A32,F20.6)')' Ref. Asymptotic energy:        ',ass
 write(*,'(A32,F20.6)')' Max. E allowed above asymptote:',Max_E
 write(*,'(A35,I6)')' Total number of points in the fit:  ',count3
 write(*,*) 
!  write(*,*) order_1
!  write(*,*) order_2
!  write(*,*) order_3
!  write(*,*) order_4
!  write(*,*) maxpoints
!  write(*,*) mass
!  write(*,*) rmax
!  write(*,*) rmin
!  write(*,*) low_grid
!  write(*,*) count_seed
! write(*,*)
! write(*,*)Max_E_seed
! write(*,*)'.. more information to be included soon...'
 write(*,*)
 write(*,*)'(1) = restart'
 write(*,*)'(2) = EXIT'
 2101 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=2102)input
 if (input==1) then
   call system('clear')
   goto 100
 endif
 if (input==2) goto 666
 2102 write(*,*)'enter a valid integer: 1 or 2'
 goto 2101


 211 continue
 ! Random test to estimate the final error -----------------------------------------
 ! -------------------------------------------------------------------------------

 write(*,*)
 write(*,*)'(1) = all the volume of the configurations space'
 write(*,*)'(2) = minimal volume used to fit'
 write(*,*)'(3) = user-defined volume'
 2116 continue
 write(*,'(A45)',ADVANCE='YES')'volume of the configuration space to sample?'
 read(*,'(I2)',ADVANCE='YES',ERR=2117)input
 if(input==1)goto 2118
 if(input==2)goto 2118
 if(input==3)goto 2118
 2117 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 2116

 2118 continue
 write(*,*)
 write(*,'(A45)',ADVANCE='YES')'Number of geometries to be randomly tested ='
 read(*,*,ERR=2113)npts
 if(npts==0)goto 2113
 goto 2115
 2113 write(*,*)'enter a valid (integer) number of geometries'
 goto 2118

 2115 continue
 if (input==1) then
   xxrmin(1)=rmin(1)
   xxrmax(1)=rmax(1)
   xxrmin(2)=-0.9999d0
   xxrmax(2)=0.9999d0
   xxrmin(3)=-0.9999d0
   xxrmax(3)=0.9999d0
   xxrmin(4)=-0.9999d0*pii
   xxrmax(4)=0.9999d0*pii
 elseif (input==2) then
   xxrmin=rmin
   xxrmax=rmax
 else if (input==3) then
   write(*,*)
   write(*,*)'range of coordinates to sample:'
   write(*,'(A20)',ADVANCE='NO')'R_min (Ang.) ='
   read(*,*)xxrmin(1)
   write(*,'(A20)',ADVANCE='NO')'R_max (Ang.) ='
   read(*,*)xxrmax(1)
   write(*,'(A20)',ADVANCE='NO')'theta1_min (deg.) ='
   read(*,*)th1
   xxrmin(2)=cos(th1*pii/180.d0)
   write(*,'(A20)',ADVANCE='NO')'theta1_max (deg.) ='
   read(*,*)th1
   xxrmax(2)=cos(th1*pii/180.d0)
   write(*,'(A20)',ADVANCE='NO')'theta2_min (deg.) = '
   read(*,*)th2
   xxrmin(3)=cos(th2*pii/180.d0)
   write(*,'(A20)',ADVANCE='NO')'theta2_max (deg.) = '
   read(*,*)th2
   xxrmax(3)=cos(th2*pii/180.d0)
   write(*,'(A20)',ADVANCE='NO')'phi_min (deg.) = '
   read(*,*)phi
   xxrmin(4)=phi*pii/180.d0
   write(*,'(A20)',ADVANCE='NO')'phi_max (deg.) = '
   read(*,*)phi
   xxrmax(4)=phi*pii/180.d0
 endif

 write(*,*)
 write(*,*)'Please, be patient,'
 write(*,*)'this could take some time...'
 write(*,*)
 xrms_R1=0.d0
 xrms_R2=0.d0
 xrms_R3=0.d0
 xMAPot=0.d0
 xR_count1=0.d0
 xR_count2=0.d0
 xR_count3=0.d0
 ! prepare R-intervals for the final summary of errors 
 xdR=(xxrmax(1)-xxrmin(1))/dble(xnum_intervals_R)
 do j=1,xnum_intervals_R+1
   xR(j)=xxrmin(1)+(j-1)*xdR
 enddo
 ! select random geometries to test
 do i=1,npts
   27 call random_number(ran_vec)
   do j=1,4
     range=xxrmax(j)-xxrmin(j)
     ran_vec(j)=xxrmin(j)+range*ran_vec(j)
   enddo
   call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
   call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
   if(dist_flag==1) goto 27! exclude points with atoms too close
   xi=ran_vec(:)
   ! exclude points based on energy (seed-grid-PES estimate + low-grid cutoff)
   if(low_grid>0)then
     temp3=func_actual_seed(xi)
     if(temp3>Max_E_seed) goto 27
   endif
   ! exclude points with E > ass + E_range
   temp=func_actual_min(xi)!(use min-PES to estimate E)
!   if(subzero==0)then
     if (temp>Max_E) goto 27
!   else
!     if (temp+temp3>Max_E) goto 27
!   endif
!   if(subzero==0)then

!   call PES(xi,temp,NAME1)
     temp=func_actual(xi)
!   else
!     temp=func_actual(xi)+temp3
!   endif
   if(input==2)tampon=func(xi)  
   if(input/=2)tampon=func1(xi)  
   if((tampon==0.d0))write(*,*)'Warning, energy difference equal zero! Coordinates: ',xi,tampon,temp
!   if((tampon==0.d0).and.(temp>=Max_E))write(*,*)'Warning, energy difference equal zero! Coordinates: ',xi,tampon
   !if(tampon==0.d0)write(*,*)xi,ran_vec(4),pii,ran_vec(4)*1.998d0*pii
   if(i==500)write(*,*)'geometries:'
   if(i==500)write(*,*)'500'
   if(i==1000)write(*,*)'1000'
   if(i==2000)write(*,*)'2000'
   if(i==4000)write(*,*)'4000'
   if(i==8000)write(*,*)'8000'
   if(i==15000)write(*,*)'15000'
   if(i==30000)write(*,*)'30000'
   if(i==60000)write(*,*)'60000'
   if(i==90000)write(*,*)'90000'
   if(i==120000)write(*,*)'120000'
   if(i==150000)write(*,*)'150000'
   if(i==180000)write(*,*)'180000'
   if(i==210000)write(*,*)'210000'
   if(i==240000)write(*,*)'240000'
   if(i==270000)write(*,*)'270000'
   do j=1,xnum_intervals_R! statistics for each R-interval
     if((xR(j)<=xi(1)).and.(xi(1)<xR(j+1)))then
       xMAPot(j)=xMAPot(j)+abs(temp-ass)
       xR_count1(j)=xR_count1(j)+1.d0
       xrms_R1(j)=xrms_R1(j)+abs(tampon)
       if(temp<ass+(0.05d0))then
!       if(temp<ass+(0.05d0/hart2kcl*CONVE))then                                                      ! check units (if kcal/mol are not used...)
         xR_count2(j)=xR_count2(j)+1.d0
         xrms_R2(j)=xrms_R2(j)+abs(tampon)
       endif
!       if(temp<(0.3-236154.954-1.613))then
!write(*,*)ass,xGlob_min
!5       if(temp<(0.3+ass+xGlob_min))then
!       if(temp<((0.3/hart2kcl*CONVE)+Glob_min))then              !!! Glob_min   
!5         xR_count3(j)=xR_count3(j)+1.d0
!5         xrms_R3(j)=xrms_R3(j)+abs(tampon)
!5       endif
     endif
   enddo
 enddo
 ! global statistics
 do j=1,xnum_intervals_R
   xMAPot(xnum_intervals_R+1)=xMAPot(xnum_intervals_R+1)+xMAPot(j)
   xrms_R1(xnum_intervals_R+1)=xrms_R1(xnum_intervals_R+1)+xrms_R1(j)
   xrms_R2(xnum_intervals_R+1)=xrms_R2(xnum_intervals_R+1)+xrms_R2(j)
!5   xrms_R3(xnum_intervals_R+1)=xrms_R3(xnum_intervals_R+1)+xrms_R3(j)
   xR_count1(xnum_intervals_R+1)=xR_count1(xnum_intervals_R+1)+xR_count1(j)
   xR_count2(xnum_intervals_R+1)=xR_count2(xnum_intervals_R+1)+xR_count2(j)
!5   xR_count3(xnum_intervals_R+1)=xR_count3(xnum_intervals_R+1)+xR_count3(j)
 enddo
 do j=1,xnum_intervals_R+1
   if(xR_count1(j)==0)then
     xMAPot(j)=0.d0
     xrms_R1(j)=0.d0
   else
     xMAPot(j)=xMAPot(j)/xR_count1(j)
     xrms_R1(j)=sqrt(xrms_R1(j)/xR_count1(j))
   endif
   if(xR_count2(j)==0)then
     xrms_R2(j)=0.d0
   else
     xrms_R2(j)=sqrt(xrms_R2(j)/xR_count2(j))
   endif
!5   if(xR_count3(j)==0)then
!5     xrms_R3(j)=0.d0
!5   else
!5     xrms_R3(j)=sqrt(xrms_R3(j)/xR_count3(j))
!5   endif
 enddo
 open(unit=200,file='errores-'//NAME1(1:nline)//'.dat')
 write(200,*)'RMS1 = Global'
 write(200,*)'RMS2 = below asymptote' 
!5 write(200,*)'RMS3 = 0.3 kcal/mol (~100 cm-1) above minimum' 
 write(200,*)
!5 write(*,*)
!5 write(*,*)'RMS1 = Global'
!5 write(*,*)'RMS2 = below asymptote' 
!5 write(*,*)'RMS3 = 0.3 kcal/mol (~100 cm-1) above minimum' 
!5 write(*,*)
 write(200,'(A97)')'   Ri    Rf    <|V(R)|>      RMS1      %      pts    RMS2      %      pts                        '
!5 write(*,'(A97)')'   Ri    Rf    <|V(R)|>      RMS1      %      pts    RMS2      %      pts    RMS3      %      pts'
 do j=1,xnum_intervals_R
   x1=xrms_R1(j)*100.d0/xMAPot(j)
   x2=xrms_R2(j)*100.d0/xMAPot(j)
   if(xMAPot(j)==0)x1=0.d0
   if(xMAPot(j)==0)x2=0.d0
!5   x3=xrms_R3(j)*100.d0/xMAPot(j)
   x5=(xR(j)+xR(j+1))/2.d0
   write(200,'(2F6.2,ES13.5,2(ES11.3,F6.2,I7))')xR(j),xR(j+1),xMAPot(j),xrms_R1(j),&
    x1,int(xR_count1(j)),xrms_R2(j),x2,int(xR_count2(j))
!5   write(*,'(2F6.2,ES13.5,3(ES11.3,F6.2,I7))')xR(j),xR(j+1),xMAPot(j),xrms_R1(j),&
!5    x1,int(xR_count1(j)),xrms_R2(j),x2,int(xR_count2(j)),xrms_R3(j),x3,int(xR_count3(j))
 enddo
 write(200,*)
!5 write(*,*)
 x1=xrms_R1(xnum_intervals_R+1)*100.d0/xMAPot(xnum_intervals_R+1)
 x2=xrms_R2(xnum_intervals_R+1)*100.d0/xMAPot(xnum_intervals_R+1)
!5 x3=xrms_R3(xnum_intervals_R+1)*100.d0/xMAPot(xnum_intervals_R+1)
 write(200,'(A12,ES13.5,2(ES11.3,F6.2,I7))')'total:',xMAPot(xnum_intervals_R+1), &
 xrms_R1(xnum_intervals_R+1),x1,int(xR_count1(xnum_intervals_R+1)), &
 xrms_R2(xnum_intervals_R+1),x2,int(xR_count2(xnum_intervals_R+1))!5, &
!5 xrms_R3(xnum_intervals_R+1),x3,int(xR_count3(xnum_intervals_R+1))
!5 write(*,'(A12,ES13.5,3(ES11.3,F6.2,I7))')'total:',xMAPot(xnum_intervals_R+1), &
!5 xrms_R1(xnum_intervals_R+1),x1,int(xR_count1(xnum_intervals_R+1)), &
!5 xrms_R2(xnum_intervals_R+1),x2,int(xR_count2(xnum_intervals_R+1)), &
!5 xrms_R3(xnum_intervals_R+1),x3,int(xR_count3(xnum_intervals_R+1))
 close(200)
 write(*,*)
 write(*,*)'Errors information saved in the file: "errores-'//NAME1(1:nline)//'.dat".'
 write(*,*)
 write(*,*)'(1) = restart'
 write(*,*)'(2) = EXIT'
 2111 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=2112)input
 if (input==1) then
   call system('clear')
   goto 100
 endif
 if (input==2) goto 666
 2112 write(*,*)'enter a valid integer: 1 or 2'
 goto 2111


 212 continue 
 !  BACKUP  ------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------

 ! check if previous files containing ab initio data already exists
 inquire(file='AbINITIO.dat',exist=logica1)
 if (logica1) then
   write(*,*)'The file "AbINITIO.dat" already exists.' 
   write(*,*)
   write(*,*)'(1) = yes'
   write(*,*)'(2) = no'
   write(*,*)'(3) = restart'
   write(*,*)'(4) = EXIT'
   2121 continue
   write(*,'(A26)',ADVANCE='YES')'Do you want to replace it?' 
   read(*,'(I2)',ADVANCE='YES',ERR=2122)input
   if (input==1) then
     open(unit=200,file='AbINITIO.dat')
     goto 2124
   endif
   if (input==2) then
     2123 continue
     write(*,*)'Enter a new name for the file:'
     read(*,*)NAME2
     NAME2=trim(adjustl(NAME2))
     nline=scan(NAME2,' ')-1
     inquire(file=NAME2(1:nline),exist=logica1)
     if(logica1)then
       write(*,*)'The file "'//NAME2(1:nline)//'" already exist'
       goto 2123
     endif
     open(unit=200,file=NAME2(1:nline))
     goto 2124
   endif
   if (input==3) then
     call system('clear')
     goto 100
   endif
   if (input==4) goto 666
   2122 write(*,*)'enter a valid integer: 1, 2, 3 or 4'
   goto 2121
 endif
 open(unit=200,file='AbINITIO.dat')
 2124 continue
 IF(low_grid>0)THEN
  inquire(file='AbINITIO_low.dat',exist=logica1)
  if (logica1) then
   write(*,*)'The file "AbINITIO_low.dat" already exists.' 
   write(*,*)
   write(*,*)'(1) = yes'
   write(*,*)'(2) = no'
   write(*,*)'(3) = restart'
   write(*,*)'(4) = EXIT'
   2125 continue
   write(*,'(A26)',ADVANCE='YES')'Do you want to replace it?' 
   read(*,'(I2)',ADVANCE='YES',ERR=2126)input
   if (input==1) then
     open(unit=201,file='AbINITIO_low.dat')
     goto 2128
   endif
   if (input==2) then
     2127 continue
     write(*,*)'Enter a new name for the file:'
     read(*,*)NAME2
     NAME2=trim(adjustl(NAME2))
     nline=scan(NAME2,' ')-1
     inquire(file=NAME2(1:nline),exist=logica1)
     if(logica1)then
       write(*,*)'The file "'//NAME2(1:nline)//'" already exist'
       goto 2127
     endif
     open(unit=201,file=NAME2(1:nline))
     goto 2128
   endif
   if (input==3) then
     call system('clear')
     goto 100
   endif
   if (input==4) goto 666
   2126 write(*,*)'enter a valid integer: 1, 2, 3 or 4'
   goto 2125
  endif
  open(unit=201,file='AbINITIO_low.dat')
  2128 continue
 ENDIF
 write(*,*)
 write(*,*)'Saving high-level ab initio data...'
 ncont=0
 do i=1,count3       
   if(mod(i,(2*symparts))==1)then
    ncont=ncont+1
    jac=coords(i,:) 
!    call PES(jac,V,NAME1)
    V=func_actual(jac)
    write(200,'(I10,4f15.8,f20.8)')ncont,jac(:),V
   endif
 enddo
 close(200)
 ncont=0
 if (low_grid>0) then
  write(*,*)'Saving low-level ab initio data...'
  do i=1,count_seed
   if(mod(i,(2*symparts))==1)then
       ncont=ncont+1
       jac=coords_seed(i,:)
       V=func_actual_seed(jac)  
       write(201,'(I10,4f15.8,f20.8)')ncont,jac(:),V
   endif
  enddo
 endif
 close(201)
 write(*,*)
 write(*,*)'done!'
 write(*,*)
 write(*,*)'(1) = restart'
 write(*,*)'(2) = EXIT'
 2129 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=2130)input
 if (input==1) then
   call system('clear')
   goto 100
 endif
 if (input==2) goto 666
 2130 write(*,*)'enter a valid integer: 1 or 2'
 goto 2129

! save also the coefficients??
! to be done...
! do i=1,160 ! print coefficients 
!! do i=1,count3
!   write(652,'(4f9.3,301f15.4)')coords(i,:),b2(:,i)
! enddo

 1000 continue
 ! ---------------------------------------------------------------------------------
 ! ***  evaluate the PES  ***
 ! ---------------------------------------------------------------------------------
 write(*,*)
 write(*,*)'define coordinates:'
 write(*,'(A16)',ADVANCE='NO')'R = '
 read(*,*)jac(1)
 write(*,'(A16)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)jac(2)
 jac(2)=cos(jac(2)*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)jac(3)
 jac(3)=cos(jac(3)*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'phi (deg.) = '
 read(*,*)jac(4)
 jac(4)=jac(4)*pii/180.d0
 call PES(jac,V,NAME1)
 V=(V)*CONVE
 write(*,*)'Value of the potential:  ',V
 write(*,501)jac
 501 format('[R;cos(th1);cos(th2);phi]->',4(F6.2))
 write(*,*)
 write(*,*)'(1) = evaluate the PES again'
 write(*,*)'(2) = restart'
 write(*,*)'(3) = EXIT'
 1001 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=104)input
 if (input==1) goto 1000
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 104 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 1001


 2000 continue
 ! ---------------------------------------------------------------------------------
 ! ***  1D cut of the PES  ***
 ! ---------------------------------------------------------------------------------
 write(*,*)
 write(*,*)'(1) = R'
 write(*,*)'(2) = theta1'
 write(*,*)'(3) = theta2'
 write(*,*)'(4) = phi'
 write(*,*)'(5) = restart'
 write(*,*)'(6) = EXIT'
 2001 continue
 write(*,'(A16)',ADVANCE='YES')'U vs. variable?'
 read(*,'(I2)',ADVANCE='YES',ERR=105)input
 if (input==1) goto 2002
 if (input==2) goto 2003
 if (input==3) goto 2004
 if (input==4) goto 2005
 if (input==5) then
   call system('clear')
   goto 100
 endif
 if (input==6) goto 666
 105 write(*,*)'enter a valid integer: 1, 2, 3, 4, 5 or 6'
 goto 2001

 2002 continue
 !  U vs. R   ----------------------------------------------------------------------
 write(*,*)
 write(*,*)'Define fixed coordinates:'
 write(*,'(A16)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=cos(th1*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=cos(th2*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'phi = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,*)'----------------------'
 write(*,'(A20)',ADVANCE='NO')'min. value of R = '
 read(*,*)xrmin
 write(*,'(A20)',ADVANCE='NO')'max. value of R = '
 read(*,*)xrmax
 write(*,'(A20)',ADVANCE='NO')'number of points = '
 read(*,*)npts
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut1D-UvsR-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut1D-UvsR-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=220)
     enddo
     close(2002)
     cycle
     220 close(2002, status="delete")
     ncount=ncount-1
     goto 221
   else
     goto 221
   endif
 enddo
 221 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut1D-UvsR-'//trim(adjustl(charid))//'.dat')
 write(2002,500)acos(jac(2))*180.d0/pii,acos(jac(3))*180.d0/pii,jac(4)*180.d0/pii
 500 format('   R           Potential    (th1=',F5.1,'; th2=',F5.1,'; phi=',F5.1,')')
 ! make the corresponding 1D cut of the PES 
 dd=(xrmax-xrmin)/dble(npts-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts
   jac(1)=xrmin+dble(i-1)*dd
   call PES(jac,V,NAME1)
   poten=(V)*CONVE
   if (poten>Emax) Emax=poten
   if (poten<Emin) Emin=poten
   write(2002,'(F6.2,F20.8)') jac(1),poten
 enddo
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut1D-UvsR-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(0) = no'
 write(*,*)'(1) = yes'
 write(*,*)'Keep the file?'
 read(*,*)input
 if(input==0)then
   call system('rm cut1D-UvsR-'//trim(adjustl(charid))//'.*')
 endif

 write(*,*)
 write(*,*)'(1) = Make another graph U vs. R'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 222 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=106)input
 if (input==1) goto 2002
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 106 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 222

 2003 continue
 !  U vs. Theta1   -----------------------------------------------------------------
 write(*,*)
 write(*,*)'Define fixed coordinates:'
 write(*,'(A16)',ADVANCE='NO')'R = '
 read(*,*)jac(1)
 write(*,'(A16)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=cos(th2*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'phi = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,*)'------------------------'
 write(*,'(A20)',ADVANCE='NO')'min. value of th1 = '
 read(*,*)th1min
 write(*,'(A20)',ADVANCE='NO')'max. value of th1 = '
 read(*,*)th1max
 write(*,'(A20)',ADVANCE='NO')'number of points = '
 read(*,*)npts
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut1D-th1-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut1D-th1-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=223)
     enddo
     close(2002)
     cycle
     223 close(2002, status="delete")
     ncount=ncount-1
     goto 224
   else
     goto 224
   endif
 enddo
 224 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut1D-th1-'//trim(adjustl(charid))//'.dat')
 502 format('    th1           Potential    (R=',F5.1,'; th2=',F5.1,'; phi=',F5.1,')')
 write(2002,502)jac(1),acos(jac(3))*180.d0/pii,jac(4)*180.d0/pii
 ! make the corresponding 1D cut of the PES 
 dd=(th1max-th1min)/dble(npts-1)
 do i=1,npts
   jac(2)=cos((th1min+dble(i-1)*dd)*pii/180.d0)
   call PES(jac,V,NAME1)
   poten=(V)*CONVE
   write(2002,'(F8.2,F20.8)') acos(jac(2))*180.d0/pii,poten
 enddo
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut1D-th1-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(0) = no'
 write(*,*)'(1) = yes'
 write(*,*)'Keep the file?'
 read(*,*)input
 if(input==0)then
   call system('rm cut1D-th1-'//trim(adjustl(charid))//'.*')
 endif

 write(*,*)
 write(*,*)'(1) = Make another graph U vs. theta1'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 225 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=107)input
 if (input==1) goto 2003
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 107 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 225

 2004 continue
 !  U vs. Theta2   -----------------------------------------------------------------
 write(*,*)
 write(*,*)'Define fixed coordinates:'
 write(*,'(A16)',ADVANCE='NO')'R = '
 read(*,*)jac(1)
 write(*,'(A16)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=cos(th1*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'phi = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,*)'------------------------'
 write(*,'(A20)',ADVANCE='NO')'min. value of th2 = '
 read(*,*)th2min
 write(*,'(A20)',ADVANCE='NO')'max. value of th2 = '
 read(*,*)th2max
 write(*,'(A20)',ADVANCE='NO')'number of points = '
 read(*,*)npts
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut1D-th2-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut1D-th2-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=226)
     enddo
     close(2002)
     cycle
     226 close(2002, status="delete")
     ncount=ncount-1
     goto 227
   else
     goto 227
   endif
 enddo
 227 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut1D-th2-'//trim(adjustl(charid))//'.dat')
 503 format('    th2           Potential    (R=',F5.1,'; th1=',F5.1,'; phi=',F5.1,')')
 write(2002,503)jac(1),acos(jac(2))*180.d0/pii,jac(4)*180.d0/pii
 ! make the corresponding 1D cut of the PES 
 dd=(th2max-th2min)/dble(npts-1)
 do i=1,npts
   jac(3)=cos((th2min+dble(i-1)*dd)*pii/180.d0)
   call PES(jac,V,NAME1)
   poten=(V)*CONVE
   write(2002,'(F8.2,F20.8)') acos(jac(3))*180.d0/pii,poten
 enddo
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut1D-th2-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(0) = no'
 write(*,*)'(1) = yes'
 write(*,*)'Keep the file?'
 read(*,*)input
 if(input==0)then
   call system('rm cut1D-th2-'//trim(adjustl(charid))//'.*')
 endif

 write(*,*)
 write(*,*)'(1) = Make another graph U vs. theta2'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 228 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=108)input
 if (input==1) goto 2004
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 108 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 228

 2005 continue
 !  U vs. Phi   --------------------------------------------------------------------
 write(*,*)
 write(*,*)'Define fixed coordinates:'
 write(*,'(A16)',ADVANCE='NO')'R = '
 read(*,*)jac(1)
 write(*,'(A16)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=cos(th1*pii/180.d0)
 write(*,'(A16)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=cos(th2*pii/180.d0)
 write(*,*)'------------------------'
 write(*,'(A20)',ADVANCE='NO')'min. value of phi = '
 read(*,*)phimin
 write(*,'(A20)',ADVANCE='NO')'max. value of phi = '
 read(*,*)phimax
 write(*,'(A20)',ADVANCE='NO')'number of points = '
 read(*,*)npts
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut1D-phi-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut1D-phi-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=229)
     enddo
     close(2002)
     cycle
     229 close(2002, status="delete")
     ncount=ncount-1
     goto 230
   else
     goto 230
   endif
 enddo
 230 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut1D-phi-'//trim(adjustl(charid))//'.dat')
 504 format('    Phi           Potential    (R=',F5.1,'; th1=',F5.1,'; th2=',F5.1,')')
 write(2002,504)jac(1),acos(jac(2))*180.d0/pii,acos(jac(3))*180.d0/pii
 ! make the corresponding 1D cut of the PES 
 dd=(phimax-phimin)/dble(npts-1)
 do i=1,npts
   jac(4)=(phimin+dble(i-1)*dd)*pii/180.d0
   call PES(jac,V,NAME1)
   poten=(V)*CONVE
   write(2002,'(F8.2,F20.8)') jac(4)*180.d0/pii,poten
 enddo
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut1D-phi-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 111 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=112)input
 if(input==1)goto 113
 if(input==2)then
   call system('rm cut1D-phi-'//trim(adjustl(charid))//'.*')
   goto 113
 endif
 112 write(*,*)'enter a valid integer: 1 or 2'
 goto 111

 113 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. Phi'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 231 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=109)input
 if (input==1) goto 2005
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 109 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 231


 3000 continue 
 ! ---------------------------------------------------------------------------------
 write(*,*)
 write(*,*)'(1) = R & theta1'
 write(*,*)'(2) = R & theta2'
 write(*,*)'(3) = R & phi'
 write(*,*)'(4) = theta1 & theta2'
 write(*,*)'(5) = theta1 & phi'
 write(*,*)'(6) = theta2 & phi'
 write(*,*)'(7) = Restart'
 write(*,*)'(8) = EXIT'
 3001 continue
 write(*,'(A17)',ADVANCE='YES')'fixed variables?'
 read(*,'(I2)',ADVANCE='YES',ERR=114)input
 if (input==1) goto 3200
 if (input==2) goto 3300
 if (input==3) goto 3400
 if (input==4) goto 3500
 if (input==5) goto 3600
 if (input==6) goto 3700
 if (input==7) then
   call system('clear')
   goto 100
 endif
 if (input==8) goto 666
 114 write(*,*)'enter a valid integer: 1, 2, 3, 4, 5, 6, 7 or 8'
 goto 3001

 3200 continue
 ! ---   U vs. (R,theta1)   --------------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'phi (deg.) = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,'(A36)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=dcos(th2*pii/180.d0)
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta1 (deg.) = '
 read(*,*)th1min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta1 (deg.) = '
 read(*,*)th1max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta1-grid  = '
 read(*,*)npts2
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of R (Ang.) = '
 read(*,*)xrmin
 write(*,'(A36)',ADVANCE='NO')'max. value of R (Ang.) = '
 read(*,*)xrmax
 write(*,'(A36)',ADVANCE='NO')'number of points for R-grid  = '
 read(*,*)npts1
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-R_TH1-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-R_TH1-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3201)
     enddo
     close(2002)
     cycle
     3201 close(2002, status="delete")
     ncount=ncount-1
     goto 3202
   else
     goto 3202
   endif
 enddo
 3202 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-R_TH1-'//trim(adjustl(charid))//'.dat')
 write(2002,3203)phi,dacos(jac(3))*180.d0/pii
 3203 format('     R    Theta1         Potential    (PHI=',F9.3,'; th2=',F5.1,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(xrmax-xrmin)/dble(npts1-1)
 dd2=(th1max-th1min)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(1)=(xrmin+dble(i-1)*dd1)
   do j=1,npts2
     jac(2)=dcos((th1min+dble(j-1)*dd2)*pii/180.d0)
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')jac(1),dacos(jac(2))*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-R_TH1-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3204 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3205)input
 if(input==1)goto 3206
 if(input==2)then
   call system('rm cut2D-R_TH1-'//trim(adjustl(charid))//'.*')
   goto 3206
 endif
 3205 write(*,*)'enter a valid integer: 1 or 2'
 goto 3204

 3206 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (Theta1,Phi)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3207 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3208)input
 if (input==1) goto 3200
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3208 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3207

 3300 continue
 ! ---   U vs. (R,theta2)   --------------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'phi (deg.) = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,'(A36)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=dcos(th1*pii/180.d0)
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta2 (deg.) = '
 read(*,*)th2min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta2 (deg.) = '
 read(*,*)th2max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta2-grid  = '
 read(*,*)npts2
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of R (Ang.) = '
 read(*,*)xrmin
 write(*,'(A36)',ADVANCE='NO')'max. value of R (Ang.) = '
 read(*,*)xrmax
 write(*,'(A36)',ADVANCE='NO')'number of points for R-grid  = '
 read(*,*)npts1
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-R_TH2-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-R_TH2-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3301)
     enddo
     close(2002)
     cycle
     3301 close(2002, status="delete")
     ncount=ncount-1
     goto 3302
   else
     goto 3302
   endif
 enddo
 3302 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-R_TH2-'//trim(adjustl(charid))//'.dat')
 write(2002,3303)phi,dacos(jac(2))*180.d0/pii
 3303 format('     R    Theta2         Potential    (PHI=',F9.3,'; th1=',F5.1,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(xrmax-xrmin)/dble(npts1-1)
 dd2=(th2max-th2min)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(1)=(xrmin+dble(i-1)*dd1)
   do j=1,npts2
     jac(3)=dcos((th2min+dble(j-1)*dd2)*pii/180.d0)
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')jac(1),dacos(jac(3))*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-R_TH2-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3304 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3305)input
 if(input==1)goto 3306
 if(input==2)then
   call system('rm cut2D-R_TH21-'//trim(adjustl(charid))//'.*')
   goto 3306
 endif
 3305 write(*,*)'enter a valid integer: 1 or 2'
 goto 3304

 3306 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (Theta2,Phi)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3307 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3308)input
 if (input==1) goto 3300
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3308 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3307

 3400 continue
 ! ---   U vs. (R,phi)   -----------------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=dcos(th1*pii/180.d0)
 write(*,'(A36)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=dcos(th2*pii/180.d0)
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of phi (deg.) = '
 read(*,*)phimin
 write(*,'(A36)',ADVANCE='NO')'max. value of phi (deg.) = '
 read(*,*)phimax
 write(*,'(A36)',ADVANCE='NO')'number of points for phi-grid  = '
 read(*,*)npts2
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of R (Ang.) = '
 read(*,*)xrmin
 write(*,'(A36)',ADVANCE='NO')'max. value of R (Ang.) = '
 read(*,*)xrmax
 write(*,'(A36)',ADVANCE='NO')'number of points for R-grid  = '
 read(*,*)npts1
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-R_PHI-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-R_PHI-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3401)
     enddo
     close(2002)
     cycle
     3401 close(2002, status="delete")
     ncount=ncount-1
     goto 3402
   else
     goto 3402
   endif
 enddo
 3402 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-R_PHI-'//trim(adjustl(charid))//'.dat')
 write(2002,3403)th1,th2
 3403 format('     R      Phi          Potential    (TH1=',F9.3,'; TH2=',F9.3,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(xrmax-xrmin)/dble(npts1-1)
 dd2=(phimax-phimin)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(1)=(xrmin+dble(i-1)*dd1)
   do j=1,npts2
     jac(4)=(phimin+dble(j-1)*dd2)*pii/180.d0
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')jac(1),jac(4)*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-R_PHI-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3404 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3405)input
 if(input==1)goto 3406
 if(input==2)then
   call system('rm cut2D-R_PHI-'//trim(adjustl(charid))//'.*')
   goto 3406
 endif
 3405 write(*,*)'enter a valid integer: 1 or 2'
 goto 3404

 3406 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (R,Phi)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3407 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3408)input
 if (input==1) goto 3400
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3408 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3407

 3500 continue
 ! ---   U vs. (theta1,theta2)   ---------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'R (Ang.) = '
 read(*,*)R
 jac(1)=R
 write(*,'(A36)',ADVANCE='NO')'phi (deg.) = '
 read(*,*)phi
 jac(4)=phi*pii/180.d0
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta1 (deg.) = '
 read(*,*)th1min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta1 (deg.) = '
 read(*,*)th1max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta1-grid  = '
 read(*,*)npts1
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta2 (deg.) = '
 read(*,*)th2min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta2 (deg.) = '
 read(*,*)th2max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta2-grid  = '
 read(*,*)npts2
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-TH1_TH2-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-TH1_TH2-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3501)
     enddo
     close(2002)
     cycle
     3501 close(2002, status="delete")
     ncount=ncount-1
     goto 3502
   else
     goto 3502
   endif
 enddo
 3502 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-TH1_TH2-'//trim(adjustl(charid))//'.dat')
 write(2002,3503)R,phi
 3503 format('   Theta1  Theta2          Potential    (R=',F7.3,'; phi=',F5.1,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(th1max-th1min)/dble(npts1-1)
 dd2=(th2max-th2min)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(2)=dcos((th1min+dble(i-1)*dd1)*pii/180.d0)
   do j=1,npts2
     jac(3)=dcos((th2min+dble(j-1)*dd2)*pii/180.d0)
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')dacos(jac(2))*180.d0/pii,dacos(jac(3))*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-TH1_TH2-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3504 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3505)input
 if(input==1)goto 3506
 if(input==2)then
   call system('rm cut2D-TH1_TH2-'//trim(adjustl(charid))//'.*')
   goto 3506
 endif
 3505 write(*,*)'enter a valid integer: 1 or 2'
 goto 3504

 3506 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (Theta1,Theta2)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3507 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3508)input
 if (input==1) goto 3500
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3508 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3507

 3600 continue
 ! ---   U vs. (theta1,phi)   ------------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'R (Ang.) = '
 read(*,*)R
 jac(1)=R
 write(*,'(A36)',ADVANCE='NO')'theta2 (deg.) = '
 read(*,*)th2
 jac(3)=dcos(th2*pii/180.d0)
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta1 (deg.) = '
 read(*,*)th1min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta1 (deg.) = '
 read(*,*)th1max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta1-grid  = '
 read(*,*)npts2
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of phi (deg.) = '
 read(*,*)phimin
 write(*,'(A36)',ADVANCE='NO')'max. value of phi (deg.) = '
 read(*,*)phimax
 write(*,'(A36)',ADVANCE='NO')'number of points for phi-grid  = '
 read(*,*)npts1
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-TH1_PHI-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-TH1_PHI-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3601)
     enddo
     close(2002)
     cycle
     3601 close(2002, status="delete")
     ncount=ncount-1
     goto 3602
   else
     goto 3602
   endif
 enddo
 3602 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-TH1_PHI-'//trim(adjustl(charid))//'.dat')
 write(2002,3603)R,dacos(jac(3))*180.d0/pii
 3603 format('   Theta1   Phi           Potential    (R=',F7.3,'; th2=',F5.1,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(phimax-phimin)/dble(npts1-1)
 dd2=(th1max-th1min)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(4)=(phimin+dble(i-1)*dd1)*pii/180.d0
   do j=1,npts2
     jac(2)=dcos((th1min+dble(j-1)*dd2)*pii/180.d0)
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')dacos(jac(2))*180.d0/pii,jac(4)*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-TH1_PHI-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3604 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3605)input
 if(input==1)goto 3606
 if(input==2)then
   call system('rm cut2D-TH1_PHI-'//trim(adjustl(charid))//'.*')
   goto 3606
 endif
 3605 write(*,*)'enter a valid integer: 1 or 2'
 goto 3604

 3606 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (Theta1,Phi)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3607 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3608)input
 if (input==1) goto 3600
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3608 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3607

 3700 continue
 ! ---   U vs. (theta2,phi)   ------------------------------------------------------
 write(*,*)
 write(*,'(A26)')'Define fixed coordinates:'
 write(*,'(A36)',ADVANCE='NO')'R (Ang.) = '
 read(*,*)R
 jac(1)=R
 write(*,'(A36)',ADVANCE='NO')'theta1 (deg.) = '
 read(*,*)th1
 jac(2)=dcos(th1*pii/180.d0)
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of theta2 (deg.) = '
 read(*,*)th2min
 write(*,'(A36)',ADVANCE='NO')'max. value of theta2 (deg.) = '
 read(*,*)th2max
 write(*,'(A36)',ADVANCE='NO')'number of points for theta2-grid  = '
 read(*,*)npts2
 write(*,*)'-----------------------------------'
 write(*,'(A36)',ADVANCE='NO')'min. value of phi (deg.) = '
 read(*,*)phimin
 write(*,'(A36)',ADVANCE='NO')'max. value of phi (deg.) = '
 read(*,*)phimax
 write(*,'(A36)',ADVANCE='NO')'number of points for phi-grid  = '
 read(*,*)npts1
 ncount=1
 do j=1,1000! check if previous files exist
   write(charid,'(I4)')j
   inquire(file='cut2D-TH2_PHI-'//trim(adjustl(charid))//'.dat',exist=logica1)
   if(logica1)then
     ncount=ncount+1
     open(unit=2002,file='cut2D-TH2_PHI-'//trim(adjustl(charid))//'.dat',action='read',status='old')
     do i=1,3
       read(2002,*,end=3701)
     enddo
     close(2002)
     cycle
     3701 close(2002, status="delete")
     ncount=ncount-1
     goto 3702
   else
     goto 3702
   endif
 enddo
 3702 write(charid,'(I4)')ncount
 ! create output file
 open(unit=2002,file='cut2D-TH2_PHI-'//trim(adjustl(charid))//'.dat')
 write(2002,3703)R,dacos(jac(2))*180.d0/pii
 3703 format('   Theta2   Phi           Potential    (R=',F7.3,'; th1=',F5.1,')')
 ! make the corresponding 2D cut of the PES 
 dd1=(phimax-phimin)/dble(npts1-1)
 dd2=(th2max-th2min)/dble(npts2-1)
 Emin=200d0
 Emax=-1d6
 do i=1,npts1
   jac(4)=(phimin+dble(i-1)*dd1)*pii/180.d0
   do j=1,npts2
     jac(3)=dcos((th2min+dble(j-1)*dd2)*pii/180.d0)
     call PES(jac,V,NAME1)
     poten=(V)*CONVE
     if (poten>Emax) Emax=poten
     if (poten<Emin) Emin=poten
     write(2002,'(F8.2,F8.2,F20.8)')dacos(jac(3))*180.d0/pii,jac(4)*180.d0/pii,poten
   enddo
   write(2002,*)
 enddo
 dE=Emax-Emin
 close(2002)
 write(*,*)'----------------------'
 write(*,*)'Output data saved in the file: "cut2D-TH2_PHI-'//trim(adjustl(charid))//'.dat".'
 write(*,*)'(1) = yes'
 write(*,*)'(2) = no'
 3704 continue
 write(*,'(A16)',ADVANCE='YES')'Keep the file?'
 read(*,'(I2)',ADVANCE='YES',ERR=3705)input
 if(input==1)goto 3706
 if(input==2)then
   call system('rm cut2D-TH2_PHI-'//trim(adjustl(charid))//'.*')
   goto 3706
 endif
 3705 write(*,*)'enter a valid integer: 1 or 2'
 goto 3704

 3706 continue
 write(*,*)
 write(*,*)'(1) = Make another graph U vs. (Theta2,Phi)'
 write(*,*)'(2) = Restart'
 write(*,*)'(3) = EXIT'
 3707 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=3708)input
 if (input==1) goto 3700
 if (input==2) then
   call system('clear')
   goto 100
 endif
 if (input==3) goto 666
 3708 write(*,*)'enter a valid integer: 1, 2 or 3'
 goto 3707

! call system('clear')
! write(*,*)'To be done...'
! goto 100

 4000 continue
 ! ---------------------------------------------------------------------------------
 ! *  R-optimized 2D cut of the PES: "U vs. extended angles" (planar geometries)  **

 gnup1=0
 if(gnup==1)then! if gnuplot available...
   write(*,*)
   write(*,*)'Details for the contour plots:'
   write(*,*)'( 0 = let gnuplot choose )'
   write(*,'(A27)',ADVANCE='NO')'min. energy level = '
   read(*,*)Emin
   if(Emin==0.d0)then
      gnup1=1
      goto 4001
   endif
   write(*,'(A27)',ADVANCE='NO')'max. energy level = '
   read(*,*)Emax
   jac(2)=cos(th1*pii/180.d0)
   write(*,'(A27)',ADVANCE='NO')'energy gap between levels = '
   read(*,*)dE
 endif
 4001 write(*,*)
 write(*,*)'Please, be patient,'
 write(*,*)'this could take some time...'
 
 EGrid2=0d0
 RGrid2=0d0

 nitv=50
 ! 
 do i=1,nitv*2+1
   coord1(i)=-1.0d0*pii+dble(i-1)*(2d0*pii)/dble(nitv*2)
 enddo
 do i=1,nitv*2+1
   coord2(i)=-1.0d0*pii+dble(i-1)*(2d0*pii)/dble(nitv*2)
 enddo

 do i=nitv+1,2*nitv+1! LOOP on extended-theta1 (from zero to pi)
   jac(2)=cos(coord1(i))

   do j=1,nitv! loop1, on extended-theta2 (from -pi to zero)
     jac(4)=0d0! phi=0
     jac(3)=cos(coord2(2*nitv+2-j))
     tampon=400d0
     Rtamp=50d0
     ! first scan
     xdR=(rmax(1)-rmin(1))/dble(30)  !! cambiar a 20 !!
     do k=1,30 
       jac(1)=dble(k-1)*xdR+rmin(1)
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! second scan
     dR=xdR
     xdR=2.0d0*xdR/dble(10)
     do k=1,11 
       jac(1)=dble(k-1)*xdR+Rtamp-dR
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! third scan
     dR=xdR
     xdR=2.0d0*xdR/dble(10)
     do k=1,11 
       jac(1)=dble(k-1)*xdR+Rtamp-dR
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! last scan
     dR=xdR
     if(2.0d0*dR>0.02)then
       xdR=0.01d0
       do k=1,int((2.0d0*dR)/xdR)+1 
         jac(1)=dble(k-1)*xdR+Rtamp-dR
         call PES(jac,V,NAME1)
         tampon1=(V)*CONVE
         if(tampon1<tampon)then
           tampon=tampon1
           Rtamp=jac(1)
         endif
       enddo
     endif
     if(tampon>0d0)then
       tampon=0d0
       Rtamp=rmax(1)
     endif
     EGrid2(i,j)=tampon
     RGrid2(i,j)=Rtamp
    !write(*,*)i,j,EGrid2(i,j)
   enddo
!write(*,*)
   do j=nitv+1,2*nitv+1! loop2, on extended-theta2 (from zero to pi)
     jac(4)=pii! phi=180
     jac(3)=cos(coord2(j))
     tampon=400d0
     Rtamp=50d0
     ! first scan
     xdR=(rmax(1)-rmin(1))/dble(20)
     do k=1,20 
       jac(1)=dble(k-1)*xdR+rmin(1)
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! second scan
     dR=xdR
     xdR=2.0d0*xdR/dble(10)
     do k=1,11 
       jac(1)=dble(k-1)*xdR+Rtamp-dR
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! third scan
     dR=xdR
     xdR=2.0d0*xdR/dble(10)
     do k=1,11 
       jac(1)=dble(k-1)*xdR+Rtamp-dR
       call PES(jac,V,NAME1)
       tampon1=(V)*CONVE
       if(tampon1<tampon)then
         tampon=tampon1
         Rtamp=jac(1)
       endif
     enddo
     ! last scan
     dR=xdR
     if(2.0d0*dR>0.02)then
       xdR=0.01d0
       do k=1,int((2.0d0*dR)/xdR)+1 
         jac(1)=dble(k-1)*xdR+Rtamp-dR
         call PES(jac,V,NAME1)
         tampon1=(V)*CONVE
         if(tampon1<tampon)then
           tampon=tampon1
           Rtamp=jac(1)
         endif
       enddo
     endif
     if(tampon>0d0)then
       tampon=0d0
       Rtamp=rmax(1)
     endif
     EGrid2(i,j)=tampon
     RGrid2(i,j)=Rtamp
    !write(*,*)i,j,EGrid2(i,j)
   enddo
   if(i==51)write(*,*)
   if(i==51)write(*,*)'completed:'
   if(i==51)write(*,*)' 2 %'
   if(i==55)write(*,*)'10 %'
   if(i==63)write(*,*)'25 %'
   if(i==75)write(*,*)'50 %'
   if(i==88)write(*,*)'75 %'
 enddo
 !write(*,*)
 do i=1,nitv
  do j=1,2*nitv+1
    EGrid2(i,j)=EGrid2(nitv*2+2-i,nitv*2+2-j)
    RGrid2(i,j)=RGrid2(nitv*2+2-i,nitv*2+2-j)
    !write(*,*)i,j,EGrid2(i,j)
  enddo
 enddo

 ! create output files:
 open(unit=2002,file='cut2D-ext-ang.dat')
 write(2002,*)"    Th1      Th2         U_min      R_min "
 ! print the corresponding 2D cut of the PES 
 do i=1,2*nitv+1
  do j=1,2*nitv+1
!    write(2002,'(2F8.3,X,F20.4)') coord1(i),coord2(j),EGrid2(i,j),RGrid2(i,j)
!    write(2002,'(2F9.2,F16.6,F8.2)') coord1(i)*180.d0/pii,coord2(j)*180.d0/pii-180.d0,EGrid2(i,j),RGrid2(i,j)
    write(2002,'(2F9.2,F16.6,F8.2)') coord1(i)*180.d0/pii,coord2(j)*180.d0/pii,EGrid2(i,j),RGrid2(i,j)
  enddo
  write(2002,*)
 enddo
 close(2002)
 open(unit=2003,file='cut2D-ext-ang-Ematrix.dat')
 ! print the corresponding 2D cut of the PES 
 do i=1,2*nitv+1
   write(2003,'(101F14.4)')EGrid2(i,:)
 enddo
 close(2003)

 write(*,*)'Done!'
 write(*,*)'Output data saved in the files:'
 write(*,*)'"cut2D-ext-ang.dat" and "cut2D-ext-ang-Ematrix.dat".'


!write(*,*)
!write(*,*)
!write(*,*)EGrid2(1:5,1:5)


 write(*,*)
 write(*,*)'(1) = Restart'
 write(*,*)'(2) = EXIT'
 232 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=110)input
 if (input==1) then
   call system('clear')
   goto 100
 endif
 if (input==2) goto 666
 110 write(*,*)'enter a valid integer: 1 or 2'
 goto 232


 5000 continue
 ! ---------------------------------------------------------------------------------
 ! **                        Prepare "AUTOSURF_XGW.dat"                           **
 ! ---------------------------------------------------------------------------------
 open(unit=200,file='input-AUTOSURF-PES.dat',status='old',action='read')
 236 read(200,'(A100)')line
 ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
 if (ncont1==0) goto 236
 do i=1,6+natom
   read(200,*)
 enddo
 do i=1,natom
   read(200,*) mass(i)! masses of all atoms
 enddo  
 close(200)

 goto 235! TEST... to use, comment this line
 deallocate(ref1,mass,ref2)
 open(unit=200,file='input-AUTOSURF-PES-test.dat',status='old',action='read')
 read(200,*) natom1! number of atoms in frag1
 read(200,*) natom2! number of atoms in frag2
 natom=natom1+natom2
 allocate(ref1(3*natom1),mass(natom),ref2(3*natom2))
 do i=1,natom
   read(200,*) mass(i)! masses of all atoms
 enddo  
 do i=1,3*natom1
   read(200,*) ref1(i)! Cartesian positions for fragment 1 atoms
 enddo
 do i=1,3*natom2
   read(200,*) ref2(i)! Cartesian positions for fragment 1 atoms
 enddo
 close(200)
 235 continue

 ! open the "AUTOSURF_XGW.dat" file
 NAME2='AUTOSURF_XGW.dat'
 open(unit=200,file=NAME2)

 x1=0.d0
 do i=1,natom1
   x1=x1+mass(i)! mass of fragment 1
 enddo  
 xmas1=x1
 write(6,*)'mass_fragment1 = ',xmas1
 write(200,'(F20.11)')x1
 x1=0.d0
 do i=1,natom2
   x1=x1+mass(natom1+i)! mass of fragment 2
 enddo  
 write(200,'(F20.11)')x1
 xmas2=x1
 write(6,*)'mass_fragment2 = ',xmas2

 Emin=200.d2! find ab initio geometry with lower energy &
 xass=200.d2! energy of the ab initio geometry with largest R
 x1=0.d0
 do i=1,count3       
   if(mod(i,(2*symparts))==1)then
     jac=coords(i,:) 
     V=func_actual(jac)
     if (V<Emin) then
       Emin=V
       xi=coords(i,:) 
     endif
     if (x1<coords(i,1)) then
       xass=V
       x1=coords(i,1)
     endif
   endif
 enddo
 write(200,'(F20.11)')xi(1)
 write(200,'(F20.11)')xi(2)
 write(200,'(F20.11)')xi(3)
 write(200,'(F20.11)')xi(4)

 !!   FRAGMENT 1   !!
 write(6,*)
 write(6,*)'*** FRAGMENT 1 ****'
 allocate(xcor(natom1))
 ! check if C.M.(x) is zero 
 do i=1,natom1
  xcor(i)=ref1(1+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom1
  x1=x1+mass(i)*xcor(i)
 enddo
 if((x1/xmas1)>1.d-6)write(6,*)'error 1'
 ! check if C.M.(y) is zero 
 do i=1,natom1
  xcor(i)=ref1(2+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom1
  x1=x1+mass(i)*xcor(i)
 enddo
 if((x1/xmas1)>1.d-6)write(6,*)'error 1'
 ! check if C.M.(z) is zero 
 do i=1,natom1
  xcor(i)=ref1(3+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom1
  x1=x1+mass(i)*xcor(i)
 enddo
 x1=x1/xmas1
 if (x1<1.d-6) then
   write(6,*)'C.M. of fragment 1 is already at the origin...'
 else! shift C.M.(z)
   xcor=xcor-x1
   write(6,*)'dZ = ',x1
 endif
 ! Inertia moment
 x1=0.d0
 do i=1,natom1
  x1=mass(i)*(xcor(i)**2)+x1
 enddo
 write(6,*)'I = ',x1
 ! Rotation constant (in wave numbers)
 x1=16.857629126507d0/x1
 write(6,*)'B = ',x1
 write(200,'(F20.11)')x1

 !!   FRAGMENT 2   !!
 write(6,*)
 write(6,*)'*** FRAGMENT 2 ****'
 deallocate(xcor)
 allocate(xcor(natom2))
 ! check if C.M.(x) is zero 
 do i=1,natom2
  xcor(i)=ref2(1+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom2
  x1=x1+mass(natom1+i)*xcor(i)
 enddo
 if((x1/xmas2)>1.d-6)write(6,*)'error 1'
 ! check if C.M.(y) is zero 
 do i=1,natom2
  xcor(i)=ref2(2+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom2
  x1=x1+mass(natom1+i)*xcor(i)
 enddo
 if((x1/xmas2)>1.d-6)write(6,*)'error 1'
 ! check if C.M.(z) is zero 
 do i=1,natom2
  xcor(i)=ref2(3+3*(i-1))
 enddo
 x1=0.d0
 do i=1,natom2
  x1=x1+mass(natom1+i)*xcor(i)
 enddo
 x1=x1/xmas2
 if (x1<1.d-6) then
   write(6,*)'C.M. of fragment 1 is already at the origin...'
 else! shift C.M.(z)
   xcor=xcor-x1
   write(6,*)'dZ = ',x1
 endif
 ! Inertia moment
 x1=0.d0
 do i=1,natom2
  x1=mass(natom1+i)*(xcor(i)**2)+x1
 enddo
 write(6,*)'I = ',x1
 ! Rotation constant (in wave numbers)
 x1=16.857629126507d0/x1
 write(6,*)'B = ',x1
 write(6,*)
 write(200,'(F20.11)')x1

 ! Asymptote energy
 write(200,'(F20.11)')xass

 ! name of the PES-file
 write(200,*)trim(adjustl(NAME1))

 write(200,'(F20.11,A25)')rmin(1)*ang2bohr,'! min. R (in Bohr)'
 write(200,'(F20.11,A25)')rmax(1)*ang2bohr,'! max. R (in Bohr)'

 close(200)

 write(*,*)
 write(*,*)'(1) = Restart'
 write(*,*)'(2) = EXIT'
 5001 continue
 write(*,'(A12)',ADVANCE='YES')'What to do?'
 read(*,'(I2)',ADVANCE='YES',ERR=5002)input
 if (input==1) then
   call system('clear')
   goto 100
 endif
 if (input==2) goto 666
 5002 write(*,*)'enter a valid integer: 1 or 2'
 goto 5001




666 continue
!call system('clear')
stop

end program AUTOSURF_PLOT

