!*******************************  A U T O S U R F  *********************************!
!===================================================================================!
!-----------------------------------------------------------------------------------!
!-                                                                                 -!
!-        AUTOSURF Package: A set of programs for the automated construction       -!
!-              of Potential Energy Surfaces on van der Waals systems              -!
!-                                                                                 -!
!-----------------------------------------------------------------------------------!
!===================================================================================!
!***********************************************************************************!
!-     "AUTOSURF-ABI_rigid4D": PROGRAM for the automated construction of a         -!
!-       4D-PES on a vdW system composed of two (rigid) linear fragments.          -!
!-----------------------------------------------------------------------------------!
!      Input files: "input-AUTOSURF-ABI.dat.dat", "molproX.abi" (X = 1, 2 & 10)     !
!      Dependencies: "FUNCTIONS.f90", "MODULES.f90", "diag.f", "SUBROUTINES.f90"    !
!      Version: 1.2
!***********************************************************************************!

PROGRAM AUTOSURF_ABI_rigid4D

use nr
use nrtype
use nrutil
use dynamic_parameters
use mpi
!-----------------------------------------------------------------------------------
implicit none
 integer :: i,j,k,k2,l,ip,jp,kp,ipp,jpp,kpp,lpp,skip,jj,numpoints,numpoints_low,  &
             holepatch,term,wellfocus,nloop
 integer :: num_spec_points,count3_old,count3_new
 integer :: loop,add,kk,iter,reiter,count,restart,numprocs,numadded,numadded_act, &
             rcc,ierr,num_err_points,num_ab_err_points,new(2),code_flag,code_flag2
 integer :: pass_check,num_err_points1
 real*8,allocatable :: ind(:),ind_seed(:),ind3(:),ind5(:),ind7(:),ind9(:),       &
             ind11(:),ind13(:),seed(:,:),seed_pot(:),seed_grad(:,:),seed_low(:,:), &
             seed_pot_low(:),seed_grad_low(:,:),test_points(:,:),test_points2(:,:),&
             error_points(:,:),order_diff(:),added_points(:,:),rms_error(:),       &
             mean_error(:),dcart_stored(:,:),error_points_real(:,:),hole(:,:),     &
             hole_val(:),low_pot(:),internal(:,:),grad_int(:,:)
 integer,allocatable :: ind2(:),ind2_seed(:),ind4(:),ind6(:),ind8(:),ind10(:),  &
             ind12(:),ind14(:),indhole(:)
 real*8 :: xi(4),tampon,pii
 real*8 :: temp,temp3,temp4(4),temp2(5),ran_vec(4),jac3(4),rann,bias,dynmin,dynmax
 real*8 :: somme,somme2,increment
 real*8 :: abs_diff,range,valeur,h2wn,dist,dist2,gtol,sqrdiff,local_d,atobohr,    &
            Min_E,Max_R,E_asym,Max_E3,bot_seed,low_high_fac,seed_cut
 character(len=300) :: sys_label,line,f602,xchar3,f603
 character(len=10) :: bdate1,bdate2,bdate3
 character(len=5) :: charmyid
 character(len=8) :: xchar,xchar1
 character(len=15),dimension(4) :: charunitsE
 character(len=15),dimension(2) :: charunitsD
 character(len=40) :: NAME1
 character(len=4) :: charid
 integer :: nline,summyid,ncont,ncont1,ncont2,nid,numpoints_low_count,old_numprocs
 integer :: xpass,numpoints_high_count,n_test,xunitE,focus,only_seed_grid,status
 integer :: tot_abinitio,ntest,auto_mod,focus_onR,focus_onLR,nabi,xcontrol,xnabi
 integer,dimension(1) :: xseed,xseed1
 integer,dimension(8) :: date_time
 real*8 :: xpot,CONVE,CONVD,xdR,minR,maxR,xtemp,x1,x2,x3,x4,x5
 real*8 :: xGlob_min(5),xGlob_min1(5),xminLR
 real*8,allocatable :: xgrad(:),xR(:),xMAPot(:,:)
 real*8,allocatable :: xrms_R1(:,:),xrms_R2(:,:),xrms_R3(:,:),xrms_R4(:,:)
 real*8,allocatable :: xR_count1(:,:),xR_count2(:,:),xR_count3(:,:),xR_count4(:,:)
 logical :: logica1, logica2
 real*8,parameter :: hart2wn=219474.6313702d0, hart2meV=27211.38602d0
 real*8,parameter :: ang2bohr=1.0d0/0.529177249d0, bohr2ang=0.529177249d0
 real*8,parameter :: hart2kcl=4.359744650D-18/4184*6.022140857D23
 integer,parameter :: xunitD=2, xnum_intervals_R=20
 data charunitsE/'hartree','kcal/mol','wave numbers','meV'/
 data charunitsD/'Bohr','Angstroms'/
 real*8, dimension(400000) :: xx,eee
 real*8 :: th1,th2,phi,xrr,th1max,th2max,phimax,xrmax,th1min,th2min,phimin,xrmin,dd
 integer,parameter :: send_data_tag = 2001, return_data_tag = 2002
 integer,parameter :: midbond = 0

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
!***********************************************************************************
!***********************************************************************************

call date_and_time(bdate1,bdate2,bdate3,date_time)
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)!
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)!
write(*,*) 'Process ', myid, ' of ', numprocs, ' is alive'
inquire(file='input-AUTOSURF-ABI.dat',exist=logica1)
if(.not.logica1)then
 if(myid==0)write(*,*)
 if(myid==0)write(*,*)'ERROR: The file: input-AUTOSURF-ABI.dat does not exist !! '
 stop 
endif
open(unit=myid+1000,file='input-AUTOSURF-ABI.dat',status='old',action='read')
111 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'GENERAL INFORMATION:')
if(ncont1==0)goto 111
read(myid+1000,*)
read(myid+1000,*)sys_label !
nline=scan(sys_label,' ')-1
if (nline>15) nline=15
!read(myid+1000,*)restart
!read(myid+1000,*)xseed
read(myid+1000,*)code_flag2
read(myid+1000,*)ab_flag2
read(myid+1000,*)xnabi
close(myid+1000)
restart=0  !
xseed=0    !
xunitE=2 ! 
nabi=1
call sobseq(ran_vec,4)
do i=1,11
  call sobseq(ran_vec)
enddo
!xseed1=10
!CALL RANDOM_SEED(PUT=xseed1) 
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0)then
  open(unit=110,file='output.dat')
  ncont1=1
  do j=1,1000!--
    write(charid,'(I4)')j
    inquire(file='outputABI_'//sys_label(1:nline)//'_'//trim(adjustl(charid))//'.dat',exist=logica1)
    if(logica1)then
      ncont1=ncont1+1
      open(unit=2002,file='outputABI_'//sys_label(1:nline)//'_'//trim(adjustl(charid))//'.dat',action='read',status='old')
      do i=1,3
        read(2002,*,end=220)
      enddo
      close(2002)
      cycle
      220 close(2002,status="delete")
      ncont1=ncont1-1
      goto 221
    else
      goto 221
    endif
  enddo
  221 write(charid,'(I4)')ncont1
  open(unit=100,file='outputABI_'//sys_label(1:nline)//'_'//trim(adjustl(charid))//'.dat',action='write')
  write(100,*)'System:  '//sys_label(1:nline)//' '
  write(100,*)
  write(100,'(A74)')' *************************************************************************'
  write(100,'(A74)')' *         The automated calculation of ab initio data has started       *'
  write(100,'(A74)')' *         System date (dd/mm/yyyy): '//bdate1(7:8)//'/'//bdate1(5:6)//'/' &
                     //bdate1(1:4)//'                          *'
  write(100,'(A74)')' *         System time (hh:min:sec): '//bdate2(1:2)//':'//bdate2(3:4)//':' &
                     //bdate2(5:6)//'                            *'
  write(100,'(A74)')' *************************************************************************'
  write(100,*)
  write(100,*)
    do j=1,2000
      write(charmyid,'(I5)')j-1
      inquire(file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//'.dat', &
               exist=logica1)
      if(logica1)then  
       open(unit=200,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//  &
                            '.dat')
       close(200, status="delete")
      endif
    enddo
  write(100,'(A74)')' -------------------------------------------------------------------------'
  write(100,*)'                        STARTING MPI ENVIROMENT                          '
  write(100,*)
  write(100,*)'Energies are in kcal/mol, distances are in Angstroms and angles in degrees.'
  write(100,'(A30,I6)')' Total number of processes:   ',numprocs
endif

!-----------------------------------------------------------------------------------
open(unit=myid+1000,file='input-AUTOSURF-ABI.dat',status='old',action='read')
112 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
if (ncont1==0) goto 112
read(myid+1000,*)
read(myid+1000,*) exch  !
read(myid+1000,*) flip1 ! 
read(myid+1000,*) flip2 !
read(myid+1000,*) natom1! 
read(myid+1000,*) natom2! 
natom=natom1+natom2 ! 
nbdist=natom*(natom-1)/2! 
allocate(ref1(3*natom1),ref2(3*natom2))
allocate(symb(natom),mass(natom))
do i=1,natom
 read(myid+1000,*) symb(i)!
enddo
do i=1,natom
 read(myid+1000,*) mass(i)!
enddo  
do i=1,3*natom1
 read(myid+1000,*) ref1(i)! 
enddo
do i=1,3*natom2
 read(myid+1000,*) ref2(i)! 
enddo
close(myid+1000)
nloop=100! 
focus=0! 
rmin(1)=0.5d0! 
rmax(1)=100.d0! 
focus_onR=0!
minR=rmin(1)
maxR=rmax(1)
numpoints=20000! 
ab_flag=1! 
x1=dble(numprocs)
numadded=int(x1)
x1=10
x1=x1*dble(numprocs)
num_ab_err_points=int(x1)
maxpoints=200000
order_1=6! 
order_2=6! 
order_3=6! 
order_4=6!
dist_tol=1.25d0
num_err_points=150000
num_err_points=num_err_points/numprocs! 
low_grid=2
numpoints_low=11
subzero=0
pii=dacos(-1d0)
term=0! 
xpass=0!
alpha=-1d0! 
wellfocus=0!  
symparts=((exch+1)*(flip1+1)*(flip2+1))! 
num_err_points1=25000/numprocs! 
!num_err_points1=2000! 
if (num_err_points<2000) num_err_points1=num_err_points
only_seed_grid=0! 
if (nloop==0) then
  only_seed_grid=1
  nloop=2
endif
W_a=0.432d0*3/max(natom1,natom2)! 
epss=1d-14
count3=numpoints
numpoints_low_count=0! 
numpoints_high_count=0! 
rmax(2)=0.9999d0
if(flip1>0)then
  rmin(2)=0.0001d0
else
  rmin(2)=-0.9999d0
endif
rmax(3)=0.9999d0
if(flip2>0)then
  rmin(3)=0.0001d0
else
  rmin(3)=-0.9999d0
endif
rmax(4)=3.14059265d0
rmin(4)=0.001d0!  
allocate(xR(xnum_intervals_R+1))
xdR=(rmax(1)-rmin(1))/dble(xnum_intervals_R)
do j=1,xnum_intervals_R+1
 xR(j)=rmin(1)+(j-1)*xdR
enddo
if(xunitE==1) CONVE=1.d0
if(xunitE==2) CONVE=hart2kcl
if(xunitE==3) CONVE=hart2wn
if(xunitE==4) CONVE=hart2meV
if(xunitD==1) CONVD=1.d0
if(xunitD==2) CONVD=bohr2ang
ugrad=CONVE/CONVD
order_1_min=3   
order_2_min=3
order_3_min=3
order_4_min=3
call basis_size(4,order_1,order_2,order_3,order_4,basis_1)  
call basis_size(4,order_1-1,order_2-1,order_3-1,order_4-1,basis_2) 
call basis_size(4,order_1_min,order_2_min,order_3_min,order_4_min,basis_3)  
if(ab_flag==2)then
  zz=6
else
  zz=4
endif
if(ab_flag2==2)then
  zz_low=6
else
  zz_low=4
endif
zz4=20 
allocate(internal(symparts,4),grad_int(symparts,4))
allocate(jac(4),jac2(4),cart(3*(natom)))
allocate(seed(maxpoints,4),seed_pot(maxpoints),coords(2*symparts*maxpoints,4),pot(2*symparts*maxpoints))
allocate(order_diff(400000),b2_minimal(basis_3,2*symparts*maxpoints))
allocate(dcart(3*(natom)),ind13(num_ab_err_points),ind14(num_ab_err_points),bdist(nbdist))
allocate(d(2*symparts*maxpoints),b2(basis_1,2*symparts*maxpoints),b2_lower(basis_2,2*symparts*maxpoints))
allocate(rms_error(numprocs),mean_error(numprocs),error_points(num_err_points*numprocs,4+1))
allocate(error_points_real(num_ab_err_points,4+1),test_points(400000,4),test_points2(400000,4))    !!!!
allocate(hole(numprocs,4),hole_val(numprocs),indhole(numprocs),ind6(num_err_points*numprocs))
allocate(xrms_R1(numprocs,xnum_intervals_R+1),xrms_R2(numprocs,xnum_intervals_R+1))
allocate(xR_count1(numprocs,xnum_intervals_R+1),xR_count2(numprocs,xnum_intervals_R+1))
allocate(xrms_R3(numprocs,xnum_intervals_R+1),xrms_R4(numprocs,xnum_intervals_R+1))
allocate(xR_count3(numprocs,xnum_intervals_R+1),xR_count4(numprocs,xnum_intervals_R+1))
allocate(xMAPot(numprocs,xnum_intervals_R+1))
allocate(seed_grad(maxpoints,3*(natom)),grad(2*symparts*maxpoints,4),xgrad(3*natom))
allocate(low_pot(maxpoints),seed_low(maxpoints,4),seed_pot_low(maxpoints))
allocate(seed_grad_low(maxpoints,3*(natom)))
dcart=0d0
601 FORMAT(I10,4f20.15,f20.8)              ! 
write(f602,'( "(I10,4f20.14,f20.8,",I3,"f20.8)" )')natom*3 !!
write(f603,'( "(I10,f20.14,f20.8,",I3,"f20.8)" )')natom*3 !! 
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0)then
  write(100,*)'Initialization completed.'
  write(100,*)
endif
 14 CONTINUE
 IF (restart==1) THEN
   inquire(file='AbINITIO_low.dat',exist=logica1)
   if(.not.logica1)then! 
     if(myid==0)then
       open(unit=334,file='AbINITIO_low.dat')
       do j=1,2000! 
         write(charmyid,'(I5)')j-1
         inquire(file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//&
                        '.dat',exist=logica2)
         if(logica2)then  
           open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//&
                                '.dat', action='read',status='old')
           do i=1,maxpoints
             if(ab_flag2==2)then
              read(230,*,end=15)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
              write(334,f602)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
             else
              read(230,*,end=15)ncont,seed_low(i,:),low_pot(i)
              write(334,601)ncont,seed_low(i,:),low_pot(i)
             endif
           enddo
           15 close(230, status="delete")
         endif
       enddo
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call MPI_BCAST(numpoints_low_count,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if(numpoints_low_count==0)then
       restart=0
       if(myid==0)write(100,*)'RESTART POINT corrupted.'
       if(myid==0)write(100,*)'Restarting calculation from scratch... '
       goto 14
     endif
     do i=1,numpoints_low_count
       call MPI_BCAST(seed_low(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(low_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       if(ab_flag2==2)then
        call MPI_BCAST(seed_grad_low(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       endif
     enddo
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if(numpoints_low_count>=numpoints_low)then! 
       if(myid==0)close(334)
       if(myid==0)write(100,*)'All low-level ab initio energies have been computed.'
       if(myid==0)write(100,*)'RESTART POINT successfully created!'
       if(myid==0)write(100,*)'Ab initio data saved in output file: "AbINITIO_low.dat".'
       xpass=1
       if(numpoints_low_count>numpoints_low)then
         if(myid==0)write(110,*)'Number of points increased from ',numpoints_low,' to ',numpoints_low_count
         numpoints_low=numpoints_low_count
       endif
     endif
     if(myid==0)write(100,*)'Low-level calculations started, but not completed.'
     if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_low_count
     if(myid==0)write(100,*)'Restarting calculations to complete the Low-level ab initio grid.'
     if(myid==0)write(100,*)
   else! 
     if(myid==0)then
       open(unit=334,file='AbINITIO_low.dat',status='old')
       do i=1,maxpoints
         if(ab_flag2==2)then
          read(334,*,end=16)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
         else
          read(334,*,end=16)ncont,seed_low(i,:),low_pot(i)
         endif
       enddo
       16 backspace(unit=334)
       if(numpoints_low_count/=numpoints_low)then
         do j=1,2000! 
           write(charmyid,'(I5)')j-1
           inquire(file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                          '.dat',exist=logica2)
           if(logica2)then  
             open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                                  '.dat', action='read',status='old')
             do i=1,maxpoints
               if(ab_flag2==2)then
                read(230,*,end=17)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
                write(334,f602)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
               else
                read(230,*,end=17)ncont,seed_low(i,:),low_pot(i)
                write(334,601)ncont,seed_low(i,:),low_pot(i)
               endif
             enddo
17           close(230, status="delete")
           endif
         enddo
       endif
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call MPI_BCAST(numpoints_low_count,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     do i=1,numpoints_low_count
       call MPI_BCAST(seed_low(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(low_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       if(ab_flag2==2)then
        call MPI_BCAST(seed_grad_low(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       endif
     enddo
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if(numpoints_low_count>=numpoints_low)then
      if(myid==0)close(334)
      if(myid==0)write(100,*)'All low-level ab initio energies have been computed.'
      xpass=1
      if(numpoints_low_count>numpoints_low)then
        if(myid==0)write(110,*)'Number of points increased from ',numpoints_low,' to ',numpoints_low_count
        numpoints_low=numpoints_low_count
      endif
     endif
     if(myid==0)write(100,*)'Low-level calculations started, but not completed.'
     if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_low_count
     if(myid==0)write(100,*)'Restarting calculations to complete the Low-level ab initio grid.'
     if(myid==0)write(100,*)
   endif
 ENDIF! 
 155 continue
 if(myid==0)then
  write(xchar3,'(I10)')nabi
  open(unit=1000,file='input-AUTOSURF-ABI.dat',status='old',action='read')
  114 read(1000,'(A100)')line
  ncont1=INDEX(line,'CODE CONTROL:')
  if (ncont1==0) goto 114
  read(1000,*)
  read(1000,*)xcontrol
  if (xcontrol==0) goto 1000
  if (xcontrol==1) goto 1001
  if (xcontrol==21)goto 1021
  stop
  1000 continue
  numpoints_low=1
  read(1000,*)th1
  jac(2)=dcos(th1*pii/180.d0)
  read(1000,*)th2
  jac(3)=dcos(th2*pii/180.d0)
  read(1000,*)phi
  jac(4)=phi*pii/180.d0
  read(1000,*)xrr
  jac(1)=xrr
  close(1000)
  if (nabi==1) then
    write(100,500)xrr,th1,th2,phi
    500 format(' Fixed coordinates:  (xrr=',F6.2,';th1=',F5.1,'; th2=',F5.1,'; phi=',F5.1,')')
    write(100,*)
  endif
  seed_low(1,:)=jac(:)
  goto 666
  1001 continue
  read(1000,*)numpoints_low
  close(1000) 
  write(100,*)'*** Header used for the ab initio calculation: "molpro'//trim(adjustl(xchar3))//'.abi"'
  write(100,*)'Geometries will be obtained from the external file: "GEOMETRIES.dat"'
  write(100,*)'Number of geometries to be computed:',numpoints_low
  write(100,*)
  open(unit=6660,file='GEOMETRIES.dat',status='old')
  do i=1,numpoints_low
    if(ab_flag2==1)read(6660,*)ncont1,ran_vec(1:4)
    if(ab_flag2==2)read(6660,*)ncont1,ran_vec(1:4)
    seed_low(i,:)=ran_vec(:)
  enddo
  close(6660)
  goto 666
  1021 continue
  read(1000,*)numpoints_low
  read(1000,*)th1
  jac(2)=dcos(th1*pii/180.d0)
  read(1000,*)th2
  jac(3)=dcos(th2*pii/180.d0)
  read(1000,*)phi
  jac(4)=phi*pii/180.d0
  read(1000,*)xrmin
  read(1000,*)xrmax
  close(1000)
  write(100,*)'*** Header used for the ab initio calculation: "molpro'//trim(adjustl(xchar3))//'.abi"'
  write(100,*)'1D cut of the PES (U vs. R)'
  write(100,*)'Number of geometries to be computed:',numpoints_low
  write(100,501)th1,th2,phi
  501 format(' Fixed coordinates:  (th1=',F5.1,'; th2=',F5.1,'; phi=',F5.1,')')
  write(100,*)
  dd=(xrmax-xrmin)/dble(numpoints_low-1)
  do i=1,numpoints_low
    jac(1)=xrmin+dble(i-1)*dd
    seed_low(i,:)=jac(:)
  enddo
  goto 666
 666 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_BCAST(numpoints_low,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do i=1,numpoints_low
   call MPI_BCAST(seed_low(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 ncont=0
 ncont1=0! 
 write(charmyid,'(I5)')myid
 if(myid/=0)open(unit=1000+myid,file='LGRID-'//sys_label(1:nline)//'-Proc'// &
                                       trim(adjustl(charmyid))//'.dat', action='write')
 61 continue! 
 IF(myid/=0)THEN! 
   call MPI_SEND(myid,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,ierr)
   call MPI_RECV(ncont,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,status,ierr)
   if(ncont==0)then! 
     close(1000+myid)
     write(*,*) 'Process:',myid,', calculated energies:',ncont1
     goto 63
   else! 
     ncont1=ncont1+1
     jac(:)=seed_low(ncont,:)
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     call ab_initio(cart,symb,mass,poten,dcart,code_flag2,ab_flag2,natom,myid,nabi,midbond,jac(1))
     poten=poten*CONVE
     low_pot(ncont)=poten
     if(ab_flag2==2)then
       dcart=dcart*ugrad  
       seed_grad_low(ncont,:)=dcart
       write(1000+myid,*)ncont,jac(:),low_pot(ncont),seed_grad_low(ncont,:)
     else
       write(1000+myid,*)ncont,jac(:),low_pot(ncont)
     endif
     goto 61
   endif
 ELSE! 
   call MPI_RECV(nid,1,MPI_INT,MPI_ANY_SOURCE,send_data_tag,MPI_COMM_WORLD,status,ierr)
   if(ncont<numpoints_low)then! 
     ncont=ncont+1
     call MPI_SEND(ncont,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
   else
     ncont1=ncont1+1
     write(*,*) 'MASTER, sending final message to process:', nid
     call MPI_SEND(0,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
     if (ncont1==numprocs-1) goto 63! 
   endif
   goto 61
 ENDIF
 63 continue
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0)then
   do j=2,numprocs
     write(charmyid,'(I5)')j-1
     open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                          '.dat',action='read',status='old')
     do i=1,numpoints_low
       if(ab_flag2==2)then
        read(230,*,end=13)ncont,jac(:),xpot,xgrad(:)
        low_pot(ncont)=xpot
        seed_grad_low(ncont,:)=xgrad(:)
       else
        read(230,*,end=13)ncont,jac(:),xpot
        low_pot(ncont)=xpot
       endif
     enddo
     13 close(230, status="delete")
   enddo 
   IF (xcontrol==0) THEN
     if(ab_flag2==2)then
       write(100,f602)i,seed_low(1,:),low_pot(1),seed_grad_low(1,:)
     else
       write(xchar3,'(I10)')nabi
       write(100,'(F12.6,F20.8,A16)')seed_low(1,1),low_pot(1),'molpro'//trim(adjustl(xchar3))//'.abi'
     endif
   ELSEIF (xcontrol==1) THEN
     do i=1,numpoints_low
       if(ab_flag2==2)then
        write(100,f602)i,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
       else
        write(100,601)i,seed_low(i,:),low_pot(i)
       endif
     enddo
     write(100,*)
   ELSEIF (xcontrol==21) THEN
     do i=1,numpoints_low
       if(ab_flag2==2)then
        write(100,f603)i,seed_low(i,1),low_pot(i),seed_grad_low(i,:)
       else
        write(100,'(I4,F12.6,F20.8)')i,seed_low(i,1),low_pot(i)
       endif
     enddo
     write(100,*)
   ENDIF
   ncont=0
   do i=1,numpoints_low
     do j=i+1,numpoints_low
       x1=low_pot(i)-low_pot(j)
       if((seed_low(i,1)==seed_low(j,1)).and.(dabs(x1)<=1.d-7))then
         ncont=ncont+1
       endif
     enddo
   enddo
   if(ncont/=0)write(100,*)
   if(ncont/=0)write(100,*)'***********'
   if(ncont/=0)write(100,*)'* Warning * Geometries with the same "R" and "E" in the computed data:',ncont
   if(ncont/=0)write(100,*)'***********'
   if(ncont/=0)write(100,*)
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if (nabi<xnabi) then
   nabi=nabi+1
   goto 155
 endif
 if(myid==0)then
   write(100,*)'All the ab initio energies have been computed.'
   call date_and_time(bdate1,bdate2,bdate3,date_time)
   write(100,*)
   write(100,*)
   write(100,'(A74)')'                        MPI ENVIROMENT CLOSED                             '
   write(100,'(A74)')' -------------------------------------------------------------------------'
   write(100,'(A74)')' *           System date (dd/mm/yyyy): '//bdate1(7:8)//'/'//bdate1(5:6)//'/'&
                      //bdate1(1:4)//'                        *'
   write(100,'(A74)')' *           System time (hh:min:sec): '//bdate2(1:2)//':'//bdate2(3:4)//':'&
                      //bdate2(5:6)//'                          *'
   write(100,'(A74)')' -------------------------------------------------------------------------'
   close(100)
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_FINALIZE(rcc)
 stop

stop
END PROGRAM AUTOSURF_ABI_rigid4D
