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
!-     "AUTOSURF-PES_rigid4D": PROGRAM for the automated construction of a         -!
!-       4D-PES on a vdW system composed of two (rigid) linear fragments.          -!
!-----------------------------------------------------------------------------------!
!      Input files: "input-AUTOSURF-PES.dat.dat", "molproX.abi" (X = 1, 2 & 10)     !
!      Dependencies: "FUNCTIONS_rigid4D.f90", "MODULES_rigid4D.f90", "diag.f",      !
!                    "ab_initio.F90", and "SUBROUTINES_rigid4D.f90"                 !
!      Version: 1.3                                                                 !
!***********************************************************************************!

PROGRAM AUTOSURF_PES_rigid4D

use nr
use nrtype
use nrutil
use dynamic_parameters
use mpi
!-----------------------------------------------------------------------------------
implicit none
 integer :: i,j,k,k2,l,ip,jp,kp,ipp,jpp,kpp,lpp,skip,jj,numpoints,numpoints_low,  &
             holepatch,term,wellfocus,nloop,nloop1
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
!5
 character(len=300) :: sys_label,line,f602
 character(len=10) :: bdate1,bdate2,bdate3
 character(len=5) :: charmyid
 character(len=8) :: xchar,xchar1
 character(len=15),dimension(4) :: charunitsE
 character(len=15),dimension(2) :: charunitsD
 character(len=40) :: NAME1
 integer :: nline,summyid,ncont,ncont1,ncont2,nid,numpoints_low_count,old_numprocs
 integer :: xpass,numpoints_high_count,n_test,xunitE,focus,only_seed_grid,status
 integer :: tot_abinitio,ntest,auto_mod,focus_onR,focus_onLR,xseed
 integer,dimension(:),allocatable :: xseed1,xseed2
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
 real*8, dimension(400000) :: xx,eee,th1,th2,phi
 integer,parameter :: send_data_tag = 2001, return_data_tag = 2002
 integer,parameter :: midbond = 0 ! 

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
call RANDOM_SEED(size=i) 
allocate(xseed1(i),xseed2(i))
inquire(file='input-AUTOSURF-PES.dat',exist=logica1)
if(.not.logica1)then
 if(myid==0)write(*,*)
 if(myid==0)write(*,*)'ERROR: The file: input-AUTOSURF-PES.dat does not exist !! '
 stop 
endif
open(unit=myid+1000,file='input-AUTOSURF-PES.dat',status='old',action='read')
111 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'GENERAL INFORMATION:')
if(ncont1==0)goto 111
read(myid+1000,*)
read(myid+1000,*)sys_label ! 
nline=scan(sys_label,' ')-1
if (nline>15) nline=15
read(myid+1000,*) restart  !
read(myid+1000,*) xseed    ! 
xunitE=2 ! kcal/mol
close(myid+1000)
call sobseq(ran_vec,4)
!do i=1,11
do i=1,2000
  call sobseq(ran_vec)
enddo
if(xseed==0)then!
  xseed1=date_time(5)+date_time(6)+date_time(7)+date_time(8)
  if(restart==0)xseed1=10
  CALL RANDOM_SEED(PUT=xseed1) 
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0)then
  open(unit=110,file='output.dat')
  write(110,*) 'All processes are alive!'
  if(restart==1)then
    inquire(file='GENERAL-'//sys_label(1:nline)//'.dat',exist=logica1)
    if(.not.logica1)then
     write(110,*)
     write(110,*)'Warning: The file: GENERAL-'//sys_label(1:nline)//'.dat does not exist !! '
     call system('echo "'//sys_label(1:nline)//'" > GENERAL-'//sys_label(1:nline)//'.dat')
    endif
    open(unit=100,file='GENERAL-'//sys_label(1:nline)//'.dat',status='old')
    11 read(100,*,end=12)
    goto 11
    12 backspace(unit=100)
    if (.not.logica1) then
      write(100,*)
      write(100,*)'Warning: the "GENERAL.dat" file was not present at the time this'
      write(100,*)'         calculation started, and had to be created!'
    endif
    write(100,*)
    write(100,*)
    write(100,'(A74)')' ========================================================================='
    write(100,'(A74)')' *           The construction of the PES has been restarted...           *'
    write(100,'(A74)')' *           System date (dd/mm/yyyy): '//bdate1(7:8)//'/'// &
                       bdate1(5:6)//'/'//bdate1(1:4)//'                        *'
    write(100,'(A74)')' *           System time (hh:min:sec): '//bdate2(1:2)//':'// &
                       bdate2(3:4)//':'//bdate2(5:6)//'                          *'
    write(100,'(A74)')' ========================================================================='
    write(100,*)
    if(xseed/=0)then! 
      xseed2=date_time(5)+date_time(6)+date_time(7)+date_time(8)!
      call RANDOM_SEED(PUT=xseed2) 
    endif
  else!
    inquire(file='AbINITIO_low.dat',exist=logica1)
    if(logica1)then
     call system('mv AbINITIO_low.dat AbINITIO_low-temp.dat')
     write(110,*)
     write(110,*)'Previous (LOW-LEVEL-PES) data found. A copy was saved in: "AbINITIO_low-temp.dat"'
    endif
    inquire(file='AbINITIO.dat',exist=logica1)
    if(logica1)then
     call system('mv AbINITIO.dat AbINITIO-temp.dat')
     write(110,*)
     write(110,*) 'Previous (HIGH-LEVEL-PES) data found. A copy was saved in: "AbINITIO-temp.dat"'
    endif
    inquire(file='GENERAL-'//sys_label(1:nline)//'.dat',exist=logica1)
    if(logica1)then
     call system('mv GENERAL-'//sys_label(1:nline)//'.dat GENERAL-temp.dat')
     write(110,*) 'Previous GENERAL file found. A copy was saved in: "GENERAL-temp.dat"'
    endif
    open(unit=100,file='GENERAL-'//sys_label(1:nline)//'.dat', action='write')
    write(100,*)'System:  '//sys_label(1:nline)//' '
    write(100,*)
    if(logica1)then
     write(100,*) 'Previous GENERAL file found. A copy was saved in: "GENERAL-temp.dat"'
    endif
    write(100,*)
    write(100,'(A74)')' *************************************************************************'
    write(100,'(A74)')' *           The automated construction of the PES has started           *'
    write(100,'(A74)')' *           System date (dd/mm/yyyy): '//bdate1(7:8)//'/'//bdate1(5:6)//'/' &
                       //bdate1(1:4)//'                        *'
    write(100,'(A74)')' *           System time (hh:min:sec): '//bdate2(1:2)//':'//bdate2(3:4)//':' &
                       //bdate2(5:6)//'                          *'
    write(100,'(A74)')' *************************************************************************'
    write(100,*)
    write(100,*)
    if(xseed/=0)then
       xseed2=xseed
       call RANDOM_SEED(PUT=xseed2)
    endif
    inquire(file='AbINITIO_low.dat',exist=logica1)
    if(logica1)then
     open(unit=200,file='AbINITIO_low.dat')
     close(200, status="delete")
    endif
    inquire(file='AbINITIO.dat',exist=logica1)
    if(logica1)then
     open(unit=200,file='AbINITIO.dat')
     close(200, status="delete")
    endif
    do j=1,2000
      write(charmyid,'(I5)')j-1
      inquire(file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//'.dat', &
               exist=logica1)
      if(logica1)then  
       open(unit=200,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//  &
                            '.dat')
       close(200, status="delete")
      endif
      inquire(file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//'.dat', &
               exist=logica1)
      if(logica1)then  
       open(unit=200,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))//  &
                            '.dat')
       close(200, status="delete")
      endif
    enddo
  endif
  write(100,'(A74)')' -------------------------------------------------------------------------'
  write(100,*)'                        STARTING MPI ENVIROMENT                          '
  write(100,*)
  write(100,*)'Unless explicitly specified, energies are '
  write(100,*)'in '//trim(adjustl(charunitsE(xunitE)))// &
               ', distances are in '//trim(adjustl(charunitsD(xunitD)))//'.'
  write(100,'(A30,I6)')' Total number of processes:   ',numprocs
endif
open(unit=myid+1000,file='input-AUTOSURF-PES.dat',status='old',action='read')
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
113 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'CODE CONTROL:')
if (ncont1==0) goto 113
read(myid+1000,*)
read(myid+1000,*) acc! 
read(myid+1000,*) nloop1!  
nloop=nloop1
acc=acc/hart2kcl! 
read(myid+1000,*) E_range!
E_range=E_range/hart2kcl! 
read(myid+1000,*) focus! 
read(myid+1000,*) increment!
increment=increment/hart2kcl! 
read(myid+1000,*) rmin(1)!
read(myid+1000,*) rmax(1)!
if(rmin(1)>rmax(1))then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "rmax(1)" should be larger than "rmin(1)".'
  call MPI_FINALIZE(rcc)
  stop
endif
read(myid+1000,*) focus_onR! 
read(myid+1000,*) minR! 
read(myid+1000,*) maxR! 
if(focus_onR==0)then
  minR=rmin(1)
  maxR=rmax(1)
else
 if(minR>maxR)then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "maxR" should be larger than "minR".'
  call MPI_FINALIZE(rcc)
  stop
 endif
endif
read(myid+1000,*) focus_onLR! 
read(myid+1000,*) xminLR! 
if(focus_onLR==1)then
  if((xminLR>=rmax(1)).or.(xminLR<rmin(1)))then
    if(myid==0)write(110,*)'Error found in the input file!' 
    if(myid==0)write(110,*)'The value of "xminLR" should be: rmin(1) < xminLR < rmax(1)'
    if(myid==0)write(110,*)'xminLR = ',xminLR    
    if(myid==0)write(110,*)'rmin(1) = ', rmin(1)
    if(myid==0)write(110,*)'rmax(1) = ', rmax(1)
    call MPI_FINALIZE(rcc)
    stop
  endif
endif
read(myid+1000,*) numpoints! 
if (numpoints==0) numpoints=(2000*natom1*natom2)/(9*flip1*flip2)
read(myid+1000,*) ab_flag! 
read(myid+1000,*) code_flag! 
read(myid+1000,*) x1! 
if(x1<=0.d0)then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "numadded" should be a positive number'
  call MPI_FINALIZE(rcc)
  stop
endif
x1=x1*dble(numprocs)
numadded=int(x1)
n_test=0
ntest=n_test
x1=10
x1=x1*dble(numprocs)
num_ab_err_points=int(x1)
maxpoints=100000
read(myid+1000,*) order_1! 
read(myid+1000,*) order_2! 
read(myid+1000,*) order_3! 
read(myid+1000,*) order_4!
dist_tol=1.25d0
read(myid+1000,*) num_err_points! 
if(num_err_points<numprocs)then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "num_err_points" should be greater than the number of processors:',numprocs
  call MPI_FINALIZE(rcc)
  stop
endif
num_err_points=num_err_points/numprocs! 
read(myid+1000,*) low_grid! 
if ((code_flag==4).and.(low_grid/=2)) then
  inquire(file='molpro1.abi',exist=logica1)
  if(.not.logica1)then
    if(myid==0)write(*,*)
    if(myid==0)write(*,*)'ERROR: The file: molpro1.abi does not exist !! '
    stop 
  else
    inquire(file='molpro2.abi',exist=logica2)
    if(.not.logica2)then
      if(myid==0)call system('cp molpro1.abi molpro2.abi')
    endif
  endif
elseif ((code_flag==1).and.(low_grid/=2)) then
  inquire(file='gaussi1.abi',exist=logica1)
  if(.not.logica1)then
    if(myid==0)write(*,*)
    if(myid==0)write(*,*)'ERROR: The file: gaussi1.abi does not exist !! '
    stop 
  else
    inquire(file='gaussi2.abi',exist=logica2)
    if(.not.logica2)then
      if(myid==0)call system('cp gaussi1.abi gaussi2.abi')
    endif
  endif
endif
if(low_grid>0)then
 read(myid+1000,*)
 read(myid+1000,*) numpoints_low! 
 if (numpoints_low==0) numpoints_low=(10000*natom1*natom2)/(36*flip1*flip2)! 
 read(myid+1000,*) code_flag2! 
 if (code_flag2==4) then
   inquire(file='molpro10.abi',exist=logica1)
   if(.not.logica1)then
     if(myid==0)write(*,*)
     if(myid==0)write(*,*)'ERROR: The file: molpro10.abi does not exist !! '
     stop 
   endif
 elseif (code_flag2==1) then
   inquire(file='gaussi10.abi',exist=logica1)
   if(.not.logica1)then
     if(myid==0)write(*,*)
     if(myid==0)write(*,*)'ERROR: The file: gaussi10.abi does not exist !! '
     stop 
   endif
 endif
 read(myid+1000,*) ab_flag2! 
 low_high_fac=1.2d0
 seed_cut=3.d0
endif
subzero=0
close(myid+1000)
pii=dacos(-1d0)
term=0! 
xpass=0! 
alpha=-1.0d0! 
xbeta=1.0d0
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
acc=acc*CONVE
E_range=E_range*CONVE
increment=increment*CONVE
order_1_min=3   
order_2_min=3
order_3_min=3
order_4_min=3
call basis_size(4,order_1,order_2,order_3,order_4,basis_1)  
if(numpoints*2*symparts<basis_1)then
  if(myid==0)write(110,*)'Basis too big for number of points. Increase value of "numpoints" on input file' 
  if(myid==0)write(110,*)'numpoints =',numpoints,'     basis_1 =',basis_1
  call MPI_FINALIZE(rcc)
  stop
endif
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
allocate(xMAPot(numprocs,xnum_intervals_R+1),xgrad(3*natom))
if(ab_flag==2)then
  allocate(seed_grad(maxpoints,3*(natom)),grad(2*symparts*maxpoints,4))
endif
if(low_grid>0)then
  allocate(low_pot(maxpoints),seed_low(maxpoints,4),seed_pot_low(maxpoints))
  if(ab_flag2==2)then
   allocate(seed_grad_low(maxpoints,3*(natom)))
  endif
endif
dcart=0d0
601 FORMAT(I10,4f20.15,f20.8)              !
write(f602,'( "(I10,4f20.14,f20.8,",I3,"f20.8)" )')natom*3 !! 
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0)then
  write(100,'(A30,I6)')' High-degree basis size:      ',basis_1
  write(100,'(A30,I6)')' Lower-degree basis size:     ',basis_2
  write(100,'(A30,I6)')' Size of minimal basis:       ',basis_3
  write(100,*)'Initialization completed.'
  write(100,*)
endif
if (low_grid>0) then
 if(myid==0)write(110,*)
 if(myid==0)write(110,*)'*** Low-level-PES'
 if(myid==0)write(110,*)
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
           do i=numpoints_low_count+1,maxpoints
             if(ab_flag2==2)then
              read(230,*,end=15)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
              write(334,f602)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
             else
              read(230,*,end=15)ncont,seed_low(i,:),low_pot(i)
              write(334,601)ncont,seed_low(i,:),low_pot(i)
             endif
             numpoints_low_count=numpoints_low_count+1
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
       goto 18! 
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
         numpoints_low_count=numpoints_low_count+1
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
             do i=numpoints_low_count+1,maxpoints
               if(ab_flag2==2)then
                read(230,*,end=17)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
                write(334,f602)ncont1,seed_low(i,:),low_pot(i),seed_grad_low(i,:)
               else
                read(230,*,end=17)ncont,seed_low(i,:),low_pot(i)
                write(334,601)ncont,seed_low(i,:),low_pot(i)
               endif
               numpoints_low_count=numpoints_low_count+1
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
      goto 18! 
     endif
     if(myid==0)write(100,*)'Low-level calculations started, but not completed.'
     if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_low_count
     if(myid==0)write(100,*)'Restarting calculations to complete the Low-level ab initio grid.'
     if(myid==0)write(100,*)
   endif
 ENDIF!
 if(myid==0)then
   if(xseed==0)then! 
     do i=1,numpoints_low_count*150! 
      call sobseq(ran_vec)
     enddo
   endif
   write(100,*)'*** Low-level grid computation'
   write(100,'(A18,I6)')' Number of points:',numpoints_low
   do i=numpoints_low_count+1,numpoints_low
874  if(xseed==0)then!
       call sobseq(ran_vec)
     else! 
       call random_number(ran_vec)
     endif
     ran_vec(1)=dexp(dlog(rmin(1))+ran_vec(1)*(log(rmax(1))-log(rmin(1))))!
     do j=2,4  
       range = rmax(j) - rmin(j)
       ran_vec(j)=rmin(j)+range*ran_vec(j)
     enddo
     seed_low(i,:)=ran_vec(:)
     call INT_Cart(cart,seed_low(i,:),mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if (dist_flag==1) goto 874
   enddo
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do i=numpoints_low_count+1,numpoints_low
   call MPI_BCAST(seed_low(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 ncont=numpoints_low_count
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
     call ab_initio(cart,symb,mass,poten,dcart,code_flag2,ab_flag2,natom,myid,10,midbond,jac(1))
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
     call MPI_SEND(0,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
     if (ncont1==numprocs-1) goto 63!
   endif
   goto 61
 ENDIF
 63 continue
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0)then
   if(restart==0)open(unit=334,file='AbINITIO_low.dat')
   do j=2,numprocs
     write(charmyid,'(I5)')j-1
     open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                          '.dat',action='read',status='old')
     do i=1,numpoints_low
       if(ab_flag2==2)then
        read(230,*,end=13)ncont,jac(:),xpot,xgrad(:)
        low_pot(ncont)=xpot
        seed_grad_low(ncont,:)=xgrad(:)
        write(334,f602)ncont,jac(:),xpot,xgrad(:)
       else
        read(230,*,end=13)ncont,jac(:),xpot
        low_pot(ncont)=xpot
        write(334,601)ncont,jac(:),xpot
       endif
     enddo
13   close(230, status="delete")
   enddo 
   close(334)
   write(100,*)'All low-level ab initio energies have been computed.'
   write(100,*)'RESTART POINT successfully created!'
   write(100,*)'Ab initio data saved in the file: "AbINITIO_low.dat"'
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do i=numpoints_low_count+1,numpoints_low
   call MPI_BCAST(low_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   if(ab_flag2==2)call MPI_BCAST(seed_grad_low(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 enddo
 18 continue! 
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0)then
   ncont=0
   do i=1,numpoints_low
     do j=i+1,numpoints_low
       x1=low_pot(i)-low_pot(j)
       if((seed_low(i,1)==seed_low(j,1)).and.(dabs(x1)<=1.d-7))then
         ncont=ncont+1
       endif
     enddo
   enddo
   if(ncont/=0)write(110,*)'***********'
   if(ncont/=0)write(110,*)'* Warning * there are ',ncont,' geometries with the same "R" and "E" in "AbINITIO_low.dat" !!'
   if(ncont/=0)write(110,*)'***********'
 endif
 count_seed=0
 Max_E3=2d2
 Min_E=-1d9
 Max_R=0d0
 do i=1,numpoints_low
   if(seed_low(i,1)>Max_R)then
    Max_R=seed_low(i,1)! 
    E_asym=low_pot(i)! 
   endif
   if(low_pot(i)<Max_E3)then
    Max_E3=low_pot(i)!
   endif     
   if(low_pot(i)>Min_E) then
    Min_E=low_pot(i)! 
   endif
 enddo
 bot_seed=Max_E3!
 Max_E_seed=E_asym+E_range*low_high_fac! 
 do i=1,numpoints_low
   if(low_pot(i)<E_asym+E_range*seed_cut)then
    count_seed=count_seed+1! 
   endif
 enddo
 allocate(coords_seed(2*symparts*maxpoints,4),ind_seed(maxpoints*2*symparts))
 allocate(ind2_seed(maxpoints*2*symparts),d_seed(maxpoints*2*symparts))
 allocate(pot_seed(maxpoints*2*symparts),b2_seed(basis_3,maxpoints*2*symparts))
 if(ab_flag2==2)then
  allocate(grad_seed(maxpoints*2*symparts,4))
 endif
 if(myid==0)then
   if(xpass/=1)write(100,*) '*** Low-level-PES fit'
   if(xpass/=1)write(100,'(A40,F15.3)')' Lowest energy in the low-level grid:   ',bot_seed 
   if(xpass/=1)write(100,'(A40,F15.3)')' Highest energy in the low-level grid:  ',Min_E
   if(xpass/=1)write(100,'(A11,F7.3,A23,F14.3)')' Largest R: ',Max_R,' Ref. asymptote energy:',E_asym
   if(xpass/=1)write(100,'(A40,F15.3)')' Maximum energy allowed above asymptote:',E_asym+E_range*seed_cut
   write(110,'(A40,F15.3)')' Lowest energy in the low-level grid:  ',bot_seed 
   write(110,'(A40,F15.3)')' Highest energy in the low-level grid: ',Min_E
   write(110,'(A11,F7.3,A23,F14.3)')' Largest R: ',Max_R,' Ref. asymptote energy:',E_asym
   write(110,'(A40,F15.3)')' Maximum energy allowed above asymptote:',E_asym+E_range*seed_cut
 endif

 if(low_grid==2)then! 
   if(myid==0)then
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
    write(100,'(A74)')'                 Only the low-level-grid was computed.                    '
    close(100)
    call system('rm fort*')
    call system('rm Molscript*.x')
    call system('rm pun*')
   endif
   call MPI_FINALIZE(rcc)
   stop
 endif
 count=0
 do i=1,numpoints_low
   if(low_pot(i)<E_asym+E_range*seed_cut)then
     if(ab_flag2==2)then
      dcart=seed_grad_low(i,:)
     endif
     call symmetry2(seed_low(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,ab_flag2,exch,flip1,flip2)
     do k2=0,1!
       do k=1,symparts
         count=count+1
         coords_seed(count,:)=internal(k,:)
         coords_seed(count,4)=internal(k,4)*(-1d0)**k2
         pot_seed(count)=low_pot(i)
         if(ab_flag2==2)then
          grad_seed(count,:)=grad_int(k,:)
          grad_seed(count,4)=grad_int(k,4)*(-1d0)**k2
         endif
       enddo
     enddo
   endif
 enddo
 count_seed=count
 if(myid==0)then
   write(100,'(A36,I19)')' Total number of points in the fit: ',count_seed
   write(110,'(A36,I19)')' Total number of points in the fit: ',count_seed
   write(110,*)
 endif
 do i=1,count_seed
   count=0 
   do ip=1,count_seed         
    count=count+1
    call dist_metric(coords_seed(i,:),coords_seed(ip,:),W_a,ind_seed(count))
   enddo
   call indexxy(count_seed,ind_seed,ind2_seed)
   d_seed(i)=ind_seed(ind2_seed(zz4))
 enddo
 support=min(count_seed,8*basis_3/(4*(ab_flag2-1)+1))
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do i=myid+1,count_seed,numprocs
   count=0
   Jac=coords_seed(i,:)
   do ip=1,count_seed
     count=count+1
     Jac2=coords_seed(ip,:)
     call dist_metric(jac,jac2,W_a,somme)
     somme=somme**2
     ind_seed(count)=dsqrt(dexp(-((somme)/d_seed(ip)**2))/ &
                          (((somme)/d_seed(ip)**2)**(zz_low/2)+epss))      
   enddo
   call indexxy(count_seed,ind_seed,ind2_seed)
   call basis_calc_seed(i,ind_seed,ind2_seed)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do j=0,numprocs-1
   do i=j+1,count_seed,numprocs
    call MPI_BCAST(b2_seed(:,i),basis_3,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
   enddo
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0) then
   write(100,*)'Low-level-PES completed.'
   write(100,*)
 endif
else!
 if(myid==0)write(100,*)'No low-level ab initio calculations requested.'
endif !
IF(restart==1)THEN
  inquire(file='AbINITIO.dat',exist=logica1)
  if(.not.logica1)then!
    if(myid==0)then
      open(unit=222,file='AbINITIO.dat')
      do j=1,1000!
        write(charmyid,'(I5)')j-1
        inquire(file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                      '.dat',exist=logica2)
        if(logica2)then  
          open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                               '.dat', action='read',status='old')
          do i=numpoints_high_count+1,maxpoints
            if(ab_flag==2)then
             read(230,*,end=21)ncont1,seed(i,:),seed_pot(i),seed_grad(i,:)
             write(222,f602)ncont1,seed(i,:),seed_pot(i),seed_grad(i,:)
            else
             read(230,*,end=21)ncont,seed(i,:),seed_pot(i)
             write(222,601)ncont,seed(i,:),seed_pot(i)
            endif
            numpoints_high_count=numpoints_high_count+1
          enddo
21        close(230, status="delete")
        endif
      enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(numpoints_high_count,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=1,numpoints_high_count
      call MPI_BCAST(seed(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seed_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(ab_flag==2)then
       call MPI_BCAST(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      endif
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(numpoints_high_count/=0)then
      if(numpoints_high_count>=numpoints)then
       if(myid==0)close(222)
       if(myid==0)write(110,*)
       if(myid==0)write(110,*)'*** High-level ab initio seed-grid have been computed.'
       if(myid==0)write(100,*)'High-level ab initio seed-grid have been computed.'
       if(myid==0)write(100,*)'RESTART POINT successfully created!'
       if(myid==0)write(100,*)'Ab initio data saved in output file: "AbINITIO.dat".'
       xpass=2
       goto 24!
      endif
      if(myid==0)write(100,*)'High-level seed-grid calculation started, but not completed.'
      if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_high_count
      if(myid==0)write(100,*)'Restarting calculations to complete high-level ab initio seed-grid.'
      if(myid==0)write(100,*)
    endif
  else! 
    if(myid==0)then
      open(unit=222,file='AbINITIO.dat',status='old')
      do i=1,maxpoints
        if(ab_flag==2)then
         read(222,*,end=22)ncont1,seed(i,:),seed_pot(i),seed_grad(i,:)
        else
         read(222,*,end=22)ncont,seed(i,:),seed_pot(i)
        endif
        numpoints_high_count=numpoints_high_count+1
      enddo
22    backspace(unit=222)
      do j=1,1000! 
        write(charmyid,'(I5)')j-1
        inquire(file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                       '.dat',exist=logica2)
        if(logica2)then  
          open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                               '.dat', action='read',status='old')
          do i=numpoints_high_count+1,maxpoints
            if(ab_flag==2)then
             read(230,*,end=25)ncont1,seed(i,:),seed_pot(i),seed_grad(i,:)
             write(222,f602)ncont1,seed(i,:),seed_pot(i),seed_grad(i,:)
            else
             read(230,*,end=25)ncont,seed(i,:),seed_pot(i)
             write(222,601)ncont,seed(i,:),seed_pot(i)
            endif
            numpoints_high_count=numpoints_high_count+1
          enddo
25        close(230, status="delete")
        endif
      enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(numpoints_high_count,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=1,numpoints_high_count
      call MPI_BCAST(seed(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seed_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(ab_flag==2)then
       call MPI_BCAST(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      endif
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(numpoints_high_count<numpoints)then
      if(myid==0)write(100,*)'High-level seed-grid calculation started, but not completed.'
      if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_high_count
      if(myid==0)write(100,*)'Restarting calculations to complete high-level ab initio seed-grid.'
      if(myid==0)write(100,*)
    else
      if(myid==0)close(222)
      if(myid==0)write(110,*)'*** High-level ab initio seed-grid have been computed.'
      if(myid==0)write(100,*)'High-level ab initio seed-grid have been computed.'
      xpass=2
      goto 24!
    endif
  endif
ENDIF!
if(myid==0)then    
  if(xseed==0)then!
    ncont=numpoints_high_count
    if(numpoints_low_count>=numpoints_low)ncont=ncont+numpoints_low_count
    do i=1,ncont*150
     call sobseq(ran_vec)
    enddo
  endif
  write(110,*) '*** High-level seed grid computation'
  write(110,*)
  write(100,*) '*** High-level seed grid computation'
  write(100,'(A18,I6)')' Number of points:',numpoints
  do i=numpoints_high_count+1,numpoints
873 if(xseed==0)then!
      call sobseq(ran_vec)
    else! 
      call random_number(ran_vec)
    endif
    ran_vec(1)=dexp(dlog(rmin(1))+ran_vec(1)*(dlog(rmax(1))-dlog(rmin(1))))!
    do j=2,4  
      range = rmax(j) - rmin(j)
      ran_vec(j)=rmin(j)+range*ran_vec(j)
    enddo
    seed(i,:)=ran_vec(:)
    call INT_Cart(cart,seed(i,:),mass,natom1,natom2,ref1,ref2)
    call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
    if (dist_flag==1) goto 873
    if(low_grid>0)then! 
      xi=seed(i,:)
      temp=func_actual_seed(xi)
      if(temp>Max_E_seed) goto 873
    endif
  enddo
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do i=numpoints_high_count+1,numpoints
  call MPI_BCAST(seed(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
ncont=numpoints_high_count
ncont1=0!
write(charmyid,'(I5)')myid
if(myid/=0)open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'// &
                                     trim(adjustl(charmyid))//'.dat', action='write')
64 continue! 
IF(myid/=0)THEN!
  call MPI_SEND(myid,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,ierr)
  call MPI_RECV(ncont,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,status,ierr)
  if(ncont==0)then!
    close(1000+myid)
    write(*,*) 'Process:',myid,', calculated energies:',ncont1
    goto 65
  else! 
    ncont1=ncont1+1
    jac(:)=seed(ncont,:)
    pass_check=1 
    call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
    801 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,pass_check,midbond,jac(1))
    if(poten>199d0)then
      if(pass_check<2) then
       pass_check=pass_check+1
       goto 801
      endif
    endif
    poten=poten*CONVE
    seed_pot(ncont)=poten
    if(ab_flag==2)then
      dcart=dcart*ugrad 
      seed_grad(ncont,:)=dcart
      write(1000+myid,*)ncont,jac(:),seed_pot(ncont),seed_grad(ncont,:)
    else
      write(1000+myid,*)ncont,jac(:),seed_pot(ncont)
    endif
    goto 64
  endif
ELSE!
  call MPI_RECV(nid,1,MPI_INT,MPI_ANY_SOURCE,send_data_tag,MPI_COMM_WORLD,status,ierr)
  if(ncont<count3)then!
    ncont=ncont+1
    call MPI_SEND(ncont,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
  else
    ncont1=ncont1+1
    call MPI_SEND(0,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
    if (ncont1==numprocs-1) goto 65! 
  endif
  goto 64
ENDIF
65 continue
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0)then
  if(restart==0)open(unit=222,file='AbINITIO.dat')
  do j=2,numprocs
    write(charmyid,'(I5)')j-1
    open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                            '.dat', action='read',status='old')
    do i=1,count3
      if(ab_flag==2)then
       read(230,*,end=23)ncont,jac(:),xpot,xgrad(:)
       seed_pot(ncont)=xpot
       seed_grad(ncont,:)=xgrad(:)
       write(222,f602)ncont,jac(:),xpot,xgrad(:)
      else
       read(230,*,end=23)ncont,jac(:),xpot
       seed_pot(ncont)=xpot
       write(222,601)ncont,jac(:),xpot
      endif
    enddo
23  close(230, status="delete")
  enddo 
  close(222)
  write(100,*)'High-level ab initio seed-grid have been computed.'
  write(100,*)'RESTART POINT successfully created!'
  write(100,*)'Ab initio data saved in the file: "AbINITIO.dat".' 
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do i=1,count3
  call MPI_BCAST(seed_pot(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(ab_flag==2)call MPI_BCAST(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
enddo
24 continue! 
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
count3=max(numpoints_high_count,numpoints)
if(myid==0)then
  ncont=0
  do i=1,count3
    do j=i+1,count3
      x1=seed_pot(i)-seed_pot(j)
      if((seed(i,1)==seed(j,1)).and.(dabs(x1)<=1.d-7))then
        ncont=ncont+1
      endif
    enddo
  enddo
  if(ncont/=0)write(110,*)'***********'
  if(ncont/=0)write(110,*)'* Warning *  there are ',ncont,' geometries with the same "R" and "E" in "AbINITIO.dat" !!'
  if(ncont/=0)write(110,*)'***********'
  tot_abinitio=count3!
endif
Max_E=2d2
Min_E=-1d9
Max_R=0d0
do i=1,count3
  if(seed(i,1)>Max_R) then
   Max_R=seed(i,1)! 
   E_asym=seed_pot(i)!
  endif
  if(seed_pot(i)<Max_E) then! 
   Max_E=seed_pot(i)!
  endif
  if(seed_pot(i)>Min_E) then! 
   Min_E=seed_pot(i)! 
  endif
enddo
Glob_min=Max_E!
Max_E=E_asym+E_range! 
E_limit=E_asym+1.2d0*E_range
if(myid==0) then
  write(100,'(A40,F15.3)')' Lowest energy in the high-level grid:  ',Glob_min
  write(100,'(A40,F15.3)')' Highest energy in the high-level grid: ',Min_E
  write(100,'(A11,F7.3,A23,F14.3)')' Largest R: ',Max_R,' Ref. asymptote energy:',E_asym
  write(100,'(A40,F15.3)')' Maximum energy allowed above asymptote:',Max_E+0.2d0*E_range
  write(110,'(A29,F15.3)')' Lowest energy in the grid:  ',Glob_min
  write(110,'(A29,F15.3)')' Highest energy in the grid: ',Min_E
  write(110,'(A11,F7.3,A23,F14.3)')' Largest R: ',Max_R,' Ref. asymptote energy:',E_asym
  write(110,'(A40,F15.3)')' Maximum energy allowed above asymptote:',Max_E+0.2d0*E_range
endif
count=0
do i=1,count3 
  if(seed_pot(i)<Max_E+0.2d0*E_range)then
    if(ab_flag==2)then
     dcart=seed_grad(i,:)
    endif
    call symmetry2(seed(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,  &
                    ab_flag,exch,flip1,flip2)
    do k2=0,1!
      do k=1,symparts
        count=count+1
        coords(count,:)=internal(k,:)
        coords(count,4)=internal(k,4)*(-1d0)**k2
        pot(count)=seed_pot(i)
        if(ab_flag==2.and.subzero==0)then
          grad(count,:)=grad_int(k,:)
          grad(count,4)=grad_int(k,4)*(-1d0)**k2
        endif
        if(subzero==1)then
          xi=coords(count,:)
          temp2=dfunc_actual_seed(xi)
          pot(count)=seed_pot(i)-temp2(1)
          if(ab_flag==2)then
            do ip=1,4
             grad(count,ip)=grad_int(k,ip)-temp2(1+ip)
            enddo
          endif
        endif
      enddo
    enddo
  else
    if(myid==0)write(110,*)'energy excluded from HGRID: ',i,seed_pot(i)
  endif
enddo
count3=count
if(myid==0)then
  write(100,'(A36,I19)')' Total number of points in the fit: ',count3
  write(110,'(A36,I19)')' Total number of points in the fit: ',count3
endif
allocate(ind(count3),ind2(count3))
do i=1,count3  
  count=0 
  do ip=1,count3         
    count=count+1
    call dist_metric(coords(i,:),coords(ip,:),W_a,ind(count))
  enddo
  call indexxy(count3,ind,ind2)
  d(i)=ind(ind2(zz4))
enddo
support=min(count3,4*basis_1/(4*(ab_flag-1)+1)) 
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do i=myid+1,count3,numprocs  
  count=0
  Jac=coords(i,:)
  do ip=1,count3
    count=count+1
    Jac2=coords(ip,:)
    call dist_metric(jac,jac2,W_a,somme)
    somme=somme**2
    ind(count)=dsqrt(dexp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss))      
  enddo
  call indexxy(count3,ind,ind2)
  call basis_calc(i,ind,ind2)
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do j=0,numprocs-1
   do i=j+1,count3,numprocs
      call MPI_BCAST(b2(:,i),basis_1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(b2_lower(:,i),basis_2,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(b2_minimal(:,i),basis_3,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
   enddo
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(myid==0) then
 write(100,*)'High-level-seed-PES completed.'
 write(110,*)'High-level-seed-PES completed.'
endif
goto 26!
if(myid==0)then
 open(unit=333,file='DEBUGGING.dat')
 jac3(:)=coords(12,:)+0.01d0
 write(333,*) jac3
 xi=jac3
 temp=func_actual(xi)
 write(333,*) temp,pot(12)
 write(333,*)
 temp2=dfunc_actual_anal1(xi)
 write(333,*) temp2
 write(333,*)
 temp2=dfunc_actual_anal2(xi)
 write(333,*) temp2
 write(333,*)
 temp=func(xi)
 write(333,*) temp
 temp4=dfunc(xi)
 write(333,*) temp4
 write(333,*)
 xi=jac3
 xi(1)=jac3(1)+0.0001d0
 temp=func_actual(xi)
 xi(1)=jac3(1)-0.0001d0
 temp3=func_actual(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(2)=jac3(2)+0.0001d0
 temp=func_actual(xi)
 xi(2)=jac3(2)-0.0001d0
 temp3=func_actual(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(3)=jac3(3)+0.0001d0
 temp=func_actual(xi)
 xi(3)=jac3(3)-0.0001d0
 temp3=func_actual(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(4)=jac3(4)+0.0001d0
 temp=func_actual(xi)
 xi(4)=jac3(4)-0.0001d0
 temp3=func_actual(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 temp4=dfunc(xi)
 write(333,*) 'temp4'
 write(333,*) temp4
 write(333,*)
 xi=jac3
 xi=jac3
 xi(1)=jac3(1)+0.0001d0
 temp=func(xi)
 xi(1)=jac3(1)-0.0001d0
 temp3=func(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(2)=jac3(2)+0.0001d0
 temp=func(xi)
 xi(2)=jac3(2)-0.0001d0
 temp3=func(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(3)=jac3(3)+0.0001d0
 temp=func(xi)
 xi(3)=jac3(3)-0.0001d0
 temp3=func(xi)
 write(333,*) (temp-temp3)/0.0002d0
 xi=jac3
 xi(4)=jac3(4)+0.0001d0
 temp=func(xi)
 xi(4)=jac3(4)-0.0001d0
 temp3=func(xi)
 write(333,*) (temp-temp3)/0.0002d0
 close(333)
endif! 
26 continue
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
add=1+(maxpoints-(count3/symparts/2))/numadded!
if(add<1)then
  if(myid==0)write(110,*)'Error: Number of loops lower than 1...'
  if(myid==0)write(110,*)'Increase the value of "maxpoints" (in the input file) & restart...'
  call MPI_FINALIZE(rcc)
  stop
endif
if(myid==0)then
 write(100,*)
 write(100,*)'***  START AUTOMATIC DATA POINT GENERATION  ***'
 write(100,'(A29,I6)')'      max. number of cycles: ',nloop1
 write(100,'(A23,ES12.5,A16)')'      accuracy target: ',acc,' '//adjustl(charunitsE(xunitE))
 write(100,'(A54,I6)')'      points used (per processor) to estimate errors: ',num_err_points
 if(focus_onLR==1)then
   write(100,'(A40,F5.2,A4,F5.2)')'      focus only in the long-range: R > ',xminLR
 else
   if(focus==1)write(100,'(A26,F7.3,A24)') &
     '      focus on energies of ',increment,' '//trim(adjustl(charunitsE(xunitE)))//' below E_asym'
   if(focus_onR==1)write(100,'(A50,F5.2,A4,F5.2)') &
     '      focus only on intermolecular distances from ',minR,' to ',maxR
 endif
 write(100,*)
 write(100,'(A48)')' LOOP      RMS          MEAN DEV.      # points  '
 write(110,*)
 write(110,*)'***  START AUTOMATIC DATA POINT GENERATION  ***'
endif
if(restart==0)then
  xseed=1000+xseed+(myid+1)! 
  xseed2=xseed
  CALL RANDOM_SEED(PUT=xseed2) 
else
  xseed=date_time(7)+date_time(8)+(myid+1)!
  xseed2=xseed
  CALL RANDOM_SEED(PUT=xseed2) 
endif
if(restart==1)then! 
  ncont=0
  if(numpoints_high_count>=numpoints)ncont=ncont+numpoints_high_count
  if(numpoints_low_count>=numpoints_low)ncont=ncont+numpoints_low_count
  do i=1,ncont*150
   call sobseq(ran_vec)
  enddo
endif
DO loop=1,add !
 numpoints_high_count=count3/(2*symparts)
 if(myid==0)write(110,*)
 if(myid==0)write(110,*)'Loop: ',loop-1
 ncont=0
 rms_error=0d0
 mean_error=0d0
 error_points=0d0 
 hole_val=Glob_min
 hole(:,1)=rmax(1)
 hole(:,2)=rmax(2)
 hole(:,3)=rmax(3)
 hole(:,4)=rmax(4)
 do i=1,num_err_points!
   19 call random_number(ran_vec)
   ncont=ncont+1
   if(ncont==1000000)then!
     write(*,*)'After a million iterations, processor ',myid
     write(*,*)'have only found ',i,' of the ',num_err_points,'points requested.'
     write(*,*)
     write(*,*)'Please, increase the volume of the available configuration space'
     write(*,*)'by loosen the "energy-focus" (if any);'
     write(*,*)'or reduce the number of points used to randomly estimate the error'
     write(*,*)'at the beginning of each loop ("num_err_points", in the input file)'
     write(*,*)'and restart the calculation.'
     stop
   endif
   if(focus_onLR==1)then
     range=rmax(1)-xminLR!
     ran_vec(1)=xminLR+range*ran_vec(1)
   else
     range=maxR-minR!
     ran_vec(1)=minR+range*ran_vec(1)
   endif
   do j=2,4
     range=rmax(j)-rmin(j)
     ran_vec(j)=rmin(j)+range*ran_vec(j)
   enddo
   call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
   call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
   if(dist_flag==1) goto 19! 
   xi=ran_vec(:)
   if(low_grid>0)then! 
     temp3=func_actual_seed(xi)!
     if(temp3>Max_E_seed) goto 19
   endif
   if(focus_onLR==1)goto 191!
   if(focus>0)then! 
     if (func_actual(xi)>E_asym+(0.05d0/hart2kcl*CONVE)-increment) goto 19
   else 
     if (wellfocus>0)then! 
       if (func_actual(xi)>E_asym+(0.05d0/hart2kcl*CONVE)) goto 19
     endif
   endif
   191 continue
   if(subzero==0)then
     temp=func_actual_min(xi)
     if(temp>Max_E) goto 19!
     temp=func_actual(xi)
   else
     temp=func_actual_min(xi)
     if(temp+temp3>Max_E) goto 19
     temp=func_actual(xi)+temp3
   endif
   if(temp<hole_val(myid+1))then
     hole_val(myid+1)=temp
     hole(myid+1,1:4)=xi(1:4)
   endif
   tampon=func(xi)
   if(i<=num_err_points1)then 
     error_points(myid*num_err_points1+i,1:4)=xi(1:4)
     error_points(myid*num_err_points1+i,4+1)=dabs(tampon)
   endif
   rms_error(myid+1)=rms_error(myid+1)+dabs(tampon)
   mean_error(myid+1)=mean_error(myid+1)+dsqrt(dabs(tampon))
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do j=0,numprocs-1   
   do i=1,num_err_points1
    call MPI_BCAST(error_points(j*num_err_points1+i,:),5,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
   enddo
 enddo
 do i=0,numprocs-1
   call MPI_BCAST(rms_error(i+1),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(mean_error(i+1),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(hole_val(i+1),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(hole(i+1,:),4,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 holepatch=0
 call indexxy(numprocs,hole_val,indhole)
 if(hole_val(indhole(1))<Glob_min-0.1d0/hart2kcl*CONVE)then! 
   holepatch=1
   if(myid==0) then
     write(110,*)'   Holes found...'
     do i=1,numprocs
      write(110,'(1f20.3,4f14.6)') hole_val(indhole(i)),hole(indhole(i),:)
     enddo
   endif
 endif
 do i=1,numprocs-1
   rms_error(1)=rms_error(1)+rms_error(1+i)
   mean_error(1)=mean_error(1)+mean_error(1+i)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0)then
   rms_error(1)=dsqrt(rms_error(1)/dble(numprocs*num_err_points))
   mean_error(1)=mean_error(1)/dble(numprocs*num_err_points)
   somme=rms_error(1)
   somme2=mean_error(1)
   if(holepatch==0)write(100,'(I4,2ES15.6,I10)')loop-1,somme,somme2,(count3/(2*symparts))     
   if(holepatch==1)write(100,'(I4,2ES15.6,I10,A21)')loop-1,somme,somme2,  &
                               (count3/(2*symparts)),'       Holes found...'
   if((somme<acc).or.(loop>nloop))then
     term=1!
   endif
   if(somme<3d0*acc.and.wellfocus<1.and.focus==0.and.restart==1) then
     wellfocus=1
     if(focus_onLR==0)then
       if (term==0) write(100,*)
       if (term==0) write(100,*) '--- switching focus to negative energies'
     endif
   endif
   OPEN (UNIT=652,FILE='PES-'//sys_label(1:nline),FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL',STATUS='REPLACE',POSITION='REWIND')
   write(652) count3
   write(652) order_1
   write(652) order_2
   write(652) order_3
   write(652) order_4
   write(652) maxpoints
   write(652) mass
   write(652) rmax
   write(652) rmin
   write(652) Max_E
   write(652) low_grid
   write(652) count_seed
   do i=1,count3
    write(652) b2(:,i)
   enddo
   do i=1,count3
    write(652) b2_lower(:,i)
   enddo
   do i=1,count3
    write(652) b2_minimal(:,i)
   enddo
   do i=1,count3
    write(652) d(i)
   enddo
   do i=1,count3
    write(652) coords(i,:)
   enddo
   if(low_grid>0)then!
     write(652) Max_E_seed
     do i=1,count_seed
      write(652) b2_seed(:,i)
     enddo
     do i=1,count_seed
      write(652) d_seed(i)
     enddo
     do i=1,count_seed
      write(652) coords_seed(i,:)      
     enddo
   endif
   close(652)
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_BCAST(term,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(wellfocus,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 IF((term>0).or.(only_seed_grid==1))THEN! 
   if(myid==0)then
     write(110,*)'The PES is complete! Preparing summary...'
     call date_and_time(bdate1,bdate2,bdate3,date_time)
     write(100,*)
     write(100,*)
     write(100,*)'                        MPI ENVIROMENT CLOSED                          '
     write(100,'(A74)')' -------------------------------------------------------------------------'
     write(100,'(A74)')' *           System date (dd/mm/yyyy): '//bdate1(7:8)//'/'//bdate1(5:6)// &
                        '/'//bdate1(1:4)//'                        *'
     write(100,'(A74)')' *           System time (hh:min:sec): '//bdate2(1:2)//':'//bdate2(3:4)// &
                        ':'//bdate2(5:6)//'                          *'
     write(100,'(A74)')' -------------------------------------------------------------------------'
     if(term>0)then
       write(100,*)'  Congratulations! A new Potential Energy Surface have been constructed  '
     else
       write(100,*)'                  Only the seed-grid-PES was computed.                   '
     endif
     if(loop>1)write(100,*)
     if(loop>1)write(100,*)'Total number of calculated ab initio points: ',tot_abinitio
     if(loop>1)write(100,*)'Lowest ab initio energy calculated:',xGlob_min(5)
     if(loop>1)write(100,*)'(intermolecular distance is in Angstroms, angles are in degrees) '
     if(loop>1)write(100,*)'  R  =',xGlob_min(1)
     if(loop>1)write(100,*)' th1 =',dacos(xGlob_min(2))*180.d0/pii
     if(loop>1)write(100,*)' th2 =',dacos(xGlob_min(3))*180.d0/pii
     if(loop>1)write(100,*)' phi =',xGlob_min(4)*180.d0/pii
   endif
   xrms_R1=0.d0
   xrms_R2=0.d0
   xrms_R3=0.d0
   xrms_R4=0.d0
   xMAPot=0.d0
   xR_count1=0.d0
   xR_count2=0.d0
   xR_count3=0.d0
   xR_count4=0.d0
   xGlob_min1(5)=2.d2
   do i=1,200000/numprocs
     27 call random_number(ran_vec)
     ran_vec(1)=rmin(1)+ran_vec(1)*(rmax(1)-rmin(1))
     ran_vec(2)=-0.999d0+ran_vec(2)*1.998d0
     ran_vec(3)=-0.999d0+ran_vec(3)*1.998d0
     ran_vec(4)=(-0.999d0+ran_vec(4)*1.998d0)*pii
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if(dist_flag==1) goto 27!
     xi=ran_vec(:)
     if(low_grid>0)then!
       temp3=func_actual_seed(xi)!
       if(temp3>Max_E_seed) goto 27
     endif
     if(subzero==0)then
       temp=func_actual(xi)
       if (temp>Max_E) goto 27
     else
       temp=func_actual(xi)+temp3
       if (temp>Max_E) goto 27
     endif
     if(temp<xGlob_min1(5))then
       xGlob_min1(1)=xi(1)
       xGlob_min1(2)=dacos(xi(2))*180.d0/pii
       xGlob_min1(3)=dacos(xi(3))*180.d0/pii
       xGlob_min1(4)=xi(4)*180.d0/pii
       xGlob_min1(5)=temp
     endif
     tampon=func1(xi)     
     do j=1,xnum_intervals_R
       if((xR(j)<=xi(1)).and.(xi(1)<xR(j+1)))then
         xMAPot(myid+1,j)=xMAPot(myid+1,j)+dabs(temp-E_asym)
         xR_count1(myid+1,j)=xR_count1(myid+1,j)+1.d0
         xrms_R1(myid+1,j)=xrms_R1(myid+1,j)+dabs(tampon)
         if(temp<E_asym+(0.05d0/hart2kcl*CONVE))then
           xR_count2(myid+1,j)=xR_count2(myid+1,j)+1.d0
           xrms_R2(myid+1,j)=xrms_R2(myid+1,j)+dabs(tampon)
         endif
         if(temp<((0.3/hart2kcl*CONVE)+Glob_min))then
           xR_count3(myid+1,j)=xR_count3(myid+1,j)+1.d0
           xrms_R3(myid+1,j)=xrms_R3(myid+1,j)+dabs(tampon)
         endif
         if(focus>0)then
           if(temp<(E_asym+(0.05d0/hart2kcl*CONVE)-increment))then
            xR_count4(myid+1,j)=xR_count4(myid+1,j)+1.d0
            xrms_R4(myid+1,j)=xrms_R4(myid+1,j)+dabs(tampon)
           endif
         endif
       endif
     enddo
   enddo
   do j=1,xnum_intervals_R
     xMAPot(myid+1,xnum_intervals_R+1)=xMAPot(myid+1,xnum_intervals_R+1)+xMAPot(myid+1,j)
     xrms_R1(myid+1,xnum_intervals_R+1)=xrms_R1(myid+1,xnum_intervals_R+1)+xrms_R1(myid+1,j)
     xrms_R2(myid+1,xnum_intervals_R+1)=xrms_R2(myid+1,xnum_intervals_R+1)+xrms_R2(myid+1,j)
     xrms_R3(myid+1,xnum_intervals_R+1)=xrms_R3(myid+1,xnum_intervals_R+1)+xrms_R3(myid+1,j)
     if(focus>0)xrms_R4(myid+1,xnum_intervals_R+1)=xrms_R4(myid+1,xnum_intervals_R+1)+xrms_R4(myid+1,j)
     xR_count1(myid+1,xnum_intervals_R+1)=xR_count1(myid+1,xnum_intervals_R+1)+xR_count1(myid+1,j)
     xR_count2(myid+1,xnum_intervals_R+1)=xR_count2(myid+1,xnum_intervals_R+1)+xR_count2(myid+1,j)
     xR_count3(myid+1,xnum_intervals_R+1)=xR_count3(myid+1,xnum_intervals_R+1)+xR_count3(myid+1,j)
     if(focus>0)xR_count4(myid+1,xnum_intervals_R+1)=xR_count4(myid+1,xnum_intervals_R+1)+xR_count4(myid+1,j)
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do i=0,numprocs-1
    do j=1,xnum_intervals_R+1
      call MPI_BCAST(xMAPot(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xrms_R1(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xrms_R2(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xrms_R3(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      if(focus>0)call MPI_BCAST(xrms_R4(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xR_count1(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xR_count2(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xR_count3(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
      if(focus>0)call MPI_BCAST(xR_count4(i+1,j),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
    enddo
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   if(myid==0)then
     write(xchar1,'(I8)')count3/(2*symparts)
     do j=1,xnum_intervals_R+1
       do i=1,numprocs-1
        xMAPot(1,j)=xMAPot(1,j)+xMAPot(1+i,j)
        xrms_R1(1,j)=xrms_R1(1,j)+xrms_R1(1+i,j)
        xrms_R2(1,j)=xrms_R2(1,j)+xrms_R2(1+i,j)
        xrms_R3(1,j)=xrms_R3(1,j)+xrms_R3(1+i,j)
        if(focus>0)xrms_R4(1,j)=xrms_R4(1,j)+xrms_R4(1+i,j)
        xR_count1(1,j)=xR_count1(1,j)+xR_count1(1+i,j)
        xR_count2(1,j)=xR_count2(1,j)+xR_count2(1+i,j)
        xR_count3(1,j)=xR_count3(1,j)+xR_count3(1+i,j)
        if(focus>0)xR_count4(1,j)=xR_count4(1,j)+xR_count4(1+i,j)
       enddo
       if(xR_count1(1,j)==0)then
         xMAPot(1,j)=0.d0
         xrms_R1(1,j)=0.d0
       else
         xMAPot(1,j)=xMAPot(1,j)/xR_count1(1,j)
         xrms_R1(1,j)=dsqrt(xrms_R1(1,j)/xR_count1(1,j))
       endif
       if(xR_count2(1,j)==0)then
         xrms_R2(1,j)=0.d0
       else
         xrms_R2(1,j)=dsqrt(xrms_R2(1,j)/xR_count2(1,j))
       endif
       if(xR_count3(1,j)==0)then
         xrms_R3(1,j)=0.d0
       else
         xrms_R3(1,j)=dsqrt(xrms_R3(1,j)/xR_count3(1,j))
       endif
       if(focus>0)then
         if(xR_count4(1,j)==0)then
           xrms_R4(1,j)=0.d0
         else
           xrms_R4(1,j)=dsqrt(xrms_R4(1,j)/xR_count4(1,j))
         endif
       endif
     enddo
     write(100,*)
     write(100,*)'Summary of errors:'
     write(100,*)'RMS1 = Global'
     write(100,*)'RMS2 = below asymptote' 
     write(100,*)'RMS3 = 0.3 kcal/mol (~100 cm-1) above minimum' 
     if(focus>0.and.focus_onLR==0)write(100,*)'RMS4 = focus energy-range'
     write(100,*)
     if(focus>0.and.focus_onLR==0)then
       write(100,'(A125)')'   Ri    Rf    <|V(R)|>     RMS1      %       pts    RMS2      %       &
                            pts    RMS3      %       pts    RMS4      %       pts'
     else
       write(100,'(A100)')'   Ri    Rf    <|V(R)|>      RMS1      %       pts    RMS2      %       &
                            pts    RMS3      %       pts'
     endif
     do j=1,xnum_intervals_R
       if(xMAPot(1,j)==0)then
         x1=0.d0
         x2=0.d0
         x3=0.d0
         x4=0.d0
       else
         x1=xrms_R1(1,j)*100.d0/xMAPot(1,j)
         x2=xrms_R2(1,j)*100.d0/xMAPot(1,j)
         x3=xrms_R3(1,j)*100.d0/xMAPot(1,j)
         x4=xrms_R4(1,j)*100.d0/xMAPot(1,j)
       endif
       x5=(xR(j)+xR(j+1))/2.d0
       IF(focus>0.and.focus_onLR==0)THEN
         if((xR(j)<=minR).and.(minR<xR(j+1)))write(100,'(A125)')'-----------------------------------&
         ------------------------------------------------------------------------------------------'
         write(100,'(2F6.2,ES13.5,4(ES11.3,F6.2,I8))')xR(j),xR(j+1),xMAPot(1,j),&
          xrms_R1(1,j),x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j)),&
          xrms_R4(1,j),x4,int(xR_count4(1,j))
         if((xR(j)<maxR).and.(maxR<=xR(j+1)))write(100,'(A125)')'-----------------------------------&
         ------------------------------------------------------------------------------------------'
       ELSE
         if((xR(j)<=minR).and.(minR<xR(j+1)))write(100,'(A101)')'-------------------------&
         ---------------------------------------------------------------------------'
         write(100,'(2F6.2,ES13.5,3(ES11.3,F6.2,I8))')xR(j),xR(j+1),xMAPot(1,j),xrms_R1(1,j),&
          x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j))
         if((xR(j)<maxR).and.(maxR<=xR(j+1)))write(100,'(A101)')'-------------------------&
         ---------------------------------------------------------------------------'
       ENDIF
     enddo
     write(100,*)
     x1=xrms_R1(1,xnum_intervals_R+1)*100.d0/xMAPot(1,xnum_intervals_R+1)
     x2=xrms_R2(1,xnum_intervals_R+1)*100.d0/xMAPot(1,xnum_intervals_R+1)
     x3=xrms_R3(1,xnum_intervals_R+1)*100.d0/xMAPot(1,xnum_intervals_R+1)
     x4=xrms_R4(1,xnum_intervals_R+1)*100.d0/xMAPot(1,xnum_intervals_R+1)
     if(focus>0.and.focus_onLR==0)then
       write(100,'(A12,ES13.5,4(ES11.3,F6.2,I8))')'total:',xMAPot(1,xnum_intervals_R+1), &
       xrms_R1(1,xnum_intervals_R+1),x1,int(xR_count1(1,xnum_intervals_R+1)), &
       xrms_R2(1,xnum_intervals_R+1),x2,int(xR_count2(1,xnum_intervals_R+1)), &
       xrms_R3(1,xnum_intervals_R+1),x3,int(xR_count3(1,xnum_intervals_R+1)), &
       xrms_R4(1,xnum_intervals_R+1),x4,int(xR_count4(1,xnum_intervals_R+1))
     else
       write(100,'(A12,ES13.5,3(ES11.3,F6.2,I8))')'total:',xMAPot(1,xnum_intervals_R+1), &
       xrms_R1(1,xnum_intervals_R+1),x1,int(xR_count1(1,xnum_intervals_R+1)), &
       xrms_R2(1,xnum_intervals_R+1),x2,int(xR_count2(1,xnum_intervals_R+1)), &
       xrms_R3(1,xnum_intervals_R+1),x3,int(xR_count3(1,xnum_intervals_R+1))
     endif
     write(100,*)
     write(100,*)'Lowest energy found in the random sampling:',xGlob_min1(5)
     write(100,*)' (intermolecular distance is in Angstroms, angles are in degrees) '
     write(100,*)'  R  =',xGlob_min1(1)
     write(100,*)' th1 =',xGlob_min1(2)
     write(100,*)' th2 =',xGlob_min1(3)
     write(100,*)' phi =',xGlob_min1(4)
     close(100)
     call system('rm fort*')
     call system('rm pun*')
     write(xchar1,'(I8)')count3/(2*symparts)
     call system('mv PES-'//sys_label(1:nline)//' PES-'//sys_label(1:nline)//'-'//trim(adjustl(xchar1)))
     inquire(file='output.dat',exist=logica1)
     if(logica1)call system('mv output.dat output-'//trim(adjustl(xchar1))//'.dat')
   endif
   call MPI_FINALIZE(rcc)
   stop
 ENDIF
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 IF(focus_onLR==1)THEN! 
   if(myid==0)then! 
    ncont=0
    do i=1,numadded
     28 call sobseq(ran_vec)
     range=rmax(1)-xminLR!
     ran_vec(1)=xminLR+range*ran_vec(1)!
     do j=2,4
       range = rmax(j) - rmin(j)
       ran_vec(j)=rmin(j)+range*ran_vec(j)
     enddo
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if (dist_flag==1) ncont=ncont+1
     if (dist_flag==1) goto 28!
     jac=ran_vec(:)
     if(low_grid>0)then!
       temp3=func_actual_seed(jac)!
       if(temp3>Max_E_seed) ncont=ncont+1
       if(temp3>Max_E_seed) goto 28
     endif
     if(subzero==0)then
       temp=func_actual_min(jac)!
       if (temp>Max_E) ncont=ncont+1
       if (temp>Max_E) goto 28
       temp=func_actual(jac)
       if (temp>Max_E) ncont=ncont+1
       if (temp>Max_E) goto 28
     else
       temp=func_actual_min(jac)!
       if (temp+temp3>Max_E) ncont=ncont+1
       if (temp+temp3>Max_E) goto 28
       temp=func_actual(jac)
       if (temp+temp3>Max_E) ncont=ncont+1
       if (temp+temp3>Max_E) goto 28
     endif
     seed(i,:)=jac
    enddo
    write(110,*)'*** LONG RANGE...'
    if(ncont/=0)write(110,*)'Skipped geometries: ',ncont
    write(110,*)'Selected geometries (and estimated energy)'
    do kk=1,numadded
      write(110,'(I10,4f14.8,3x,f20.8)') kk,seed(kk,:),func_actual(seed(kk,:))
    enddo
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do i=1,numadded
     call MPI_BCAST(seed(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   write(charmyid,'(I5)')myid
   open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                             '.dat',action='write')
   do i=myid+1,numadded,numprocs
     jac=seed(i,:)
     pass_check=1
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     804 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,pass_check,midbond,jac(1))
     poten=poten*CONVE
     seed_pot(i)=poten
     if(ab_flag==2)then
      dcart=dcart*ugrad 
      seed_grad(i,:)=dcart 
      write(1000+myid,f602)numpoints_high_count+i,jac(:),seed_pot(i),seed_grad(i,:)
     else
      write(1000+myid,601)numpoints_high_count+i,jac(:),seed_pot(i)
     endif
   enddo
   close(1000+myid)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do j=0,numprocs-1
     do i=j+1,numadded,numprocs
       call MPI_BCAST(seed_pot(i),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
       if(ab_flag==2)then
        call MPI_BCAST(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
       endif
     enddo
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   numadded_act=numadded
   if(holepatch==1)holepatch=0
   GOTO 189!
 ENDIF!
 if(ntest/=0)then!
 IF(mod(loop,ntest)==ntest-1)THEN!
   rms_error=0d0
   mean_error=0d0
   if(myid==0)then!
    do i=1,num_ab_err_points
     20 call random_number(ran_vec)
     if(focus_onLR==1)then
       range=rmax(1)-xminLR!
       ran_vec(1)=xminLR+range*ran_vec(1)
     else
       range=maxR-minR!
       ran_vec(1)=minR+range*ran_vec(1)
     endif
     do j=2,4
       range=rmax(j)-rmin(j)
       ran_vec(j)=rmin(j)+range*ran_vec(j)
     enddo
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if (dist_flag==1) goto 20!
     jac=ran_vec(:)
     if(low_grid>0)then!
       temp3=func_actual_seed(jac)
       if (temp3>Max_E_seed) goto 20
     endif
     if(focus_onLR==1)goto 192!
     if(focus>0)then! 
       if (func_actual(jac)>E_asym+(0.05d0/hart2kcl*CONVE)-increment) goto 20
     else
       if(wellfocus>0)then!
        if (func_actual(jac)>E_asym+(0.05d0/hart2kcl*CONVE)) goto 20
       endif
     endif
     192 temp=func_actual_min(jac)!
     if(subzero==0)then
       if (temp>Max_E) goto 20!
     else
       if (temp+temp3>Max_E) goto 20
     endif
     seed(i,:)=jac
    enddo
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do i=1,num_ab_err_points
     call MPI_BCAST(seed(i,:),4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   write(charmyid,'(I5)')myid
   open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                             '.dat',action='write')
   do i=myid+1,num_ab_err_points,numprocs
     jac=seed(i,:)
     pass_check=1
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     803 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,pass_check,midbond,jac(1))
     poten=poten*CONVE
     seed_pot(i)=poten
     if(ab_flag==2)then
      dcart=dcart*ugrad 
      seed_grad(i,:)=dcart 
      write(1000+myid,f602)numpoints_high_count+i,jac(:),seed_pot(i),seed_grad(i,:)
     else
      write(1000+myid,601)numpoints_high_count+i,jac(:),seed_pot(i)
     endif
     temp3=func_actual_seed(jac)
     if (subzero==0) tampon=func_actual(jac)
     if (subzero==1) tampon=func_actual(jac)+temp3
     rms_error(myid+1)=rms_error(myid+1)+(tampon-poten)**2
     mean_error(myid+1)=mean_error(myid+1)+dabs(tampon-poten)
   enddo
   close(1000+myid)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do i=0,numprocs-1
     call MPI_BCAST(rms_error(i+1),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(mean_error(i+1),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   do i=1,numprocs-1
     rms_error(1)=rms_error(1)+rms_error(1+i)
     mean_error(1)=mean_error(1)+mean_error(1+i)
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   if(myid==0) then
     rms_error(1)=dsqrt(rms_error(1)/dble(num_ab_err_points))
     mean_error(1)=mean_error(1)/dble(num_ab_err_points)
     somme=rms_error(1)
     somme2=mean_error(1)
     write(xchar,'(I8)')count3
     write(xchar1,'(I8)')count3/(2*symparts)
     write(100,'(A4,2ES15.6,I10,A21)')'REAL',somme,somme2,num_ab_err_points,  &
           '   '//trim(adjustl(xchar1))//'(+sym. '//trim(adjustl(xchar))//')'
   endif
   do j=0,numprocs-1
     do i=j+1,num_ab_err_points,numprocs
       call MPI_BCAST(seed_pot(i),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
       if(ab_flag==2)then
        call MPI_BCAST(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
       endif
     enddo
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   numadded_act=num_ab_err_points
   if(holepatch==1)holepatch=0
   GOTO 189!
 ENDIF! 
 endif
 allocate(ind7(count3),ind8(count3))
 count7=0
 do i=1,count3,2*symparts
   do ip=1,count3
     call dist_metric(coords(i,:),coords(ip,:),W_a,ind7(ip))
   enddo
   call indexxy(count3,ind7,ind8)
   do j=2,5!
     if(ind8(j)>i)then
       xtemp=0.5d0*(coords(i,1)+coords(ind8(j),1))
       if((focus_onLR==1).and.(xtemp>=xminLR))then!
         count7=count7+1
         do k=1,4
          test_points(count7,k)=0.5d0*(coords(i,k)+coords(ind8(j),k))
         enddo
       endif
       if(focus_onLR==0)then
         if((xtemp>=minR).and.(xtemp<=maxR))then!
           do k=1,4
             xi(k)=0.5d0*(coords(i,k)+coords(ind8(j),k))
           enddo
           if(focus>0)then
             if(func_actual(xi)<E_asym+(0.05d0/hart2kcl*CONVE)-increment)then
               count7=count7+1
               do k=1,4
                 test_points(count7,k)=xi(k)
               enddo
             endif
           else
             if(wellfocus>0)then
               if(func_actual(xi)<E_asym+(0.05d0/hart2kcl*CONVE))then
                 count7=count7+1
                 do k=1,4
                   test_points(count7,k)=xi(k)
                 enddo
               endif
             else
               count7=count7+1
               do k=1,4
                 test_points(count7,k)=xi(k)
               enddo
             endif
           endif
         endif
       endif
     endif
   enddo
 enddo
 test_points2=test_points
 count=0
 do kk=1,count7
   count=count+1
   test_points2(count,:)=test_points(kk,:)
   do jj=1,count-1
     if(dsqrt((test_points2(count,1)-test_points2(count-jj,1))**2+(test_points2(count,2)- &
       test_points2(count-jj,2))**2+(test_points2(count,3)-test_points2(count-jj,3))**2+ &
       (test_points2(count,4)-test_points2(count-jj,4))**2)<3d-3)then
       count=count-1
       goto 1777
     endif
   enddo
 1777 enddo
 count7=count
 test_points=test_points2
 call indexxy(num_err_points1*numprocs,error_points(:,5),ind6)
 goto 1778!
 if(myid==0)then
   ncont1=100*numprocs
   if (num_err_points1<100) ncont1=num_err_points1*numprocs
   do i=1,ncont1
     write(820,'(4f10.3,1f20.3)') error_points(ind6(num_err_points1*numprocs+1-i),:)
   enddo
   write(820,*)
 endif
 1778 continue
 ncont1=15*numadded
 if (ncont1>num_err_points1*numprocs) ncont1=num_err_points1*numprocs
 do i=1,ncont1
   count7=count7+1
   test_points(count7,1:4)=error_points(ind6(num_err_points1*numprocs+1-i),1:4)
 enddo
 deallocate(ind7,ind8)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if(myid==0)write(110,*) 'Num. minimization starting points:' , count7
 write(3000+myid,*)
 write(3000+myid,*)'Loop: ',loop-1
 gtol=1d-4
 do i=myid+1,count7,numprocs      
   reiter=5
   iter=5
   xi=test_points(i,:)
   sqrdiff=dabs(func(xi))
   write(3000+myid,'(4f10.4)') xi
   write(3000+myid,*) i,sqrdiff
   if(dabs(func(xi))>1d-4*acc/hart2kcl*CONVE/hart2kcl*CONVE)then
     call frprmn(xi,gtol,iter,valeur)
     sqrdiff=dabs(valeur)
   else
     sqrdiff=0d0
   endif
   test_points(i,:)=xi
   order_diff(i)=sqrdiff      
   write(3000+myid,'(4f10.4)') xi      
   write(3000+myid,*) i,sqrdiff
   write(3000+myid,*)      
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do j=0,numprocs-1   
    do i=j+1,count7,numprocs
       call MPI_BCAST(test_points(i,:),4,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(order_diff(i),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
    enddo
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 allocate(ind3(count7),ind4(count7))
 ind3(1:count7)=order_diff(1:count7)
 call indexxy(count7,ind3,ind4)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 numadded_act=numadded
 skip=0
 count=0
 do kk=1,numadded     
   count=count+1
   177 continue
   if(skip>count7-numadded) then
     if(myid==0)write(110,*) 'Warning! Not enough unique-points...',kk
     count=count-1
     numadded_act=kk-1
     goto 188
   endif
   seed(count,:)=test_points(ind4(count7+1-skip-kk),:)
   if(holepatch==1.and.kk==numadded)then
     seed(count,:)=hole(indhole(1),:)
     holepatch=0
   endif
   call INT_Cart(cart,seed(count,:),mass,natom1,natom2,ref1,ref2)
   call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
   if(dist_flag==1) then
     skip=skip+1
     goto 177
   endif
   do j=1,4
     if(seed(count,j)>rmax(j).or.seed(count,j)<rmin(j))then
       skip=skip+1
       goto 177
     endif
   enddo
   IF(focus_onLR==0)THEN
     if(focus>0)then
       if(func_actual(seed(count,:))>E_asym+(0.05d0/hart2kcl*CONVE)-increment)then
         skip=skip+1
         goto 177
       endif
     else
       if(wellfocus>0)then
         if(func_actual(seed(count,:))>E_asym+(0.05d0/hart2kcl*CONVE))then
           skip=skip+1
           goto 177
         endif
       endif
     endif
   ENDIF
   deallocate(ind,ind2)
   allocate(ind(count),ind2(count))
   local_d=0d0
   if(count>1)then
     do ip=1,count!
       call dist_metric(seed(count,:),seed(ip,:),W_a,ind(ip))
     enddo
     call indexxy(count,ind,ind2)
     local_d=ind(ind2(2)) !
     if(local_d<3d-1) then
       skip=skip+1
       goto 177
     endif
   endif
   NAME1='AbINITIO.dat'
   call search_abinitio_dat(seed(count,:),maxpoints,4,NAME1,ab_flag,natom,dist_flag)
   if(dist_flag==1) then
     skip=skip+1
     goto 177
   endif
 enddo
 188 continue
 if (numadded_act<numadded) write(110,*)'Warning! the number of points added is lower than "numadded"'
 if(myid==0) then! 
   write(110,*)'Selected geometries (and estimated energy)'
   do kk=1,numadded_act
    write(110,'(I10,4f14.8,3x,f20.8)') kk,seed(kk,:),func_actual(seed(kk,:))
   enddo
 endif
 write(charmyid,'(I5)')myid
 open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                           '.dat',action='write')
 do kk=myid+1,numadded_act,numprocs!
   jac=seed(count+1-kk,:)
   pass_check=1
   call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
   802 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,pass_check,midbond,jac(1))
   if(poten>199d0)then
     if(pass_check<2) then
       pass_check=pass_check+1
       goto 802
     endif
   endif
   poten=poten*CONVE
   seed_pot(count+1-kk)=poten
   if(ab_flag==2)then
     dcart=dcart*ugrad 
     seed_grad(count+1-kk,:)=dcart
     write(1000+myid,*)numpoints_high_count+kk,jac(:),seed_pot(count+1-kk),seed_grad(count+1-kk,:)
   else
     write(1000+myid,*)numpoints_high_count+kk,jac(:),seed_pot(count+1-kk)
   endif
 enddo
 close(1000+myid)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do j=0,numprocs-1
   do i=j+1,numadded_act,numprocs
     if(ab_flag==2)then
       call MPI_BCAST(seed_grad(count+1-i,:),3*(natom),MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     endif
     call MPI_BCAST(seed_pot(count+1-i),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(seed(count+1-i,:),4,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
   enddo
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 deallocate(ind3,ind4)
 189 CONTINUE! 
 if(myid==0)then
   open(unit=222,file='AbINITIO.dat',status='old')
   do i=1,maxpoints
      read(222,*,end=31)
   enddo
   31 backspace(unit=222)
   do j=1,numprocs
     write(charmyid,'(I5)')j-1
     open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                         '.dat', action='read',status='old')
     do i=1,numadded_act
       if(ab_flag==2)then
         read(230,*,end=32)ncont1,jac(:),xpot,xgrad(:)
         write(222,f602)ncont1,jac(:),xpot,xgrad(:)
       else
         read(230,*,end=32)ncont1,jac(:),xpot
         write(222,601)ncont1,jac(:),xpot
       endif
     enddo
     32 close(230, status="delete")
   enddo
   close(222)
 endif
 count=count3
 do i=1,numadded_act
   if(ab_flag==2)then
      dcart=seed_grad(i,:)
   endif
   call symmetry2(seed(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,  &
                   ab_flag,exch,flip1,flip2)
   if(seed_pot(i)<E_limit)then! 
     do k2=0,1!
      do k=1,symparts
        count=count+1
        coords(count,:)=internal(k,:)
        coords(count,4)=internal(k,4)*(-1d0)**k2
        pot(count)=seed_pot(i)
        if(ab_flag==2.and.subzero==0)then
          grad(count,:)=grad_int(k,:)
          grad(count,4)=grad_int(k,4)*(-1d0)**k2
        endif
        if(subzero==1)then
          xi=coords(count,:)
          temp2=dfunc_actual_seed(xi)
          pot(count)=seed_pot(i)-temp2(1)
          if(ab_flag==2)then
            do ip=1,4
             grad(count,ip)=grad_int(k,ip)-temp2(1+ip)
            enddo
          endif
        endif
      enddo
     enddo
   endif
 enddo
 count3=count
 if(myid==0)then
   ncont=0
   open(unit=222,file='AbINITIO.dat',status='old')
   do i=1,100000
     read(222,*,end=33)j,xx(i),th1(i),th2(i),phi(i),eee(i)
     ncont=ncont+1
   enddo
   33 continue
   tot_abinitio=ncont! 
   ncont1=0
   do i=1,ncont
     do j=i+1,ncont
       x1=eee(i)-eee(j)
       if((xx(i)==xx(j)).and.(dabs(x1)<=1.d-7))then
         ncont1=ncont1+1
       endif
     enddo
   enddo
   if(ncont1/=0)write(110,*)'***********'
   if(ncont1/=0)write(110,*)'* Warning * there are ',ncont1,' geometries with the same "R" and "E" in "AbINITIO.dat" !!'
   if(ncont1/=0)write(110,*)'***********'
 close(222)
 endif
 dynmax=1d2
 dynmin=-1d9
 if(loop==1)Min_E=-1d9!
 Max_E=2d0
 ncont=0
 ncont1=0
 ncont2=0
 do i=1,count3
   if(subzero==0)then
     if(pot(i)<Max_E)then
      Max_E=pot(i)
      xGlob_min(1)=coords(i,1)
      xGlob_min(2)=coords(i,2)
      xGlob_min(3)=coords(i,3)
      xGlob_min(4)=coords(i,4)
      xGlob_min(5)=Max_E
     endif
     if(pot(i)>Min_E)then
       ncont1=ncont1+1
       Min_E=pot(i)
     endif
     if(coords(i,1)>Max_R)then
      ncont2=ncont2+1
      Max_R=coords(i,1)! 
      E_asym=pot(i)! 
     endif
   else
     xi=coords(i,:)
     temp=func_actual_seed(xi)
     if(pot(i)+temp<Max_E)then
      Max_E=pot(i)+temp
      xGlob_min(1)=coords(i,1)
      xGlob_min(2)=coords(i,2)
      xGlob_min(3)=coords(i,3)
      xGlob_min(4)=coords(i,4)
      xGlob_min(5)=Max_E
     endif
     if(pot(i)+temp>Min_E)then
       ncont1=ncont1+1
       Min_E=pot(i)+temp
     endif
     if (pot(i)<dynmax) dynmax=pot(i)
     if (pot(i)>dynmin) dynmin=pot(i)
     if(coords(i,1)>Max_R) then
      ncont2=ncont2+1
      Max_R=coords(i,1)!
      E_asym=pot(i)+temp!
     endif
   endif
 enddo
 if(Max_E<Glob_min)then
   ncont=ncont+1
   Glob_min=Max_E
 endif
 if((myid==0).and.(ncont>0))write(110,*)' ** new lowest E:',Glob_min
 if((myid==0).and.(ncont>0))write(110,'(A18,4F10.4)')'    coordinates:',xGlob_min(1:4)
 if((myid==0).and.(ncont1>0).and.(loop==1))write(110,*)'        highest E:',Min_E
 if((myid==0).and.(ncont1>0).and.(loop/=1))write(110,*)' ** new highest E:',Min_E
 if((myid==0).and.(ncont2>0))write(110,*)' ** new asymptote:',E_asym
 Max_E=E_asym+E_range
 E_limit=E_asym+1.2d0*E_range!
 if(myid==0) then
   if(ncont2>0)write(110,*)' ** max. E above asymptote:',Max_E
   if(subzero==1)then
     write(110,*) 'dynamic range',dynmax,dynmin
   endif
 endif
 if(myid==0)write(110,*)' ** Points used to fit:',count3/2/symparts
 if(myid==0)write(110,*)'   + symmetry partners:',count3
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 deallocate(ind,ind2)
 allocate(ind(count3),ind2(count3))
 do i=1,count3  
   count=0 
   do ip=1,count3         
     count=count+1
     call dist_metric(coords(i,:),coords(ip,:),W_a,ind(count))
   enddo
   call indexxy(count3,ind,ind2)
   d(i)=ind(ind2(zz4))
 enddo
 deallocate(ind,ind2)
 allocate(ind(count3),ind2(count3))
 support=min(count3,4*basis_1/(4*(ab_flag-1)+1)) 
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do i=myid+1,count3,numprocs  
   count=0
   Jac=coords(i,:) 
   do ip=1,count3
     count=count+1
     Jac2=coords(ip,:)
     call dist_metric(jac,jac2,W_a,somme)
     somme=somme**2
     ind(count)=dsqrt(exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss))      
   enddo
   call indexxy(count3,ind,ind2)
   call basis_calc(i,ind,ind2)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 do j=0,numprocs-1
   do i=j+1,count3,numprocs
     call MPI_BCAST(b2(:,i),basis_1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(b2_lower(:,i),basis_2,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(b2_minimal(:,i),basis_3,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
   enddo
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
ENDDO! 
if(myid==0)then
  call date_and_time(bdate1,bdate2,bdate3,date_time)
  write(100,*)
  write(100,*)
  write(100,*)'                        MPI ENVIROMENT CLOSED                          '
  write(100,'(A74)')' -------------------------------------------------------------------------'
  write(100,'(A74)')' *           System date (dd/mm/yyyy): '//bdate1(7:8)//'/'//bdate1(5:6)// &
                     '/'//bdate1(1:4)//'                        *'
  write(100,'(A74)')' *           System time (hh:min:sec): '//bdate2(1:2)//':'//bdate2(3:4)// &
                     ':'//bdate2(5:6)//'                          *'
  write(100,'(A74)')' -------------------------------------------------------------------------'
  write(100,'(A74)')'             The maximum number of cycles have been reached !             '
  write(100,'(A74)')'To continue, increase the value of "maxpoints" in the input file & restart'
  write(100,*)
  write(100,*)'Total number of calculated ab initio points: ',tot_abinitio
  write(100,*)'Lowest ab initio energy calculated:',xGlob_min(5)
  write(100,*)'(intermolecular distance is in Angstroms, angles are in degrees) '
  write(100,*)'  R  =',xGlob_min(1)
  write(100,*)' th1 =',dacos(xGlob_min(2))*180.d0/pii
  write(100,*)' th2 =',dacos(xGlob_min(3))*180.d0/pii
  write(100,*)' phi =',xGlob_min(4)*180.d0/pii
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(rcc)

stop
END PROGRAM AUTOSURF_PES_rigid4D
