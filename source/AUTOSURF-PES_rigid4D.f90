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
!      Dependencies: "FUNCTIONS.f90", "MODULES.f90", "diag.f", "SUBROUTINES.f90"    !
!      Version: 1.05
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
             holepatch,term,wellfocus
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
 integer,allocatable :: ind2(:),ind2_seed(:),ind4(:),ind6(:),ind8(:),ind10(:),   &
             ind12(:),ind14(:),indhole(:)
 real*8 :: xi(4),tampon,pii
 real*8 :: temp,temp3,temp4(4),temp2(5),ran_vec(4),jac3(4),rann,bias,dynmin,dynmax
 real*8 :: somme,somme2,increment
 real*8 :: abs_diff,range,valeur,h2wn,dist,dist2,gtol,sqrdiff,local_d,atobohr,     &
            Min_E,Max_R,E_asym,Max_E3,bot_seed,low_high_fac,seed_cut
!5
 character(len=300) :: sys_label,line,f602
 character(len=10) :: bdate1,bdate2,bdate3
 character(len=5) :: charmyid
 character(len=8) :: xchar,xchar1
 character(len=15),dimension(4) :: charunitsE
 character(len=15),dimension(2) :: charunitsD
 integer :: nline,summyid,ncont,ncont1,ncont2,nid,numpoints_low_count,old_numprocs,focus_onR
 integer :: xpass,numpoints_high_count,n_test,xunitE,focus,only_seed_grid,status,focus_onLR
 integer :: tot_abinitio,ntest,auto_mod
 integer,dimension(1) :: xseed,xseed1
 integer,dimension(8) :: date_time
 real*8 :: xpot,CONVE,CONVD,xdR,minR,maxR,xtemp,x1,x2,x3,x4,x5,xGlob_min(5),xGlob_min1(5),xminLR
 real*8,allocatable :: xgrad(:),xR(:),xMAPot(:,:)
 real*8,allocatable :: xrms_R1(:,:),xrms_R2(:,:),xrms_R3(:,:),xrms_R4(:,:)
 real*8,allocatable :: xR_count1(:,:),xR_count2(:,:),xR_count3(:,:),xR_count4(:,:)
 logical :: logica1, logica2
 ! conversion factors
 real*8,parameter :: hart2wn=219474.6313702d0, hart2meV=27211.38602d0
 real*8,parameter :: ang2bohr=1.0d0/0.529177249d0, bohr2ang=0.529177249d0
 real*8,parameter :: hart2kcl=4.359744650D-18/4184*6.022140857D23
! real*8,parameter :: hart2kcl=627.5095d0
 integer,parameter :: xunitD=2, xnum_intervals_R=20
 data charunitsE/'hartree','kcal/mol','wave numbers','meV'/
 data charunitsD/'Bohr','Angstroms'/
 real*8, dimension(400000) :: xx,eee,th1,th2,phi
 integer,parameter :: send_data_tag = 2001, return_data_tag = 2002

!-----------------------------------------------------------------------------------
! Interface blocks
!-----------------------------------------------------------------------------------
INTERFACE! Negative of squared difference surface
  FUNCTION func(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi     
    REAL(SP) :: func
  END FUNCTION func
end interface
INTERFACE! Ibidem. "func", but defined in ALL configuration space
  FUNCTION func1(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi     
    REAL(SP) :: func1
  END FUNCTION func1
end interface
INTERFACE! Energy & analytic gradient of largest basis and high-level ab initio
  FUNCTION dfunc_actual_anal1(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_anal1
  END FUNCTION dfunc_actual_anal1
end interface
INTERFACE! Energy & analytic gradient of minimal basis and low-level ab initio
  FUNCTION dfunc_actual_seed(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_seed
  END FUNCTION dfunc_actual_seed
end interface
INTERFACE! Energy & analytic gradient of secondary basis and high-level ab initio
  FUNCTION dfunc_actual_anal2(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP),dimension(size(xi)+1) :: dfunc_actual_anal2
  END FUNCTION dfunc_actual_anal2
end interface
INTERFACE! Energy of largest basis and high-level ab initio
  FUNCTION func_actual(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual
  END FUNCTION func_actual
end interface
INTERFACE! Energy of minimal basis and high-level ab initio
  FUNCTION func_actual_min(xi)   
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_min
  END FUNCTION func_actual_min
end interface
INTERFACE! Energy of minimal basis and low-level ab initio
  FUNCTION func_actual_seed(xi)  
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_seed
  END FUNCTION func_actual_seed
end interface
INTERFACE! Analytic gradient of negative of squared difference surface
          ! 'func' and 'dfunc' must be what is used by canned CJ-minimization code
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

! check if input file exist...
inquire(file='input-AUTOSURF-PES.dat',exist=logica1)
if(.not.logica1)then
 write(*,*)
 write(*,*)'ERROR: The file: input-AUTOSURF-PES.dat does not exist !! '
 stop 
endif

open(unit=10,file='input-AUTOSURF-PES.dat',status='old',action='read')
111 read(10,'(A100)')line
ncont1=INDEX(line,'GENERAL INFORMATION:')
if(ncont1==0)goto 111
read(10,*)
read(10,*)sys_label ! System LABEL 
nline=scan(sys_label,' ')-1
if (nline>15) nline=15
read(10,*) restart  ! 1=yes, 0=no
read(10,*) xseed    ! seed for the random-number-generator algorithm
read(10,*) xunitE   ! PES units, Energies
close(10)

! Initialize the Sobol sequence
call sobseq(ran_vec,4)
do i=1,11
  call sobseq(ran_vec)
enddo
if(xseed(1)==0)then! init. random sequence using the time the calculation started
  xseed1=date_time(5)+date_time(6)+date_time(7)+date_time(8)
  if(restart==0)xseed1=10
  CALL RANDOM_SEED(PUT=xseed1) 
endif

!***********************************************************************************
!                       INITIALIZATION OF THE  MPI ENVIRONMENT                       
!***********************************************************************************
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)! find out MY process ID..
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)! ..and total number of processes
write(*,*) 'Process ', myid, ' of ', numprocs, ' is alive'

if(myid==0)then
  if(restart==1)then
    ! check if the GENERAL file exist  ---------------------------------------------
    inquire(file='GENERAL-'//sys_label(1:nline)//'.dat',exist=logica1)
    if(.not.logica1)then
     write(*,*)
     write(*,*)'ERROR: The file: GENERAL-'//sys_label(1:nline)//'.dat does not exist !! '
     write(*,*)'Value for variable "restart" should be changed to zero in the input file...'
     stop 
    endif
    open(unit=100,file='GENERAL-'//sys_label(1:nline)//'.dat',status='old')
    11 read(100,*,end=12)
    goto 11
    12 backspace(unit=100)
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
    if(xseed(1)/=0)then! reset the random sequence using the time the calculation was restarted
      xseed=date_time(5)+date_time(6)+date_time(7)+date_time(8)!seed = hour + min + sec + mili-sec
      CALL RANDOM_SEED(PUT=xseed) 
    endif
  else! (in the case of a fresh-start calculation)
    ! Open the GENERAL output file  ------------------------------------------------
    open(unit=100,file='GENERAL-'//sys_label(1:nline)//'.dat', action='write')
    write(100,*)'System:  '//sys_label(1:nline)//' '
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
    ! Set the generated random-number-sequence  ------------------------------------
    if(xseed(1)/=0)CALL RANDOM_SEED(PUT=xseed)
    ! Remove output files (if any) from previous calculations  ---------------------
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
  write(100,*)'Energies are in '//trim(adjustl(charunitsE(xunitE)))// &
               ', distances are in '//trim(adjustl(charunitsD(xunitD)))//'.'
  write(100,'(A30,I6)')' Total number of processes:   ',numprocs
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if (myid==0) then
  open(unit=110,file='output.dat')
  write(110,*) 'All processes are alive!'
endif

!-----------------------------------------------------------------------------------
!     Read input file
!-----------------------------------------------------------------------------------
open(unit=myid+1000,file='input-AUTOSURF-PES.dat',status='old',action='read')
112 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
if (ncont1==0) goto 112
read(myid+1000,*)
!# FRAGMENTS INFORMATION:
read(myid+1000,*) exch  ! are the two fragments identical? 1=yes, 0=no
read(myid+1000,*) flip1 ! frag 1 symmetric upon 180 degree flip (e.g., CO2)? 1=yes, 0=no
read(myid+1000,*) flip2 ! frag 2 symmetric upon 180 degree flip (e.g., CO2)? 1=yes, 0=no
read(myid+1000,*) natom1! number of atoms in frag1
read(myid+1000,*) natom2! number of atoms in frag2
natom=natom1+natom2 ! total number of atoms
nbdist=natom*(natom-1)/2! number of internuclear distances
allocate(ref1(3*natom1),ref2(3*natom2))
allocate(symb(natom),mass(natom))
do i=1,natom
 read(myid+1000,*) symb(i)! element labels of all atoms
enddo
do i=1,natom
 read(myid+1000,*) mass(i)! masses of all atoms
enddo  
do i=1,3*natom1
 read(myid+1000,*) ref1(i)! Cartesian positions for fragment 1 atoms
enddo
do i=1,3*natom2
 read(myid+1000,*) ref2(i)! Cartesian positions for fragment 2 atoms
enddo

113 read(myid+1000,'(A100)')line
ncont1=INDEX(line,'CODE CONTROL:')
if (ncont1==0) goto 113
read(myid+1000,*)
!# CODE CONTROL:
read(myid+1000,*) acc! accuracy target in kcal/mol
acc=acc/hart2kcl! convert to hartree
read(myid+1000,*) E_range! energy range in kcal/mol above asymptote
E_range=E_range/hart2kcl! convert to hartree
read(myid+1000,*) focus! 0=no, 1=yes; if "focus=0", all energy range is considered
read(myid+1000,*) increment!(kcal/mol) considered energies= asymptotic energy - "increment"; if "focus=1"
increment=increment/hart2kcl! convert to hartree
read(myid+1000,*) rmin(1)! minimum value of R (distance between centers of mass)
read(myid+1000,*) rmax(1)! maximum value of R
if(rmin(1)>rmax(1))then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "rmax(1)" should be larger than "rmin(1)".'
  call MPI_FINALIZE(rcc)
  stop
endif
read(myid+1000,*) focus_onR! 0=no, 1=yes; if "focus_onR=0", all R-range is considered
read(myid+1000,*) minR! minimum R considered; if "focus_onR=1"
read(myid+1000,*) maxR! maximum R considered; if "focus_onR=1"
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
read(myid+1000,*) focus_onLR! 0=no, 1=yes; if "focus_onLR=1", only long-range is considered
read(myid+1000,*) xminLR! minimum R considered; if "focus_onLR=1" (after low- and seed-grid are computed)
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
read(myid+1000,*) numpoints! # high-level seed points
if (numpoints==0) numpoints=(2000*natom1*natom2)/(9*flip1*flip2)
read(myid+1000,*) ab_flag! 1=single point energies, 2=gradients
read(myid+1000,*) code_flag! 1=Gaussian, 2=Molpro, 3=Aces II
read(myid+1000,*) x1! # of points to add per cycle (usually numprocs)
if(x1<=0.d0)then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "numadded" should be a positive number'
  call MPI_FINALIZE(rcc)
  stop
endif
x1=x1*dble(numprocs)
numadded=int(x1)
read(myid+1000,*) n_test! test accuracy via random high-level ab initio calc. every "n_test" cycles
ntest=n_test
read(myid+1000,*) x1! # of ab initio test points, tested every "n_test" iterations
x1=x1*dble(numprocs)
num_ab_err_points=int(x1)
read(myid+1000,*) maxpoints! maximum number of high level ab intio data points
read(myid+1000,*) order_1! maximum power of R=exp(alpha*r)
read(myid+1000,*) order_2! L1 max in angular basis
read(myid+1000,*) order_3! L2 max in angular basis
read(myid+1000,*) order_4! L=L1+L2 max in angular basis
read(myid+1000,*) dist_tol! minimum internuclear distance (between atoms in different fragments)
read(myid+1000,*) num_err_points! total # of points to randomly estimate errors at the beginning of each loop
if(num_err_points<numprocs)then
  if(myid==0)write(110,*)'Error found in the input file!' 
  if(myid==0)write(110,*)'The value of "num_err_points" should be greater than the number of processors:',numprocs
  call MPI_FINALIZE(rcc)
  stop
endif
num_err_points=num_err_points/numprocs! # of points to randomly estimate errors per processor
read(myid+1000,*) low_grid! 0=no, 1=yes, 2=compute only the low-level-grid
!read(myid+1000,*) extra_D! include extra (previously computed) ab initio points?  0=no, 1=yes
!0               ! [extra_D] include extra (high-level) ab initio points?  0=no, 1=yes
!extra-data.dat  ! [dat_label] name of the input file containing extra data, only needed if: "extra_D=1"
if(low_grid>0)then
 read(myid+1000,*)
 read(myid+1000,*) numpoints_low! number of random low-level seed points
 if (numpoints_low==0) numpoints_low=(10000*natom1*natom2)/(36*flip1*flip2)! = 1.25*numpoints
 read(myid+1000,*) code_flag2! 1=Gaussian, 2=Molpro, 3=Aces II, for low level
 read(myid+1000,*) ab_flag2! 1=single point energies, 2=gradients, for low level
 read(myid+1000,*) low_high_fac! multiple of E_range from asymptote on low grid, to guide high config space (e.g. 1.2d0)
 read(myid+1000,*)seed_cut! seed_cut: multiple of E_range from asymptote on low grid kept in
                           ! low grid (optimize to minimize, but exclude holes in high grid) 
                           ! (e.g. test 2.0,3.0,4.0...until no holes found).
! read(myid+1000,*) subzero! 1-yes, 0=no; subtract zeroth-order fit (low_grid fit) from high fit, not currently recommended
endif
subzero=0
close(myid+1000)

!write(*,*)'lei todo bien',num_err_points,auto_mod
!stop
!-----------------------------------------------------------------------------------
!     Miscellaneous definitions  -------------------------------------------------
!-------------------------------------------------------------------------------
pii=dacos(-1d0)
term=0! Termination flag. The random generation of new points will continue until "term=1"
xpass=0! is switched to 1 after Low-level-PES fit is done, and 2 after seed-grid-PES is computed
alpha=-1d0! coefficient in R=exp(alpha*r) coordinate
wellfocus=0! is switched automatically to 1 when error is below 3-times the accuracy target  
symparts=((exch+1)*(flip1+1)*(flip2+1))! total number of symmetry partners
num_err_points1=25000/numprocs! include only 50000 in the minimization process
!num_err_points1=2000! include only the first 2000 points per processor in the minimization process
if (num_err_points<2000) num_err_points1=num_err_points
only_seed_grid=0! compute only the high-level seed-grid? 0=no, 1=yes
W_a=0.432d0*3/max(natom1,natom2)! scaling for R in dist metric, suggested: 1/fraglength
epss=1d-14
count3=numpoints
! calculated ab initio points
numpoints_low_count=0! will be modified (if necessary) during a restart
numpoints_high_count=0! will be modified (if necessary) during a restart
! adjust (reduce) angular range depending on the system's symmetry
rmax(2)=0.999d0
if(flip1>0)then
  rmin(2)=0.001d0
else
  rmin(2)=-0.999d0
endif
rmax(3)=0.999d0
if(flip2>0)then
  rmin(3)=0.001d0
else
  rmin(3)=-0.999d0
endif
rmax(4)=3.14059265d0
rmin(4)=0.001d0! for two linear fragments, the system is always symmetric on phi   
! prepare R-intervals for the final summary of errors 
allocate(xR(xnum_intervals_R+1))
xdR=(rmax(1)-rmin(1))/dble(xnum_intervals_R)
do j=1,xnum_intervals_R+1
 xR(j)=rmin(1)+(j-1)*xdR
enddo
! conversion factors
if(xunitE==1) CONVE=1.d0
if(xunitE==2) CONVE=hart2kcl
if(xunitE==3) CONVE=hart2wn
if(xunitE==4) CONVE=hart2meV
if(xunitD==1) CONVD=1.d0
if(xunitD==2) CONVD=bohr2ang
! convert units
ugrad=CONVE/CONVD
acc=acc*CONVE
E_range=E_range*CONVE
increment=increment*CONVE
! basis for "minimal fit" to high level data, also used to fit low-level-grid. 
order_1_min=3   
order_2_min=3
order_3_min=3
order_4_min=3

! calculate the size of high-degree basis:
call basis_size(4,order_1,order_2,order_3,order_4,basis_1)  
! check for minimum basis support (including sym. partners)
if(numpoints*2*symparts<basis_1)then
  if(myid==0)write(110,*)'Basis too big for number of points. Increase value of "numpoints" on input file' 
  if(myid==0)write(110,*)'numpoints =',numpoints,'     basis_1 =',basis_1
  call MPI_FINALIZE(rcc)
  stop
endif
! calculate the size of lower-degree basis:
call basis_size(4,order_1-1,order_2-1,order_3-1,order_4-1,basis_2) 
! calculate the size of minimal basis
call basis_size(4,order_1_min,order_2_min,order_3_min,order_4_min,basis_3)  
! set interpolative weight function depending on whether or not gradient data
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
!** ALLOCATE **! 
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
if(ab_flag==2)then
  allocate(seed_grad(maxpoints,3*(natom)),grad(2*symparts*maxpoints,4))
endif
if(low_grid>0)then
  allocate(low_pot(maxpoints),seed_low(maxpoints,4),seed_pot_low(maxpoints))
  if(ab_flag2==2)then
   allocate(seed_grad_low(maxpoints,3*(natom)), xgrad(3*natom))
  endif
endif
dcart=0d0
! format for the output files
601 FORMAT(I10,4f20.15,f20.8)              ! AbINITIO.dat & AbINITIO_low.dat
write(f602,'( "(I10,4f20.14,f20.8,",I3,"f20.8)" )')natom*3 !! with gradients                        !increase precision
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

if(myid==0)then
  write(100,'(A30,I6)')' High-degree basis size:      ',basis_1
  write(100,'(A30,I6)')' Lower-degree basis size:     ',basis_2
  write(100,'(A30,I6)')' Size of minimal basis:       ',basis_3
  write(100,*)'Initialization completed.'
  write(100,*)
endif

if (low_grid>0) then
!-----------------------------------------------------------------------------------
!     COMPUTE LOW-LEVEL AB INITIO GRID
!-----------------------------------------------------------------------------------
 if(myid==0)write(110,*)
 if(myid==0)write(110,*)'*** Low-level-PES'
 if(myid==0)write(110,*)

 14 CONTINUE
 IF (restart==1) THEN
!  check if Low-level ab initio grid was computed completely; if not, restart calculations
   inquire(file='AbINITIO_low.dat',exist=logica1)
   if(.not.logica1)then! if the file "AbINITIO_low.dat" does not exist
     if(myid==0)then
       open(unit=334,file='AbINITIO_low.dat')
       do j=1,2000! check if "LGRID-..." files exist
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
     if(numpoints_low_count>=numpoints_low)then! low level grid is already computed
       if(myid==0)close(334)
       if(myid==0)write(100,*)'All low-level ab initio energies have been computed.'
       if(myid==0)write(100,*)'RESTART POINT successfully created!'
       if(myid==0)write(100,*)'Ab initio data saved in output file: "AbINITIO_low.dat".'
       xpass=1
       if(numpoints_low_count>numpoints_low)then
         if(myid==0)write(110,*)'Number of points increased from ',numpoints_low,' to ',numpoints_low_count
         numpoints_low=numpoints_low_count
       endif
       goto 18! skip low-level-grid calculation
     endif
     if(myid==0)write(100,*)'Low-level calculations started, but not completed.'
     if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_low_count
     if(myid==0)write(100,*)'Restarting calculations to complete the Low-level ab initio grid.'
     if(myid==0)write(100,*)
   else! if the file "AbINITIO_low.dat" exist
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
         do j=1,2000! check if "LGRID-..." files exist
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
      goto 18! skip low-level-grid calculation
     endif
     if(myid==0)write(100,*)'Low-level calculations started, but not completed.'
     if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_low_count
     if(myid==0)write(100,*)'Restarting calculations to complete the Low-level ab initio grid.'
     if(myid==0)write(100,*)
   endif
 ENDIF! IF (restart=1)

! select the geometries for the low-level ab initio grid  --------------------------
 if(myid==0)then
   if(xseed(1)==0)then! advance the Sobol sequence in case of a restart
     do i=1,numpoints_low_count*150! if (restart = 0) -> (numpoints_low_count = 0)
      call sobseq(ran_vec)
     enddo
   endif
   write(100,*)'*** Low-level grid computation'
   write(100,'(A18,I6)')' Number of points:',numpoints_low
   do i=numpoints_low_count+1,numpoints_low
874  if(xseed(1)==0)then! Sobol sequence
       call sobseq(ran_vec)
     else! random sequence
       call random_number(ran_vec)
     endif
     ran_vec(1)=dexp(dlog(rmin(1))+ran_vec(1)*(log(rmax(1))-log(rmin(1))))! bias on R
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

! compute the low-level ab initio energies  ----------------------------------------
 ncont=numpoints_low_count
 ncont1=0! count the number of points calculated on each processor
 write(charmyid,'(I5)')myid
 if(myid/=0)open(unit=1000+myid,file='LGRID-'//sys_label(1:nline)//'-Proc'// &
                                       trim(adjustl(charmyid))//'.dat', action='write')
 61 continue! LOAD BALANCING IMPLEMENTED
 IF(myid/=0)THEN! I'm a slave process...
   ! write(*,*) 'Process ', myid, ', asking the master for a geometry'
   call MPI_SEND(myid,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,ierr)
   ! receive answer from the master
   call MPI_RECV(ncont,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,status,ierr)
   ! write(*,*) 'Process ', myid, ', the master answered: ',ncont
   if(ncont==0)then! all geometries have been computed
     close(1000+myid)
     write(110,*) 'Process:',myid,', calculated energies:',ncont1
     goto 63
   else! compute ab initio energy
     ncont1=ncont1+1
     jac(:)=seed_low(ncont,:)
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     call ab_initio(cart,symb,mass,poten,dcart,code_flag2,ab_flag2,natom,myid,0,0,10)
     poten=poten*CONVE
     low_pot(ncont)=poten
     if(ab_flag2==2)then
       dcart=dcart*ugrad  
       seed_grad_low(ncont,:)=dcart
       ! save geometries, energies and gradients (if computed) 
       write(1000+myid,*)ncont,jac(:),low_pot(ncont),seed_grad_low(ncont,:)
     else
       write(1000+myid,*)ncont,jac(:),low_pot(ncont)
     endif
     goto 61
   endif
 ELSE! I'm the master...
   ! receive request from slave
   call MPI_RECV(nid,1,MPI_INT,MPI_ANY_SOURCE,send_data_tag,MPI_COMM_WORLD,status,ierr)
   ! write(*,*) 'MASTER, new request received from process:', nid
   if(ncont<numpoints_low)then! send a new geometry to the slave
     ncont=ncont+1
     call MPI_SEND(ncont,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
     ! write(*,*) 'MASTER, new geometry (',ncont,'), sent to process:', nid
   else
     ncont1=ncont1+1
     ! write(*,*) 'MASTER, sending final message to process:', nid
     call MPI_SEND(0,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
     if (ncont1==numprocs-1) goto 63! all the slaves finished
   endif
   goto 61
 ENDIF
 63 continue
 ! write(*,*) '*** Process ', myid, ', all done!'
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

 do j=1,numprocs-1
   write(charmyid,'(I5)')j 
   open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                        '.dat', action='read')
   do i=1,numpoints_low
     read(230,*,end=2301)ncont
     call MPI_BCAST(low_pot(ncont),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     if(ab_flag2==2)then
      call MPI_BCAST(seed_grad_low(ncont,:),3*(natom),MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
     endif
   enddo
   2301 close(230)
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 IF(myid==2)THEN
   call MPI_SEND(low_pot,numpoints_low,MPI_DOUBLE_PRECISION,0,send_data_tag,MPI_COMM_WORLD,ierr)
   if(ab_flag2==2)then
     do i=1,numpoints_low
      call MPI_SEND(seed_grad_low(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,send_data_tag,MPI_COMM_WORLD,ierr)
     enddo
   endif
 ENDIF
 IF(myid==0)THEN
   call MPI_RECV(low_pot,numpoints_low,MPI_DOUBLE_PRECISION,2,send_data_tag,MPI_COMM_WORLD,status,ierr)
   if(ab_flag2==2)then
     do i=1,numpoints_low
      call MPI_RECV(seed_grad_low(i,:),3*(natom),MPI_DOUBLE_PRECISION,2,send_data_tag,MPI_COMM_WORLD,status,ierr)
     enddo
   endif
 ENDIF
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!***  RESTART POINT # 1  ***********************************************************
 if(myid==0)then
   if(restart==0)open(unit=334,file='AbINITIO_low.dat')
   do j=2,numprocs
     write(charmyid,'(I5)')j-1
     open(unit=230,file='LGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                          '.dat',action='read',status='old')
     do i=1,numpoints_low
       if(ab_flag2==2)then
        read(230,*,end=13)ncont,jac(:),xpot,xgrad(:)
        write(334,f602)ncont,jac(:),xpot,xgrad(:)
       else
        read(230,*,end=13)ncont,jac(:),xpot
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
 18 continue! skip low-level-grid calculation in a restart...
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

 ! check for repeated geometries on "AbINITIO_low.dat"
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
 ! check range of data, set asymptote, max. and min. energy in low-level-grid
 count_seed=0
 Max_E3=2d2
 Min_E=-1d9
 Max_R=0d0
 do i=1,numpoints_low
   if(seed_low(i,1)>Max_R)then
    Max_R=seed_low(i,1)! find the point with largest R
    E_asym=low_pot(i)! energy for point with largest R is set as asymptote energy
   endif
   if(low_pot(i)<Max_E3)then
    Max_E3=low_pot(i)! Lowest energy in the low-level ab initio seed grid
   endif     
   if(low_pot(i)>Min_E) then
    Min_E=low_pot(i)! Highest energy in the low-level ab initio seed grid
   endif
 enddo
 bot_seed=Max_E3! (Lowest energy...) Allows looking for holes
 Max_E_seed=E_asym+E_range*low_high_fac! sets limit for low-level guide to high-level config. space
 do i=1,numpoints_low
   if(low_pot(i)<E_asym+E_range*seed_cut)then
    count_seed=count_seed+1! qualifying points from LGRID for the guide fit
   endif
 enddo
 allocate(coords_seed(2*symparts*count_seed,4),ind_seed(count_seed*2*symparts))
 allocate(ind2_seed(count_seed*2*symparts),d_seed(count_seed*2*symparts))
 allocate(pot_seed(count_seed*2*symparts),b2_seed(basis_3,count_seed*2*symparts))
 if(ab_flag2==2)then
  allocate(grad_seed(count_seed*2*symparts,4))
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

 if(low_grid==2)then! TERMINATE: no high-level-PES needed.
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

 ! Prepare to fit  -----------------------------------------------------------------
 count=0
 do i=1,numpoints_low
   ! select qualifying points from LGRID for the guide fit
   if(low_pot(i)<E_asym+E_range*seed_cut)then
     if(ab_flag2==2)then
      dcart=seed_grad_low(i,:)
     endif
     ! Include symmetry permutations 
     call symmetry2(seed_low(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,ab_flag2,exch,flip1,flip2)
     do k2=0,1! reflection to the other side of torsion
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
 ! Compute local density ("d") for each point in low-level-fit grid  ---------------
 do i=1,count_seed
   count=0 
   do ip=1,count_seed         
    count=count+1
    call dist_metric(coords_seed(i,:),coords_seed(ip,:),W_a,ind_seed(count))
   enddo
   call indexxy(count_seed,ind_seed,ind2_seed)
   d_seed(i)=ind_seed(ind2_seed(zz4))
 enddo
 ! set number of highest weighted neighboring points included in each 
 ! local fit for low-level data (less if gradient data)
 support=min(count_seed,8*basis_3/(4*(ab_flag2-1)+1))
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!***  DO THE FITS  *****************************************************************
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

else! if (low_grid<=0)
 if(myid==0)write(100,*)'No low-level ab initio calculations requested.'
endif ! low-level-PES finished

!-----------------------------------------------------------------------------------
!     COMPUTE HIGH-LEVEL AB INITIO SEED GRID
!-----------------------------------------------------------------------------------

IF(restart==1)THEN
  ! check if high-level ab initio seed-grid was computed completely; if not, restart calculations
  inquire(file='AbINITIO.dat',exist=logica1)
  if(.not.logica1)then! if the file "AbINITIO.dat" does not exist
    if(myid==0)then
      open(unit=222,file='AbINITIO.dat')
      do j=1,1000! check if "HGRID-..." files exist
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
       goto 24! skip high-level-grid calculation
      endif
      if(myid==0)write(100,*)'High-level seed-grid calculation started, but not completed.'
      if(myid==0)write(100,*)'Calculated geometries so far: ',numpoints_high_count
      if(myid==0)write(100,*)'Restarting calculations to complete high-level ab initio seed-grid.'
      if(myid==0)write(100,*)
    endif
  else! if the file "AbINITIO.dat" exist
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
      do j=1,1000! check if "HGRID-..." files exist
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
      goto 24! skip high-level-grid calculation
    endif
  endif
ENDIF! IF (restart==1)

! select random geometries for the high-level seed grid  ---------------------------
if(myid==0)then    
  if(xseed(1)==0)then! advance the Sobol sequence in case of a restart
    ! if (restart = 0) -> (numpoints_low_count + numpoints_high_count = 0)
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
873 if(xseed(1)==0)then! Sobol sequence
      call sobseq(ran_vec)
    else! random sequence
      call random_number(ran_vec)
    endif
    ran_vec(1)=dexp(dlog(rmin(1))+ran_vec(1)*(dlog(rmax(1))-dlog(rmin(1))))! bias on R
    do j=2,4  
      range = rmax(j) - rmin(j)
      ran_vec(j)=rmin(j)+range*ran_vec(j)
    enddo
    seed(i,:)=ran_vec(:)
    call INT_Cart(cart,seed(i,:),mass,natom1,natom2,ref1,ref2)
    call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
    if (dist_flag==1) goto 873
    if(low_grid>0)then! filter high energies using the low-grid-PES as reference
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

! compute the high-level ab initio energies  ----------------------------------------
ncont=numpoints_high_count
ncont1=0! count the number of points calculated on each processor
write(charmyid,'(I5)')myid
if(myid/=0)open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'// &
                                     trim(adjustl(charmyid))//'.dat', action='write')
64 continue! LOAD BALANCING IMPLEMENTED
IF(myid/=0)THEN! I'm a slave process...
  ! write(*,*) 'Process ', myid, ', asking the master for a geometry'
  call MPI_SEND(myid,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,ierr)
  ! receive answer from the master
  call MPI_RECV(ncont,1,MPI_INT,0,send_data_tag,MPI_COMM_WORLD,status,ierr)
  ! write(*,*) 'Process ', myid, ', the master answered: ',ncont
  if(ncont==0)then! all geometries have been computed
    close(1000+myid)
    write(110,*) 'Process:',myid,', calculated energies:',ncont1
    goto 65
  else! compute ab initio energy
    ncont1=ncont1+1
    jac(:)=seed(ncont,:)
    pass_check=1 
    call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
    801 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,0,0,pass_check)
    ! try molpro1.abi, and then molpro2.abi if no convergence 
    ! ("pass_check" limits the number of attempts)
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
      ! save geometries, energies and gradients (if computed) 
      write(1000+myid,*)ncont,jac(:),seed_pot(ncont),seed_grad(ncont,:)
    else
      write(1000+myid,*)ncont,jac(:),seed_pot(ncont)
    endif
    goto 64
  endif
ELSE! I'm the master...
  ! receive request from slave
  call MPI_RECV(nid,1,MPI_INT,MPI_ANY_SOURCE,send_data_tag,MPI_COMM_WORLD,status,ierr)
  ! write(*,*) 'MASTER, new request from process:', nid
  if(ncont<count3)then! send a new geometry to the slave
    ncont=ncont+1
    call MPI_SEND(ncont,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
    ! write(*,*) 'MASTER, new geometry (',ncont,'), sent to process:', nid
  else
    ncont1=ncont1+1
    ! write(*,*) 'MASTER, sending final message to process:', nid
    call MPI_SEND(0,1,MPI_INT,nid,send_data_tag,MPI_COMM_WORLD,ierr)
    if (ncont1==numprocs-1) goto 65! all the slaves finished
  endif
  goto 64
ENDIF
65 continue
! write(*,*) '*** Process ', myid, ', ya terminé...'
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

do j=1,numprocs-1
  write(charmyid,'(I5)')j
  open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                       '.dat', action='read')
  do i=1,count3
    read(230,*,end=2300)ncont
    call MPI_BCAST(seed_pot(ncont),1,MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
    if(ab_flag==2)then
     call MPI_BCAST(seed_grad(ncont,:),3*(natom),MPI_DOUBLE_PRECISION,j,MPI_COMM_WORLD,ierr)
    endif
  enddo
  2300 close(230)
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF(myid==2)THEN
  call MPI_SEND(seed_pot,count3,MPI_DOUBLE_PRECISION,0,send_data_tag,MPI_COMM_WORLD,ierr)
  if(ab_flag==2)then
    do i=1,count3
     call MPI_SEND(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,0,send_data_tag,MPI_COMM_WORLD,ierr)
    enddo
  endif
ENDIF
IF(myid==0)THEN
  call MPI_RECV(seed_pot,count3,MPI_DOUBLE_PRECISION,2,send_data_tag,MPI_COMM_WORLD,status,ierr)
  if(ab_flag==2)then
    do i=1,count3
     call MPI_RECV(seed_grad(i,:),3*(natom),MPI_DOUBLE_PRECISION,2,send_data_tag,MPI_COMM_WORLD,status,ierr)
    enddo
  endif
ENDIF
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!***  RESTART POINT # 2  ***********************************************************
if(myid==0)then
  if(restart==0)open(unit=222,file='AbINITIO.dat')
  do j=2,numprocs
    write(charmyid,'(I5)')j-1
    open(unit=230,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                            '.dat', action='read',status='old')
    do i=1,count3
      if(ab_flag==2)then
       read(230,*,end=23)ncont,jac(:),xpot,xgrad(:)
       write(222,f602)ncont,jac(:),xpot,xgrad(:)
      else
       read(230,*,end=23)ncont,jac(:),xpot
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
24 continue! skip high-level-grid calculation in a restart...
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

count3=max(numpoints_high_count,numpoints)
! check that no repeated geometries exist on "AbINITIO.dat"
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
  tot_abinitio=count3! save the number of ab initio points calculated so far...
endif
! check range of data, set asymptote, max. and min. energy on high-level seed grid
Max_E=2d2
Min_E=-1d9
Max_R=0d0
do i=1,count3
  if(seed(i,1)>Max_R) then
   Max_R=seed(i,1)! find the point with largest R
   E_asym=seed_pot(i)! energy for point with largest R is set as asymptote energy
  endif
  if(seed_pot(i)<Max_E) then! Lowest energy in the high-level ab initio seed grid
   Max_E=seed_pot(i)! Allows looking for holes
  endif
  if(seed_pot(i)>Min_E) then! Highest energy in the high-level ab initio seed grid 
   Min_E=seed_pot(i)! Allows looking for holes
  endif
enddo
Glob_min=Max_E! Global minimum
Max_E=E_asym+E_range! maximum energy above asymptote
! sets a warning, that will be triggered if 
! data added later is more than 50% above range of interest
E_limit=E_asym+1.7d0*E_range
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
! Prepare to fit  ------------------------------------------------------------------
count=0
do i=1,count3 
  ! Enters qualifying points from HGRID into high-level guide fit
  ! Include all energies within 120% of range of interest,
  ! This also excludes problematic points whose energies are set to 1d2 by the ab initio script.
  if(seed_pot(i)<Max_E+0.2d0*E_range)then
    if(ab_flag==2)then
     dcart=seed_grad(i,:)
    endif
    ! Include symmetry permutations
    call symmetry2(seed(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,  &
                    ab_flag,exch,flip1,flip2)
    do k2=0,1! reflected to the other side of torsion
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

! compute local density ("d") for each point in high-level-fit grid  ---------------
do i=1,count3  
  count=0 
  do ip=1,count3         
    count=count+1
    call dist_metric(coords(i,:),coords(ip,:),W_a,ind(count))
  enddo
  call indexxy(count3,ind,ind2)
  d(i)=ind(ind2(zz4))
enddo
! implement local expansions for high-level data
support=min(count3,4*basis_1/(4*(ab_flag-1)+1)) 
!if(myid==0) then! print local density if necessary 
!  do i=1,count3
!   write(680,*) d(i)
!  enddo
!endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     DO THE FITS  *****************************************************************
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
!!!!subroutine, stored expansions (b2, and b2_lower) given in module "dynamic parameters"           ! ???
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

goto 26! comment this line to use...
! Debugging test functions   -------------------------------------------------------
if(myid==0)then
 open(unit=333,file='DEBUGGING.dat')
 jac3(:)=coords(12,:)+0.01d0
 write(333,*) jac3
 ! Example of calling various functions to test, analytic vs. numerical grads, etc..
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
endif! debugging test funcs.
26 continue
call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!***********************************************************************************
!===================================================================================
!     START LOOPING OVER, 
!     ADDING AUTOMATICALLY DETERMINED POINTS UNTIL THE DESIRED ACCURACY IS REACHED
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

add=1+(maxpoints-(count3/symparts/2))/numadded! maximum # of cycles
if(add<1)then
  if(myid==0)write(110,*)'Error: Number of loops lower than 1...'
  if(myid==0)write(110,*)'Increase the value of "maxpoints" (in the input file) & restart...'
  call MPI_FINALIZE(rcc)
  stop
endif

if(myid==0)then
 write(100,*)
 write(100,*)'***  START AUTOMATIC DATA POINT GENERATION  ***'
 write(100,'(A29,I6)')'      max. number of cycles: ',add
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

! reset the random number seed differently for each process 
if(restart==0)then
  xseed=1000+xseed+(myid+1)! seed = original seed + (myid+1): allows exact reproduction if needed
  CALL RANDOM_SEED(PUT=xseed) 
else
  xseed=date_time(7)+date_time(8)+(myid+1)! seed = sec + mili-sec + (myid+1)
  CALL RANDOM_SEED(PUT=xseed) 
endif

! a Sobol sequence will be always used to improve the long range in case "focus_onLR=1"
if(restart==1)then! advance the Sobol sequence if neccesary
  ncont=0
  if(numpoints_high_count>=numpoints)ncont=ncont+numpoints_high_count
  if(numpoints_low_count>=numpoints_low)ncont=ncont+numpoints_low_count
  do i=1,ncont*150
   call sobseq(ran_vec)
  enddo
endif

!***********************************************************************************
DO loop=1,add ! Loop over, adding automatically determined points  
!***********************************************************************************
 numpoints_high_count=count3/(2*symparts)
 if(myid==0)write(110,*)
 if(myid==0)write(110,*)'Loop: ',loop-1

 !----------------------------------------------------------------------------------
 ! Random test to estimate the error, at the beginning of each loop  -------------
 !------------------------------------------------------------------------------
 ! Also finds holes ...
 ncont=0
 rms_error=0d0
 mean_error=0d0
 error_points=0d0 
 hole_val=Glob_min
 hole(:,1)=rmax(1)
 hole(:,2)=rmax(2)
 hole(:,3)=rmax(3)
 hole(:,4)=rmax(4)
 do i=1,num_err_points! select random geometries to test (per processor)  ------
   19 call random_number(ran_vec)
   ncont=ncont+1
   if(ncont==1000000)then! avoid and endless loop
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
     range=rmax(1)-xminLR! only test the long-range, as indicated in the input file
     ran_vec(1)=xminLR+range*ran_vec(1)
   else
     range=maxR-minR! only test the R-range indicated in the input file
     ran_vec(1)=minR+range*ran_vec(1)
   endif
   do j=2,4
     range=rmax(j)-rmin(j)
     ran_vec(j)=rmin(j)+range*ran_vec(j)
   enddo
   call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
   call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
   if(dist_flag==1) goto 19! exclude points if the atoms are too close
   xi=ran_vec(:)
   ! exclude geometries based on the estimated energy
   if(low_grid>0)then! exclude points based on energy
     temp3=func_actual_seed(xi)! (seed-grid-PES estimate + low-grid cutoff)
     if(temp3>Max_E_seed) goto 19
   endif
   if(focus_onLR==1)goto 191! avoid energy-focus if only long-range is considered
   if(focus>0)then! focus only on the energy range specified in the input file
     if (func_actual(xi)>E_asym+(0.05d0/hart2kcl*CONVE)-increment) goto 19
   else 
     if (wellfocus>0)then! exclude positive energies if error is already below 3-times acc. target
       if (func_actual(xi)>E_asym+(0.05d0/hart2kcl*CONVE)) goto 19
     endif
   endif
   191 continue
   if(subzero==0)then
     temp=func_actual_min(xi)
     if(temp>Max_E) goto 19! exclude points with E (min-PES estimate) > E_asym + E_range
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
 ! Holes (energy values below the minimum energy value computed so far)
 holepatch=0
 call indexxy(numprocs,hole_val,indhole)
 if(hole_val(indhole(1))<Glob_min-0.1d0/hart2kcl*CONVE)then! Holes found...
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
   ! writes error estimates
   if(holepatch==0)write(100,'(I4,2ES15.6,I10)')loop-1,somme,somme2,(count3/(2*symparts))     
   if(holepatch==1)write(100,'(I4,2ES15.6,I10,A21)')loop-1,somme,somme2,  &
                               (count3/(2*symparts)),'       Holes found...'
   if(somme<acc) then
     term=1! TERMINATE: accuracy target have been reached, the PES is complete!
   endif
   if(somme<3d0*acc.and.wellfocus<1.and.focus==0) then
     wellfocus=1
     if(focus_onLR==0)then
       if (term==0) write(100,*)
       if (term==0) write(100,*) '--- switching focus to negative energies'
     endif
   endif

!  old "frag_data"
!  open(unit=10,file='frag_data-'//sys_label(1:nline)//'.dat')
!  do i=1,natom
!    write(10,*) symb(i)
!  enddo
!  do i=1,natom
!    write(10,*) mass(i)
!  enddo  
!  do i=1,3*natom1
!    write(10,*) ref1(i)
!  enddo
!  do i=1,3*natom2
!    write(10,*) ref2(i)
!  enddo
!  close(10)


   ! Writes out current state of the PES
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
   if(low_grid>0)then! writes out low-level-PES
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

 !----------------------------------------------------------------------------------
 IF((term>0).or.(only_seed_grid==1))THEN! TERMINATE: the PES is complete --------
 !-------------------------------------------------------------------------------
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
       write(100,*)'                Only high-level seed-grid-PES computed.                  '
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
   ! Random global test to estimate the final error ------------------------------------
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
   ! select random geometries
   do i=1,200000/numprocs
     27 call random_number(ran_vec)
     ! sample ALL the volume of the configuration space
     ran_vec(1)=rmin(1)+ran_vec(1)*(rmax(1)-rmin(1))
     ran_vec(2)=-0.999d0+ran_vec(2)*1.998d0
     ran_vec(3)=-0.999d0+ran_vec(3)*1.998d0
     ran_vec(4)=(-0.999d0+ran_vec(4)*1.998d0)*pii
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if(dist_flag==1) goto 27! exclude points with atoms too close
     xi=ran_vec(:)
     if(low_grid>0)then! exclude points based on energy
       temp3=func_actual_seed(xi)! (seed-grid-PES estimate + low-grid cutoff)
       if(temp3>Max_E_seed) goto 27
     endif
     ! exclude points with E > E_asym + E_range
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
     !*** statistics for each R-interval ***
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
   ! global statistics
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
     open(unit=200,file='errores_'//trim(adjustl(xchar1))//'.dat')
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
       x1=xrms_R1(1,j)*100.d0/xMAPot(1,j)
       x2=xrms_R2(1,j)*100.d0/xMAPot(1,j)
       x3=xrms_R3(1,j)*100.d0/xMAPot(1,j)
       x4=xrms_R4(1,j)*100.d0/xMAPot(1,j)
       x5=(xR(j)+xR(j+1))/2.d0
       IF(focus>0.and.focus_onLR==0)THEN
         if((xR(j)<=minR).and.(minR<xR(j+1)))write(100,'(A125)')'-----------------------------------&
         ------------------------------------------------------------------------------------------'
         write(200,'(F6.2,ES13.5,4(ES11.3,F6.2,I8))')x5,xMAPot(1,j),xrms_R1(1,j), &
          x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j)),&
          xrms_R4(1,j),x4,int(xR_count4(1,j))
         write(100,'(2F6.2,ES13.5,4(ES11.3,F6.2,I8))')xR(j),xR(j+1),xMAPot(1,j),&
          xrms_R1(1,j),x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j)),&
          xrms_R4(1,j),x4,int(xR_count4(1,j))
         if((xR(j)<maxR).and.(maxR<=xR(j+1)))write(100,'(A125)')'-----------------------------------&
         ------------------------------------------------------------------------------------------'
       ELSE
         if((xR(j)<=minR).and.(minR<xR(j+1)))write(100,'(A101)')'-------------------------&
         ---------------------------------------------------------------------------'
         write(200,'(F6.2,ES13.5,3(ES11.3,F6.2,I8))')x5,xMAPot(1,j),xrms_R1(1,j), &
          x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j))
         write(100,'(2F6.2,ES13.5,3(ES11.3,F6.2,I8))')xR(j),xR(j+1),xMAPot(1,j),xrms_R1(1,j),&
          x1,int(xR_count1(1,j)),xrms_R2(1,j),x2,int(xR_count2(1,j)),xrms_R3(1,j),x3,int(xR_count3(1,j))
         if((xR(j)<maxR).and.(maxR<=xR(j+1)))write(100,'(A101)')'-------------------------&
         ---------------------------------------------------------------------------'
       ENDIF
     enddo
     close(200)
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
     ! rename the final PES-data-file, including the number of ab initio points
     write(xchar1,'(I8)')count3/(2*symparts)
     call system('mv PES-'//sys_label(1:nline)//' PES-'//sys_label(1:nline)//'-'//trim(adjustl(xchar1)))
     ! rename (if exist) the 'output' file, including the number of ab initio points
     inquire(file='output.dat',exist=logica1)
     if(logica1)call system('mv output.dat output-'//trim(adjustl(xchar1))//'.dat')
   endif
   call MPI_FINALIZE(rcc)
   stop
 ENDIF
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

 !----------------------------------------------------------------------------------
 IF(focus_onLR==1)THEN! focus only in the long range   --------------------------
 !-------------------------------------------------------------------------------
   if(myid==0)then! select random (Sobol) geometries for the test   -----------
    ncont=0
    do i=1,numadded
     28 call sobseq(ran_vec)
     range=rmax(1)-xminLR! only test the long-range, as indicated in the input file
     ran_vec(1)=xminLR+range*ran_vec(1)! now unbiased sampling on R...
     do j=2,4
       range = rmax(j) - rmin(j)
       ran_vec(j)=rmin(j)+range*ran_vec(j)
     enddo
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if (dist_flag==1) ncont=ncont+1
     if (dist_flag==1) goto 28! exclude points if the atoms are too close
     jac=ran_vec(:)
     if(low_grid>0)then! exclude points based on energy
       temp3=func_actual_seed(jac)! (seed-grid-PES estimate + low-grid cutoff)
       if(temp3>Max_E_seed) ncont=ncont+1
       if(temp3>Max_E_seed) goto 28
     endif
     ! exclude points with E > E_asym + E_range
     if(subzero==0)then
       temp=func_actual_min(jac)! (min-PES estimate)
       if (temp>Max_E) ncont=ncont+1
       if (temp>Max_E) goto 28
       temp=func_actual(jac)
       if (temp>Max_E) ncont=ncont+1
       if (temp>Max_E) goto 28
     else
       temp=func_actual_min(jac)! (min-PES estimate)
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
   ! Compute the ab initio energies   -------------------------------------------
   write(charmyid,'(I5)')myid
   open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                             '.dat',action='write')
   do i=myid+1,numadded,numprocs
     jac=seed(i,:)
     pass_check=1
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     804 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,0,0,pass_check)
     ! try molpro1.abi, and then molpro2.abi if no convergence 
     ! ("pass_check" limits the number of attempts)
!     if(poten>199d0)then
!       if(pass_check<2) then                                                                       !!!!!!!!!! revisar  
!        pass_check=pass_check+1
!        goto 803 
!       endif
!     endif
!     if (poten>199d0) goto xxx28 !!!!!!!!! revisar  
     poten=poten*CONVE
     seed_pot(i)=poten
     if(ab_flag==2)then
      dcart=dcart*ugrad 
      seed_grad(i,:)=dcart 
      ! save geometries, energies and gradients (if computed) 
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
   GOTO 189! skip the automatic data-point-location & calculation
 ENDIF!IF(focus_onLR==1)

 if(ntest/=0)then! only if "ntest" non equal zero:
 !==================================================================================
 ! TEST ACCURACY VIA RANDOM HIGH-LEVEL AB INITIO CALCULATIONS   ===================
 !--------------------------------------------------------------------------------
 IF(mod(loop,ntest)==ntest-1)THEN! test every "n_test" cycles   ----------------
 !------------------------------------------------------------------------------
   rms_error=0d0
   mean_error=0d0
   if(myid==0)then! select random geometries for the test   -------------------
    do i=1,num_ab_err_points
     20 call random_number(ran_vec)
     if(focus_onLR==1)then
       range=rmax(1)-xminLR! only test the long-range, as indicated in the input file
       ran_vec(1)=xminLR+range*ran_vec(1)
     else
       range=maxR-minR! only test the R-range indicated in the input file
       ran_vec(1)=minR+range*ran_vec(1)
     endif
     do j=2,4
       range=rmax(j)-rmin(j)
       ran_vec(j)=rmin(j)+range*ran_vec(j)
     enddo
     call INT_Cart(cart,ran_vec,mass,natom1,natom2,ref1,ref2)
     call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
     if (dist_flag==1) goto 20! exclude points if the atoms are too close
     jac=ran_vec(:)
     !*** exclude geometries based on the estimated energy
     if(low_grid>0)then! seed-grid-PES: estimate + low-level-grid: cutoff
       temp3=func_actual_seed(jac)
       if (temp3>Max_E_seed) goto 20
     endif
     if(focus_onLR==1)goto 192! avoid energy-focus if only long-range is considered
     if(focus>0)then! focus only on the energy range specified in the input file
       if (func_actual(jac)>E_asym+(0.05d0/hart2kcl*CONVE)-increment) goto 20
     else
       if(wellfocus>0)then! exclude positive energies if error is already below 3-times acc. target
        if (func_actual(jac)>E_asym+(0.05d0/hart2kcl*CONVE)) goto 20
       endif
     endif
     192 temp=func_actual_min(jac)! (min-PES estimate)
     if(subzero==0)then
       if (temp>Max_E) goto 20! exclude points with E > E_asym + E_range
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
   ! Compute the ab initio energies   -------------------------------------------
   write(charmyid,'(I5)')myid
   open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                             '.dat',action='write')
   do i=myid+1,num_ab_err_points,numprocs
     jac=seed(i,:)
     pass_check=1
     call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
     803 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,0,0,pass_check)
     ! try molpro1.abi, and then molpro2.abi if no convergence 
     ! ("pass_check" limits the number of attempts)
!     if(poten>199d0)then
!       if(pass_check<2) then                                                                       !!!!!!!!!! revisar  
!        pass_check=pass_check+1
!        goto 803 
!       endif
!     endif
!     if (poten>199d0) goto 20 !!!!!!!!! revisar  
     poten=poten*CONVE
     seed_pot(i)=poten
     if(ab_flag==2)then
      dcart=dcart*ugrad 
      seed_grad(i,:)=dcart 
      ! save geometries, energies and gradients (if computed) 
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
     ! writes error estimates
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
   GOTO 189! skip the automatic data-point-location & calculation
 ENDIF! test accuracy every "n_test" cycles
 endif


 !----------------------------------------------------------------------------------
 ! AUTOMATIC DATA POINT LOCATION FOR THE NEXT LOOP  ------------------------------
 !------------------------------------------------------------------------------

 ! COMPUTE MINIMIZATION STARTING POINTS   **************************************
 allocate(ind7(count3),ind8(count3))
 count7=0
 ! find midpoints between computed data points in HGRID   ----------------------
 do i=1,count3,2*symparts
   do ip=1,count3
     call dist_metric(coords(i,:),coords(ip,:),W_a,ind7(ip))
   enddo
   call indexxy(count3,ind7,ind8)
   do j=2,5!zz4/5
     if(ind8(j)>i)then
       xtemp=0.5d0*(coords(i,1)+coords(ind8(j),1))
       ! if only "long-range" is considered (focus_onLR=1)
       if((focus_onLR==1).and.(xtemp>=xminLR))then! select only midpoints in the adequate R-range
         count7=count7+1
         do k=1,4
          test_points(count7,k)=0.5d0*(coords(i,k)+coords(ind8(j),k))
         enddo
       endif
       if(focus_onLR==0)then
         if((xtemp>=minR).and.(xtemp<=maxR))then! select midpoints in the adequate R-range
           do k=1,4
             xi(k)=0.5d0*(coords(i,k)+coords(ind8(j),k))
           enddo
           ! and in the adequate energy range
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
 ! remove duplicates   ---------------------------------------------------------
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
 ! include worst points in the initial random test   ---------------------------
 call indexxy(num_err_points1*numprocs,error_points(:,5),ind6)
 goto 1778! print worst points if needed (comment this line to use)  
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

 ! LOOP OVER MINIMIZATIONS   ***************************************************
 write(3000+myid,*)
 write(3000+myid,*)'Loop: ',loop-1
!write(*,*)'proc:',myid,'1'
 gtol=1d-4
 do i=myid+1,count7,numprocs      
   reiter=5
   iter=5
   ! report initial geometry and squared diff. value   
   xi=test_points(i,:)
   sqrdiff=dabs(func(xi))
   write(3000+myid,'(4f10.4)') xi
   write(3000+myid,*) i,sqrdiff
   ! find geometry that maximize square diff. value
   if(dabs(func(xi))>1d-4*acc/hart2kcl*CONVE/hart2kcl*CONVE)then
     call frprmn(xi,gtol,iter,valeur)
     sqrdiff=dabs(valeur)
   else
     sqrdiff=0d0
   endif
   test_points(i,:)=xi
   order_diff(i)=sqrdiff      
   ! report final geometry and maximized square diff. value
   write(3000+myid,'(4f10.4)') xi      
   write(3000+myid,*) i,sqrdiff
   write(3000+myid,*)      
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!write(*,*)'proc:',myid,'2'
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
!write(*,*)'proc:',myid,'3'

 ! ADD NEW POINTS   ************************************************************
 numadded_act=numadded
 skip=0
 count=0
 do kk=1,numadded     
   count=count+1
   177 seed(count,:)=test_points(ind4(count7+1-skip-kk),:)
   ! include (as last geometry) the lowest 'hole' found (if any)
   if(holepatch==1.and.kk==numadded)then
     seed(count,:)=hole(indhole(1),:)
     holepatch=0
   endif
   ! avoid geometries with atoms too close
   call INT_Cart(cart,seed(count,:),mass,natom1,natom2,ref1,ref2)
   call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
   if(dist_flag==1) then
     skip=skip+1
     goto 177
   endif
   ! avoid geometries outside the symmetry subspace
   do j=1,4
     if(seed(count,j)>rmax(j).or.seed(count,j)<rmin(j))then
       skip=skip+1
       goto 177
     endif
   enddo
   ! avoid geometries outside the adequate energy range
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
   ! avoid two new geometries in the same region (similar configurations)
   deallocate(ind,ind2)
   allocate(ind(count),ind2(count))
   local_d=0d0
   if(count>1)then
     do ip=1,count! compute "d"
       call dist_metric(seed(count,:),seed(ip,:),W_a,ind(ip))
     enddo
     call indexxy(count,ind,ind2)
     local_d=ind(ind2(2)) ! ind(ind2(zz4/2))
     if(local_d<3d-1) then
       skip=skip+1
       if(skip>count7-numadded) then
         write(110,*) 'not enough unique points',kk
         count=count-1
         numadded_act=kk-1
         goto 188
       endif
       goto 177
     endif
   endif
 enddo
 188 continue
!write(*,*)'proc:',myid,'4'
 if (numadded_act<numadded) write(110,*)'Warning! the number of points added is lower than "numadded"'
 if(myid==0) then! print out the new geometries to be computed each cycle
   write(110,*)'Selected geometries (and estimated energy)'
   do kk=1,numadded_act
    write(110,'(I10,4f14.8,3x,f20.8)') kk,seed(kk,:),func_actual(seed(kk,:))
   enddo
 endif

 ! COMPUTE AB INITIO ENERGIES AT NEWLY DETERMINED POINTS   *********************
 write(charmyid,'(I5)')myid
 open(unit=1000+myid,file='HGRID-'//sys_label(1:nline)//'-Proc'//trim(adjustl(charmyid))// &
                           '.dat',action='write')
 do kk=myid+1,numadded_act,numprocs! compute energies    
   jac=seed(count+1-kk,:)
   pass_check=1
   call INT_Cart(cart,jac,mass,natom1,natom2,ref1,ref2)
   802 call ab_initio(cart,symb,mass,poten,dcart,code_flag,ab_flag,natom,myid,0,0,pass_check)
   ! try molpro1.abi, and then molpro2.abi if no convergence 
   ! ("pass_check" limits the number of attempts)
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
     ! save geometries, energies and gradients (if computed) 
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

 !----------------------------------------------------------------------------------
 189 CONTINUE! skip the automatic data-point-location & calculation  -----------
 !-----------------------------------------------------------------------------

 ! save new geometries, energies and gradients (if computed) into "AbINITIO.dat"
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

 ! Prepare to fit  -----------------------------------------------------------------
 count=count3
 do i=1,numadded_act
   if(ab_flag==2)then
      dcart=seed_grad(i,:)
   endif
   ! Include symmetry permutations
   call symmetry2(seed(i,:),dcart,internal,grad_int,mass,ref1,ref2,natom1,natom2,  &
                   ab_flag,exch,flip1,flip2)
   if(seed_pot(i)<E_limit)then! this exclude failed or very high-energy points
     do k2=0,1! reflected to the other side of torsion
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

 ! check that no repeated geometries exist on "AbINITIO.dat" after adding new data
 if(myid==0)then
   ncont=0
   open(unit=222,file='AbINITIO.dat',status='old')
   do i=1,100000
     read(222,*,end=33)j,xx(i),th1(i),th2(i),phi(i),eee(i)
     ncont=ncont+1
   enddo
   33 continue
   tot_abinitio=ncont! save the number of ab initio points calculated so far...
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

 ! check range of data, set asymptote, max. & min. energy on high-level grid
 dynmax=1d2
 dynmin=-1d9
 if(loop==1)Min_E=-1d9! in case Emax was removed from the seed-grid to fit
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
      Max_R=coords(i,1)! find the point with largest R  
      E_asym=pot(i)! energy for point with largest R is set as asymptote energy
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
      Max_R=coords(i,1)! find the point with largest R  
      E_asym=pot(i)+temp! energy for point with largest R is set as asymptote energy
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
 E_limit=E_asym+1.7d0*E_range! limits the max. value of E used in the fit
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

 ! compute local density ("d") for each point in high-level-fit grid
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
 ! implement local expansions for high-level data
 support=min(count3,4*basis_1/(4*(ab_flag-1)+1)) 
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!***  DO THE FITS  *****************************************************************
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

ENDDO! loop over, adding automatically determined points

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
