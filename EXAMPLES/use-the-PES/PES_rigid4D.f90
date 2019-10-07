!*******************************  A U T O S U R F  *********************************
!===================================================================================
!-----------------------------------------------------------------------------------
!-                                                                                 -
!-        AUTOSURF Package: A set of programs for the automated construction       -
!-              of Potential Energy Surfaces on van der Waals systems              -
!-                                                                                 -
!-----------------------------------------------------------------------------------
!===================================================================================
!***********************************************************************************
!-       "PES_rigid4D": SUBROUTINES for ...                                      -
!-        ver. 2.0                                                                 -
!-----------------------------------------------------------------------------------
!-        Input files: "input-AUTOSURF-PES.dat" & "PES-file"                       -
!-                                                                                 -


MODULE dynamic_parameters

implicit none
save
public

real*8,allocatable :: b2(:,:),b2_lower(:,:),b2_minimal(:,:),b2_seed(:,:),d(:)
real*8,allocatable :: d_seed(:),Jac(:),Jac2(:),coords(:,:),coords_seed(:,:)
real*8,allocatable :: cart(:),bdist(:),ref1(:),ref2(:)
real*8,allocatable :: pot(:),pot_seed(:),grad(:,:),grad_seed(:,:),mass(:),dcart(:)
real*8 :: epss,acc,rmax(4),rmin(4),poten,ugrad,Max_E,Max_E_seed,E_range,dist_tol
real*8 :: E_limit,W_a,alpha,xbeta,Glob_min,ass
integer :: basis_1,basis_2,basis_3,order_1,order_2,order_3,order_4,count3,zz,zz_low,zz4
integer :: support,count7,ddd,myid,natom,ab_flag,ab_flag2,lab,permfac,maxpoints,natom1
integer :: natom2,nbdist,order_1_min,order_2_min,order_3_min,order_4_min,count_seed
integer :: low_grid,subzero,dist_flag,exch,flip1,flip2,symparts
character(len=3),allocatable :: symb(:)

END MODULE dynamic_parameters


MODULE nrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0D0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	TYPE sprs2_sp
		INTEGER(I4B) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
	TYPE sprs2_dp
		INTEGER(I4B) :: n,len
		REAL(DP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE nrtype


MODULE nr
	INTERFACE
		SUBROUTINE airy(x,ai,bi,aip,bip)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: ai,bi,aip,bip
		END SUBROUTINE airy
	END INTERFACE
	INTERFACE
		SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
		USE nrtype
		INTEGER(I4B), INTENT(INOUT) :: iter
		REAL(SP), INTENT(INOUT) :: yb
		REAL(SP), INTENT(IN) :: ftol,temptr
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE amebsa
	END INTERFACE
	INTERFACE
		SUBROUTINE amoeba(p,y,ftol,func,iter)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE amoeba
	END INTERFACE
	INTERFACE
		SUBROUTINE anneal(x,y,iorder)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		END SUBROUTINE anneal
	END INTERFACE
	INTERFACE
		SUBROUTINE asolve(b,x,itrnsp)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: b
		REAL(DP), DIMENSION(:), INTENT(OUT) :: x
		INTEGER(I4B), INTENT(IN) :: itrnsp
		END SUBROUTINE asolve
	END INTERFACE
	INTERFACE
		SUBROUTINE atimes(x,r,itrnsp)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: r
		INTEGER(I4B), INTENT(IN) :: itrnsp
		END SUBROUTINE atimes
	END INTERFACE
	INTERFACE
		SUBROUTINE avevar(data,ave,var)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), INTENT(OUT) :: ave,var
		END SUBROUTINE avevar
	END INTERFACE
	INTERFACE
		SUBROUTINE balanc(a)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE balanc
	END INTERFACE
	INTERFACE
		SUBROUTINE banbks(a,m1,m2,al,indx,b)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE banbks
	END INTERFACE
	INTERFACE
		SUBROUTINE bandec(a,m1,m2,al,indx,d)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
		REAL(SP), INTENT(OUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
		END SUBROUTINE bandec
	END INTERFACE
	INTERFACE
		SUBROUTINE banmul(a,m1,m2,x,b)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: b
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		END SUBROUTINE banmul
	END INTERFACE
	INTERFACE
		SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
		USE nrtype
		REAL(SP), INTENT(IN) :: d1,d2
		REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
		REAL(SP), DIMENSION(4,4), INTENT(OUT) :: c
		END SUBROUTINE bcucof
	END INTERFACE
	INTERFACE
		SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,&
			ansy1,ansy2)
		USE nrtype
		REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
		REAL(SP), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
		REAL(SP), INTENT(OUT) :: ansy,ansy1,ansy2
		END SUBROUTINE bcuint
	END INTERFACE
	INTERFACE beschb
		SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
		USE nrtype
		REAL(DP), INTENT(IN) :: x
		REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
		END SUBROUTINE beschb_s
!BL
		SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
		END SUBROUTINE beschb_v
	END INTERFACE
	INTERFACE bessi
		FUNCTION bessi_s(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi_s
		END FUNCTION bessi_s
!BL
		FUNCTION bessi_v(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi_v
		END FUNCTION bessi_v
	END INTERFACE
	INTERFACE bessi0
		FUNCTION bessi0_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi0_s
		END FUNCTION bessi0_s
!BL
		FUNCTION bessi0_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi0_v
		END FUNCTION bessi0_v
	END INTERFACE
	INTERFACE bessi1
		FUNCTION bessi1_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi1_s
		END FUNCTION bessi1_s
!BL
		FUNCTION bessi1_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi1_v
		END FUNCTION bessi1_v
	END INTERFACE
	INTERFACE
		SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,xnu
		REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
		END SUBROUTINE bessik
	END INTERFACE
	INTERFACE bessj
		FUNCTION bessj_s(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj_s
		END FUNCTION bessj_s
!BL
		FUNCTION bessj_v(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj_v
		END FUNCTION bessj_v
	END INTERFACE
	INTERFACE bessj0
		FUNCTION bessj0_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj0_s
		END FUNCTION bessj0_s
!BL
		FUNCTION bessj0_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj0_v
		END FUNCTION bessj0_v
	END INTERFACE
	INTERFACE bessj1
		FUNCTION bessj1_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj1_s
		END FUNCTION bessj1_s
!BL
		FUNCTION bessj1_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj1_v
		END FUNCTION bessj1_v
	END INTERFACE
	INTERFACE bessjy
		SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,xnu
		REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
		END SUBROUTINE bessjy_s
!BL
		SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
		USE nrtype
		REAL(SP), INTENT(IN) :: xnu
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
		END SUBROUTINE bessjy_v
	END INTERFACE
	INTERFACE bessk
		FUNCTION bessk_s(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk_s
		END FUNCTION bessk_s
!BL
		FUNCTION bessk_v(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk_v
		END FUNCTION bessk_v
	END INTERFACE
	INTERFACE bessk0
		FUNCTION bessk0_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk0_s
		END FUNCTION bessk0_s
!BL
		FUNCTION bessk0_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk0_v
		END FUNCTION bessk0_v
	END INTERFACE
	INTERFACE bessk1
		FUNCTION bessk1_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk1_s
		END FUNCTION bessk1_s
!BL
		FUNCTION bessk1_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk1_v
		END FUNCTION bessk1_v
	END INTERFACE
	INTERFACE bessy
		FUNCTION bessy_s(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy_s
		END FUNCTION bessy_s
!BL
		FUNCTION bessy_v(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy_v
		END FUNCTION bessy_v
	END INTERFACE
	INTERFACE bessy0
		FUNCTION bessy0_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy0_s
		END FUNCTION bessy0_s
!BL
		FUNCTION bessy0_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy0_v
		END FUNCTION bessy0_v
	END INTERFACE
	INTERFACE bessy1
		FUNCTION bessy1_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy1_s
		END FUNCTION bessy1_s
!BL
		FUNCTION bessy1_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy1_v
		END FUNCTION bessy1_v
	END INTERFACE
	INTERFACE beta
		FUNCTION beta_s(z,w)
		USE nrtype
		REAL(SP), INTENT(IN) :: z,w
		REAL(SP) :: beta_s
		END FUNCTION beta_s
!BL
		FUNCTION beta_v(z,w)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
		REAL(SP), DIMENSION(size(z)) :: beta_v
		END FUNCTION beta_v
	END INTERFACE
	INTERFACE betacf
		FUNCTION betacf_s(a,b,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP) :: betacf_s
		END FUNCTION betacf_s
!BL
		FUNCTION betacf_v(a,b,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(size(x)) :: betacf_v
		END FUNCTION betacf_v
	END INTERFACE
	INTERFACE betai
		FUNCTION betai_s(a,b,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP) :: betai_s
		END FUNCTION betai_s
!BL
		FUNCTION betai_v(a,b,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(size(a)) :: betai_v
		END FUNCTION betai_v
	END INTERFACE
	INTERFACE bico
		FUNCTION bico_s(n,k)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n,k
		REAL(SP) :: bico_s
		END FUNCTION bico_s
!BL
		FUNCTION bico_v(n,k)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n,k
		REAL(SP), DIMENSION(size(n)) :: bico_v
		END FUNCTION bico_v
	END INTERFACE
	INTERFACE
		FUNCTION bnldev(pp,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: pp
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: bnldev
		END FUNCTION bnldev
	END INTERFACE
	INTERFACE
		FUNCTION brent(ax,bx,cx,func,tol,xmin)
		USE nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: brent
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION brent
	END INTERFACE
	INTERFACE
		SUBROUTINE broydn(x,check)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		LOGICAL(LGT), INTENT(OUT) :: check
		END SUBROUTINE broydn
	END INTERFACE
	INTERFACE
		SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE bsstep
	END INTERFACE
	INTERFACE
		SUBROUTINE caldat(julian,mm,id,iyyy)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: julian
		INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
		END SUBROUTINE caldat
	END INTERFACE
	INTERFACE
		FUNCTION chder(a,b,c)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chder
		END FUNCTION chder
	END INTERFACE
	INTERFACE chebev
		FUNCTION chebev_s(a,b,c,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP) :: chebev_s
		END FUNCTION chebev_s
!BL
		FUNCTION chebev_v(a,b,c,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
		REAL(SP), DIMENSION(size(x)) :: chebev_v
		END FUNCTION chebev_v
	END INTERFACE
	INTERFACE
		FUNCTION chebft(a,b,n,func)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: chebft
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION chebft
	END INTERFACE
	INTERFACE
		FUNCTION chebpc(c)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chebpc
		END FUNCTION chebpc
	END INTERFACE
	INTERFACE
		FUNCTION chint(a,b,c)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chint
		END FUNCTION chint
	END INTERFACE
	INTERFACE
		SUBROUTINE choldc(a,p)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: p
		END SUBROUTINE choldc
	END INTERFACE
	INTERFACE
		SUBROUTINE cholsl(a,p,b,x)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		END SUBROUTINE cholsl
	END INTERFACE
	INTERFACE
		SUBROUTINE chsone(bins,ebins,knstrn,df,chsq,prob)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: knstrn
		REAL(SP), INTENT(OUT) :: df,chsq,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
		END SUBROUTINE chsone
	END INTERFACE
	INTERFACE
		SUBROUTINE chstwo(bins1,bins2,knstrn,df,chsq,prob)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: knstrn
		REAL(SP), INTENT(OUT) :: df,chsq,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: bins1,bins2
		END SUBROUTINE chstwo
	END INTERFACE
	INTERFACE
		SUBROUTINE cisi(x,ci,si)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: ci,si
		END SUBROUTINE cisi
	END INTERFACE
	INTERFACE
		SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
		USE nrtype
		INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
		REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
		END SUBROUTINE cntab1
	END INTERFACE
	INTERFACE
		SUBROUTINE cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
		USE nrtype
		INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
		REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
		END SUBROUTINE cntab2
	END INTERFACE
	INTERFACE
		FUNCTION convlv(data,respns,isign)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), DIMENSION(:), INTENT(IN) :: respns
		INTEGER(I4B), INTENT(IN) :: isign
		REAL(SP), DIMENSION(size(data)) :: convlv
		END FUNCTION convlv
	END INTERFACE
	INTERFACE
		FUNCTION correl(data1,data2)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), DIMENSION(size(data1)) :: correl
		END FUNCTION correl
	END INTERFACE
	INTERFACE
		SUBROUTINE cosft1(y)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		END SUBROUTINE cosft1
	END INTERFACE
	INTERFACE
		SUBROUTINE cosft2(y,isign)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE cosft2
	END INTERFACE
	INTERFACE
		SUBROUTINE covsrt(covar,maska)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		END SUBROUTINE covsrt
	END INTERFACE
	INTERFACE
		SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
		REAL(SP), INTENT(IN) :: alpha,beta
		REAL(SP), DIMENSION(:), INTENT(OUT):: x
		END SUBROUTINE cyclic
	END INTERFACE
	INTERFACE
		SUBROUTINE daub4(a,isign)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE daub4
	END INTERFACE
	INTERFACE dawson
		FUNCTION dawson_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: dawson_s
		END FUNCTION dawson_s
!BL
		FUNCTION dawson_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dawson_v
		END FUNCTION dawson_v
	END INTERFACE
	INTERFACE
		FUNCTION dbrent(ax,bx,cx,func,dbrent_dfunc,tol,xmin)
		USE nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: dbrent
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
!BL
			FUNCTION dbrent_dfunc(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: dbrent_dfunc
			END FUNCTION dbrent_dfunc
		END INTERFACE
		END FUNCTION dbrent
	END INTERFACE
	INTERFACE
		SUBROUTINE ddpoly(c,x,pd)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
		END SUBROUTINE ddpoly
	END INTERFACE
	INTERFACE
		FUNCTION decchk(string,ch)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
		CHARACTER(1), INTENT(OUT) :: ch
		LOGICAL(LGT) :: decchk
		END FUNCTION decchk
	END INTERFACE
	INTERFACE
		SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: gtol
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(p)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: p
			REAL(SP) :: func
			END FUNCTION func
!BL
			FUNCTION dfunc(p)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: p
			REAL(SP), DIMENSION(size(p)) :: dfunc
			END FUNCTION dfunc
		END INTERFACE
		END SUBROUTINE dfpmin
	END INTERFACE
	INTERFACE
		FUNCTION dfridr(func,x,h,err)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), INTENT(OUT) :: err
		REAL(SP) :: dfridr
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION dfridr
	END INTERFACE
	INTERFACE
		SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
		USE nrtype
		REAL(SP), INTENT(IN) :: w,delta,a,b
		REAL(SP), INTENT(OUT) :: corre,corim,corfac
		REAL(SP), DIMENSION(:), INTENT(IN) :: endpts
		END SUBROUTINE dftcor
	END INTERFACE
	INTERFACE
		SUBROUTINE dftint(func,a,b,w,cosint,sinint)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b,w
		REAL(SP), INTENT(OUT) :: cosint,sinint
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE dftint
	END INTERFACE
	INTERFACE
		SUBROUTINE difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: is1,isf,jsf,k,k1,k2
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: s
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: y
		END SUBROUTINE difeq
	END INTERFACE
	INTERFACE
		FUNCTION eclass(lista,listb,n)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: lista,listb
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), DIMENSION(n) :: eclass
		END FUNCTION eclass
	END INTERFACE
	INTERFACE
		FUNCTION eclazz(equiv,n)
		USE nrtype
		INTERFACE
			FUNCTION equiv(i,j)
			USE nrtype
			LOGICAL(LGT) :: equiv
			INTEGER(I4B), INTENT(IN) :: i,j
			END FUNCTION equiv
		END INTERFACE
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), DIMENSION(n) :: eclazz
		END FUNCTION eclazz
	END INTERFACE
	INTERFACE
		FUNCTION ei(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: ei
		END FUNCTION ei
	END INTERFACE
	INTERFACE
		SUBROUTINE eigsrt(d,v)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
		END SUBROUTINE eigsrt
	END INTERFACE
	INTERFACE elle
		FUNCTION elle_s(phi,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: phi,ak
		REAL(SP) :: elle_s
		END FUNCTION elle_s
!BL
		FUNCTION elle_v(phi,ak)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
		REAL(SP), DIMENSION(size(phi)) :: elle_v
		END FUNCTION elle_v
	END INTERFACE
	INTERFACE ellf
		FUNCTION ellf_s(phi,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: phi,ak
		REAL(SP) :: ellf_s
		END FUNCTION ellf_s
!BL
		FUNCTION ellf_v(phi,ak)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
		REAL(SP), DIMENSION(size(phi)) :: ellf_v
		END FUNCTION ellf_v
	END INTERFACE
	INTERFACE ellpi
		FUNCTION ellpi_s(phi,en,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: phi,en,ak
		REAL(SP) :: ellpi_s
		END FUNCTION ellpi_s
!BL
		FUNCTION ellpi_v(phi,en,ak)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
		REAL(SP), DIMENSION(size(phi)) :: ellpi_v
		END FUNCTION ellpi_v
	END INTERFACE
	INTERFACE
		SUBROUTINE elmhes(a)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE elmhes
	END INTERFACE
	INTERFACE erf
		FUNCTION erf_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erf_s
		END FUNCTION erf_s
!BL
		FUNCTION erf_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erf_v
		END FUNCTION erf_v
	END INTERFACE
	INTERFACE erfc
		FUNCTION erfc_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erfc_s
		END FUNCTION erfc_s
!BL
		FUNCTION erfc_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erfc_v
		END FUNCTION erfc_v
	END INTERFACE
	INTERFACE erfcc
		FUNCTION erfcc_s(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erfcc_s
		END FUNCTION erfcc_s
!BL
		FUNCTION erfcc_v(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erfcc_v
		END FUNCTION erfcc_v
	END INTERFACE
	INTERFACE
		SUBROUTINE eulsum(sum,term,jterm)
		USE nrtype
		REAL(SP), INTENT(INOUT) :: sum
		REAL(SP), INTENT(IN) :: term
		INTEGER(I4B), INTENT(IN) :: jterm
		END SUBROUTINE eulsum
	END INTERFACE
	INTERFACE
		FUNCTION evlmem(fdt,d,xms)
		USE nrtype
		REAL(SP), INTENT(IN) :: fdt,xms
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP) :: evlmem
		END FUNCTION evlmem
	END INTERFACE
	INTERFACE expdev
		SUBROUTINE expdev_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE expdev_s
!BL
		SUBROUTINE expdev_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE expdev_v
	END INTERFACE
	INTERFACE
		FUNCTION expint(n,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: expint
		END FUNCTION expint
	END INTERFACE
	INTERFACE factln
		FUNCTION factln_s(n)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: factln_s
		END FUNCTION factln_s
!BL
		FUNCTION factln_v(n)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
		REAL(SP), DIMENSION(size(n)) :: factln_v
		END FUNCTION factln_v
	END INTERFACE
	INTERFACE factrl
		FUNCTION factrl_s(n)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: factrl_s
		END FUNCTION factrl_s
!BL
		FUNCTION factrl_v(n)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
		REAL(SP), DIMENSION(size(n)) :: factrl_v
		END FUNCTION factrl_v
	END INTERFACE
	INTERFACE
		SUBROUTINE fasper(x,y,ofac,hifac,px,py,jmax,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(IN) :: ofac,hifac
		INTEGER(I4B), INTENT(OUT) :: jmax
		REAL(SP), INTENT(OUT) :: prob
		REAL(SP), DIMENSION(:), POINTER :: px,py
		END SUBROUTINE fasper
	END INTERFACE
	INTERFACE
		SUBROUTINE fdjac(x,fvec,df)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
		END SUBROUTINE fdjac
	END INTERFACE
	INTERFACE
		SUBROUTINE fgauss(x,a,y,dyda)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
		END SUBROUTINE fgauss
	END INTERFACE
	INTERFACE
		SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
		END SUBROUTINE fit
	END INTERFACE
	INTERFACE
		SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
		REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
		END SUBROUTINE fitexy
	END INTERFACE
	INTERFACE
		SUBROUTINE fixrts(d)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		END SUBROUTINE fixrts
	END INTERFACE
	INTERFACE
		FUNCTION fleg(x,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: fleg
		END FUNCTION fleg
	END INTERFACE
	INTERFACE
		SUBROUTINE flmoon(n,nph,jd,frac)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n,nph
		INTEGER(I4B), INTENT(OUT) :: jd
		REAL(SP), INTENT(OUT) :: frac
		END SUBROUTINE flmoon
	END INTERFACE
	INTERFACE four1
!BL
		SUBROUTINE four1_sp(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE four1_alt(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE four1_gather(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_gather
	END INTERFACE
	INTERFACE
		SUBROUTINE four2(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B),INTENT(IN) :: isign
		END SUBROUTINE four2
	END INTERFACE
	INTERFACE
		SUBROUTINE four2_alt(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four2_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE four3(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B),INTENT(IN) :: isign
		END SUBROUTINE four3
	END INTERFACE
	INTERFACE
		SUBROUTINE four3_alt(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four3_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE fourcol(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourcol
	END INTERFACE
	INTERFACE
		SUBROUTINE fourcol_3d(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourcol_3d
	END INTERFACE
	INTERFACE
		SUBROUTINE fourn_gather(data,nn,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourn_gather
	END INTERFACE
	INTERFACE
!BL
		SUBROUTINE fourrow_sp(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourrow_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE fourrow_3d(data,isign)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourrow_3d
	END INTERFACE
	INTERFACE
		FUNCTION fpoly(x,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: fpoly
		END FUNCTION fpoly
	END INTERFACE
	INTERFACE
		SUBROUTINE fred2(a,b,t,f,w,g,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
		INTERFACE
			FUNCTION g(t)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t
			REAL(SP), DIMENSION(size(t)) :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(size(t),size(s)) :: ak
			END FUNCTION ak
		END INTERFACE
		END SUBROUTINE fred2
	END INTERFACE
	INTERFACE
		FUNCTION fredin(x,a,b,t,f,w,g,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
		REAL(SP), DIMENSION(size(x)) :: fredin
		INTERFACE
			FUNCTION g(t)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t
			REAL(SP), DIMENSION(size(t)) :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(size(t),size(s)) :: ak
			END FUNCTION ak
		END INTERFACE
		END FUNCTION fredin
	END INTERFACE
	INTERFACE
		SUBROUTINE frenel(x,s,c)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: s,c
		END SUBROUTINE frenel
	END INTERFACE
	INTERFACE
		SUBROUTINE frprmn(p,ftol,iter,fret)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		END SUBROUTINE frprmn
	END INTERFACE
	INTERFACE
		SUBROUTINE ftest(data1,data2,f,prob)
		USE nrtype
		REAL(SP), INTENT(OUT) :: f,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE ftest
	END INTERFACE
	INTERFACE
		FUNCTION gamdev(ia)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: ia
		REAL(SP) :: gamdev
		END FUNCTION gamdev
	END INTERFACE
	INTERFACE gammln
		FUNCTION gammln_s(xx)
		USE nrtype
		REAL(SP), INTENT(IN) :: xx
		REAL(SP) :: gammln_s
		END FUNCTION gammln_s
!BL
		FUNCTION gammln_v(xx)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		REAL(SP), DIMENSION(size(xx)) :: gammln_v
		END FUNCTION gammln_v
	END INTERFACE
	INTERFACE gammp
		FUNCTION gammp_s(a,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP) :: gammp_s
		END FUNCTION gammp_s
!BL
		FUNCTION gammp_v(a,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(size(a)) :: gammp_v
		END FUNCTION gammp_v
	END INTERFACE
	INTERFACE gammq
		FUNCTION gammq_s(a,x)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP) :: gammq_s
		END FUNCTION gammq_s
!BL
		FUNCTION gammq_v(a,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(size(a)) :: gammq_v
		END FUNCTION gammq_v
	END INTERFACE
	INTERFACE gasdev
		SUBROUTINE gasdev_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE gasdev_s
!BL
		SUBROUTINE gasdev_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE gasdev_v
	END INTERFACE
	INTERFACE
		SUBROUTINE gaucof(a,b,amu0,x,w)
		USE nrtype
		REAL(SP), INTENT(IN) :: amu0
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaucof
	END INTERFACE
	INTERFACE
		SUBROUTINE gauher(x,w)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gauher
	END INTERFACE
	INTERFACE
		SUBROUTINE gaujac(x,w,alf,bet)
		USE nrtype
		REAL(SP), INTENT(IN) :: alf,bet
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaujac
	END INTERFACE
	INTERFACE
		SUBROUTINE gaulag(x,w,alf)
		USE nrtype
		REAL(SP), INTENT(IN) :: alf
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaulag
	END INTERFACE
	INTERFACE
		SUBROUTINE gauleg(x1,x2,x,w)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gauleg
	END INTERFACE
	INTERFACE
		SUBROUTINE gaussj(a,b)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
		END SUBROUTINE gaussj
	END INTERFACE
	INTERFACE gcf
		FUNCTION gcf_s(a,x,gln)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP) :: gcf_s
		END FUNCTION gcf_s
!BL
		FUNCTION gcf_v(a,x,gln)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP), DIMENSION(size(a)) :: gcf_v
		END FUNCTION gcf_v
	END INTERFACE
	INTERFACE
		FUNCTION golden(ax,bx,cx,func,tol,xmin)
		USE nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: golden
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION golden
	END INTERFACE
	INTERFACE gser
		FUNCTION gser_s(a,x,gln)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP) :: gser_s
		END FUNCTION gser_s
!BL
		FUNCTION gser_v(a,x,gln)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP), DIMENSION(size(a)) :: gser_v
		END FUNCTION gser_v
	END INTERFACE
	INTERFACE
		SUBROUTINE hqr(a,wr,wi)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE hqr
	END INTERFACE
	INTERFACE
		SUBROUTINE hunt(xx,x,jlo)
		USE nrtype
		INTEGER(I4B), INTENT(INOUT) :: jlo
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		END SUBROUTINE hunt
	END INTERFACE
	INTERFACE
		SUBROUTINE hypdrv(s,ry,rdyds)
		USE nrtype
		REAL(SP), INTENT(IN) :: s
		REAL(SP), DIMENSION(:), INTENT(IN) :: ry
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
		END SUBROUTINE hypdrv
	END INTERFACE
	INTERFACE
		FUNCTION hypgeo(a,b,c,z)
		USE nrtype
		COMPLEX(SPC), INTENT(IN) :: a,b,c,z
		COMPLEX(SPC) :: hypgeo
		END FUNCTION hypgeo
	END INTERFACE
	INTERFACE
		SUBROUTINE hypser(a,b,c,z,series,deriv)
		USE nrtype
		COMPLEX(SPC), INTENT(IN) :: a,b,c,z
		COMPLEX(SPC), INTENT(OUT) :: series,deriv
		END SUBROUTINE hypser
	END INTERFACE
	INTERFACE
		FUNCTION icrc(crc,buf,jinit,jrev)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
		INTEGER(I2B), INTENT(IN) :: crc,jinit
		INTEGER(I4B), INTENT(IN) :: jrev
		INTEGER(I2B) :: icrc
		END FUNCTION icrc
	END INTERFACE
	INTERFACE
		FUNCTION igray(n,is)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n,is
		INTEGER(I4B) :: igray
		END FUNCTION igray
	END INTERFACE
	INTERFACE
		RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: index
		INTEGER, OPTIONAL, INTENT(IN) :: partial
		END SUBROUTINE index_bypack
	END INTERFACE
	INTERFACE indexx
		SUBROUTINE indexx_sp(arr,index)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_sp
		SUBROUTINE indexx_i4b(iarr,index)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_i4b
	END INTERFACE
	INTERFACE
		FUNCTION interp(uc)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
		REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
		END FUNCTION interp
	END INTERFACE
	INTERFACE
		FUNCTION rank(indx)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		INTEGER(I4B), DIMENSION(size(indx)) :: rank
		END FUNCTION rank
	END INTERFACE
	INTERFACE
		FUNCTION irbit1(iseed)
		USE nrtype
		INTEGER(I4B), INTENT(INOUT) :: iseed
		INTEGER(I4B) :: irbit1
		END FUNCTION irbit1
	END INTERFACE
	INTERFACE
		FUNCTION irbit2(iseed)
		USE nrtype
		INTEGER(I4B), INTENT(INOUT) :: iseed
		INTEGER(I4B) :: irbit2
		END FUNCTION irbit2
	END INTERFACE
	INTERFACE
		SUBROUTINE jacobi(a,d,v,nrot)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: nrot
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		END SUBROUTINE jacobi
	END INTERFACE
	INTERFACE
		SUBROUTINE jacobn(x,y,dfdx,dfdy)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
		END SUBROUTINE jacobn
	END INTERFACE
	INTERFACE
		FUNCTION julday(mm,id,iyyy)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: mm,id,iyyy
		INTEGER(I4B) :: julday
		END FUNCTION julday
	END INTERFACE
	INTERFACE
		SUBROUTINE kendl1(data1,data2,tau,z,prob)
		USE nrtype
		REAL(SP), INTENT(OUT) :: tau,z,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE kendl1
	END INTERFACE
	INTERFACE
		SUBROUTINE kendl2(tab,tau,z,prob)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
		REAL(SP), INTENT(OUT) :: tau,z,prob
		END SUBROUTINE kendl2
	END INTERFACE
	INTERFACE
		FUNCTION kermom(y,m)
		USE nrtype
		REAL(DP), INTENT(IN) :: y
		INTEGER(I4B), INTENT(IN) :: m
		REAL(DP), DIMENSION(m) :: kermom
		END FUNCTION kermom
	END INTERFACE
	INTERFACE
		SUBROUTINE ks2d1s(x1,y1,quadvl,d1,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
		REAL(SP), INTENT(OUT) :: d1,prob
		INTERFACE
			SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
			USE nrtype
			REAL(SP), INTENT(IN) :: x,y
			REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
			END SUBROUTINE quadvl
		END INTERFACE
		END SUBROUTINE ks2d1s
	END INTERFACE
	INTERFACE
		SUBROUTINE ks2d2s(x1,y1,x2,y2,d,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
		REAL(SP), INTENT(OUT) :: d,prob
		END SUBROUTINE ks2d2s
	END INTERFACE
	INTERFACE
		SUBROUTINE ksone(data,func,d,prob)
		USE nrtype
		REAL(SP), INTENT(OUT) :: d,prob
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE ksone
	END INTERFACE
	INTERFACE
		SUBROUTINE kstwo(data1,data2,d,prob)
		USE nrtype
		REAL(SP), INTENT(OUT) :: d,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE kstwo
	END INTERFACE
	INTERFACE
		SUBROUTINE laguer(a,x,its)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: its
		COMPLEX(SPC), INTENT(INOUT) :: x
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
		END SUBROUTINE laguer
	END INTERFACE
	INTERFACE
		SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
		REAL(SP), INTENT(OUT) :: chisq
		INTERFACE
			SUBROUTINE funcs(x,arr)
			USE nrtype
			REAL(SP),INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE lfit
	END INTERFACE
	INTERFACE
		SUBROUTINE linbcg(b,x,itol,tol,itmax,iter,err)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: b
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
		INTEGER(I4B), INTENT(IN) :: itol,itmax
		REAL(DP), INTENT(IN) :: tol
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(DP), INTENT(OUT) :: err
		END SUBROUTINE linbcg
	END INTERFACE
	INTERFACE
		SUBROUTINE dlinmin(p,xi,fret)
		USE nrtype
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
		END SUBROUTINE dlinmin
	END INTERFACE
	INTERFACE
		SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		REAL(SP), INTENT(IN) :: fold,stpmax
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		REAL(SP), INTENT(OUT) :: f
		LOGICAL(LGT), INTENT(OUT) :: check
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP) :: func
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE lnsrch
	END INTERFACE
	INTERFACE
		FUNCTION locate(xx,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B) :: locate
		END FUNCTION locate
	END INTERFACE
	INTERFACE
		FUNCTION lop(u)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
		REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
		END FUNCTION lop
	END INTERFACE
	INTERFACE
		SUBROUTINE lubksb(a,indx,b)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE lubksb
	END INTERFACE
	INTERFACE
		SUBROUTINE ludcmp(a,indx,d)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
		REAL(SP), INTENT(OUT) :: d
		END SUBROUTINE ludcmp
	END INTERFACE
	INTERFACE
		SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
			maxexp,eps,epsneg,xmin,xmax)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp,&
			minexp,negep,ngrd
		REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
		END SUBROUTINE machar
	END INTERFACE
	INTERFACE
		SUBROUTINE medfit(x,y,a,b,abdev)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: a,b,abdev
		END SUBROUTINE medfit
	END INTERFACE
	INTERFACE
		SUBROUTINE memcof(data,xms,d)
		USE nrtype
		REAL(SP), INTENT(OUT) :: xms
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d
		END SUBROUTINE memcof
	END INTERFACE
	INTERFACE
		SUBROUTINE mgfas(u,maxcyc)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		INTEGER(I4B), INTENT(IN) :: maxcyc
		END SUBROUTINE mgfas
	END INTERFACE
	INTERFACE
		SUBROUTINE mglin(u,ncycle)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		INTEGER(I4B), INTENT(IN) :: ncycle
		END SUBROUTINE mglin
	END INTERFACE
	INTERFACE
		SUBROUTINE midexp(funk,aa,bb,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midexp
	END INTERFACE
	INTERFACE
		SUBROUTINE midinf(funk,aa,bb,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midinf
	END INTERFACE
	INTERFACE
		SUBROUTINE midpnt(func,a,b,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE midpnt
	END INTERFACE
	INTERFACE
		SUBROUTINE midsql(funk,aa,bb,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midsql
	END INTERFACE
	INTERFACE
		SUBROUTINE midsqu(funk,aa,bb,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midsqu
	END INTERFACE
	INTERFACE
		RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
		USE nrtype
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP) :: func
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			END FUNCTION func
		END INTERFACE
		REAL(SP), DIMENSION(:), INTENT(IN) :: regn
		INTEGER(I4B), INTENT(IN) :: ndim,npts
		REAL(SP), INTENT(IN) :: dith
		REAL(SP), INTENT(OUT) :: ave,var
		END SUBROUTINE miser
	END INTERFACE
	INTERFACE
		SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), INTENT(IN) :: xs,htot
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE mmid
	END INTERFACE
	INTERFACE
		SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
		USE nrtype
		REAL(SP), INTENT(INOUT) :: ax,bx
		REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE mnbrak
	END INTERFACE
	INTERFACE
		SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(SP), INTENT(IN) :: tolx,tolf
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
			END SUBROUTINE usrfun
		END INTERFACE
		END SUBROUTINE mnewt
	END INTERFACE
	INTERFACE
		SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
		USE nrtype
		REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		END SUBROUTINE moment
	END INTERFACE
	INTERFACE
		SUBROUTINE mp2dfr(a,s,n,m)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), INTENT(OUT) :: m
		CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
		END SUBROUTINE mp2dfr
	END INTERFACE
	INTERFACE
		SUBROUTINE mpdiv(q,r,u,v,n,m)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpdiv
	END INTERFACE
	INTERFACE
		SUBROUTINE mpinv(u,v,n,m)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpinv
	END INTERFACE
	INTERFACE
		SUBROUTINE mpmul(w,u,v,n,m)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpmul
	END INTERFACE
	INTERFACE
		SUBROUTINE mppi(n)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		END SUBROUTINE mppi
	END INTERFACE
	INTERFACE
		SUBROUTINE mprove(a,alud,indx,b,x)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:), INTENT(IN) :: b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		END SUBROUTINE mprove
	END INTERFACE
	INTERFACE
		SUBROUTINE mpsqrt(w,u,v,n,m)
		USE nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpsqrt
	END INTERFACE
	INTERFACE
		SUBROUTINE mrqcof(x,y,sig,a,maska,alpha,beta,chisq,funcs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,a,sig
		REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
		REAL(SP), INTENT(OUT) :: chisq
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		INTERFACE
			SUBROUTINE funcs(x,a,yfit,dyda)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
			REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE mrqcof
	END INTERFACE
	INTERFACE
		SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
		REAL(SP), INTENT(OUT) :: chisq
		REAL(SP), INTENT(INOUT) :: alamda
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		INTERFACE
			SUBROUTINE funcs(x,a,yfit,dyda)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
			REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE mrqmin
	END INTERFACE
	INTERFACE
		SUBROUTINE newt(x,check)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		LOGICAL(LGT), INTENT(OUT) :: check
		END SUBROUTINE newt
	END INTERFACE
	INTERFACE
		SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
		REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
!BL
			SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
			REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
			REAL(SP), INTENT(INOUT) :: x
			REAL(SP), INTENT(IN) :: htry,eps
			REAL(SP), INTENT(OUT) :: hdid,hnext
				INTERFACE
				SUBROUTINE derivs(x,y,dydx)
					USE nrtype
					REAL(SP), INTENT(IN) :: x
					REAL(SP), DIMENSION(:), INTENT(IN) :: y
					REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
					END SUBROUTINE derivs
				END INTERFACE
			END SUBROUTINE rkqs
		END INTERFACE
		END SUBROUTINE odeint
	END INTERFACE
	INTERFACE
		SUBROUTINE orthog(anu,alpha,beta,a,b)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
		REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
		END SUBROUTINE orthog
	END INTERFACE
	INTERFACE
		SUBROUTINE pade(cof,resid)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
		REAL(SP), INTENT(OUT) :: resid
		END SUBROUTINE pade
	END INTERFACE
	INTERFACE
		FUNCTION pccheb(d)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP), DIMENSION(size(d)) :: pccheb
		END FUNCTION pccheb
	END INTERFACE
	INTERFACE
		SUBROUTINE pcshft(a,b,d)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		END SUBROUTINE pcshft
	END INTERFACE
	INTERFACE
		SUBROUTINE pearsn(x,y,r,prob,z)
		USE nrtype
		REAL(SP), INTENT(OUT) :: r,prob,z
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		END SUBROUTINE pearsn
	END INTERFACE
	INTERFACE
		SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob)
		USE nrtype
		INTEGER(I4B), INTENT(OUT) :: jmax
		REAL(SP), INTENT(IN) :: ofac,hifac
		REAL(SP), INTENT(OUT) :: prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), POINTER :: px,py
		END SUBROUTINE period
	END INTERFACE
	INTERFACE plgndr
		FUNCTION plgndr_s(l,m,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: l,m
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: plgndr_s
		END FUNCTION plgndr_s
!BL
		FUNCTION plgndr_v(l,m,x)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: l,m
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: plgndr_v
		END FUNCTION plgndr_v
	END INTERFACE
	INTERFACE
		FUNCTION poidev(xm)
		USE nrtype
		REAL(SP), INTENT(IN) :: xm
		REAL(SP) :: poidev
		END FUNCTION poidev
	END INTERFACE
	INTERFACE
		FUNCTION polcoe(x,y)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(size(x)) :: polcoe
		END FUNCTION polcoe
	END INTERFACE
	INTERFACE
		FUNCTION polcof(xa,ya)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), DIMENSION(size(xa)) :: polcof
		END FUNCTION polcof
	END INTERFACE
	INTERFACE
		SUBROUTINE poldiv(u,v,q,r)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
		REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
		END SUBROUTINE poldiv
	END INTERFACE
	INTERFACE
		SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE polin2
	END INTERFACE
	INTERFACE
		SUBROUTINE polint(xa,ya,x,y,dy)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE polint
	END INTERFACE
	INTERFACE
		SUBROUTINE powell(p,xi,ftol,iter,fret)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), INTENT(OUT) :: fret
		END SUBROUTINE powell
	END INTERFACE
	INTERFACE
		FUNCTION predic(data,d,nfut)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data,d
		INTEGER(I4B), INTENT(IN) :: nfut
		REAL(SP), DIMENSION(nfut) :: predic
		END FUNCTION predic
	END INTERFACE
	INTERFACE
		FUNCTION probks(alam)
		USE nrtype
		REAL(SP), INTENT(IN) :: alam
		REAL(SP) :: probks
		END FUNCTION probks
	END INTERFACE
	INTERFACE psdes
		SUBROUTINE psdes_s(lword,rword)
		USE nrtype
		INTEGER(I4B), INTENT(INOUT) :: lword,rword
		END SUBROUTINE psdes_s
!BL
		SUBROUTINE psdes_v(lword,rword)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
		END SUBROUTINE psdes_v
	END INTERFACE
	INTERFACE
		SUBROUTINE pwt(a,isign)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE pwt
	END INTERFACE
	INTERFACE
		SUBROUTINE pwtset(n)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		END SUBROUTINE pwtset
	END INTERFACE
	INTERFACE pythag
!BL
		FUNCTION pythag_sp(a,b)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: pythag_sp
		END FUNCTION pythag_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE pzextr(iest,xest,yest,yz,dy)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: iest
		REAL(SP), INTENT(IN) :: xest
		REAL(SP), DIMENSION(:), INTENT(IN) :: yest
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
		END SUBROUTINE pzextr
	END INTERFACE
	INTERFACE
		SUBROUTINE qrdcmp(a,c,d,sing)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
		LOGICAL(LGT), INTENT(OUT) :: sing
		END SUBROUTINE qrdcmp
	END INTERFACE
	INTERFACE
		FUNCTION qromb(func,a,b)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qromb
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qromb
	END INTERFACE
	INTERFACE
		FUNCTION qromo(func,a,b,choose)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qromo
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		INTERFACE
			SUBROUTINE choose(funk,aa,bb,s,n)
			USE nrtype
			REAL(SP), INTENT(IN) :: aa,bb
			REAL(SP), INTENT(INOUT) :: s
			INTEGER(I4B), INTENT(IN) :: n
			INTERFACE
				FUNCTION funk(x)
				USE nrtype
				REAL(SP), DIMENSION(:), INTENT(IN) :: x
				REAL(SP), DIMENSION(size(x)) :: funk
				END FUNCTION funk
			END INTERFACE
			END SUBROUTINE choose
		END INTERFACE
		END FUNCTION qromo
	END INTERFACE
	INTERFACE
		SUBROUTINE qroot(p,b,c,eps)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), INTENT(INOUT) :: b,c
		REAL(SP), INTENT(IN) :: eps
		END SUBROUTINE qroot
	END INTERFACE
	INTERFACE
		SUBROUTINE qrsolv(a,c,d,b)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE qrsolv
	END INTERFACE
	INTERFACE
		SUBROUTINE qrupdt(r,qt,u,v)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
		REAL(SP), DIMENSION(:), INTENT(IN) :: v
		END SUBROUTINE qrupdt
	END INTERFACE
	INTERFACE
		FUNCTION qsimp(func,a,b)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qsimp
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qsimp
	END INTERFACE
	INTERFACE
		FUNCTION qtrap(func,a,b)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qtrap
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qtrap
	END INTERFACE
	INTERFACE
		SUBROUTINE quadct(x,y,xx,yy,fa,fb,fc,fd)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx,yy
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadct
	END INTERFACE
	INTERFACE
		SUBROUTINE quadmx(a)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
		END SUBROUTINE quadmx
	END INTERFACE
	INTERFACE
		SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadvl
	END INTERFACE
	INTERFACE
		FUNCTION ran(idum)
		INTEGER(selected_int_kind(9)), INTENT(INOUT) :: idum
		REAL :: ran
		END FUNCTION ran
	END INTERFACE
	INTERFACE ran0
		SUBROUTINE ran0_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran0_s
!BL
		SUBROUTINE ran0_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran0_v
	END INTERFACE
	INTERFACE ran1
		SUBROUTINE ran1_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran1_s
!BL
		SUBROUTINE ran1_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran1_v
	END INTERFACE
	INTERFACE ran2
		SUBROUTINE ran2_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran2_s
!BL
		SUBROUTINE ran2_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran2_v
	END INTERFACE
	INTERFACE ran3
		SUBROUTINE ran3_s(harvest)
		USE nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran3_s
!BL
		SUBROUTINE ran3_v(harvest)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran3_v
	END INTERFACE
	INTERFACE
		SUBROUTINE ratint(xa,ya,x,y,dy)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE ratint
	END INTERFACE
	INTERFACE
		SUBROUTINE ratlsq(func,a,b,mm,kk,cof,dev)
		USE nrtype
		REAL(DP), INTENT(IN) :: a,b
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
		REAL(DP), INTENT(OUT) :: dev
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(DP), DIMENSION(:), INTENT(IN) :: x
			REAL(DP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE ratlsq
	END INTERFACE
	INTERFACE ratval
		FUNCTION ratval_s(x,cof,mm,kk)
		USE nrtype
		REAL(DP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
		REAL(DP) :: ratval_s
		END FUNCTION ratval_s
!BL
		FUNCTION ratval_v(x,cof,mm,kk)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
		REAL(DP), DIMENSION(size(x)) :: ratval_v
		END FUNCTION ratval_v
	END INTERFACE
	INTERFACE rc
		FUNCTION rc_s(x,y)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: rc_s
		END FUNCTION rc_s
!BL
		FUNCTION rc_v(x,y)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(size(x)) :: rc_v
		END FUNCTION rc_v
	END INTERFACE
	INTERFACE rd
		FUNCTION rd_s(x,y,z)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y,z
		REAL(SP) :: rd_s
		END FUNCTION rd_s
!BL
		FUNCTION rd_v(x,y,z)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
		REAL(SP), DIMENSION(size(x)) :: rd_v
		END FUNCTION rd_v
	END INTERFACE
	INTERFACE realft
!BL
		SUBROUTINE realft_sp(data,isign,zdata)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
		END SUBROUTINE realft_sp
	END INTERFACE
	INTERFACE
		RECURSIVE FUNCTION recur1(a,b) RESULT(u)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a)) :: u
		END FUNCTION recur1
	END INTERFACE
	INTERFACE
		FUNCTION recur2(a,b,c)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
		REAL(SP), DIMENSION(size(a)) :: recur2
		END FUNCTION recur2
	END INTERFACE
	INTERFACE
		SUBROUTINE relax(u,rhs)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
		END SUBROUTINE relax
	END INTERFACE
	INTERFACE
		SUBROUTINE relax2(u,rhs)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
		END SUBROUTINE relax2
	END INTERFACE
	INTERFACE
	FUNCTION resid(u,rhs)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
		REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
		END FUNCTION resid
	END INTERFACE
	INTERFACE rf
		FUNCTION rf_s(x,y,z)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y,z
		REAL(SP) :: rf_s
		END FUNCTION rf_s
!BL
		FUNCTION rf_v(x,y,z)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
		REAL(SP), DIMENSION(size(x)) :: rf_v
		END FUNCTION rf_v
	END INTERFACE
	INTERFACE rj
		FUNCTION rj_s(x,y,z,p)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y,z,p
		REAL(SP) :: rj_s
		END FUNCTION rj_s
!BL
		FUNCTION rj_v(x,y,z,p)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
		REAL(SP), DIMENSION(size(x)) :: rj_v
		END FUNCTION rj_v
	END INTERFACE
	INTERFACE
		SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rk4
	END INTERFACE
	INTERFACE
		SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkck
	END INTERFACE
	INTERFACE
		SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: vstart
		REAL(SP), INTENT(IN) :: x1,x2
		INTEGER(I4B), INTENT(IN) :: nstep
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkdumb
	END INTERFACE
	INTERFACE
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkqs
	END INTERFACE
	INTERFACE
		SUBROUTINE rlft2(data,spec,speq,isign)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
		COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: spec
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: speq
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE rlft2
	END INTERFACE
	INTERFACE
		SUBROUTINE rlft3(data,spec,speq,isign)
		USE nrtype
		REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(OUT) :: spec
		COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: speq
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE rlft3
	END INTERFACE
	INTERFACE
		SUBROUTINE rotate(r,qt,i,a,b)
		USE nrtype
		REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
		INTEGER(I4B), INTENT(IN) :: i
		REAL(SP), INTENT(IN) :: a,b
		END SUBROUTINE rotate
	END INTERFACE
	INTERFACE
		SUBROUTINE rsolv(a,d,b)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE rsolv
	END INTERFACE
	INTERFACE
		FUNCTION rstrct(uf)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
		REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
		END FUNCTION rstrct
	END INTERFACE
	INTERFACE
		FUNCTION rtbis(func,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtbis
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtbis
	END INTERFACE
	INTERFACE
		FUNCTION rtflsp(func,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtflsp
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtflsp
	END INTERFACE
	INTERFACE
		FUNCTION rtnewt(funcd,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtnewt
		INTERFACE
			SUBROUTINE funcd(x,fval,fderiv)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), INTENT(OUT) :: fval,fderiv
			END SUBROUTINE funcd
		END INTERFACE
		END FUNCTION rtnewt
	END INTERFACE
	INTERFACE
		FUNCTION rtsafe(funcd,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtsafe
		INTERFACE
			SUBROUTINE funcd(x,fval,fderiv)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), INTENT(OUT) :: fval,fderiv
			END SUBROUTINE funcd
		END INTERFACE
		END FUNCTION rtsafe
	END INTERFACE
	INTERFACE
		FUNCTION rtsec(func,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtsec
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtsec
	END INTERFACE
	INTERFACE
		SUBROUTINE rzextr(iest,xest,yest,yz,dy)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: iest
		REAL(SP), INTENT(IN) :: xest
		REAL(SP), DIMENSION(:), INTENT(IN) :: yest
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
		END SUBROUTINE rzextr
	END INTERFACE
	INTERFACE
		FUNCTION savgol(nl,nrr,ld,m)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
		REAL(SP), DIMENSION(nl+nrr+1) :: savgol
		END FUNCTION savgol
	END INTERFACE
	INTERFACE
		SUBROUTINE scrsho(func)
		USE nrtype
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE scrsho
	END INTERFACE
	INTERFACE
		FUNCTION select(k,arr)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		REAL(SP) :: select
		END FUNCTION select
	END INTERFACE
	INTERFACE
		FUNCTION select_bypack(k,arr)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		REAL(SP) :: select_bypack
		END FUNCTION select_bypack
	END INTERFACE
	INTERFACE
		SUBROUTINE select_heap(arr,heap)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
		END SUBROUTINE select_heap
	END INTERFACE
	INTERFACE
		FUNCTION select_inplace(k,arr)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP) :: select_inplace
		END FUNCTION select_inplace
	END INTERFACE
	INTERFACE
		SUBROUTINE simplx(a,m1,m2,m3,icase,izrov,iposv)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: m1,m2,m3
		INTEGER(I4B), INTENT(OUT) :: icase
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: izrov,iposv
		END SUBROUTINE simplx
	END INTERFACE
	INTERFACE
		SUBROUTINE simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
		USE nrtype
		REAL(SP), INTENT(IN) :: xs,htot
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE simpr
	END INTERFACE
	INTERFACE
		SUBROUTINE sinft(y)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		END SUBROUTINE sinft
	END INTERFACE
	INTERFACE
		SUBROUTINE slvsm2(u,rhs)
		USE nrtype
		REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
		REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
		END SUBROUTINE slvsm2
	END INTERFACE
	INTERFACE
		SUBROUTINE slvsml(u,rhs)
		USE nrtype
		REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
		REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
		END SUBROUTINE slvsml
	END INTERFACE
	INTERFACE
		SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
		USE nrtype
		REAL(SP), INTENT(IN) :: uu,emmc
		REAL(SP), INTENT(OUT) :: sn,cn,dn
		END SUBROUTINE sncndn
	END INTERFACE
	INTERFACE
		FUNCTION snrm(sx,itol)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: sx
		INTEGER(I4B), INTENT(IN) :: itol
		REAL(DP) :: snrm
		END FUNCTION snrm
	END INTERFACE
	INTERFACE
		SUBROUTINE sobseq(x,init)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
		END SUBROUTINE sobseq
	END INTERFACE
	INTERFACE
		SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,nb,y)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: itmax,nb
		REAL(SP), INTENT(IN) :: conv,slowc
		REAL(SP), DIMENSION(:), INTENT(IN) :: scalv
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: y
		END SUBROUTINE solvde
	END INTERFACE
	INTERFACE
		SUBROUTINE sor(a,b,c,d,e,f,u,rjac)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), INTENT(IN) :: rjac
		END SUBROUTINE sor
	END INTERFACE
	INTERFACE
		SUBROUTINE sort(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort
	END INTERFACE
	INTERFACE
		SUBROUTINE sort2(arr,slave)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
		END SUBROUTINE sort2
	END INTERFACE
	INTERFACE
		SUBROUTINE sort3(arr,slave1,slave2)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
		END SUBROUTINE sort3
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_bypack(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_bypack
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_byreshape(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_byreshape
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_heap(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_heap
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_pick(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_pick
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_radix(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_radix
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_shell(arr)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_shell
	END INTERFACE
	INTERFACE
		SUBROUTINE spctrm(p,k,ovrlap,unit,n_window)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: p
		INTEGER(I4B), INTENT(IN) :: k
		LOGICAL(LGT), INTENT(IN) :: ovrlap
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: n_window,unit
		END SUBROUTINE spctrm
	END INTERFACE
	INTERFACE
		SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
		END SUBROUTINE spear
	END INTERFACE
	INTERFACE sphbes
		SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
		END SUBROUTINE sphbes_s
!BL
		SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
		END SUBROUTINE sphbes_v
	END INTERFACE
	INTERFACE
		SUBROUTINE splie2(x1a,x2a,ya,y2a)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
		END SUBROUTINE splie2
	END INTERFACE
	INTERFACE
		FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP) :: splin2
		END FUNCTION splin2
	END INTERFACE
	INTERFACE
		SUBROUTINE spline(x,y,yp1,ypn,y2)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(IN) :: yp1,ypn
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
		END SUBROUTINE spline
	END INTERFACE
	INTERFACE
		FUNCTION splint(xa,ya,y2a,x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: splint
		END FUNCTION splint
	END INTERFACE
	INTERFACE sprsax
		SUBROUTINE sprsax_dp(sa,x,b)
		USE nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION (:), INTENT(IN) :: x
		REAL(DP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprsax_dp
!BL
		SUBROUTINE sprsax_sp(sa,x,b)
		USE nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION (:), INTENT(IN) :: x
		REAL(SP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprsax_sp
	END INTERFACE
	INTERFACE sprsdiag
		SUBROUTINE sprsdiag_dp(sa,b)
		USE nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION(:), INTENT(OUT) :: b
		END SUBROUTINE sprsdiag_dp
!BL
		SUBROUTINE sprsdiag_sp(sa,b)
		USE nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION(:), INTENT(OUT) :: b
		END SUBROUTINE sprsdiag_sp
	END INTERFACE
	INTERFACE sprsin
		SUBROUTINE sprsin_sp(a,thresh,sa)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), INTENT(IN) :: thresh
		TYPE(sprs2_sp), INTENT(OUT) :: sa
		END SUBROUTINE sprsin_sp
!BL
		SUBROUTINE sprsin_dp(a,thresh,sa)
		USE nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(DP), INTENT(IN) :: thresh
		TYPE(sprs2_dp), INTENT(OUT) :: sa
		END SUBROUTINE sprsin_dp
	END INTERFACE
	INTERFACE
		SUBROUTINE sprstp(sa)
		USE nrtype
		TYPE(sprs2_sp), INTENT(INOUT) :: sa
		END SUBROUTINE sprstp
	END INTERFACE
	INTERFACE sprstx
		SUBROUTINE sprstx_dp(sa,x,b)
		USE nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION (:), INTENT(IN) :: x
		REAL(DP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprstx_dp
!BL
		SUBROUTINE sprstx_sp(sa,x,b)
		USE nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION (:), INTENT(IN) :: x
		REAL(SP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprstx_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stifbs
	END INTERFACE
	INTERFACE
		SUBROUTINE stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stiff
	END INTERFACE
	INTERFACE
		SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
		REAL(SP), INTENT(IN) :: xs,htot
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stoerm
	END INTERFACE
	INTERFACE svbksb
!BL
		SUBROUTINE svbksb_sp(u,w,v,b,x)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
		REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		END SUBROUTINE svbksb_sp
	END INTERFACE
	INTERFACE svdcmp
!BL
		SUBROUTINE svdcmp_sp(a,w,v)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		END SUBROUTINE svdcmp_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		REAL(SP), INTENT(OUT) :: chisq
		INTERFACE
			FUNCTION funcs(x,n)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			INTEGER(I4B), INTENT(IN) :: n
			REAL(SP), DIMENSION(n) :: funcs
			END FUNCTION funcs
		END INTERFACE
		END SUBROUTINE svdfit
	END INTERFACE
	INTERFACE
		SUBROUTINE svdvar(v,w,cvm)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
		REAL(SP), DIMENSION(:), INTENT(IN) :: w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
		END SUBROUTINE svdvar
	END INTERFACE
	INTERFACE
		FUNCTION toeplz(r,y)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
		REAL(SP), DIMENSION(size(y)) :: toeplz
		END FUNCTION toeplz
	END INTERFACE
	INTERFACE
		SUBROUTINE tptest(data1,data2,t,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE tptest
	END INTERFACE
	INTERFACE
		SUBROUTINE tqli(d,e,z)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
		REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
		END SUBROUTINE tqli
	END INTERFACE
	INTERFACE
		SUBROUTINE trapzd(func,a,b,s,n)
		USE nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE trapzd
	END INTERFACE
	INTERFACE
		SUBROUTINE tred2(a,d,e,novectors)
		USE nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
		LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
		END SUBROUTINE tred2
	END INTERFACE
!	On a purely serial machine, for greater efficiency, remove
!	the generic name tridag from the following interface,
!	and put it on the next one after that.
	INTERFACE tridag
		RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
		REAL(SP), DIMENSION(:), INTENT(OUT) :: u
		END SUBROUTINE tridag_par
	END INTERFACE
	INTERFACE
		SUBROUTINE tridag_ser(a,b,c,r,u)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
		REAL(SP), DIMENSION(:), INTENT(OUT) :: u
		END SUBROUTINE tridag_ser
	END INTERFACE
	INTERFACE
		SUBROUTINE ttest(data1,data2,t,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE ttest
	END INTERFACE
	INTERFACE
		SUBROUTINE tutest(data1,data2,t,prob)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE tutest
	END INTERFACE
	INTERFACE
		SUBROUTINE twofft(data1,data2,fft1,fft2)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
		END SUBROUTINE twofft
	END INTERFACE
	INTERFACE
		FUNCTION vander(x,q)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
		REAL(DP), DIMENSION(size(x)) :: vander
		END FUNCTION vander
	END INTERFACE
	INTERFACE
		SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: region
		INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
		REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
		INTERFACE
			FUNCTION func(pt,wgt)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: pt
			REAL(SP), INTENT(IN) :: wgt
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE vegas
	END INTERFACE
	INTERFACE
		SUBROUTINE voltra(t0,h,t,f,g,ak)
		USE nrtype
		REAL(SP), INTENT(IN) :: t0,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: t
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: f
		INTERFACE
			FUNCTION g(t)
			USE nrtype
			REAL(SP), INTENT(IN) :: t
			REAL(SP), DIMENSION(:), POINTER :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			USE nrtype
			REAL(SP), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(:,:), POINTER :: ak
			END FUNCTION ak
		END INTERFACE
		END SUBROUTINE voltra
	END INTERFACE
	INTERFACE
		SUBROUTINE wt1(a,isign,wtstep)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		INTERFACE
			SUBROUTINE wtstep(a,isign)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
			INTEGER(I4B), INTENT(IN) :: isign
			END SUBROUTINE wtstep
		END INTERFACE
		END SUBROUTINE wt1
	END INTERFACE
	INTERFACE
		SUBROUTINE wtn(a,nn,isign,wtstep)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
		INTEGER(I4B), INTENT(IN) :: isign
		INTERFACE
			SUBROUTINE wtstep(a,isign)
			USE nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
			INTEGER(I4B), INTENT(IN) :: isign
			END SUBROUTINE wtstep
		END INTERFACE
		END SUBROUTINE wtn
	END INTERFACE
	INTERFACE
		FUNCTION wwghts(n,h,kermom)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: h
		REAL(SP), DIMENSION(n) :: wwghts
		INTERFACE
			FUNCTION kermom(y,m)
			USE nrtype
			REAL(DP), INTENT(IN) :: y
			INTEGER(I4B), INTENT(IN) :: m
			REAL(DP), DIMENSION(m) :: kermom
			END FUNCTION kermom
		END INTERFACE
		END FUNCTION wwghts
	END INTERFACE
	INTERFACE
		SUBROUTINE zbrac(func,x1,x2,succes)
		USE nrtype
		REAL(SP), INTENT(INOUT) :: x1,x2
		LOGICAL(LGT), INTENT(OUT) :: succes
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE zbrac
	END INTERFACE
	INTERFACE
		SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
		USE nrtype
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), INTENT(OUT) :: nb
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE zbrak
	END INTERFACE
	INTERFACE
		FUNCTION zbrent(func,x1,x2,tol)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,tol
		REAL(SP) :: zbrent
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION zbrent
	END INTERFACE
	INTERFACE
		SUBROUTINE zrhqr(a,rtr,rti)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
		END SUBROUTINE zrhqr
	END INTERFACE
	INTERFACE
		FUNCTION zriddr(func,x1,x2,xacc)
		USE nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: zriddr
		INTERFACE
			FUNCTION func(x)
			USE nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION zriddr
	END INTERFACE
	INTERFACE
		SUBROUTINE zroots(a,roots,polish)
		USE nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
		LOGICAL(LGT), INTENT(IN) :: polish
		END SUBROUTINE zroots
	END INTERFACE
END MODULE nr




MODULE nrutil
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
	INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
	INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
	INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
	INTEGER(I4B), PARAMETER :: NPAR_POLY=8
	INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
	INTERFACE array_copy
		MODULE PROCEDURE array_copy_r, array_copy_i
	END INTERFACE
	INTERFACE swap
		MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
			swap_cv,swap_cm, &
			masked_swap_rs,masked_swap_rv,masked_swap_rm
	END INTERFACE
	INTERFACE reallocate
		MODULE PROCEDURE reallocate_rv,reallocate_rm,&
			reallocate_iv,reallocate_im,reallocate_hv
	END INTERFACE
	INTERFACE imaxloc
		MODULE PROCEDURE imaxloc_r,imaxloc_i
	END INTERFACE
	INTERFACE assert
		MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
	END INTERFACE
	INTERFACE assert_eq
		MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
	END INTERFACE
	INTERFACE arth
		MODULE PROCEDURE arth_r, arth_i
	END INTERFACE
	INTERFACE geop
		MODULE PROCEDURE geop_r, geop_i, geop_c, geop_dv
	END INTERFACE
	INTERFACE cumsum
		MODULE PROCEDURE cumsum_r,cumsum_i
	END INTERFACE
	INTERFACE poly
		MODULE PROCEDURE poly_rr,poly_rrv,&
			poly_rc,poly_cc,poly_msk_rrv
	END INTERFACE
	INTERFACE poly_term
		MODULE PROCEDURE poly_term_rr,poly_term_cc
	END INTERFACE
	INTERFACE outerprod
		MODULE PROCEDURE outerprod_r
	END INTERFACE
	INTERFACE outerdiff
		MODULE PROCEDURE outerdiff_r,outerdiff_i
	END INTERFACE
	INTERFACE scatter_add
		MODULE PROCEDURE scatter_add_r
	END INTERFACE
	INTERFACE scatter_max
		MODULE PROCEDURE scatter_max_r
	END INTERFACE
	INTERFACE diagadd
		MODULE PROCEDURE diagadd_rv,diagadd_r
	END INTERFACE
	INTERFACE diagmult
		MODULE PROCEDURE diagmult_rv,diagmult_r
	END INTERFACE
	INTERFACE get_diag
		MODULE PROCEDURE get_diag_rv
	END INTERFACE
	INTERFACE put_diag
		MODULE PROCEDURE put_diag_rv, put_diag_r
	END INTERFACE
CONTAINS
!BL
	SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
	REAL(SP), DIMENSION(:), INTENT(IN) :: src
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
	INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
	n_copied=min(size(src),size(dest))
	n_not_copied=size(src)-n_copied
	dest(1:n_copied)=src(1:n_copied)
	END SUBROUTINE array_copy_r
!BL
	SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
	INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
	n_copied=min(size(src),size(dest))
	n_not_copied=size(src)-n_copied
	dest(1:n_copied)=src(1:n_copied)
	END SUBROUTINE array_copy_i
!BL
!BL
	SUBROUTINE swap_i(a,b)
	INTEGER(I4B), INTENT(INOUT) :: a,b
	INTEGER(I4B) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_i
!BL
	SUBROUTINE swap_r(a,b)
	REAL(SP), INTENT(INOUT) :: a,b
	REAL(SP) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_r
!BL
	SUBROUTINE swap_rv(a,b)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(SP), DIMENSION(SIZE(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_rv
!BL
	SUBROUTINE swap_c(a,b)
	COMPLEX(SPC), INTENT(INOUT) :: a,b
	COMPLEX(SPC) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_c
!BL
	SUBROUTINE swap_cv(a,b)
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
	COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cv
!BL
	SUBROUTINE swap_cm(a,b)
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
	COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cm
!BL
	SUBROUTINE masked_swap_rs(a,b,mask)
	REAL(SP), INTENT(INOUT) :: a,b
	LOGICAL(LGT), INTENT(IN) :: mask
	REAL(SP) :: swp
	if (mask) then
		swp=a
		a=b
		b=swp
	end if
	END SUBROUTINE masked_swap_rs
!BL
	SUBROUTINE masked_swap_rv(a,b,mask)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(a)) :: swp
	where (mask)
		swp=a
		a=b
		b=swp
	end where
	END SUBROUTINE masked_swap_rv
!BL
	SUBROUTINE masked_swap_rm(a,b,mask)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
	LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
	where (mask)
		swp=a
		a=b
		b=swp
	end where
	END SUBROUTINE masked_swap_rm
!BL
!BL
	FUNCTION reallocate_rv(p,n)
	REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: nold,ierr
	allocate(reallocate_rv(n),stat=ierr)
	if (ierr /= 0) call &
		nrerror('reallocate_rv: problem in attempt to allocate memory')
	if (.not. associated(p)) RETURN
	nold=size(p)
	reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
	deallocate(p)
	END FUNCTION reallocate_rv
!BL
	FUNCTION reallocate_iv(p,n)
	INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: nold,ierr
	allocate(reallocate_iv(n),stat=ierr)
	if (ierr /= 0) call &
		nrerror('reallocate_iv: problem in attempt to allocate memory')
	if (.not. associated(p)) RETURN
	nold=size(p)
	reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
	deallocate(p)
	END FUNCTION reallocate_iv
!BL
	FUNCTION reallocate_hv(p,n)
	CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: nold,ierr
	allocate(reallocate_hv(n),stat=ierr)
	if (ierr /= 0) call &
		nrerror('reallocate_hv: problem in attempt to allocate memory')
	if (.not. associated(p)) RETURN
	nold=size(p)
	reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
	deallocate(p)
	END FUNCTION reallocate_hv
!BL
	FUNCTION reallocate_rm(p,n,m)
	REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B) :: nold,mold,ierr
	allocate(reallocate_rm(n,m),stat=ierr)
	if (ierr /= 0) call &
		nrerror('reallocate_rm: problem in attempt to allocate memory')
	if (.not. associated(p)) RETURN
	nold=size(p,1)
	mold=size(p,2)
	reallocate_rm(1:min(nold,n),1:min(mold,m))=&
		p(1:min(nold,n),1:min(mold,m))
	deallocate(p)
	END FUNCTION reallocate_rm
!BL
	FUNCTION reallocate_im(p,n,m)
	INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B) :: nold,mold,ierr
	allocate(reallocate_im(n,m),stat=ierr)
	if (ierr /= 0) call &
		nrerror('reallocate_im: problem in attempt to allocate memory')
	if (.not. associated(p)) RETURN
	nold=size(p,1)
	mold=size(p,2)
	reallocate_im(1:min(nold,n),1:min(mold,m))=&
		p(1:min(nold,n),1:min(mold,m))
	deallocate(p)
	END FUNCTION reallocate_im
!BL
	FUNCTION ifirstloc(mask)
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	INTEGER(I4B) :: ifirstloc
	INTEGER(I4B), DIMENSION(1) :: loc
	loc=maxloc(merge(1,0,mask))
	ifirstloc=loc(1)
	if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
	END FUNCTION ifirstloc
!BL
	FUNCTION imaxloc_r(arr)
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B) :: imaxloc_r
	INTEGER(I4B), DIMENSION(1) :: imax
	imax=maxloc(arr(:))
	imaxloc_r=imax(1)
	END FUNCTION imaxloc_r
!BL
	FUNCTION imaxloc_i(iarr)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
	INTEGER(I4B), DIMENSION(1) :: imax
	INTEGER(I4B) :: imaxloc_i
	imax=maxloc(iarr(:))
	imaxloc_i=imax(1)
	END FUNCTION imaxloc_i
!BL
	FUNCTION iminloc(arr)
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(1) :: imin
	INTEGER(I4B) :: iminloc
	imin=minloc(arr(:))
	iminloc=imin(1)
	END FUNCTION iminloc
!BL
	SUBROUTINE assert1(n1,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1
	if (.not. n1) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert1'
	end if
	END SUBROUTINE assert1
!BL
	SUBROUTINE assert2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2
	if (.not. (n1 .and. n2)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert2'
	end if
	END SUBROUTINE assert2
!BL
	SUBROUTINE assert3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2,n3
	if (.not. (n1 .and. n2 .and. n3)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert3'
	end if
	END SUBROUTINE assert3
!BL
	SUBROUTINE assert4(n1,n2,n3,n4,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2,n3,n4
	if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert4'
	end if
	END SUBROUTINE assert4
!BL
	SUBROUTINE assert_v(n,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, DIMENSION(:), INTENT(IN) :: n
	if (.not. all(n)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert_v'
	end if
	END SUBROUTINE assert_v
!BL
	FUNCTION assert_eq2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2
	INTEGER :: assert_eq2
	if (n1 == n2) then
		assert_eq2=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq2'
	end if
	END FUNCTION assert_eq2
!BL
	FUNCTION assert_eq3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2,n3
	INTEGER :: assert_eq3
	if (n1 == n2 .and. n2 == n3) then
		assert_eq3=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq3'
	end if
	END FUNCTION assert_eq3
!BL
	FUNCTION assert_eq4(n1,n2,n3,n4,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2,n3,n4
	INTEGER :: assert_eq4
	if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
		assert_eq4=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq4'
	end if
	END FUNCTION assert_eq4
!BL
	FUNCTION assert_eqn(nn,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, DIMENSION(:), INTENT(IN) :: nn
	INTEGER :: assert_eqn
	if (all(nn(2:) == nn(1))) then
		assert_eqn=nn(1)
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eqn'
	end if
	END FUNCTION assert_eqn
!BL
	SUBROUTINE nrerror(string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	write (*,*) 'nrerror: ',string
!	STOP 'program terminated by nrerror'
	END SUBROUTINE nrerror
!BL
	FUNCTION arth_r(first,increment,n)
	REAL(SP), INTENT(IN) :: first,increment
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: arth_r
	INTEGER(I4B) :: k,k2
	REAL(SP) :: temp
	if (n > 0) arth_r(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_r(k)=arth_r(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_r(k)=arth_r(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth_r
!BL
	FUNCTION arth_i(first,increment,n)
	INTEGER(I4B), INTENT(IN) :: first,increment,n
	INTEGER(I4B), DIMENSION(n) :: arth_i
	INTEGER(I4B) :: k,k2,temp
	if (n > 0) arth_i(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_i(k)=arth_i(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_i(k)=arth_i(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth_i
!BL
!BL
	FUNCTION geop_r(first,factor,n)
	REAL(SP), INTENT(IN) :: first,factor
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: geop_r
	INTEGER(I4B) :: k,k2
	REAL(SP) :: temp
	if (n > 0) geop_r(1)=first
	if (n <= NPAR_GEOP) then
		do k=2,n
			geop_r(k)=geop_r(k-1)*factor
		end do
	else
		do k=2,NPAR2_GEOP
			geop_r(k)=geop_r(k-1)*factor
		end do
		temp=factor**NPAR2_GEOP
		k=NPAR2_GEOP
		do
			if (k >= n) exit
			k2=k+k
			geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
			temp=temp*temp
			k=k2
		end do
	end if
	END FUNCTION geop_r
!BL
	FUNCTION geop_i(first,factor,n)
	INTEGER(I4B), INTENT(IN) :: first,factor,n
	INTEGER(I4B), DIMENSION(n) :: geop_i
	INTEGER(I4B) :: k,k2,temp
	if (n > 0) geop_i(1)=first
	if (n <= NPAR_GEOP) then
		do k=2,n
			geop_i(k)=geop_i(k-1)*factor
		end do
	else
		do k=2,NPAR2_GEOP
			geop_i(k)=geop_i(k-1)*factor
		end do
		temp=factor**NPAR2_GEOP
		k=NPAR2_GEOP
		do
			if (k >= n) exit
			k2=k+k
			geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
			temp=temp*temp
			k=k2
		end do
	end if
	END FUNCTION geop_i
!BL
	FUNCTION geop_c(first,factor,n)
	COMPLEX(SP), INTENT(IN) :: first,factor
	INTEGER(I4B), INTENT(IN) :: n
	COMPLEX(SP), DIMENSION(n) :: geop_c
	INTEGER(I4B) :: k,k2
	COMPLEX(SP) :: temp
	if (n > 0) geop_c(1)=first
	if (n <= NPAR_GEOP) then
		do k=2,n
			geop_c(k)=geop_c(k-1)*factor
		end do
	else
		do k=2,NPAR2_GEOP
			geop_c(k)=geop_c(k-1)*factor
		end do
		temp=factor**NPAR2_GEOP
		k=NPAR2_GEOP
		do
			if (k >= n) exit
			k2=k+k
			geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
			temp=temp*temp
			k=k2
		end do
	end if
	END FUNCTION geop_c
!BL
	FUNCTION geop_dv(first,factor,n)
	REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), DIMENSION(size(first),n) :: geop_dv
	INTEGER(I4B) :: k,k2
	REAL(DP), DIMENSION(size(first)) :: temp
	if (n > 0) geop_dv(:,1)=first(:)
	if (n <= NPAR_GEOP) then
		do k=2,n
			geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
		end do
	else
		do k=2,NPAR2_GEOP
			geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
		end do
		temp=factor**NPAR2_GEOP
		k=NPAR2_GEOP
		do
			if (k >= n) exit
			k2=k+k
			geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
				spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
			temp=temp*temp
			k=k2
		end do
	end if
	END FUNCTION geop_dv
!BL
!BL
	RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP), OPTIONAL, INTENT(IN) :: seed
	REAL(SP), DIMENSION(size(arr)) :: ans
	INTEGER(I4B) :: n,j
	REAL(SP) :: sd
	n=size(arr)
	if (n == 0_i4b) RETURN
	sd=0.0_sp
	if (present(seed)) sd=seed
	ans(1)=arr(1)+sd
	if (n < NPAR_CUMSUM) then
		do j=2,n
			ans(j)=ans(j-1)+arr(j)
		end do
	else
		ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
		ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
	end if
	END FUNCTION cumsum_r
!BL
	RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
	INTEGER(I4B), DIMENSION(size(arr)) :: ans
	INTEGER(I4B) :: n,j,sd
	n=size(arr)
	if (n == 0_i4b) RETURN
	sd=0_i4b
	if (present(seed)) sd=seed
	ans(1)=arr(1)+sd
	if (n < NPAR_CUMSUM) then
		do j=2,n
			ans(j)=ans(j-1)+arr(j)
		end do
	else
		ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
		ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
	end if
	END FUNCTION cumsum_i
!BL
!BL
	RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP), OPTIONAL, INTENT(IN) :: seed
	REAL(SP), DIMENSION(size(arr)) :: ans
	INTEGER(I4B) :: n,j
	REAL(SP) :: sd
	n=size(arr)
	if (n == 0_i4b) RETURN
	sd=1.0_sp
	if (present(seed)) sd=seed
	ans(1)=arr(1)*sd
	if (n < NPAR_CUMPROD) then
		do j=2,n
			ans(j)=ans(j-1)*arr(j)
		end do
	else
		ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
		ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
	end if
	END FUNCTION cumprod
!BL
!BL
	FUNCTION poly_rr(x,coeffs)
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
	REAL(SP) :: poly_rr
	REAL(SP) :: pow
	REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
	INTEGER(I4B) :: i,n,nn
	n=size(coeffs)
	if (n <= 0) then
		poly_rr=0.0_sp
	else if (n < NPAR_POLY) then
		poly_rr=coeffs(n)
		do i=n-1,1,-1
			poly_rr=x*poly_rr+coeffs(i)
		end do
	else
		allocate(vec(n+1))
		pow=x
		vec(1:n)=coeffs
		do
			vec(n+1)=0.0_sp
			nn=ishft(n+1,-1)
			vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
			if (nn == 1) exit
			pow=pow*pow
			n=nn
		end do
		poly_rr=vec(1)
		deallocate(vec)
	end if
	END FUNCTION poly_rr
!BL
	FUNCTION poly_rc(x,coeffs)
	COMPLEX(SPC), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
	COMPLEX(SPC) :: poly_rc
	COMPLEX(SPC) :: pow
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
	INTEGER(I4B) :: i,n,nn
	n=size(coeffs)
	if (n <= 0) then
		poly_rc=0.0_sp
	else if (n < NPAR_POLY) then
		poly_rc=coeffs(n)
		do i=n-1,1,-1
			poly_rc=x*poly_rc+coeffs(i)
		end do
	else
		allocate(vec(n+1))
		pow=x
		vec(1:n)=coeffs
		do
			vec(n+1)=0.0_sp
			nn=ishft(n+1,-1)
			vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
			if (nn == 1) exit
			pow=pow*pow
			n=nn
		end do
		poly_rc=vec(1)
		deallocate(vec)
	end if
	END FUNCTION poly_rc
!BL
	FUNCTION poly_cc(x,coeffs)
	COMPLEX(SPC), INTENT(IN) :: x
	COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
	COMPLEX(SPC) :: poly_cc
	COMPLEX(SPC) :: pow
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
	INTEGER(I4B) :: i,n,nn
	n=size(coeffs)
	if (n <= 0) then
		poly_cc=0.0_sp
	else if (n < NPAR_POLY) then
		poly_cc=coeffs(n)
		do i=n-1,1,-1
			poly_cc=x*poly_cc+coeffs(i)
		end do
	else
		allocate(vec(n+1))
		pow=x
		vec(1:n)=coeffs
		do
			vec(n+1)=0.0_sp
			nn=ishft(n+1,-1)
			vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
			if (nn == 1) exit
			pow=pow*pow
			n=nn
		end do
		poly_cc=vec(1)
		deallocate(vec)
	end if
	END FUNCTION poly_cc
!BL
	FUNCTION poly_rrv(x,coeffs)
	REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
	REAL(SP), DIMENSION(size(x)) :: poly_rrv
	INTEGER(I4B) :: i,n,m
	m=size(coeffs)
	n=size(x)
	if (m <= 0) then
		poly_rrv=0.0_sp
	else if (m < n .or. m < NPAR_POLY) then
		poly_rrv=coeffs(m)
		do i=m-1,1,-1
			poly_rrv=x*poly_rrv+coeffs(i)
		end do
	else
		do i=1,n
			poly_rrv(i)=poly_rr(x(i),coeffs)
		end do
	end if
	END FUNCTION poly_rrv
!BL
	FUNCTION poly_msk_rrv(x,coeffs,mask)
	REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
	poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
	END FUNCTION poly_msk_rrv
!BL
!BL
!BL
	RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), INTENT(IN) :: b
	REAL(SP), DIMENSION(size(a)) :: u
	INTEGER(I4B) :: n,j
	n=size(a)
	if (n <= 0) RETURN
	u(1)=a(1)
	if (n < NPAR_POLYTERM) then
		do j=2,n
			u(j)=a(j)+b*u(j-1)
		end do
	else
		u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
		u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
	end if
	END FUNCTION poly_term_rr
!BL
	RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
	COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
	COMPLEX(SPC), INTENT(IN) :: b
	COMPLEX(SPC), DIMENSION(size(a)) :: u
	INTEGER(I4B) :: n,j
	n=size(a)
	if (n <= 0) RETURN
	u(1)=a(1)
	if (n < NPAR_POLYTERM) then
		do j=2,n
			u(j)=a(j)+b*u(j-1)
		end do
	else
		u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
		u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
	end if
	END FUNCTION poly_term_cc
!BL
!BL
	FUNCTION zroots_unity(n,nn)
	INTEGER(I4B), INTENT(IN) :: n,nn
	COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
	INTEGER(I4B) :: k
	REAL(SP) :: theta
	zroots_unity(1)=1.0
	theta=TWOPI/n
	k=1
	do
		if (k >= nn) exit
		zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
		zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
			zroots_unity(2:min(k,nn-k))
		k=2*k
	end do
	END FUNCTION zroots_unity
!BL
	FUNCTION outerprod_r(a,b)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
	outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerprod_r
!BL
!BL
	FUNCTION outerdiv(a,b)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
	outerdiv = spread(a,dim=2,ncopies=size(b)) / &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiv
!BL
	FUNCTION outersum(a,b)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a),size(b)) :: outersum
	outersum = spread(a,dim=2,ncopies=size(b)) + &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outersum
!BL
	FUNCTION outerdiff_r(a,b)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
	outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_r
!BL
	FUNCTION outerdiff_d(a,b)
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
	outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_d
!BL
	FUNCTION outerdiff_i(a,b)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
	INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
	outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_i
!BL
	FUNCTION outerand(a,b)
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
	LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
	outerand = spread(a,dim=2,ncopies=size(b)) .and. &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerand
!BL
	SUBROUTINE scatter_add_r(dest,source,dest_index)
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
	REAL(SP), DIMENSION(:), INTENT(IN) :: source
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
	INTEGER(I4B) :: m,n,j,i
	n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
	m=size(dest)
	do j=1,n
		i=dest_index(j)
		if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
	end do
	END SUBROUTINE scatter_add_r
	SUBROUTINE scatter_max_r(dest,source,dest_index)
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
	REAL(SP), DIMENSION(:), INTENT(IN) :: source
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
	INTEGER(I4B) :: m,n,j,i
	n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
	m=size(dest)
	do j=1,n
		i=dest_index(j)
		if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
	end do
	END SUBROUTINE scatter_max_r
!BL
	SUBROUTINE diagadd_rv(mat,diag)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	REAL(SP), DIMENSION(:), INTENT(IN) :: diag
	INTEGER(I4B) :: j,n
	n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
	do j=1,n
		mat(j,j)=mat(j,j)+diag(j)
	end do
	END SUBROUTINE diagadd_rv
!BL
	SUBROUTINE diagadd_r(mat,diag)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	REAL(SP), INTENT(IN) :: diag
	INTEGER(I4B) :: j,n
	n = min(size(mat,1),size(mat,2))
	do j=1,n
		mat(j,j)=mat(j,j)+diag
	end do
	END SUBROUTINE diagadd_r
!BL
	SUBROUTINE diagmult_rv(mat,diag)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	REAL(SP), DIMENSION(:), INTENT(IN) :: diag
	INTEGER(I4B) :: j,n
	n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
	do j=1,n
		mat(j,j)=mat(j,j)*diag(j)
	end do
	END SUBROUTINE diagmult_rv
!BL
	SUBROUTINE diagmult_r(mat,diag)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	REAL(SP), INTENT(IN) :: diag
	INTEGER(I4B) :: j,n
	n = min(size(mat,1),size(mat,2))
	do j=1,n
		mat(j,j)=mat(j,j)*diag
	end do
	END SUBROUTINE diagmult_r
!BL
	FUNCTION get_diag_rv(mat)
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
	REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
	INTEGER(I4B) :: j
	j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
	do j=1,size(mat,1)
		get_diag_rv(j)=mat(j,j)
	end do
	END FUNCTION get_diag_rv
!BL
!BL
	SUBROUTINE put_diag_rv(diagv,mat)
	REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	INTEGER(I4B) :: j,n
	n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
	do j=1,n
		mat(j,j)=diagv(j)
	end do
	END SUBROUTINE put_diag_rv
!BL
	SUBROUTINE put_diag_r(scal,mat)
	REAL(SP), INTENT(IN) :: scal
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
	INTEGER(I4B) :: j,n
	n = min(size(mat,1),size(mat,2))
	do j=1,n
		mat(j,j)=scal
	end do
	END SUBROUTINE put_diag_r
!BL
	SUBROUTINE unit_matrix(mat)
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
	INTEGER(I4B) :: i,n
	n=min(size(mat,1),size(mat,2))
	mat(:,:)=0.0_sp
	do i=1,n
		mat(i,i)=1.0_sp
	end do
	END SUBROUTINE unit_matrix
!BL
	FUNCTION upper_triangle(j,k,extra)
	INTEGER(I4B), INTENT(IN) :: j,k
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
	LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
	INTEGER(I4B) :: n
	n=0
	if (present(extra)) n=extra
	upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
	END FUNCTION upper_triangle
!BL
	FUNCTION lower_triangle(j,k,extra)
	INTEGER(I4B), INTENT(IN) :: j,k
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
	LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
	INTEGER(I4B) :: n
	n=0
	if (present(extra)) n=extra
	lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
	END FUNCTION lower_triangle
!BL
	FUNCTION vabs(v)
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP) :: vabs
	vabs=sqrt(dot_product(v,v))
	END FUNCTION vabs
!BL
END MODULE nrutil





!***********************************************************************************
!!                                                                                !!
!! Fitted range: rmin(1) < R < rmax(1), as specified in the input file.           !!
!! Jac3(1) is R, the distance between centers of mass (in Angstroms).             !!
!! Jac3(2) and Jac3(3) are cos(theta1) and cos(theta2) and range from (-1,1).     !!
!! Jac3(4) is the dihedral angle, in radians, with range: (0,2pi).                !!
!! NAME1 is the name of the PES-file generated by AUTOSURF-PES.                   !!
!! Subroutine PES(jac3,V,NAME1) returns the potential "V".                        !!
!! Output energy is in kcal/mol.                                                  !!
!!                                                                                !!
!***********************************************************************************
SUBROUTINE PES(jac3,V,NAME1)  

use dynamic_parameters
!-----------------------------------------------------------------------------------
implicit none
 integer,parameter :: XDIM=4
 character (len=40) :: NAME1 
 character (len=300) :: line
 integer :: i,initflag,nline,ncont1
 real*8 :: xi(XDIM),jac3(XDIM)
 real*8 :: temp3,v,pii
 real*8,allocatable :: cart3(:)
 logical :: logica1
 save initflag
 data initflag /1/
!-----------------------------------------------------------------------------------
! Interface blocks
!-----------------------------------------------------------------------------------
INTERFACE! Energy of largest basis and high-level ab initio
  FUNCTION func_actual(xi)    
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: xi
      REAL(SP) :: func_actual
  END FUNCTION func_actual
end interface
INTERFACE! Energy of minimal basis and low-level ab initio
  FUNCTION func_actual_seed(xi)  
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: xi
      REAL(SP) :: func_actual_seed
  END FUNCTION func_actual_seed
END INTERFACE
!-----------------------------------------------------------------------------------

pii=dacos(-1d0)
nline=scan(NAME1,' ')-1

! check if PES file exists...
inquire(file=NAME1(1:nline),exist=logica1)
if(.not.logica1)then
 write(*,*)
 write(*,*)'ERROR: The file: ',NAME1(1:nline),' does not exist !! '
 stop 
endif
! check if input file exist...
inquire(file='input-AUTOSURF-PES.dat',exist=logica1)
if(.not.logica1)then
 write(*,*)
 write(*,*)'ERROR: The file: input-AUTOSURF-PES.dat does not exist !! '
 stop 
endif

!***********************************************************************************
!                            INITIALIZATION
!***********************************************************************************
IF(initflag==1)THEN! initialize   

  ddd=4
  alpha=-1.0d0! coefficient in R=exp(alpha*r**xbeta) coordinate
  xbeta=1.0
  epss=1d-14
  zz=4
  zz_low=6
  zz4=20
  dist_tol=1.25d0
  ! basis for "minimal fit" to high level data, also used to fit low level grid.
  order_1_min=3  
  order_2_min=3
  order_3_min=3
  order_4_min=3

  allocate(jac(XDIM),jac2(XDIM))
  ! Read from the #FRAGMENTS INFORMATION section of the input file
  open(unit=10,file='input-AUTOSURF-PES.dat',status='old',action='read')
  112 read(10,'(A100)')line
  ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
  if (ncont1==0) goto 112
  read(10,*)
  read(10,*)exch
  read(10,*)flip1
  read(10,*)flip2
  symparts=((exch+1)*(flip1+1)*(flip2+1))
  read(10,*) natom1! number of atoms in frag1
  read(10,*) natom2! number of atoms in frag2
  natom=natom1+natom2
  nbdist=natom*(natom-1)/2
  W_a=0.432d0*3/max(natom1,natom2)! scaling for R in dist metric, suggested: 1/fraglength
  allocate(ref1(3*natom1),ref2(3*natom2),bdist(nbdist),cart(3*(natom)))
  allocate(symb(natom),mass(natom))
  do i=1,natom
    read(10,*) symb(i)! element labels of all atoms
  enddo
  do i=1,natom
    read(10,*) mass(i)! masses of all atoms
  enddo  
  do i=1,3*natom1
    read(10,*) ref1(i)! Cartesian positions for fragment 1 atoms
  enddo
  do i=1,3*natom2
    read(10,*) ref2(i)! Cartesian positions for fragment 2 atoms
  enddo
  close(10)
  !# PES INFORMATION:
  OPEN (UNIT=652, FILE=NAME1(1:nline), FORM='UNFORMATTED', ACCESS='SEQUENTIAL') 
  read(652) count3
  read(652) order_1
  read(652) order_2
  read(652) order_3
  read(652) order_4
  read(652) maxpoints
  read(652) mass
  read(652) rmax
  read(652) rmin
  read(652) Max_E
  read(652) low_grid
  read(652) count_seed
  call basis_size(XDIM,order_1,order_2,order_3,order_4,basis_1)
  call basis_size(XDIM,order_1-1,order_2-1,order_3-1,order_4-1,basis_2)  
  call basis_size(XDIM,order_1_min,order_2_min,order_3_min,order_4_min,basis_3) 
  allocate(b2(basis_1,count3),b2_lower(basis_2,count3),b2_minimal(basis_3,count3),d(count3),coords(count3,XDIM)) 
  allocate(b2_seed(basis_3,count_seed),d_seed(count_seed),coords_seed(count_seed,XDIM))
  b2=0d0
  b2_lower=0d0
  b2_minimal=0d0
  b2_seed=0d0
  d=0d0
  d_seed=0d0
  coords=0d0
  coords_seed=0d0
  do i=1,count3 
     read(652) b2(:,i)       
  enddo
  do i=1,count3      
     read(652) b2_lower(:,i)       
  enddo
  do i=1,count3      
     read(652) b2_minimal(:,i)       
  enddo  
  do i=1,count3       
     read(652) d(i) 
  enddo
  do i=1,count3        
     read(652) coords(i,:)
  enddo
  if(low_grid>0)then
     read(652) Max_E_seed
     do i=1,count_seed      
        read(652) b2_seed(:,i)       
     enddo
     do i=1,count_seed 
        read(652) d_seed(i)       
     enddo
     do i=1,count_seed        
        read(652) coords_seed(i,:)       
     enddo
  endif
  close(652)
  initflag=2

  xi(1)=rmax(1)
  xi(2)=1.d0  
  xi(3)=1.d0  
  xi(4)=0.d0  
  ass=func_actual(xi)

ENDIF
!***********************************************************************************

allocate(cart3(3*(natom)))
xi=jac3

if(xi(1)<rmin(1))then
  v=Max_E
  !write(*,*) 'coords outside fitted range'
  !write(*,*)  xi(1),rmax(1),rmin(1)
  return
endif
if(xi(4)>pii)then
  xi(4)=xi(4)-2d0*pii
endif
if(xi(4)<-pii)then
  xi(4)=xi(4)+2d0*pii
endif
if(xi(2)>1d0)then
  xi(2)=2d0-xi(2)
endif
if(xi(3)>1d0)then
  xi(3)=2d0-xi(3)
endif
if(xi(2)<-1d0)then
  xi(2)=-2d0-xi(2)
endif
if(xi(3)<-1d0)then
  xi(3)=-2d0-xi(3)
endif

call INT_Cart(cart3,xi,mass,natom1,natom2,ref1,ref2)
!write(*,*)
!write(*,*)xi
!write(*,*)
!write(*,*)cart3
!write(*,*)
call cart_to_bdist_inter(cart3,natom1,natom2,dist_tol,dist_flag)
if(dist_flag==1) then
  v=Max_E-ass
  !write(*,*) 'bdist less than distol'
  return
endif

if(low_grid>0)then
  temp3=func_actual_seed(xi)
  if(temp3>Max_E_seed)then
    v=Max_E-ass
    !write(*,*)xi(1),'hit ceiling on low grid'
    return
  endif
endif

temp3=func_actual(xi)
if(temp3>Max_E)then
  temp3=Max_E
  !write(*,*)xi(1),'hit ceiling (func_actual)'
endif

V=temp3-ass
return

END SUBROUTINE PES


!***********************************************************************************

!                    S U B R O U T I N E S   &   F U N C T I O N S                           

!***********************************************************************************

subroutine cart_to_bdist_inter(x,natom1,natom2,dist_tol,flag)
  implicit none
  integer :: i,j,k,flag,natom1,natom2
  real*8 :: x(3*(natom1+natom2)),summ,dist_tol
  flag=0
  do i=1,natom1
     do j=natom1+1,natom1+natom2
        summ=0d0
        do k=1,3
           summ=summ+(x(3*(i-1)+k)-x(3*(j-1)+k))**2
        enddo
        if(sqrt(summ)<dist_tol)then
           flag=1
        endif
     enddo
  enddo
  return
end subroutine cart_to_bdist_inter

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      B A S I S _ S I Z E
! ----------------------------------------------------------------------------------
! Calculate the size of the basis set

! *** Input ***
! d         <-- 
! order_1   <-- maximum power of R = exp(alpha*r)
! order_2   <-- maximum value of L1
! order_3   <-- maximum value of L2
! order_4   <-- maximum value of L = L1 + L2

subroutine basis_size(d,order_1,order_2,order_3,order_4,basis)
  
 implicit none
 integer :: d,order_1,order_2,order_3,order_4,count,basis,l1,l2,m
   
!!!!!!!!basis calc
 count=0
 !  count=d*order_1
  
 do l1=0,order_2
   do l2=0,order_3
     if((l1+l2)<order_4+1)then
       do m=0,min(l1,l2)
        count=count+1
       enddo
     endif
   enddo
 enddo
 basis=count*(order_1)+1
! basis=count+1

return
end subroutine basis_size

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      I N T _ C a r t
! ----------------------------------------------------------------------------------
! Known the internal coordinates for a given configuration:
! internal2(1) -> R
! internal2(2) -> cos(theta1)
! internal2(3) -> cos(theta2)
! internal2(4) -> phi
! the Cartesian coordinates for all atoms in the system (cart) are calculated

! *** Input ***
! internal2 <-- vector containing internal coordinates
! mass      <-- masses of all atoms
! natom1    <-- number of atoms in fragment 1
! natom2    <-- number of atoms in fragment 2
! ref1      <-- Cartesian coord. of atoms in frag. 1, placed along z-axis, 1st atom at origin
! ref2      <-- Cartesian coord. of atoms in frag. 2, placed along z-axis, 1st atom at origin

subroutine INT_Cart(cart,internal2,mass,natom1,natom2,ref1,ref2)

 implicit none
 integer :: i,k,kp,natom1,natom2
 real*8 :: internal(6),internal2(4),cart((natom1+natom2)*3),mass(natom1+natom2),   &
            ref1(natom1*3),ref1_temp(natom1*3),ref2(natom2*3),ref2_temp(natom2*3)
 real*8 :: cart_ref1(3,natom1),cm(3),cart_ref2(3,natom2),cart_ref1t(3,natom1),     &
            cart_ref2t(3,natom2),U_rot(3,3)
 real*8 :: pii,gamma1,gamma2,beta1,beta2,alpha1,alpha2
 
 pii=acos(-1d0)
 internal(1)=internal2(1)
 internal(2)=0d0
 internal(3)=internal2(2)
 internal(4)=0d0
 internal(5)=internal2(3)
 internal(6)=internal2(4)
 ref1_temp=ref1
 ref2_temp=ref2

 ! set c.m. of fragment 1 at origin
 call rm_cmass(ref1_temp,mass(1:natom1),natom1,natom1)
 ! set c.m. of fragment 2 at origin
 call rm_cmass(ref2_temp,mass(natom1+1:natom1+natom2),natom2,natom2)

 ! Cartesian coordinates of c.m. for fragment 2
 cm(1)=0d0
 cm(2)=0d0
 cm(3)=internal(1)

 alpha1=0d0
 gamma1=internal(2)
 beta1=acos(internal(3))

 ! U = Z1(alpha1) Y2(beta1) Z3(gamma1)
 !ZYZ for proper Euler angles
 U_rot(1,1)=cos(alpha1)*cos(beta1)*cos(gamma1)-sin(alpha1)*sin(gamma1)
 U_rot(1,2)=-cos(alpha1)*cos(beta1)*sin(gamma1)-sin(alpha1)*cos(gamma1)
 U_rot(1,3)=cos(alpha1)*sin(beta1)

 U_rot(2,1)=sin(alpha1)*cos(beta1)*cos(gamma1)+cos(alpha1)*sin(gamma1)
 U_rot(2,2)=-sin(alpha1)*cos(beta1)*sin(gamma1)+cos(alpha1)*cos(gamma1)
 U_rot(2,3)=sin(alpha1)*sin(beta1)

 U_rot(3,1)=-sin(beta1)*cos(gamma1)
 U_rot(3,2)=sin(beta1)*sin(gamma1)
 U_rot(3,3)=cos(beta1)

 call vec_to_mat2(ref1_temp,cart_ref1,natom1)
 call rotmol(natom1,cart_ref1,cart_ref1t,U_rot)
 call mat_to_vec2(cart_ref1t,ref1_temp,natom1)

 gamma2=internal(4)
 beta2=acos(internal(5))
 alpha2=-internal(6)

 U_rot(1,1)=cos(alpha2)*cos(beta2)*cos(gamma2)-sin(alpha2)*sin(gamma2)
 U_rot(1,2)=-cos(alpha2)*cos(beta2)*sin(gamma2)-sin(alpha2)*cos(gamma2)
 U_rot(1,3)=cos(alpha2)*sin(beta2)

 U_rot(2,1)=sin(alpha2)*cos(beta2)*cos(gamma2)+cos(alpha2)*sin(gamma2)
 U_rot(2,2)=-sin(alpha2)*cos(beta2)*sin(gamma2)+cos(alpha2)*cos(gamma2)
 U_rot(2,3)=sin(alpha2)*sin(beta2)

 U_rot(3,1)=-sin(beta2)*cos(gamma2)
 U_rot(3,2)=sin(beta2)*sin(gamma2)
 U_rot(3,3)=cos(beta2)

 call vec_to_mat2(ref2_temp,cart_ref2,natom2)
 call rotmol(natom2,cart_ref2,cart_ref2t,U_rot)
 call mat_to_vec2(cart_ref2t,ref2_temp,natom2)

 do k=1,natom2
   do kp=1,3
    ref2_temp((k-1)*3+kp)=ref2_temp((k-1)*3+kp)+cm(kp)      
   enddo
 enddo

 do i=1,3*natom2
   cart(3*natom1+i)=ref2_temp(i)
 enddo

 do i=1,3*natom1
   cart(i)=ref1_temp(i)
 enddo

return
end subroutine INT_Cart

subroutine rm_cmass(cart,mass,natom,natom1)
integer :: k,kp,natom,natom1
real*8 :: mass(natom),cart(natom*3),mtot,cmass1(3)
mtot=0d0
do k=1,natom1
   mtot=mtot+mass(k)
enddo
cmass1=0d0
do k=1,natom1
   do kp=1,3
      cmass1(kp)=cmass1(kp)+cart((k-1)*3+kp)*mass(k)
   enddo
enddo
cmass1=cmass1/mtot

do k=1,natom
   do kp=1,3
      cart((k-1)*3+kp)=cart((k-1)*3+kp)-cmass1(kp)      
   enddo
enddo
return
end subroutine rm_cmass


subroutine vec_to_mat2(cart_perms,cart_mat,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_mat(kp,k)=cart_perms((k-1)*3+kp)
   enddo
enddo
return
end subroutine vec_to_mat2

subroutine mat_to_vec2(cart_mat,cart_perms,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_perms((k-1)*3+kp)=cart_mat(kp,k)
   enddo
enddo
return
end subroutine mat_to_vec2

subroutine dist_metric(jac,jac2,scale,dist)
integer :: i
real*8 :: jac(4),jac2(4),scale,dist,temp(4),pii
pii=acos(-1d0)
temp(1)=((jac(1)-jac2(1))*scale)**2
temp(2)=(acos(jac(2))-acos(jac2(2)))**2
temp(3)=(acos(jac(3))-acos(jac2(3)))**2
temp(4)=jac(4)-jac2(4)
if(temp(4)>pii)then
   temp(4)=temp(4)-2d0*pii
endif
if(temp(4)<-pii)then
   temp(4)=temp(4)+2d0*pii
endif
temp(4)=(temp(4)**2)*sqrt((1d0-jac(2)**2)*(1d0-jac2(2)**2)*(1d0-jac(3)**2)*(1d0-jac2(3)**2))
dist=0d0
do i=1,4
   dist=dist+temp(i)
enddo
dist=sqrt(dist)
return
end subroutine dist_metric


!!! index.f90
!c-----------------------------------------------------------------------

SUBROUTINE indexxy(n,arr,indx)
  INTEGER :: n
  integer,parameter :: nstack=50, m=7
  INTEGER ::indx(n),istack(nstack)
  REAL*8 :: arr(n)
  INTEGER :: i,indxt,ir,itemp,j,jstack,k,l
  REAL*8 :: a
 
  do j=1,n
     indx(j)=j
  enddo
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,1,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        enddo
        i=0
2       indx(i+1)=indxt
     enddo
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l+1)).gt.arr(indx(l)))then
        itemp=indx(l+1)
        indx(l+1)=indx(l)
        indx(l)=itemp
     endif
     i=l+1
     j=ir
     indxt=indx(l)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
!     if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
     if(jstack.gt.NSTACK)write(6,*)'NSTACK too small in indexx'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     endif
  endif
  goto 1
  return

end subroutine indexxy


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D B R E N T (ax,bx,cx,func,dfunc,tol,xmin)
! ----------------------------------------------------------------------------------
!  

! *** Input ***
!         <-- 

FUNCTION dbrent(ax,bx,cx,func,dfunc,tol,xmin)

USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(SP), INTENT(IN) :: ax,bx,cx,tol
REAL(SP), INTENT(OUT) :: xmin
REAL(SP) :: dbrent
INTERFACE
  FUNCTION func(x)
  USE nrtype
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: func
  END FUNCTION func
  FUNCTION dfunc(x)
  USE nrtype
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: dfunc
  END FUNCTION dfunc
END INTERFACE

INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: ZEPS=1.0e-3_sp*epsilon(ax)
INTEGER(I4B) :: iter
REAL(SP) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
	u,u1,u2,v,w,x,xm
LOGICAL :: ok1,ok2
a=min(ax,cx)
b=max(ax,cx)
v=bx
w=v
x=v
e=0.0
fx=func(x)
fv=fx
fw=fx
dx=dfunc(x)
dv=dx
dw=dx
do iter=1,ITMAX
	xm=0.5_sp*(a+b)
	tol1=tol*dabs(x)+ZEPS
	tol2=2.0_sp*tol1
	if (dabs(x-xm) <= (tol2-0.5_sp*(b-a))) exit
	if (dabs(e) > tol1) then
		d1=2.0_sp*(b-a)
		d2=d1
		if (dw /= dx) d1=(w-x)*dx/(dx-dw)
		if (dv /= dx) d2=(v-x)*dx/(dx-dv)
		u1=x+d1
		u2=x+d2
		ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
		ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
		olde=e
		e=d
		if (ok1 .or. ok2) then
			if (ok1 .and. ok2) then
				d=merge(d1,d2, dabs(d1) < dabs(d2))
			else
				d=merge(d1,d2,ok1)
			end if
			if (dabs(d) <= dabs(0.5_sp*olde)) then
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) &
					d=dsign(tol1,xm-x)
			else
				e=merge(a,b, dx >= 0.0)-x
				d=0.5_sp*e
			end if
		else
			e=merge(a,b, dx >= 0.0)-x
			d=0.5_sp*e
		end if
	else
		e=merge(a,b, dx >= 0.0)-x
		d=0.5_sp*e
	end if
	if (dabs(d) >= tol1) then
		u=x+d
		fu=func(u)
	else
		u=x+sign(tol1,d)
		fu=func(u)
		if (fu > fx) exit
	end if
	du=dfunc(u)
	if (fu <= fx) then
		if (u >= x) then
			a=x
		else
			b=x
		end if
		call mov3(v,fv,dv,w,fw,dw)
		call mov3(w,fw,dw,x,fx,dx)
		call mov3(x,fx,dx,u,fu,du)
	else
		if (u < x) then
			a=u
		else
			b=u
		end if
		if (fu <= fw .or. w == x) then
			call mov3(v,fv,dv,w,fw,dw)
			call mov3(w,fw,dw,u,fu,du)
		else if (fu <= fv .or. v == x .or. v == w) then
			call mov3(v,fv,dv,u,fu,du)
		end if
	end if
end do
if (iter > ITMAX) call nrerror('dbrent: exceeded maximum iterations')
xmin=x
dbrent=fx
CONTAINS
SUBROUTINE mov3(a,b,c,d,e,f)
REAL(SP), INTENT(IN) :: d,e,f
REAL(SP), INTENT(OUT) :: a,b,c
 a=d
 b=e
 c=f
END SUBROUTINE mov3
END FUNCTION dbrent


!                                F U N C T I O N S                             

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      f u n c _ a c t u a l (xi)
! ----------------------------------------------------------------------------------
! Energy of largest basis and high-level ab initio

! *** Input ***
! xi        <-- vector containing the internal coordinates

function func_actual(xi)

use nrtype
USE dynamic_parameters
implicit none

integer,parameter :: XDIM=4
REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual
integer :: i,ip,quitt,l1,l2,jj,R,M,count
real*8 :: temp,weight,norm,somme,jac3(XDIM),jac4(XDIM),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=dacos(-1d0)
jac3=xi
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1))
allocate(PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))

count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)
quitt=0! number of expansions included in the interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo

12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1)**xbeta)
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)
norm=0d0
temp=0d0
do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count3+1-i)) 
   !write(665,*) weight
   norm=norm+weight
   temp=temp+weight*b2(1,jj)
   count=1
   do R=1,order_1
     do L1=0,order_2
       do L2=0,order_3
         if((L1+L2)<order_4+1)then
           do M=0,min(L1,L2)
            count=count+1
            temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*dcos(dble(M)*jac4(4))
           enddo
         endif
       enddo
     enddo
   enddo
!  endif
enddo

func_actual=temp/norm

return
end function func_actual

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      f u n c _ a c t u a l _ s e e d
! ----------------------------------------------------------------------------------
!  Energy of minimal basis and low-level ab initio

! *** Input ***
! xi        <-- vector containing the internal coordinates

function func_actual_seed(xi)

use nrtype
USE dynamic_parameters
implicit none

integer,parameter :: XDIM=4
REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_seed
integer :: i,ip,quitt,l1,l2,jj,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(XDIM),jac4(XDIM),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=dacos(-1d0)

jac3=xi
allocate(ind7(count_seed),ind8(count_seed),PM1(0:order_2_min+1,0:order_2_min+1))
allocate(PM2(0:order_3_min+1,0:order_3_min+1),PD1(0:order_2_min+1,0:order_2_min+1))
allocate(PD2(0:order_3_min+1,0:order_3_min+1))

count=0
do ip=1,count_seed 
  count=count+1
  Jac4=coords_seed(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo
call indexxy(count_seed,ind7,ind8)
quitt=0! number of expansions are included in interpolation
do ip=1,count_seed
   if(ind7(ind8(count_seed))/ind7(ind8(count_seed+1-ip))>1d11) goto 12
   quitt=quitt+1
enddo

12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1)**xbeta)
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
call LPMN(order_3_min+1,order_3_min,order_3_min,jac4(3),PM2,PD2)
norm=0d0
temp=0d0
do i=1,quitt
  jj=ind8(count_seed+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count_seed+1-i)) 
   ! write(665,*) weight
   norm=norm+weight
   temp=temp+weight*b2_seed(1,jj)
   count=1
   do R=1,order_1_min
     do L1=0,order_2_min
       do L2=0,order_3_min
         if((L1+L2)<order_4_min+1)then
           do M=0,min(L1,L2)
             count=count+1
             temp=temp+weight*b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)* &
                  dcos(dble(M)*jac4(4))
           enddo
         endif
       enddo
     enddo
   enddo
!  endif
enddo

func_actual_seed=temp/norm

return
end function func_actual_seed

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      F U N C (xi)
! ----------------------------------------------------------------------------------
! Negative of squared difference surface

! *** Input ***
! xi        <-- 

function func(xi)

use nrtype
USE dynamic_parameters
implicit none

INTERFACE
  FUNCTION func_actual_min(xi)
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_min
  END FUNCTION func_actual_min
end interface
INTERFACE
  FUNCTION func_actual_seed(xi)
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_seed
  END FUNCTION func_actual_seed
end interface

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),tampon,h2wn,dist,temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=acos(-1d0)

do j=1,4
  if(xi(j)>rmax(j).or.xi(j)<rmin(j))then
    func=0d0
    return
  endif
enddo

call INT_Cart(cart,xi,mass,natom1,natom2,ref1,ref2)
call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
if(dist_flag==1) then
  func=0d0
  return
endif

if(low_grid>0) then
  temp=func_actual_seed(xi)
  if(temp>Max_E_seed)then
    func=0d0
    return
  endif
endif

temp1=func_actual_min(xi)
if(subzero==0)then
  if(temp1>Max_E)then
    func=0d0
    return
  endif
endif
if(subzero==1)then
  if(temp1+temp>Max_E)then
    func=0d0
    return
  endif
endif

allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1))
allocate(PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))

jac3=xi
count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)

quitt=0! number of expansions included in interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo
!write(701,*) quitt
12 norm=0d0

temp=0d0
jac3=xi
jac4=jac3
jac4(1)=exp(alpha*jac4(1)**xbeta)
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)
do i=1,quitt
  jj=ind8(count3+1-i)
  weight=ind7(ind8(count3+1-i)) 
  norm=norm+weight
  temp=temp+weight*b2(1,jj)
  count=1
  do R=1,order_1
    do L1=0,order_2
      do L2=0,order_3
        if((L1+L2)<order_4+1)then
          do M=0,min(L1,L2)
           count=count+1
           temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
          enddo
        endif
      enddo
    enddo
  enddo
enddo
func=temp/norm
   
! lower basis
order_1=order_1-1
order_2=order_2-1
order_3=order_3-1
order_4=order_4-1
norm=0d0
temp=0d0
do i=1,quitt  
  jj=ind8(count3+1-i)   
  weight=ind7(ind8(count3+1-i)) 
  norm=norm+weight   
  temp=temp+weight*b2_lower(1,jj)
  count=1
  do R=1,order_1
    do L1=0,order_2
      do L2=0,order_3
        if((L1+L2)<order_4+1)then
          do M=0,min(L1,L2)
           count=count+1
           temp=temp+weight*b2_lower(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
          enddo
        endif
      enddo
    enddo
  enddo
enddo
deallocate(ind7,ind8)
func=-abs(func-temp/norm)**2
!write(550+myid,*) xi
!write(550+myid,*) func

! restore order of basis   
order_1=order_1+1
order_2=order_2+1
order_3=order_3+1
order_4=order_4+1

return
end function func

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      F U N C 1 (xi)
! ----------------------------------------------------------------------------------
! Negative of squared difference surface. Ibidem. func(xi), but modified to be used 
! instead in all configuration space, and not only in the region where new geometries 
! are located in the case of symmetric molecules. If no symmetry exist: func1 = func

! *** Input ***
! xi        <-- 

function func1(xi)

use nrtype
USE dynamic_parameters
implicit none

INTERFACE
  FUNCTION func_actual_min(xi)
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_min
  END FUNCTION func_actual_min
end interface
INTERFACE
  FUNCTION func_actual_seed(xi)
    use nrtype
    USE dynamic_parameters
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP) :: func_actual_seed
  END FUNCTION func_actual_seed
end interface

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func1
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),tampon,h2wn,dist,temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=acos(-1d0)

 if(xi(1)>rmax(1).or.xi(1)<rmin(1))then
   func1=0d0
   return
 endif
do j=2,3
 if(xi(j)>1d0.or.xi(j)<-1d0)then
   func1=0d0
   return
 endif
enddo
 if(xi(4)>pii.or.xi(4)<(-1.d0*pii))then
   func1=0d0
   return
 endif


!do j=1,4
!  if(xi(j)>rmax(j).or.xi(j)<rmin(j))then
!    func1=0d0
!    return
!  endif
!enddo


call INT_Cart(cart,xi,mass,natom1,natom2,ref1,ref2)
call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
if(dist_flag==1) then
  func1=0d0
  return
endif

if(low_grid>0) then
  temp=func_actual_seed(xi)
  if(temp>Max_E_seed)then
    func1=0d0
    return
  endif
endif

temp1=func_actual_min(xi)
if(subzero==0)then
  if(temp1>Max_E)then
    func1=0d0
    return
  endif
endif
if(subzero==1)then
  if(temp1+temp>Max_E)then
    func1=0d0
    return
  endif
endif

allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1))
allocate(PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))

jac3=xi
count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)

quitt=0! number of expansions included in interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo
!write(701,*) quitt
12 norm=0d0

temp=0d0
jac3=xi
jac4=jac3
jac4(1)=exp(alpha*jac4(1)**xbeta)
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)
do i=1,quitt
  jj=ind8(count3+1-i)
  weight=ind7(ind8(count3+1-i)) 
  norm=norm+weight
  temp=temp+weight*b2(1,jj)
  count=1
  do R=1,order_1
    do L1=0,order_2
      do L2=0,order_3
        if((L1+L2)<order_4+1)then
          do M=0,min(L1,L2)
           count=count+1
           temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
          enddo
        endif
      enddo
    enddo
  enddo
enddo
func1=temp/norm
   
! lower basis
order_1=order_1-1
order_2=order_2-1
order_3=order_3-1
order_4=order_4-1
norm=0d0
temp=0d0
do i=1,quitt  
  jj=ind8(count3+1-i)   
  weight=ind7(ind8(count3+1-i)) 
  norm=norm+weight   
  temp=temp+weight*b2_lower(1,jj)
  count=1
  do R=1,order_1
    do L1=0,order_2
      do L2=0,order_3
        if((L1+L2)<order_4+1)then
          do M=0,min(L1,L2)
           count=count+1
           temp=temp+weight*b2_lower(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
          enddo
        endif
      enddo
    enddo
  enddo
enddo
deallocate(ind7,ind8)
func1=-abs(func1-temp/norm)**2
!write(550+myid,*) xi
!write(550+myid,*) func1

! restore order of basis   
order_1=order_1+1
order_2=order_2+1
order_3=order_3+1
order_4=order_4+1

return
end function func1

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      F U N C _ A C T U A L _ M I N (xi)
! ----------------------------------------------------------------------------------
! Energy of minimal basis and high-level ab initio

! *** Input ***
! xi        <-- 

function func_actual_min(xi)

use nrtype
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_min
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=acos(-1d0)
jac3=xi
allocate(ind7(count3),ind8(count3),PM1(0:order_2_min+1,0:order_2_min+1))
allocate(PM2(0:order_3_min+1,0:order_3_min+1),PD1(0:order_2_min+1,0:order_2_min+1))
allocate(PD2(0:order_3_min+1,0:order_3_min+1))

count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)

quitt=0! number of expansions included in interpolation
do ip=1,count3
   if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
   quitt=quitt+1
enddo
!write(701,*) quitt 

12 jac3=xi
Jac4=jac3
jac4(1)=exp(alpha*jac4(1)**xbeta)
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
call LPMN(order_3_min+1,order_3_min,order_3_min,jac4(3),PM2,PD2)
norm=0d0
temp=0d0
do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count3+1-i)) 
   ! write(665,*) weight
   norm=norm+weight
   temp=temp+weight*b2_minimal(1,jj)
   count=1
   do R=1,order_1_min
     do L1=0,order_2_min
       do L2=0,order_3_min
         if((L1+L2)<order_4_min+1)then
           do M=0,min(L1,L2)
             count=count+1
             temp=temp+weight*b2_minimal(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)* &
                  cos(dble(M)*jac4(4))
           enddo
         endif
       enddo
     enddo
   enddo
!  endif
enddo

func_actual_min=temp/norm

return
end function func_actual_min

!C ROTMOL
!C rotate a molecule

      SUBROUTINE rotmol (n, x, y, u)
!      IMPLICIT CHARACTER (A-Z)
implicit none
 integer :: n,i
 real*8 :: x(3, n),y(3, n),u(3, 3)
!      DOUBLEPRECISION x (3, n)
!      DOUBLEPRECISION y (3, n)
!      DOUBLEPRECISION u (3, 3)
 real*8 :: yx, yy, yz
!      DOUBLEPRECISION yx, yy, yz

      DO 10000 i = 1, n
        yx = u(1, 1) * x(1, i) + u(1, 2) * x(2, i) + u(1, 3) * x(3, i)
        yy = u(2, 1) * x(1, i) + u(2, 2) * x(2, i) + u(2, 3) * x(3, i)
        yz = u(3, 1) * x(1, i) + u(3, 2) * x(2, i) + u(3, 3) * x(3, i)

        y (1, i) = yx
        y (2, i) = yy
        y (3, i) = yz
10000   CONTINUE

      RETURN
      END

        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
!C       =====================================================
!C       Purpose: Compute the associated Legendre functions 
!C                Pmn(x) and their derivatives Pmn'(x)
!C       Input :  x  --- Argument of Pmn(x)
!C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
!C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
!C                mm --- Physical dimension of PM and PD
!C       Output:  PM(m,n) --- Pmn(x)
!C                PD(m,n) --- Pmn'(x)
!C       =====================================================

IMPLICIT DOUBLE PRECISION (P,X)
DIMENSION ::  PM(0:MM,0:MM),PD(0:MM,0:MM)

        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-(I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)*(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END

