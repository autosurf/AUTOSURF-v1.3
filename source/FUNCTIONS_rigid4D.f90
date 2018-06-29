!*********************************  A U T O F I T  *********************************
!===================================================================================
!-----------------------------------------------------------------------------------
!-                                                                                 -
!-        AUTOFIT Package: A set of programs for the automated construction        -
!-              of Potential Energy Surfaces on van der Waals systems              -
!-                                                                                 -
!-----------------------------------------------------------------------------------
!===================================================================================
!***********************************************************************************
!-             Set of Fortran90 functions for "IMLS_rigid4D" PROGRAM              -
!***********************************************************************************



!                                F U N C T I O N S                             


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
jac4(1)=exp(alpha*jac4(1))
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
 if(xi(j)>0.999d0.or.xi(j)<-0.999d0)then
   func1=0d0
   return
 endif
enddo
 if(xi(4)>rmax(4).or.xi(4)<(-1.d0*rmax(4)))then
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
jac4(1)=exp(alpha*jac4(1))
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
!      D F U N C _ A C T U A L _ A N A L 1 (xi)
! ----------------------------------------------------------------------------------
! Energy & analytic gradient of largest basis and high-level ab initio

! *** Input ***
! xi        <-- 

function dfunc_actual_anal1(xi)

use nrtype 
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP), DIMENSION(size(xi)) :: x2,x3
REAL(SP), DIMENSION(size(xi)+1) :: dfunc_actual_anal1
real*8 :: tampon,tampon2,scale
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,count,jp,jj,kk,R,M
real*8 :: temp(4),weight,norm,h2wn,jac3(4),jac4(4),temp2(4+1),temp3,temp4,diff(4)
real*8,allocatable :: ind7(:),weight_grad(:,:),weight_grad2(:,:)
real*8,allocatable :: PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)
real*8 :: somme,somme2,norm_grad(4),grad2(4),valeur,pii

pii=acos(-1d0)
allocate(ind7(count3),ind8(count3))

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
norm=0d0
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 13
  norm=norm+ind7(ind8(count3+1-ip))
  quitt=quitt+1
enddo
13 continue

allocate(weight_grad(quitt,4),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1))
allocate(PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))
weight_grad=0d0
norm_grad=0d0
!write(*,*) 'made it15'

do i=1,quitt
  jj=ind8(count3+1-i)
  Jac4=coords(jj,:)
  scale=sqrt((1d0-jac3(2)**2)*(1d0-jac4(2)**2)*(1d0-jac3(3)**2)*(1d0-jac4(3)**2))
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  jac4(1)=(jac3(1)-jac4(1))*(W_a**2)
  jac4(2)=(acos(jac3(2))-acos(jac4(2)))/sqrt(1d0-jac3(2)**2)
  jac4(3)=(acos(jac3(3))-acos(jac4(3)))/sqrt(1d0-jac3(3)**2)
  jac4(4)=jac3(4)-jac4(4)
  if(jac4(4)>pii)then
    jac4(4)=jac4(4)-2d0*pii
  endif
  if(jac4(4)<-pii)then
    jac4(4)=jac4(4)+2d0*pii
  endif
  jac4(4)=jac4(4)*scale
  somme2=0d0
  somme2=somme/(d(jj)**2)
  temp=0d0
  if(somme>1d-10)then
    do ip=1,4
      temp(ip)=Jac4(ip)*(-2d0)*ind7(ind8(count3+1-i))*&
      ((1.0d0/(d(jj)**2))+(zz/2)*((somme2**(zz/2))/((somme2**(zz/2))+epss))*(1.0d0/(somme)))
    enddo
  else
    temp=0d0
  endif
  weight_grad(i,:)=temp
  do ip=1,4
    norm_grad(ip)=norm_grad(ip)+weight_grad(i,ip)
  enddo
enddo
!write(*,*) 'made it16'

temp2=0d0
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)
do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   temp=0d0
   grad2=0d0
   valeur=0d0
   weight=ind7(ind8(count3+1-i)) 
   valeur=valeur+weight*b2(1,jj)
   count=1
   do R=1,order_1
     do L1=0,order_2
       do L2=0,order_3
         if((L1+L2)<order_4+1)then
           do M=0,min(L1,L2)
            count=count+1
            valeur=valeur+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(1)=grad2(1)+b2(count,jj)*dble(R)*alpha*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)* &
                     cos(dble(M)*jac4(4))
            grad2(2)=grad2(2)+b2(count,jj)*(jac4(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(3)=grad2(3)+b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*jac4(4))
            grad2(4)=grad2(4)+b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)* &
                     sin(dble(M)*jac4(4)))
           enddo
         endif
       enddo
     enddo
   enddo
   temp2(1)=temp2(1)+valeur
   do k=1,4
      temp(k)=(weight/norm)*grad2(k)
   enddo
   temp2(2:4+1)=temp2(2:4+1)+temp
   do k=1,4
      temp2(k+1)=temp2(k+1)+(valeur/(weight*norm))*weight_grad(i,k)-&
                 (1.0d0/norm**2)*valeur*norm_grad(k)
   enddo
!  endif
enddo

dfunc_actual_anal1(1)=temp2(1)/norm
dfunc_actual_anal1(2:5)=temp2(2:5)

deallocate(ind7,ind8,weight_grad)

return
end function dfunc_actual_anal1


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D F U N C _ A C T U A L _ S E E D (xi)
! ----------------------------------------------------------------------------------
! Energy & analytic gradient of minimal basis and low-level ab initio

! *** Input ***
! xi        <-- 

function dfunc_actual_seed(xi)

use nrtype 
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP), DIMENSION(size(xi)) :: x2,x3
REAL(SP), DIMENSION(size(xi)+1) :: dfunc_actual_seed
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,count,jp,jj,kk,R,M
integer,allocatable :: ind8(:)
real*8 :: tampon,tampon2
real*8 :: temp(4),weight,norm,h2wn,jac3(4),jac4(4),temp2(4+1),temp3,temp4,diff(4)
real*8 :: somme,somme2,norm_grad(4),grad2(4),valeur,pii,scale
real*8,allocatable :: ind7(:),weight_grad(:,:),weight_grad2(:,:)
real*8,allocatable :: PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)

pii=acos(-1d0)
allocate(ind7(count_seed),ind8(count_seed))

jac3=xi
count=0
do ip=1,count_seed 
  count=count+1
  Jac4=coords_seed(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo
call indexxy(count_seed,ind7,ind8)

quitt=0
norm=0d0
do ip=1,count_seed
  if(ind7(ind8(count_seed))/ind7(ind8(count_seed+1-ip))>1d11) goto 13
! if(pot(ind8(count_seed+1-ip))<E_limit)then
   norm=norm+ind7(ind8(count_seed+1-ip))
! endif
  quitt=quitt+1
enddo
13 continue
allocate(weight_grad(quitt,4),PM1(0:order_2_min+1,0:order_2_min+1))
allocate(PM2(0:order_3_min+1,0:order_3_min+1),PD1(0:order_2_min+1,0:order_2_min+1))
allocate(PD2(0:order_3_min+1,0:order_3_min+1))
weight_grad=0d0
norm_grad=0d0

!write(*,*) 'made it15'
do i=1,quitt
  jj=ind8(count_seed+1-i)
  Jac4=coords(jj,:)
  scale=sqrt((1d0-jac3(2)**2)*(1d0-jac4(2)**2)*(1d0-jac3(3)**2)*(1d0-jac4(3)**2))
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  jac4(1)=(jac3(1)-jac4(1))*(W_a**2)
  jac4(2)=(acos(jac3(2))-acos(jac4(2)))/sqrt(1d0-jac3(2)**2)
  jac4(3)=(acos(jac3(3))-acos(jac4(3)))/sqrt(1d0-jac3(3)**2)
  jac4(4)=jac3(4)-jac4(4)
  if(jac4(4)>pii)then
    jac4(4)=jac4(4)-2d0*pii
  endif
  if(jac4(4)<-pii)then
    jac4(4)=jac4(4)+2d0*pii
  endif
  jac4(4)=jac4(4)*scale
  somme2=0d0
  somme2=somme/(d_seed(jj)**2)
  temp=0d0
  if(somme>1d-10)then
    do ip=1,4
     temp(ip)=Jac4(ip)*(-2d0)*ind7(ind8(count_seed+1-i))*((1.0d0/(d_seed(jj)**2))+ &
     (zz_low/2)*((somme2**(zz_low/2))/((somme2**(zz_low/2))+epss))*(1.0d0/(somme)))
    enddo
  else
    temp=0d0
  endif
  weight_grad(i,:)=temp
  do ip=1,4
    norm_grad(ip)=norm_grad(ip)+weight_grad(i,ip)
  enddo
enddo
!write(*,*) 'made it16'

temp2=0d0
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
call LPMN(order_3_min+1,order_3_min,order_3_min,jac4(3),PM2,PD2)
do i=1,quitt
  jj=ind8(count_seed+1-i)
!  if(pot(jj)<E_limit)then
   temp=0d0
   grad2=0d0
   valeur=0d0
   weight=ind7(ind8(count_seed+1-i)) 
   valeur=valeur+weight*b2_seed(1,jj)
   count=1
   do R=1,order_1_min
     do L1=0,order_2_min
       do L2=0,order_3_min
         if((L1+L2)<order_4_min+1)then
           do M=0,min(L1,L2)
            count=count+1
            valeur=valeur+weight*b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(1)=grad2(1)+b2_seed(count,jj)*dble(R)*alpha*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)* &
                     cos(dble(M)*jac4(4))
            grad2(2)=grad2(2)+b2_seed(count,jj)*(jac4(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(3)=grad2(3)+b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*jac4(4))
            grad2(4)=grad2(4)+b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)* &
                     sin(dble(M)*jac4(4)))
           enddo
         endif
       enddo
     enddo
   enddo
   temp2(1)=temp2(1)+valeur
   do k=1,4
     temp(k)=(weight/norm)*grad2(k)
   enddo
   temp2(2:4+1)=temp2(2:4+1)+temp
   do k=1,4
     temp2(k+1)=temp2(k+1)+(valeur/(weight*norm))*weight_grad(i,k)- &
                (1.0d0/norm**2)*valeur*norm_grad(k)
   enddo
!  endif
enddo
   
dfunc_actual_seed(1)=temp2(1)/norm
dfunc_actual_seed(2:5)=temp2(2:5)

deallocate(ind7,ind8,weight_grad)

return
end function dfunc_actual_seed


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D F U N C _ A C T U A L _ A N A L 2 (xi)
! ----------------------------------------------------------------------------------
! Energy & analytic gradient of secondary basis and high-level ab initio

! *** Input ***
! xi        <-- 

function dfunc_actual_anal2(xi)

use nrtype 
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP), DIMENSION(size(xi)) :: x2,x3
REAL(SP), DIMENSION(size(xi)+1) :: dfunc_actual_anal2
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,count,jp,jj,kk,R,M
integer,allocatable :: ind8(:)
real*8 :: tampon,tampon2
real*8 :: somme,somme2,norm_grad(4),grad2(4),valeur,pii,scale
real*8 :: temp(4),weight,norm,h2wn,jac3(4),jac4(4),temp2(4+1),temp3,temp4,diff(4)
real*8,allocatable :: ind7(:),weight_grad(:,:),weight_grad2(:,:)
real*8,allocatable :: PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)

pii=acos(-1d0)

! lower basis   
order_1=order_1-1
order_2=order_2-1
order_3=order_3-1
order_4=order_4-1

jac3=xi
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1))
allocate(PM2(0:order_3+1,0:order_3+1),PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))

count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)

quitt=0
norm=0d0
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 13
!  if(pot(ind8(count3+1-ip))<E_limit)then
   norm=norm+ind7(ind8(count3+1-ip))
!  endif
  quitt=quitt+1
enddo
13 allocate(weight_grad(quitt,4))

weight_grad=0d0
norm_grad=0d0
!write(*,*) 'made it15'

do i=1,quitt
  jj=ind8(count3+1-i)
  Jac4=coords(jj,:)
  scale=sqrt((1d0-jac3(2)**2)*(1d0-jac4(2)**2)*(1d0-jac3(3)**2)*(1d0-jac4(3)**2))
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  jac4(1)=(jac3(1)-jac4(1))*(W_a**2)
  jac4(2)=(acos(jac3(2))-acos(jac4(2)))/sqrt(1d0-jac3(2)**2)
  jac4(3)=(acos(jac3(3))-acos(jac4(3)))/sqrt(1d0-jac3(3)**2)
  jac4(4)=jac3(4)-jac4(4)
  if(jac4(4)>pii)then
    jac4(4)=jac4(4)-2d0*pii
  endif
  if(jac4(4)<-pii)then
    jac4(4)=jac4(4)+2d0*pii
  endif
  jac4(4)=jac4(4)*scale
  somme2=0d0
  somme2=somme/(d(jj)**2)
  temp=0d0
  if(somme>1d-10)then
    do ip=1,4
      temp(ip)=Jac4(ip)*(-2d0)*ind7(ind8(count3+1-i))*((1.0d0/(d(jj)**2))+ &
      (zz/2)*((somme2**(zz/2))/((somme2**(zz/2))+epss))*(1.0d0/somme))
    enddo
  else
    temp=0d0
  endif
  weight_grad(i,:)=temp
  do ip=1,4
    norm_grad(ip)=norm_grad(ip)+weight_grad(i,ip)
  enddo
enddo
!write(*,*) 'made it16'

temp2=0d0
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)

do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   temp=0d0
   grad2=0d0
   valeur=0d0
   weight=ind7(ind8(count3+1-i)) 
   valeur=valeur+weight*b2_lower(1,jj)
   count=1
   do R=1,order_1
     do L1=0,order_2
       do L2=0,order_3
         if((L1+L2)<order_4+1)then
           do M=0,min(L1,L2)
            count=count+1
            valeur=valeur+weight*b2_lower(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(1)=grad2(1)+b2_lower(count,jj)*dble(R)*alpha*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)* &
                     cos(dble(M)*jac4(4))
            grad2(2)=grad2(2)+b2_lower(count,jj)*(jac4(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
            grad2(3)=grad2(3)+b2_lower(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*jac4(4))
            grad2(4)=grad2(4)+b2_lower(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)* &
                     sin(dble(M)*jac4(4)))
           enddo
         endif
       enddo
     enddo
   enddo
   temp2(1)=temp2(1)+valeur
   do k=1,4
     temp(k)=(weight/norm)*grad2(k)
   enddo
   temp2(2:4+1)=temp2(2:4+1)+temp
   do k=1,4
      temp2(k+1)=temp2(k+1)+(valeur/(weight*norm))*weight_grad(i,k)- &
      (1.0d0/norm**2)*valeur*norm_grad(k)
   enddo
!  endif
enddo
   
dfunc_actual_anal2(1)=temp2(1)/norm
dfunc_actual_anal2(2:5)=temp2(2:5)

deallocate(ind7,ind8,weight_grad)

! restore order of basis   
order_1=order_1+1
order_2=order_2+1
order_3=order_3+1
order_4=order_4+1

return
end function dfunc_actual_anal2


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D F U N C (xi)
! ----------------------------------------------------------------------------------
!  Analytic gradient of negative of squared difference surface
! 'func' and 'dfunc' must be what is used by canned CJ-minimization code

! *** Input ***
! xi        <-- 

function dfunc(xi)

use nrtype 
USE dynamic_parameters
implicit none

INTERFACE
  function dfunc_actual_anal1(xi)
    use nrtype 
    USE dynamic_parameters
    implicit none
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP), DIMENSION(size(xi)) :: x2,x3
    REAL(SP), DIMENSION(size(xi)+1) :: dfunc_actual_anal1
  end function dfunc_actual_anal1
end interface
INTERFACE
  function dfunc_actual_anal2(xi)
    use nrtype 
    USE dynamic_parameters
    implicit none
    REAL(SP), DIMENSION(:), INTENT(IN) :: xi
    REAL(SP), DIMENSION(size(xi)) :: x2,x3
    REAL(SP), DIMENSION(size(xi)+1) :: dfunc_actual_anal2
  end function dfunc_actual_anal2
end interface
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
REAL(SP), DIMENSION(size(xi)) :: dfunc,x2,x3
integer :: i,j,k
real*8 :: temp(5),val1,val2,grad1(4),grad2(4),tampon,tampon1,h2wn,dist,jac3(4),pii

pii=acos(-1d0)
do j=1,4
  if(xi(j)>rmax(j).or.xi(j)<rmin(j))then
    dfunc=0d0
    return
  endif
enddo

call INT_Cart(cart,xi,mass,natom1,natom2,ref1,ref2)
call cart_to_bdist_inter(cart,natom1,natom2,dist_tol,dist_flag)
if(dist_flag==1) then
   dfunc=0d0
   return
endif

if(low_grid>0) then
  tampon=func_actual_seed(xi)
  if(tampon>Max_E_seed)then
    dfunc=0d0
    return
  endif
endif

tampon1=func_actual_min(xi)
if(subzero==0)then
  if(tampon1>Max_E)then
    dfunc=0d0
    return
  endif
endif
if(subzero==1)then
  if(tampon1+tampon>Max_E)then
    dfunc=0d0
    return
  endif
endif

temp=dfunc_actual_anal1(xi)

val1=temp(1)
grad1(:)=temp(2:5)

temp=dfunc_actual_anal2(xi)
val2=temp(1)
grad2(:)=temp(2:5)

dfunc(:)=-2d0*(val1-val2)*(grad1(:)-grad2(:))
!write(570+myid,*) xi
!write(570+myid,*) dfunc

return
end function dfunc


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D F U N C (xi)
! ----------------------------------------------------------------------------------
! Energy & analytic gradient of secondary basis and high-level ab initio

! *** Input ***
! xi        <-- 

function func_actual(xi)

use nrtype
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=acos(-1d0)
jac3=xi
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1))
allocate(PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))

count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)

quitt=0! number of expansions included in the interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo
!write(701,*) quitt 

12 jac3=xi
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
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
            temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
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
!      D F U N C (xi)
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
jac4(1)=exp(alpha*jac4(1))
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


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      D F U N C (xi)
! ----------------------------------------------------------------------------------
!  Analytic gradient of negative of squared difference surface
! 'func' and 'dfunc' must be what is used by canned CJ-minimization code

! *** Input ***
! xi        <-- 

function func_actual_seed(xi)

use nrtype
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_seed
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)

pii=acos(-1d0)

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
  ind7(count)=exp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo
call indexxy(count_seed,ind7,ind8)

quitt=0! number of expansions are included in interpolation
do ip=1,count_seed
   if(ind7(ind8(count_seed))/ind7(ind8(count_seed+1-ip))>1d11) goto 12
   quitt=quitt+1
enddo
!write(701,*) quitt 

12 jac3=xi
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
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
                  cos(dble(M)*jac4(4))
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
!      D B R E N T (ax,bx,cx,func,dfunc,tol,xmin)
! ----------------------------------------------------------------------------------
!  

! *** Input ***
! xi        <-- 

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
!BL
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
	tol1=tol*abs(x)+ZEPS
	tol2=2.0_sp*tol1
	if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) exit
	if (abs(e) > tol1) then
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
				d=merge(d1,d2, abs(d1) < abs(d2))
			else
				d=merge(d1,d2,ok1)
			end if
			if (abs(d) <= abs(0.5_sp*olde)) then
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) &
					d=sign(tol1,xm-x)
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
	if (abs(d) >= tol1) then
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
!BL
SUBROUTINE mov3(a,b,c,d,e,f)
REAL(SP), INTENT(IN) :: d,e,f
REAL(SP), INTENT(OUT) :: a,b,c
a=d
b=e
c=f
END SUBROUTINE mov3
END FUNCTION dbrent




