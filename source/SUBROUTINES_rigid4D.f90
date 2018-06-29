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
!-            Set of Fortran90 subroutines for "IMLS_rigid4D" PROGRAM              -
!***********************************************************************************



!                                S U B R O U T I N E S                              


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
 integer :: i,j,k,kp,lab,ierr,natom1,natom2
 real*8 :: internal(6),internal2(4),cart((natom1+natom2)*3),mass(natom1+natom2),   &
            ref1(natom1*3),ref1_temp(natom1*3),ref2(natom2*3),ref2_temp(natom2*3)
 real*8 :: cart_mat1(3,natom1),cart_mat(3,natom1+natom2),cart_ref1(3,natom1),cm(3),&
            cart_ref2(3,natom2),cart_ref1t(3,natom1),cart_ref2t(3,natom2),U_rot(3,3)
 real*8 :: cart_mat2(3,natom2),cart_frag2(natom2*3),quat(4),quat2(4),pii,vec1(3)
 real*8 :: gamma1,gamma2,beta1,beta2,alpha1,alpha2,vec2(3)
 
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












subroutine fact(n,p)
  integer :: n,p,i
  p=1
  do i=1,n
     p=p*i
  enddo
end subroutine fact

subroutine binomial(n,m,b)
  integer :: n,m,b,num,denom1,denom2,p
  call fact(n,p)
  num=p
  call fact(n-m,p)
  denom1=p
  call fact(m,p)
  denom2=p
  b=num/(denom1*denom2)
end subroutine binomial

subroutine beta_term(j,m,mp,x,dx,ddx)
  implicit none
  integer :: j,m,mp,n,k,cof1,cof2,cof3,lam
  real*8 :: x,px,dx,ddx,theta,ddpx,temp,alfa,beta
  k=min(j+m,j-m,j+mp,j-mp)
  if(k==j+m)then
     alfa=dble(mp-m)
     lam=mp-m
  endif
  if(k==j-m)then
     alfa=dble(m-mp)
     lam=0
  endif
  if(k==j+mp)then
     alfa=dble(m-mp)
     lam=0
  endif
  if(k==j-mp)then
     alfa=dble(mp-m)
     lam=mp-m
  endif   
  beta=dble(2*j)-dble(2*k)-alfa
  call jacobi_pol(k,alfa,beta,x,px)
  cof1=(-1)**lam
  call binomial(2*j-k,int(k+alfa),cof2)
  call binomial(int(k+beta),int(beta),cof3)
  dx=dble(cof1)*(dble(cof2)**0.5d0)*(dble(cof3)**(-0.5d0))*(sqrt((1d0-x)/2d0)**alfa)*(sqrt((1d0+x)/2d0)**beta)*px
  call jacobi_pol(k-1,alfa+1d0,beta+1d0,x,ddpx)
  ddpx=0.5d0*dble(k+alfa+beta+1)*ddpx
  temp=(beta*sqrt(0.5d0-x/2d0)**alfa*((x+1d0)/2d0)**(-1d0+dble(beta)/2d0)-alfa*((1d0-x)/2d0)**(-1d0+dble(alfa)/2d0)*sqrt((x+1d0)/2d0)**(beta))/4d0
  ddx=dble(cof1)*(dble(cof2)**0.5d0)*(dble(cof3)**(-0.5d0))* (ddpx*(sqrt((1d0-x)/2d0)**alfa)*(sqrt((1d0+x)/2d0)**beta)+ (temp)*px)  
  return
end subroutine beta_term



subroutine jacobi_pol(n,alfa,beta,x,px)
  implicit none
  integer :: i,n
  real*8 :: alfa,beta,cx(0:n),x,c1,c2,c3,c4,px
  if (n<0) then
     px=0d0
     return
  endif
  cx(0)=1d0
  if (n==0) then
     px=cx(n)
     return
  endif
  cx(1)=(1.0d0+0.5d0*(alfa+beta))*x+0.5d0*(alfa-beta)
  do i=2,n
     c1=2.0d0*dble(i)*(dble(i)+alfa+beta)*(dble(2*i-2)+alfa+beta)
     c2=(dble(2*i-1)+alfa+beta)*(dble(2*i)+alfa+beta)*(dble(2*i-2)+alfa+beta)
     c3=(dble(2*i-1)+alfa+beta)*(alfa+beta)*(alfa-beta)
     c4=-dble(2)*(dble(i-1)+alfa)*(dble(i-1)+beta)*(dble(2*i)+alfa+beta)
     cx(i)=((c3+c2*x)*cx(i-1)+c4*cx(i-2))/c1
  enddo
  px=cx(n)
  return
end subroutine jacobi_pol









subroutine Cart_INT(cart_temp,internal2,mass,natom1,natom2,ref1,ref2)
  implicit none
  integer :: i,j,k,lab,ierr,natom1,natom2
  real*8 :: internal(6),internal2(4),cart((natom1+natom2)*3),cart_temp((natom1+natom2)*3),mass(natom1+natom2),ref1(natom1*3),ref2(natom2*3)
  real*8 :: cart_mat1(3,natom1),cart_mat(3,natom1+natom2),cart_matt(3,natom1+natom2),cart_ref1(3,natom1),cart_ref2(3,natom2),U_rot(3,3),U_rot2(3,3),cm(3)
  real*8 :: cart_mat2(3,natom2),cart_frag2(natom2*3),quat(4),quat2(4),ref1_temp(natom1*3),ref2_temp(natom2*3),pii,theta,phi
  real*8 :: gamma1,gamma2,beta1,beta2,alpha1,alpha2
  pii=acos(-1d0)
  cart=cart_temp
  ref1_temp=ref1
  ref2_temp=ref2
  

!!! set center of mass of frag1 at origin
  call rm_cmass(cart,mass,natom1+natom2,natom1)
!!! set center of mass of ref1 at origin
  call rm_cmass(ref1_temp,mass(1:natom1),natom1,natom1)


!!!!!! find c.m. of frag2
  call cmass(cart,cm,mass,natom1+natom2,natom2)


!write(67,*) cm

!!!! set "R"  
  internal(1)=sqrt(cm(1)**2+cm(2)**2+cm(3)**2)
  !!!!set theta
  theta=acos(cm(3)/internal(1))

!!!!set phi
  if(abs(cos(theta))>1d0-1d-7)then
     phi=0d0
  else
     phi=atan2(cm(2),cm(1))
  endif
!!!!!!!!
!write(67,*) theta
!write(67,*) phi

call vec_to_mat2(cart,cart_mat,natom1+natom2)

U_rot=0d0
!!$U_rot(1,1)=cos(-phi)
!!$U_rot(2,2)=cos(-phi)
!!$U_rot(3,3)=1d0
!!$U_rot(1,2)=-sin(-phi)
!!$U_rot(2,1)=sin(-phi)

U_rot(1,1)=cos(phi)*cos(theta)
U_rot(2,2)=cos(phi)
U_rot(3,3)=cos(theta)
U_rot(1,2)=sin(phi)*cos(theta)
U_rot(1,3)=-sin(theta)
U_rot(2,1)=-sin(phi)
U_rot(3,1)=cos(phi)*sin(theta)
U_rot(3,2)=sin(phi)*sin(theta)
call rotmol(natom1+natom2,cart_mat,cart_matt,U_rot)
cart_mat=cart_matt



!!!!!!!!!!!!!delete after tests
!!$call mat_to_vec2(cart_mat,cart,natom1+natom2)
!!$call cmass(cart,cm,mass,natom1+natom2,natom2)
!!$write(67,*) cm
!!$stop

!!!!!!!!!!!

cart_mat1(1:3,1:natom1)=cart_mat(1:3,1:natom1)
call vec_to_mat2(ref1_temp,cart_ref1,natom1)
call qtrfit(natom1,cart_ref1,cart_mat1,mass(1:natom1),quat,U_rot,ierr)



call mat_to_vec2(cart_mat,cart,natom1+natom2)


  
  cart_frag2=cart(3*natom1+1:3*(natom1+natom2))

!  do i=1,3*natom2
!     write(66,*) cart_frag2(i)
!  enddo

!!! set center of mass of frag2 at origin
  call rm_cmass(cart_frag2,mass(natom1+1:natom1+natom2),natom2,natom2)
!!! set center of mass of ref2 at origin
  call rm_cmass(ref2_temp,mass(natom1+1:natom1+natom2),natom2,natom2)

  call vec_to_mat2(cart_frag2,cart_mat2,natom2)
  call vec_to_mat2(ref2_temp,cart_ref2,natom2)



  !!!!solve for Eulers (U_rot)
  call qtrfit(natom2,cart_ref2,cart_mat2,mass(natom1+1:natom1+natom2),quat2,U_rot2,ierr)


beta1=acos(U_rot(3,3))
if(sin(beta1)<1d-7)then
   gamma1=atan2(-U_rot(1,2),U_rot(1,1))
   alpha1=0d0
else
   alpha1=atan2(U_rot(2,3),U_rot(1,3))
   gamma1=atan2(U_rot(3,2),-U_rot(3,1))
endif
beta2=acos(U_rot2(3,3))
if(sin(beta2)<1d-7)then
   gamma2=atan2(-U_rot2(1,2),U_rot2(1,1))
   alpha2=0d0
else
   alpha2=atan2(U_rot2(2,3),U_rot2(1,3))
   gamma2=atan2(U_rot2(3,2),-U_rot2(3,1))
endif
!write(*,*) gamma1
!write(*,*) beta1
!write(*,*) gamma2
!write(*,*) beta2
!write(*,*) alpha1
!write(*,*) alpha2
internal(2)=gamma1
!internal(3)=beta1
internal(3)=cos(beta1)
internal(4)=gamma2
!internal(5)=beta2
internal(5)=cos(beta2)
internal(6)=alpha1-alpha2
if(internal(6)>pii)then
   internal(6)=internal(6)-2d0*pii
endif
if(internal(6)<-pii)then
   internal(6)=internal(6)+2d0*pii
endif

internal2(1)=internal(1)
internal2(2)=internal(3)
internal2(3)=internal(5)
!internal2(4)=cos(internal(6))
internal2(4)=abs(internal(6))!
return
end subroutine Cart_INT




subroutine cmass(cart,cm,mass,natom,natom2)
integer :: k,kp,natom,natom2
real*8 :: mass(natom),cart(natom*3),mtot,cm(3)
mtot=0d0
do k=natom-natom2+1,natom
   mtot=mtot+mass(k)
enddo
cm=0d0
do k=natom-natom2+1,natom
   do kp=1,3
      cm(kp)=cm(kp)+cart((k-1)*3+kp)*mass(k)
   enddo
enddo
cm=cm/mtot
return
end subroutine cmass
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
integer :: i,j
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
!-----------------------------------------------------------------------
subroutine dcart_dint(xtemp,mass,natom1,natom2,ref1,ref2,b)
  implicit none
  integer :: i,j,k,tabbb,i1,natom1,natom2
  real*8 :: x(4),xtemp(4),b(3*(natom1+natom2),4),x2(4),x3(4),x4(4),x5(4),d(3*(natom1+natom2)),d2(3*(natom1+natom2)),d3(3*(natom1+natom2)),d4(3*(natom1+natom2)),mass(natom1+natom2)
  real*8 :: ref1(3*natom1),ref2(3*natom2)
  x=xtemp
  do i1=1,size(x)
     x2=x
     x3=x
     x4=x
     x5=x
     x2(i1)=x2(i1)+2d-6
     x3(i1)=x3(i1)+1d-6
     x4(i1)=x4(i1)-1d-6
     x5(i1)=x5(i1)-2d-6
     call INT_Cart(d,x2,mass,natom1,natom2,ref1,ref2)
     call INT_Cart(d2,x3,mass,natom1,natom2,ref1,ref2)
     call INT_Cart(d3,x4,mass,natom1,natom2,ref1,ref2)
     call INT_Cart(d4,x5,mass,natom1,natom2,ref1,ref2)
     do j=1,3*(natom1+natom2)
        b(j,i1)=(-d(j)+8d0*d2(j)-8d0*d3(j)+d4(j))/(12d0*1d-6)
     enddo
  enddo  
  return
end subroutine dcart_dint
!-----------------------------------------------------------------------

!!$!-----------------------------------------------------------------------

subroutine basis_calc_seed(i,ind,ind2)
  use dynamic_parameters
  implicit none
  integer :: i,i2,l,j,ip,jp,k,ipp,jpp,l1,l2,l3,l4,ind2(count_seed),count,count4,info,lwork,rankk,count5,jj,kk,R,M
  real*8 :: grad_vec(4),weight,ind(count_seed),jac3(4),Jac4(4),diff(4),temp1,temp2,pii,scale,factor
  real*8,allocatable :: design3(:,:),b_seed(:),s3(:),work2(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
  pii=acos(-1d0)
  lwork=10*max((4*(ab_flag2-1)+1)*support,basis_3)
!  write(*,*) 'lwork', lwork
  allocate(work2(lwork))
  allocate(s3(basis_3))
  allocate(b_seed((4*(ab_flag2-1)+1)*support),design3((4*(ab_flag2-1)+1)*support,basis_3))
  allocate(PM1(0:order_2_min+1,0:order_2_min+1),PM2(0:order_3_min+1,0:order_3_min+1),PD1(0:order_2_min+1,0:order_2_min+1),PD2(0:order_3_min+1,0:order_3_min+1))
!!$    if(myid==0)then
!!$      write(*,*) 'support', support
!!$   endif
!!$   if(myid==1)then
!!$      write(*,*) 'support', support
!!$   endif 

  design3=0d0
 
  do i2=1,support
  
     jj=ind2(count_seed+1-i2)
     Jac4=coords_seed(jj,:)

     diff=jac4
     diff(1)=exp(alpha*diff(1))
!     diff(4)=acos(diff(4))
     call LPMN(order_2_min+1,order_2_min,order_2_min,diff(2),PM1,PD1)
     call LPMN(order_3_min+1,order_3_min,order_3_min,diff(3),PM2,PD2)
!!!!!!!!!!
     if(ab_flag2==2)then
        grad_vec=grad_seed(jj,:)
     endif
     weight=ind(ind2(count_seed+1-i2))
!!!!!!!!!!!!!!!
     
!!$     if(pot(jj)>Max_E)then
!!$        factor=1d0+(pot(jj)-Max_E)/E_range
!!$        scale=1d0/factor**2
!!$        weight=weight*scale
!!$     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     design3(i2,1)=weight
     b_seed(i2)=weight*pot_seed(jj)
     
!     write(620+myid,*) b_seed(i2)
     if(ab_flag2==2) then
        do j=1,4
           b_seed(i2+j*support)=weight*grad_vec(j)
        enddo
     endif


!!$!!!design3
     count5=1
!     write(640,*) design3(i2,1),i2,count5
     do R=1,order_1_min
        do L1=0,order_2_min
           do L2=0,order_3_min
              if((L1+L2)<order_4_min+1)then
                 do M=0,min(L1,L2)
                    count5=count5+1
                    design3(i2,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
!                    write(640,*) design3(i2,count5),i2,count5
                    if(ab_flag2==2) then
                       design3(i2+support,count5)=weight*dble(R)*alpha*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+2*support,count5)=weight*(diff(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+3*support,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+4*support,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)*sin(dble(M)*diff(4)))
                    endif
                 enddo
              endif
           enddo
        enddo
     enddo

     
     
  enddo
!!$  if(myid==0) then
!!$     write(*,*)  (4*(ab_flag2-1)+1)*support,basis_3,(4*(ab_flag2-1)+1)*support
!!$     write(*,*) 'lwork', lwork
!!$  endif
!  call dgelss(320,40,1,design3,320,b_seed,320,s3,-1d0,rankk,work2,-1,info)
  call dgelss((4*(ab_flag2-1)+1)*support,basis_3,1,design3,(4*(ab_flag2-1)+1)*support,b_seed,(4*(ab_flag2-1)+1)*support,s3,1d-12,rankk,work2,lwork,info)
!  write(600+myid,*) work2(1)
!  write(600+myid,*) info

  b2_seed(1:basis_3,i)=b_seed(1:basis_3)
  return
end subroutine basis_calc_seed



subroutine basis_calc(i,ind,ind2)
  use dynamic_parameters
  implicit none
  integer :: i,i2,l,j,ip,jp,k,ipp,jpp,l1,l2,l3,l4,ind2(count3),count,count4,info,lwork,rankk,count5,jj,kk,R,M,cond
  real*8 :: grad_vec(4),weight,ind(count3),jac3(4),Jac4(4),diff(4),temp,temp1,temp2,pii,scale,factor,somme,rcond,bench,norm
  real*8,allocatable :: design(:,:),design2(:,:),design3(:,:),b(:),b_temp(:),b_sol(:),b_lower(:),b_minimal(:),s(:),s2(:),s3(:),work2(:)
  real*8,allocatable :: PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:),design_temp(:,:),design2_temp(:,:)
  pii=acos(-1d0)
  lwork=1d7
  allocate(work2(lwork),design_temp((4*(ab_flag-1)+1)*support,basis_1),design2_temp((4*(ab_flag-1)+1)*support,basis_2))
  allocate(s(basis_1),s2(basis_2),s3(basis_3),b_temp((4*(ab_flag-1)+1)*support),b_sol((4*(ab_flag-1)+1)*support))
  allocate(design((4*(ab_flag-1)+1)*support,basis_1),b((4*(ab_flag-1)+1)*support),design3((4*(ab_flag-1)+1)*support,basis_3),b_minimal((4*(ab_flag-1)+1)*support))
  allocate(design2((4*(ab_flag-1)+1)*support,basis_2),b_lower((4*(ab_flag-1)+1)*support),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1),PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))
  
  design=0d0
  design2=0d0
  design3=0d0
  design_temp=0d0
  design2_temp=0d0
  do i2=1,support
  
     jj=ind2(count3+1-i2)
     Jac4=coords(jj,:)

     diff=jac4
     diff(1)=exp(alpha*diff(1))
!     diff(4)=acos(diff(4))
     call LPMN(order_2+1,order_2,order_2,diff(2),PM1,PD1)
     call LPMN(order_3+1,order_3,order_3,diff(3),PM2,PD2)
!!!!!!!!!!
     if(ab_flag==2)then
        grad_vec=grad(jj,:)
     endif
     weight=ind(ind2(count3+1-i2))
!!!!!!!!!!!!!!!
     
!!$     if(pot(jj)>Max_E)then
!!$        factor=1d0+(pot(jj)-Max_E)/E_range
!!$        scale=1d0/factor**2
!!$        weight=weight*scale
!!$     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     design(i2,1)=weight
     design2(i2,1)=weight
     design3(i2,1)=weight
     b(i2)=weight*pot(jj)
     if(ab_flag==2) then
        do j=1,4
           b(i2+j*support)=weight*grad_vec(j)
        enddo
     endif
!!!design    
     count=1
     do R=1,order_1
        do L1=0,order_2
           do L2=0,order_3
              if((L1+L2)<order_4+1)then
                 do M=0,min(L1,L2)
                    count=count+1
                    design(i2,count)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                    if(ab_flag==2) then
                       design(i2+support,count)=weight*dble(R)*alpha*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design(i2+2*support,count)=weight*(diff(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design(i2+3*support,count)=weight*(diff(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*diff(4))
                       design(i2+4*support,count)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)*sin(dble(M)*diff(4)))
                    endif
                 enddo
              endif
           enddo
        enddo
     enddo

     
     
!!$!!!design2 
     count4=1
     do R=1,order_1-1
        do L1=0,order_2-1
           do L2=0,order_3-1
              if((L1+L2)<order_4)then
                 do M=0,min(L1,L2)
                    count4=count4+1
                    design2(i2,count4)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                     if(ab_flag==2) then
                       design2(i2+support,count4)=weight*dble(R)*alpha*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design2(i2+2*support,count4)=weight*(diff(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design2(i2+3*support,count4)=weight*(diff(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*diff(4))
                       design2(i2+4*support,count4)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)*sin(dble(M)*diff(4)))
                    endif                    
                 enddo
              endif
           enddo
        enddo
     enddo

!!$!!!design3
     count5=1
     do R=1,order_1_min
        do L1=0,order_2_min
           do L2=0,order_3_min
              if((L1+L2)<order_4_min+1)then
                 do M=0,min(L1,L2)
                    count5=count5+1
                    design3(i2,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                    if(ab_flag==2) then
                       design3(i2+support,count5)=weight*dble(R)*alpha*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+2*support,count5)=weight*(diff(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+3*support,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*diff(4))
                       design3(i2+4*support,count5)=weight*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)*sin(dble(M)*diff(4)))
                    endif   
                 enddo
              endif
           enddo
        enddo
     enddo

     
     
  enddo
  design_temp=design
  design2_temp=design2
  b_temp=b
  b_lower=b
  b_minimal=b
  !   write(199,*) count
  call dgelss((4*(ab_flag-1)+1)*support,basis_1,1,design,(4*(ab_flag-1)+1)*support,b,(4*(ab_flag-1)+1)*support,s,-1d0,rankk,work2,lwork,info)
  b_sol=b
  norm=0d0
  somme=0d0
  do i2=2,support
     temp=0d0
     jj=ind2(count3+1-i2)
     Jac4=coords(jj,:)
     diff=jac4
     diff(1)=exp(alpha*diff(1))
     call LPMN(order_2+1,order_2,order_2,diff(2),PM1,PD1)
     call LPMN(order_3+1,order_3,order_3,diff(3),PM2,PD2)
     weight=ind(ind2(count3+1-i2))
     norm=norm+weight**2
     temp=temp+b_sol(1)
     count=1
     do R=1,order_1
        do L1=0,order_2
           do L2=0,order_3
              if((L1+L2)<order_4+1)then
                 do M=0,min(L1,L2)
                    count=count+1
                    temp=temp+b_sol(count)*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                 enddo
              endif
           enddo
        enddo
     enddo
     somme=somme+(weight*(temp-pot(jj)))**2
  enddo
  bench=sqrt(somme/norm)
  write(5000+myid,*) 'all',bench,rankk,i
  do cond=1,10
     design=design_temp
     b=b_temp
     rcond=1d-15*10d0**(cond)
     call dgelss((4*(ab_flag-1)+1)*support,basis_1,1,design,(4*(ab_flag-1)+1)*support,b,(4*(ab_flag-1)+1)*support,s,rcond,rankk,work2,lwork,info)
     norm=0d0
     somme=0d0
     do i2=2,support
        temp=0d0
        jj=ind2(count3+1-i2)
        Jac4=coords(jj,:)
        diff=jac4
        diff(1)=exp(alpha*diff(1))
        call LPMN(order_2+1,order_2,order_2,diff(2),PM1,PD1)
        call LPMN(order_3+1,order_3,order_3,diff(3),PM2,PD2)
        weight=ind(ind2(count3+1-i2))
        norm=norm+weight**2
        temp=temp+b(1)
        count=1
        do R=1,order_1
           do L1=0,order_2
              do L2=0,order_3
                 if((L1+L2)<order_4+1)then
                    do M=0,min(L1,L2)
                       count=count+1
                       temp=temp+b(count)*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                    enddo
                 endif
              enddo
           enddo
        enddo
        somme=somme+(weight*(temp-pot(jj)))**2
     enddo
     
     write(5000+myid,*) rcond,sqrt(somme/norm),rankk
     if((sqrt(somme/norm)/bench)<1.1d0)then
        b_sol=b
     else
        goto 10
     endif
  enddo
  
  
10 b2(1:basis_1,i)=b_sol(1:basis_1) 


  call dgelss((4*(ab_flag-1)+1)*support,basis_2,1,design2,(4*(ab_flag-1)+1)*support,b_lower,(4*(ab_flag-1)+1)*support,s2,-1d0,rankk,work2,lwork,info)
  b_sol=b_lower
  norm=0d0
  somme=0d0
  do i2=2,support
     temp=0d0
     jj=ind2(count3+1-i2)
     Jac4=coords(jj,:)
     diff=jac4
     diff(1)=exp(alpha*diff(1))
     call LPMN(order_2+1,order_2,order_2,diff(2),PM1,PD1)
     call LPMN(order_3+1,order_3,order_3,diff(3),PM2,PD2)
     weight=ind(ind2(count3+1-i2))
     norm=norm+weight**2
     temp=temp+b_sol(1)
     count=1
     do R=1,order_1-1
        do L1=0,order_2-1
           do L2=0,order_3-1
              if((L1+L2)<order_4)then
                 do M=0,min(L1,L2)
                    count=count+1
                    temp=temp+b_sol(count)*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                 enddo
              endif
           enddo
        enddo
     enddo
     somme=somme+(weight*(temp-pot(jj)))**2
  enddo
  bench=sqrt(somme/norm)
  write(5000+myid,*) 'all lower',bench,rankk
  do cond=1,10
     design2=design2_temp
     b_lower=b_temp
     rcond=1d-15*10d0**(cond)
     call dgelss((4*(ab_flag-1)+1)*support,basis_2,1,design2,(4*(ab_flag-1)+1)*support,b_lower,(4*(ab_flag-1)+1)*support,s2,rcond,rankk,work2,lwork,info)
     norm=0d0
     somme=0d0
     do i2=2,support
        temp=0d0
        jj=ind2(count3+1-i2)
        Jac4=coords(jj,:)
        diff=jac4
        diff(1)=exp(alpha*diff(1))
        call LPMN(order_2+1,order_2,order_2,diff(2),PM1,PD1)
        call LPMN(order_3+1,order_3,order_3,diff(3),PM2,PD2)
        weight=ind(ind2(count3+1-i2))
        norm=norm+weight**2
        temp=temp+b_lower(1)
        count=1
        do R=1,order_1-1
           do L1=0,order_2-1
              do L2=0,order_3-1
                 if((L1+L2)<order_4)then
                    do M=0,min(L1,L2)
                       count=count+1
                       temp=temp+b_lower(count)*(diff(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*diff(4))
                    enddo
                 endif
              enddo
           enddo
        enddo
        somme=somme+(weight*(temp-pot(jj)))**2
     enddo
     
     write(5000+myid,*) rcond,sqrt(somme/norm),rankk
     if((sqrt(somme/norm)/bench)<1.1d0)then
        b_sol=b_lower
     else
        goto 11
     endif
  enddo
  
  
11 b2_lower(1:basis_2,i)=b_sol(1:basis_2) 
  




  call dgelss((4*(ab_flag-1)+1)*support,basis_3,1,design3,(4*(ab_flag-1)+1)*support,b_minimal,(4*(ab_flag-1)+1)*support,s3,1d-12,rankk,work2,lwork,info)

!  b2(1:basis_1,i)=b(1:basis_1)
!  b2_lower(1:basis_2,i)=b_lower(1:basis_2)
  b2_minimal(1:basis_3,i)=b_minimal(1:basis_3)
  return
end subroutine basis_calc

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
subroutine cart_to_bdist(x,d,natomm,nbdist)
!  use sizes
  implicit none
  integer :: i,j,k,tabbb,natomm,nbdist
  real*8 :: x(3*natomm),d(nbdist),summ
  tabbb=0
  do i=1,natomm
     do j=i+1,natomm
        tabbb=tabbb+1
        summ=0d0
        do k=1,3
           summ=summ+(x(3*(i-1)+k)-x(3*(j-1)+k))**2
        enddo
        d(tabbb)=sqrt(summ)
     enddo
  enddo
  return
end subroutine cart_to_bdist
!!$!-----------------------------------------------------------------------!
subroutine map_int(inter)
 
  implicit none
  real*8 :: inter(4),internal(2,4)
  internal(1,:)=inter(:)
  internal(2,:)=inter(:)
  internal(2,2)=-internal(1,3)
  internal(2,3)=-internal(1,2)
  if(internal(1,2)<internal(2,2))then
     inter(:)=internal(1,:)
  else
     inter(:)=internal(2,:)
  endif

  return 
end subroutine map_int

subroutine map_cart(cart,mass,ref1,ref2,natom1,natom2)
 
  implicit none
  integer :: i,j,k,natom1,natom2,natom,ierr
  real*8 :: cart(3*(natom1+natom2)),cart2(3*(natom1+natom2)),internal(2,4),mass(natom1+natom2),grad(3*(natom1+natom2)),grad2(3*(natom1+natom2)),ref1(3*natom1),ref2(3*natom2)
  real*8 :: pii,cart3(3*(natom1+natom2)),grad3(3*(natom1+natom2)),grad_int(2,4),cart_perms(2,18),grad_perms(2,18),bmat(3*(natom1+natom2),4)
  real*8 :: cart_mat(3,natom1+natom2),cart_mat2(3,natom1+natom2),grad_mat(3,natom1+natom2),quat(4),U_rot(3,3)
  natom=natom1+natom2
  pii=acos(-1d0)
  call perm_cart(cart,grad,cart_perms,grad_perms)
  do j=1,2
     call Cart_INT(cart_perms(j,:),internal(j,:),mass,natom1,natom2,ref1,ref2)
  enddo
  if(internal(1,2)<internal(2,2))then
     cart(:)=cart_perms(1,:)
  else
     cart(:)=cart_perms(2,:)
  endif

  return 
end subroutine map_cart
subroutine symmetry2(int_temp,grad,internal,grad_int,mass,ref1,ref2,natom1,natom2,flag,exch,flip1,flip2)

  implicit none
  integer ::i,j,k,kk,natom1,natom2,natom,ierr,flag,exch,flip1,flip2,symparts,ind889((1+exch)*(1+flip1)*(1+flip2))
  real*8 ::cart(3*(natom1+natom2)),cart2(3*(natom1+natom2)),internal((1+exch)*(1+flip1)*(1+flip2),4),mass(natom1+natom2),grad(3*(natom1+natom2)),grad2(3*(natom1+natom2)),ref1(3*natom1),ref2(3*natom2)
  real*8 ::pii,cart3(3*(natom1+natom2)),grad3(3*(natom1+natom2)),grad_int((1+exch)*(1+flip1)*(1+flip2),4),cart_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
  real*8 ::grad_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2)),bmat(3*(natom1+natom2),4),somme
  real*8 ::cart_mat(3,natom1+natom2),cart_mat2(3,natom1+natom2),grad_mat(3,natom1+natom2),grad_matt(3,natom1+natom2),quat(4),U_rot(3,3),ind888((1+exch)*(1+flip1)*(1+flip2))
  real*8 ::grad_permst((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2)),cart_permst((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
  real*8 :: int_temp(4),grad_temp(4)
  natom=natom1+natom2
 symparts=(1+exch)*(1+flip1)*(1+flip2)
  pii=acos(-1d0)

if(flag==2)then
    call dcart_dint(int_temp,mass,natom1,natom2,ref1,ref2,bMat)
    grad_temp=matmul(transpose(bMat),grad)
endif

  call perm_int(int_temp,grad_temp,internal,grad_int,exch,flip1,flip2,natom1,natom2)
end subroutine symmetry2

subroutine perm_int(int_temp,grad_temp,internal,grad_int,exch,flip1,flip2,natom1,natom2)
  implicit none
  integer :: i,j,k,exch,flip1,flip2,natom,natom1,natom2,swap1,swap2
  real*8 :: cart(3*(natom1+natom2)),cart2(3*(natom1+natom2)),cart3(3*(natom1+natom2)),cart_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
 real*8 :: grad(3*(natom1+natom2)),grad2(3*(natom1+natom2)),grad3(3*(natom1+natom2)),grad_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
 real*8 :: internal((1+exch)*(1+flip1)*(1+flip2),4),grad_int((1+exch)*(1+flip1)*(1+flip2),4)
 real*8 :: int_temp(4),grad_temp(4),pii
 natom=natom1+natom2
pii=dacos(-1d0)
  do i=1,(1+exch)*(1+flip1)*(1+flip2)
  internal(i,:)=int_temp(:)
  grad_int(i,:)=grad_temp(:)
  enddo

if(flip1>0)then

internal(2,2)=-int_temp(2)
internal(2,4)=pii-int_temp(4)
grad_int(2,2)=-grad_temp(2)
grad_int(2,4)=-grad_temp(4)


if(flip2>0)then
!internal(3,:)=internal(2,:)
!internal(3,3)=-internal(2,3)
!grad_int(3,:)=grad_int(2,:)
!grad_int(3,3)=-grad_int(2,3)

internal(3,3)=-int_temp(3)
internal(3,4)=pii-int_temp(4)
grad_int(3,3)=-grad_temp(3)
grad_int(3,4)=-grad_temp(4)

internal(4,2)=-int_temp(2)
internal(4,3)=-int_temp(3)
grad_int(4,2)=-grad_temp(2)
grad_int(4,3)=-grad_temp(3)


if(exch>0) then

internal(5,2)=-int_temp(3)
internal(5,3)=-int_temp(2)
grad_int(5,2)=-grad_temp(3)
grad_int(5,3)=-grad_temp(2)

internal(6,2)=-int_temp(3)
internal(6,3)=int_temp(2)
internal(6,4)=pii-int_temp(4)
grad_int(6,2)=-grad_temp(3)
grad_int(6,3)=grad_temp(2)
grad_int(6,4)=-grad_temp(4)

internal(7,2)=int_temp(3)
internal(7,3)=-int_temp(2)
internal(7,4)=pii-int_temp(4)
grad_int(7,2)=grad_temp(3)
grad_int(7,3)=-grad_temp(2)
grad_int(7,4)=-grad_temp(4)

internal(8,2)=int_temp(3)
internal(8,3)=int_temp(2)
grad_int(8,2)=grad_temp(3)
grad_int(8,3)=grad_temp(2)

endif
endif
endif

if(flip1<1) then
if(flip2>0)then
!internal(2,2)=-int_temp(3)
!internal(2,4)=pii-int_temp(4)
!grad_int(2,2)=-grad_temp(3)
!grad_int(2,4)=-grad_temp(4)

internal(2,3)=-int_temp(3)
internal(2,4)=pii-int_temp(4)
grad_int(2,3)=-grad_temp(3)
grad_int(2,4)=-grad_temp(4)

endif
endif

if(flip1<1) then
if(flip2<1) then
if(exch>0) then

internal(2,2)=-int_temp(3)
internal(2,3)=-int_temp(2)

grad_int(2,2)=-grad_temp(3)
grad_int(2,3)=-grad_temp(2)

endif
endif
endif

  return
end subroutine Perm_int





subroutine symmetry(cart,grad,internal,grad_int,mass,ref1,ref2,natom1,natom2,flag,exch,flip1,flip2)
 
  implicit none
  integer :: i,j,k,kk,natom1,natom2,natom,ierr,flag,exch,flip1,flip2,symparts,ind889((1+exch)*(1+flip1)*(1+flip2))
  real*8 :: cart(3*(natom1+natom2)),cart2(3*(natom1+natom2)),internal((1+exch)*(1+flip1)*(1+flip2),4),mass(natom1+natom2),grad(3*(natom1+natom2)),grad2(3*(natom1+natom2)),ref1(3*natom1),ref2(3*natom2)
  real*8 :: pii,cart3(3*(natom1+natom2)),grad3(3*(natom1+natom2)),grad_int((1+exch)*(1+flip1)*(1+flip2),4),cart_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
  real*8 :: grad_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2)),bmat(3*(natom1+natom2),4),somme
  real*8 :: cart_mat(3,natom1+natom2),cart_mat2(3,natom1+natom2),grad_mat(3,natom1+natom2),grad_matt(3,natom1+natom2),quat(4),U_rot(3,3),ind888((1+exch)*(1+flip1)*(1+flip2))
  real*8 :: grad_permst((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2)),cart_permst((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
  natom=natom1+natom2
 symparts=(1+exch)*(1+flip1)*(1+flip2)
  pii=acos(-1d0)
  call perm_cart(cart,grad,cart_permst,grad_permst,exch,flip1,flip2,natom1,natom2)
  do j=1,symparts
     call Cart_INT(cart_permst(j,:),internal(j,:),mass,natom1,natom2,ref1,ref2)
     cart_perms(j,:)=cart_permst(j,:)
!     if(j==1)then
!        call dcart_dint(internal(j,:),mass,natom1,natom2,ref1,ref2,bMat)
!        grad_int(j,:)=matmul(transpose(bMat),grad_perms(j,:))
!     else
write(68,*) cart_permst(j,:)
write(68,*) 


  enddo

  if(flag==2)then
   do j=1,symparts  
    call INT_Cart(cart3,internal(j,:),mass,natom1,natom2,ref1,ref2)
    write(69,*) cart3(:)
    write(69,*)
    
    do k=1,symparts
    cart2=cart3
    cart_permst(k,:)=cart_perms(k,:)
    call rm_cmass(cart2,mass,natom,natom)
    call rm_cmass(cart_permst(k,:),mass,natom,natom)
    call vec_to_mat2(cart2,cart_mat,natom)
    call vec_to_mat2(cart_permst(k,:),cart_mat2,natom)
    call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
    call rotmol(natom,cart_mat2,cart_mat,U_rot)
    call mat_to_vec2(cart_mat,cart_permst(k,:),natom)

    somme=0d0
    do kk=1,3*(natom1+natom2)
    somme=somme+(cart2(kk)-cart_permst(k,kk))**2d0
    enddo
    ind888(k)=somme
    write(67,*) somme
    enddo
    write(67,*)
    call indexxy(symparts,ind888,ind889)
    cart_perms(j,:)=cart_permst(ind889(1),:)
    grad_perms(j,:)=grad_permst(ind889(1),:)
    cart2=cart3

    call rm_cmass(cart2,mass,natom,natom)
    call rm_cmass(cart_perms(j,:),mass,natom,natom)
    call vec_to_mat2(cart2,cart_mat,natom)
    call vec_to_mat2(cart_perms(j,:),cart_mat2,natom)
    call vec_to_mat2(grad_perms(j,:),grad_mat,natom)
    call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
    call rotmol(natom,grad_mat,grad_matt,U_rot)
    grad_mat=grad_matt
    call mat_to_vec2(grad_mat,grad_perms(j,:),natom)
    call dcart_dint(internal(j,:),mass,natom1,natom2,ref1,ref2,bMat)
    grad_int(j,:)=matmul(transpose(bMat),grad_perms(j,:))


!    grad_perms(j,:)=grad_permst(ind889(1),:) 
!    write(66,*) ind889(1),ind888(ind889(1)) 

!        write(*,*) cart_perms(j,:)
!        write(*,*)
!        write(*,*) grad_perms(j,:)
!        write(*,*)
!        call rm_cmass(cart2,mass,natom,natom)
!        call rm_cmass(cart_perms(j,:),mass,natom,natom)
!        call vec_to_mat2(cart2,cart_mat,natom)
!        call vec_to_mat2(cart_perms(j,:),cart_mat2,natom)
!        call vec_to_mat2(grad_perms(j,:),grad_mat,natom)
!        call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
!        write(*,*) U_rot(:,:)
!        write(*,*)

!        call rotmol(natom,grad_mat,grad_matt,U_rot)
!        grad_mat=grad_matt
!        call mat_to_vec2(grad_mat,grad_perms(j,:),natom)
!        call dcart_dint(internal(j,:),mass,natom1,natom2,ref1,ref2,bMat)
!        grad_int(j,:)=matmul(transpose(bMat),grad_perms(j,:))
        write(*,*) grad_int(j,:)
     enddo
  endif
  return 
end subroutine symmetry



subroutine perm_cart(cart,grad,cart_perms,grad_perms,exch,flip1,flip2,natom1,natom2)
  implicit none
  integer :: i,j,k,exch,flip1,flip2,natom,natom1,natom2,swap1,swap2
  real*8 :: cart(3*(natom1+natom2)),cart2(3*(natom1+natom2)),cart3(3*(natom1+natom2)),cart_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2)),grad(3*(natom1+natom2)),grad2(3*(natom1+natom2)),grad3(3*(natom1+natom2)),grad_perms((1+exch)*(1+flip1)*(1+flip2),3*(natom1+natom2))
 natom=natom1+natom2
  cart2=cart
  grad2=grad
  do i=1,(1+exch)*(1+flip1)*(1+flip2)
  cart_perms(i,:)=cart2(:)
  grad_perms(i,:)=grad2(:)
  enddo
  
swap1=floor(dble(natom1)/2d0)
swap2=floor(dble(natom2)/2d0)

if(flip1>0)then
do i=0,swap1-1
cart_perms(2,3*i+1:3*i+3)=cart_perms(1,1+natom1*3-(3*i+3):1+natom1*3-(3*i+1))
cart_perms(2,1+natom1*3-(3*i+3):1+natom1*3-(3*i+1))=cart_perms(1,3*i+1:3*i+3)
grad_perms(2,3*i+1:3*i+3)=grad_perms(1,1+natom1*3-(3*i+3):1+natom1*3-(3*i+1))
grad_perms(2,1+natom1*3-(3*i+3):1+natom1*3-(3*i+1))=grad_perms(1,3*i+1:3*i+3)
enddo


if(flip2>0)then

cart_perms(4,:)=cart_perms(2,:)
grad_perms(4,:)=grad_perms(2,:)

do i=0,swap2-1
cart_perms(3,3*natom1+3*i+1:3*natom1+3*i+3)=cart_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
cart_perms(3,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=cart_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
grad_perms(3,3*natom1+3*i+1:3*natom1+3*i+3)=grad_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
grad_perms(3,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=grad_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
enddo

do i=0,swap2-1
cart_perms(4,3*natom1+3*i+1:3*natom1+3*i+3)=cart_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
cart_perms(4,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=cart_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
grad_perms(4,3*natom1+3*i+1:3*natom1+3*i+3)=grad_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
grad_perms(4,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=grad_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
enddo


if(exch>0) then
cart_perms(5,1:3*natom1)=cart_perms(1,3*natom1+1:3*natom)
cart_perms(5,3*natom1+1:3*natom)=cart_perms(1,1:3*natom1)
cart_perms(6,1:3*natom1)=cart_perms(2,3*natom1+1:3*natom)
cart_perms(6,3*natom1+1:3*natom)=cart_perms(2,1:3*natom1)
cart_perms(7,1:3*natom1)=cart_perms(3,3*natom1+1:3*natom)
cart_perms(7,3*natom1+1:3*natom)=cart_perms(3,1:3*natom1)
cart_perms(8,1:3*natom1)=cart_perms(4,3*natom1+1:3*natom)
cart_perms(8,3*natom1+1:3*natom)=cart_perms(4,1:3*natom1)

grad_perms(5,1:3*natom1)=grad_perms(1,3*natom1+1:3*natom)
grad_perms(5,3*natom1+1:3*natom)=grad_perms(1,1:3*natom1)
grad_perms(6,1:3*natom1)=grad_perms(2,3*natom1+1:3*natom)
grad_perms(6,3*natom1+1:3*natom)=grad_perms(2,1:3*natom1)
grad_perms(7,1:3*natom1)=grad_perms(3,3*natom1+1:3*natom)
grad_perms(7,3*natom1+1:3*natom)=grad_perms(3,1:3*natom1)
grad_perms(8,1:3*natom1)=grad_perms(4,3*natom1+1:3*natom)
grad_perms(8,3*natom1+1:3*natom)=grad_perms(4,1:3*natom1)
endif
endif
endif

if(flip1<1) then
if(flip2>0)then
do i=0,swap2-1
cart_perms(2,3*natom1+3*i+1:3*natom1+3*i+3)=cart_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
cart_perms(2,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=cart_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
grad_perms(2,3*natom1+3*i+1:3*natom1+3*i+3)=grad_perms(1,1+natom*3-(3*i+3):1+natom*3-(3*i+1))
grad_perms(2,1+natom*3-(3*i+3):1+natom*3-(3*i+1))=grad_perms(1,3*natom1+3*i+1:3*natom1+3*i+3)
enddo
endif
endif

if(flip1<1) then
if(flip2<1) then
if(exch>0) then
cart_perms(2,1:3*natom1)=cart_perms(1,3*natom1+1:3*natom)
cart_perms(2,3*natom1+1:3*natom)=cart_perms(1,1:3*natom1)
grad_perms(2,1:3*natom1)=grad_perms(1,3*natom1+1:3*natom)
grad_perms(2,3*natom1+1:3*natom)=grad_perms(1,1:3*natom1)
endif
endif
endif

  return
end subroutine Perm_cart

! ******************************************************************************
  subroutine crossprod(ra,rb,rc)
! cross-product of two vectors rc = ra * rb
  implicit double precision(a-h,o-z)          
  real*8, intent(in) :: ra(3),rb(3)
  real*8, intent(out) :: rc(3)

  rc(1)=ra(2)*rb(3) - ra(3)*rb(2)
  rc(2)=ra(3)*rb(1) - ra(1)*rb(3)
  rc(3)=ra(1)*rb(2) - ra(2)*rb(1)
! normalize
  xlen=dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2)
  rc(1)=rc(1)/xlen
  rc(2)=rc(2)/xlen
  rc(3)=rc(3)/xlen  

  return
  end subroutine crossprod
!c-----------------------------------------------------------------------



!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

!!! dlinmin.f90
!c-----------------------------------------------------------------------


MODULE df1dim_mod
	USE nrtype
	INTEGER(I4B) :: ncom
	REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
!BL
	FUNCTION df1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: df1dim
	INTERFACE
		FUNCTION dfunc(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt,df
	allocate(xt(ncom),df(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	df(:)=dfunc(xt)
	df1dim=dot_product(df,xicom)
	deallocate(xt,df)
	END FUNCTION df1dim
END MODULE df1dim_mod

	SUBROUTINE dlinmin(p,xi,fret)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : mnbrak,dbrent
	USE df1dim_mod
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), TARGET :: p,xi
	REAL(SP), PARAMETER :: TOL=1.0e-8_sp
	REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'dlinmin')
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE dlinmin


!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

!!! frprmn.f90
!c-----------------------------------------------------------------------


	SUBROUTINE frprmn(p,ftol,iter,fret)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : dlinmin
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(p)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
!	INTEGER(I4B), PARAMETER :: ITMAX=10000
	INTEGER(I4B) :: ITMAX
	REAL(SP), PARAMETER :: EPS=1.0e-18_sp
	INTEGER(I4B) :: its,i
	REAL(SP) :: dgg,fp,gam,gg
	REAL(SP), DIMENSION(size(p)) :: g,h,xi,gxi
	fp=func(p)
	xi=dfunc(p)
	g=-xi
	h=g
	xi=h
	ITMAX=size(p)+1
	do its=1,ITMAX
		iter=its
		call dlinmin(p,xi,fret)
!		write(*,*)2.0_sp*abs(fret-fp),ftol*(abs(fret)+abs(fp)+EPS)
		if (2.0_sp*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
		fp=fret
		xi=dfunc(p)
!		write(*,*) its,fcalls,maxval(dabs(xi)),fp
		gg=dot_product(g,g)
!		dgg=dot_product(xi,xi)
		dgg=dot_product(xi+g,xi)
		if (gg == 0.0) RETURN
		gam=dgg/gg
		g=-xi
		h=g+gam*h
		xi=h
	end do
!	call nrerror('frprmn: maximum iterations exceeded')
	END SUBROUTINE frprmn




!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

!!! mnbrak.f90
!c-----------------------------------------------------------------------



	SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: ax,bx
	REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
	REAL(SP) :: fu,q,r,u,ulim
	fa=func(ax)
	fb=func(bx)
	if (fb > fa) then
		call swap(ax,bx)
		call swap(fa,fb)
	end if
	cx=bx+GOLD*(bx-ax)
	fc=func(cx)
	do
		if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call shft(fb,fc,fu,func(u))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		end if
		call shft(ax,bx,cx,u)
		call shft(fa,fb,fc,fu)
	end do
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak



!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

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
     if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
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
            


!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

!!! ab_parallel999.f90
!c-----------------------------------------------------------------------


!!$program ab
!!$! ab is a fake program used to made ab_initio subroutine
!!$implicit none
!!$real*8,allocatable :: q(:),dq(:)
!!$character(len=3),allocatable :: nuclei(:)
!!$character(len=3) :: ff 
!!$real*8 :: pot
!!$integer :: f1,f2,natom,j 
!!$natom=3
!!$allocate(q(natom*3),dq(natom*3),nuclei(natom))
!!$
!!$open(unit=5,file='geom')
!!$!!! read geom !!!!!!!!
!!$read(5,*)ff,f1
!!$read(5,*)ff,f2
!!$read(5,*)natom
!!$do j=1,natom
!!$read(5,*)nuclei(j),q((j-1)*3+1:j*3)
!!$enddo
!!$
!!$
!!$call ab_initio(q,nuclei,pot,dq,f1,f2,natom)
!!$write(*,*) pot
!!$end program ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  ab_initio(q,nuclei,mass,pot,dq,f1,f2,natom,myid,chk_ind,rename,inp_num)
! f1 1=Gaussian, 2=Molpro, 3= ACES II 

implicit none
real*8 :: q(natom*3),dq(natom*3),mass(natom)
real*8 :: pot
character(len=3) :: nuclei(natom)
integer :: f1,f2,f3,natom,myid,chk_ind,rename,inp_num
if(f1.eq.1) then

  call wr_gauss(q,nuclei,f2,f3,natom,myid,chk_ind,rename,inp_num) 
!  call system('./gauss.run')
  call rd_gauss(q,pot,dq,nuclei,mass,f2,f3,natom,myid)
elseif(f1.eq.2) then
  call wr_molpro(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num) 
!  call system('./molpro.run')
  call rd_molpro(q,pot,dq,nuclei,f2,f3,natom,myid,mass)
elseif(f1.eq.3) then
  call wr_aces(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num) 
!  call system('./aces.run')
  call rd_aces(q,pot,dq,nuclei,f2,f3,natom,myid,mass)
elseif(f1.eq.4) then
  call wr_molpro9(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num) 
!  call system('./molpro.run')
  call rd_molpro9(q,pot,dq,nuclei,f2,f3,natom,myid,mass)
endif 
!  write(12,*)pot
!  write(12,*)dq
end subroutine  ab_initio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rd_gauss(q,pot,dq,nuclei,mass,f2,f3,natom,myid)

  implicit none
  real*8 :: q(natom*3),dq(natom*3),q_rot(natom*3),q_temp(natom*3),mass(natom)
  real*8 :: pot,cart_mat(3,natom),cart_mat2(3,natom),grad_mat(3,natom),quat(4),U_rot(3,3)
  character(len=3) :: nuclei(natom)
  character(len=80) :: line,stdm
  character :: filename*8,str*3
  integer :: f1,f2,f3,natom,info,id,i,j,k,l,flag,myid,ierr
  write(unit=str, fmt='(I3)') myid+115
      filename= "fort." // str

  open(unit=9,file=filename)
  
  rewind(9)
  flag=0
  info=0
12 call getrec(' ITN=',5,9,line,info)
  if(info==1.and.flag==0) then
     goto 15
  endif
  if(info==0) then
     id=index(line,'=')
     read(line(id+17:id+35),'(f18.10)')pot
     flag=1
     goto 12
  endif
  if(flag==1) then
     rewind(9)
     goto 17
  endif



15  rewind(9)
  flag=0
  info=0
  do while (info==0)
     call getrec('SCF Done',8,9,line,info)
     if(info==0) then
        flag=1
        id=index(line,'=')
        read(line(id+1:id+21),'(f20.16)')pot
     endif
  enddo
  if(flag==0)then
     f3=3
     close(9)
     call  wr_gauss(q,nuclei,f2,f3,natom,myid)
     do i=1,3*natom
        write(925+myid,*) q(i)
     enddo
     open(unit=9,file=filename)
     info=0
     do while (info==0)
        call getrec('SCF Done',8,9,line,info)
        if(info==0) then
           flag=1
           id=index(line,'=')
           read(line(id+1:id+21),'(f20.16)')pot
        endif
     enddo
  endif
  if(flag==0)then
     write(925+myid,*) 'could not converge'
     do i=1,3*natom
        write(925+myid,*) q(i)
     enddo
     pot=2d2
     dq(:)=0d0
     return
  endif

  rewind(9)


!!!!  
17  if(f2.eq.2) then
     call getrec('Forces at end of L703',21,9,line,info)
     if(info==1)then
        write(925+myid,*) 'could not converge'
        do i=1,3*natom
           write(925+myid,*) q(i)
        enddo
        pot=2d2
        dq(:)=0d0
        return
     endif
     read(9,'(a80)')line
     read(line,*)stdm,i,stdm,dq(1),stdm,dq(2),stdm,dq(3)
     do j=i+1,natom
        k=(j-1)*3
        read(9,'(a80)')line
        read(line,*)stdm,l,stdm,dq(k+1),stdm,dq(k+2),stdm,dq(k+3)
     enddo
!!!  uncomment to rotate from standard orientation back to input orientation
!!$     rewind(9)
!!$     call getrec('Standard orientation:',21,9,line,info)
!!$     
!!$     do j=1,4
!!$        read(9,'(a80)')line
!!$     enddo
!!$     do j=1,natom
!!$        k=(j-1)*3
!!$        read(9,'(a80)')line
!!$        read(line,*) i,i,i,q_rot(k+1),q_rot(k+2),q_rot(k+3)
!!$     enddo
!!$     q_temp=q
!!$     call rm_cmass(q_temp,mass)
!!$     call rm_cmass(q_rot,mass)
!!$     call vec_to_mat(q_temp,cart_mat,natom)
!!$     call vec_to_mat(q_rot,cart_mat2,natom)
!!$     call vec_to_mat(dq,grad_mat,natom)
!!$     call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
!!$     call rotmol(natom,grad_mat,grad_mat,U_rot)
!!$     call mat_to_vec(grad_mat,dq,natom)
  endif
  
40 continue
  return
end subroutine rd_gauss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rd_molpro9(q,pot,dq,nuclei,f2,f3,natom,myid,mass)

implicit none
real*8 :: q(natom*3),dq(natom*3),dq_temp(natom*3),dq_temp2(natom*3),q_rot(natom*3),q_temp(natom*3),q_temp2(natom*3),mass(natom),bohr
real*8 :: pot,zz(4),z,cart_mat(3,natom),cart_mat2(3,natom),grad_mat(3,natom),quat(4),U_rot(3,3)
character(len=3) :: nuclei(natom),symbord(natom)
character(len=80) :: line,st
integer :: f1,f2,f3,natom,info,id,i,j,k,l,n,myid,skip(natom),ierr,pass
character :: filename*10,str*5,punch*8,str2*2,line2*131
logical :: var
10 pass=0
213 bohr=0.529177249d0 !value used by Molpro
write(unit=str, fmt='(I5)') myid+10000
!write(unit=str2, fmt='(I2)') myid+15
!filename= "fort" // str2 //".out"
punch="pun"//str
!open(unit=9,file=filename)

! call getrec('ATOMIC COORDINATES',18,9,line,info) 
 
!      read(9,fmt='(a80)')st
!      read(9,fmt='(a80)')st
!      read(9,fmt='(a80)')st
!      do i=1,natom
!         k=(i-1)*3
 !        read(9,*)n,symbord(i),z,q_temp2(k+1),q_temp2(k+2),q_temp2(k+3)
 !     enddo
!      rewind(9)

  inquire(file=punch,exist=var)
  if(var) then
     open(unit=10,file=punch)
     call getrec('MOLPRO_ENERGY',13,10,line,info) 
     if(info==1)then
        if(pass==0) then
           call sleep(10)
           pass=1
           close(10)
           goto 213
        endif
        write(925+myid,*) 'could not converge'
        do i=1,3*natom
           write(925+myid,*) q(i)
        enddo
        pot=2d2
        dq(:)=0d0
        return
     endif
     id=index(line,'-')
     read(line(id:id+20),'(f20.15)')pot
     close(10)
  else
     if(pass==0) then
        call sleep(10)
        pass=1
        goto 213
     endif
     write(925+myid,*) 'could not converge'
     do i=1,3*natom
        write(925+myid,*) q(i)
     enddo
     pot=2d2
     dq(:)=0d0
     return
  endif


if(f2.eq.2) then  
open(unit=10,file=punch)
rewind(10)
!write(*,*) 'made it here'
!call getrec('GRADIENT FOR STATE',18,9,line,info) 
!if(info==1)then
!   rewind(9)
do i=1,natom
   call getrec('GRADIENT',8,10,line,info)
!   if(info==1)then
!      if(pass==0) then
!         call sleep(10)
!         pass=1
!         close(9)
!         goto 213
!      endif

!      write(925+myid,*) 'could not converge'
!      do i=1,3*natom
!         write(925+myid,*) q(i)
!      enddo
!      pot=2d2
!      dq(:)=0d0
!      return
!   endif
   id=index(line,':')
   read(line(id+1:80),*) dq(3*(i-1)+1:3*(i-1)+3)
!   write(*,*) dq(3*(i-1)+1:3*(i-1)+3)
enddo
close(10)
endif
return
!   call getrec2('GRADY(',6,9,line2,info)
!   id=index(line2,'=')
!   read(line2(id+5:id+15*natom),*) dq_temp(natom+1:2*natom)
!   call getrec2('GRADZ(',6,9,line2,info)
!   id=index(line2,'=')
!   read(line2(id+5:id+15*natom),*) dq_temp(2*natom+1:3*natom)
!   do i=1,3
!      do j=1,natom
!         dq_temp2((j-1)*3+i)=dq_temp((i-1)*3+j)
!         write(*,*) dq_temp2((j-1)*3+i)
!      enddo
!   enddo
!stop
!   skip=0
!   do i=1,natom
!      do j=1,natom
!         do k=1,natom
!            if(skip(k)==j) goto 14
!         enddo
!         id=INDEX(symbord(i),nuclei(j))
!         !    write(*,*) id,i,j
!         if(id.gt.0)then
           
!            q_rot((j-1)*3+1:j*3)=q_temp2((i-1)*3+1:i*3)
!            dq((j-1)*3+1:j*3)=dq_temp2((i-1)*3+1:i*3)
!         endif
!         if(id==1) then
!            skip(i)=j
!            goto 15
!         endif
!14    enddo
!15 enddo
!   close(9)
   
   
!   goto 13 !!activate rotation
!goto 10

!endif

!!$      read(9,fmt='(a80)')st
!!$      read(9,fmt='(a80)')st
!!$      read(9,fmt='(a80)')st
!!$
!!$skip=0
!!$do i=1,natom
!!$  do j=1,natom
!!$     do k=1,natom
!!$        if(skip(k)==j) goto 12
!!$     enddo
!!$    id=INDEX(symbord(i),nuclei(j))
!!$!    write(*,*) id,i,j
!!$    if(id.gt.0)then
!!$      read(9,*)n,dq((j-1)*3+1:j*3)
!!$      q_rot((j-1)*3+1:j*3)=q_temp2((i-1)*3+1:i*3)
!!$    endif
!!$    if(id==1) then
!!$       skip(i)=j
!!$       goto 11
!!$    endif
!!$12 enddo
!!$11 enddo
!!$close(9)
!!$
!!$goto 10 !!!!!!!!! rotate if allowed to reorient
!!$
!!$13      q_rot=q_rot*bohr
!!$!      do i=1,9
!!$!         write(30+myid,*) q_rot(i)
!!$!      enddo
!!$      q_temp=q
!!$!      write(30+myid,*)
!!$!      do i=1,9
!!$! !        write(30+myid,*) q_temp(i)
!!$!      enddo
!!$      call rm_cmass2(q_temp,mass,natom)
!!$      call rm_cmass2(q_rot,mass,natom)
!!$      call vec_to_mat(q_temp,cart_mat,natom)
!!$      call vec_to_mat(q_rot,cart_mat2,natom)
!!$      call vec_to_mat(dq,grad_mat,natom)
!!$      call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
!!$      call rotmol(natom,grad_mat,grad_mat,U_rot)
!!$      call mat_to_vec(grad_mat,dq,natom)
!!$
!!$
!!$!!!!!!!!!!
!!$
!!$
!!$
!!$
!!$10 continue
 return
end subroutine rd_molpro9
subroutine rd_molpro(q,pot,dq,nuclei,f2,f3,natom,myid,mass)

implicit none
real*8 :: q(natom*3),dq(natom*3),dq_temp(natom*3),dq_temp2(natom*3),q_rot(natom*3),q_temp(natom*3),q_temp2(natom*3),mass(natom),bohr
real*8 :: pot,zz(4),z,cart_mat(3,natom),cart_mat2(3,natom),grad_mat(3,natom),quat(4),U_rot(3,3)
character(len=3) :: nuclei(natom),symbord(natom)
character(len=80) :: line,st
integer :: f1,f2,f3,natom,info,id,i,j,k,l,n,myid,skip(natom),ierr,pass
character :: filename*8,str*3,punch*6
logical :: var
pass=0
213 bohr=0.529177249d0 !value used by Molpro
write(unit=str, fmt='(I3)') myid+115
filename= "fort." // str
punch="pun"//str
open(unit=9,file=filename)

 call getrec('ATOMIC COORDINATES',18,9,line,info) 
 
      read(9,fmt='(a80)')st
      read(9,fmt='(a80)')st
      read(9,fmt='(a80)')st
      do i=1,natom
         k=(i-1)*3
         read(9,*)n,symbord(i),z,q_temp2(k+1),q_temp2(k+2),q_temp2(k+3)
      enddo
      rewind(9)

  inquire(file=punch,exist=var)
  if(var) then
     open(unit=10,file=punch)
     call getrec('MOLPRO_ENERGY',13,10,line,info) 
     if(info==1)then
        if(pass==0) then
           call sleep(10)
           pass=1
           close(10)
           goto 213
        endif
        write(925+myid,*) 'could not converge'
        do i=1,3*natom
           write(925+myid,*) q(i)
        enddo
        pot=2d2
        dq(:)=0d0
        return
     endif
     id=index(line,'-')
     read(line(id:id+20),'(f20.15)')pot
     close(10)
  else
     if(pass==0) then
        call sleep(10)
        pass=1
        goto 213
     endif
     write(925+myid,*) 'could not converge'
     do i=1,3*natom
        write(925+myid,*) q(i)
     enddo
     pot=2d2
     dq(:)=0d0
     return
  endif


if(f2.eq.1) goto 10  
rewind(9)
call getrec('GRADIENT FOR STATE',18,9,line,info) 
if(info==1)then
   rewind(9)
   call getrec('GRADX(',6,9,line,info)
   if(info==1)then
      if(pass==0) then
         call sleep(10)
         pass=1
         close(9)
         goto 213
      endif

      write(925+myid,*) 'could not converge'
      do i=1,3*natom
         write(925+myid,*) q(i)
      enddo
      pot=2d2
      dq(:)=0d0
      return
   endif
   id=index(line,'=')
   read(line(id+1:id+50),*) dq_temp(1:natom)
   call getrec('GRADY(',6,9,line,info)
   id=index(line,'=')
   read(line(id+1:id+50),*) dq_temp(natom+1:2*natom)
   call getrec('GRADZ(',6,9,line,info)
   id=index(line,'=')
   read(line(id+1:id+50),*) dq_temp(2*natom+1:3*natom)
   do i=1,3
      do j=1,natom
         dq_temp2((j-1)*3+i)=dq_temp((i-1)*3+j)
      enddo
   enddo

   skip=0
   do i=1,natom
      do j=1,natom
         do k=1,natom
            if(skip(k)==j) goto 14
         enddo
         id=INDEX(symbord(i),nuclei(j))
         !    write(*,*) id,i,j
         if(id.gt.0)then
           
            q_rot((j-1)*3+1:j*3)=q_temp2((i-1)*3+1:i*3)
            dq((j-1)*3+1:j*3)=dq_temp2((i-1)*3+1:i*3)
         endif
         if(id==1) then
            skip(i)=j
            goto 15
         endif
14    enddo
15 enddo
   close(9)
   
   
!   goto 13 !!activate rotation
goto 10

endif

      read(9,fmt='(a80)')st
      read(9,fmt='(a80)')st
      read(9,fmt='(a80)')st

skip=0
do i=1,natom
  do j=1,natom
     do k=1,natom
        if(skip(k)==j) goto 12
     enddo
    id=INDEX(symbord(i),nuclei(j))
!    write(*,*) id,i,j
    if(id.gt.0)then
      read(9,*)n,dq((j-1)*3+1:j*3)
      q_rot((j-1)*3+1:j*3)=q_temp2((i-1)*3+1:i*3)
    endif
    if(id==1) then
       skip(i)=j
       goto 11
    endif
12 enddo
11 enddo
close(9)

goto 10 !!!!!!!!! rotate if allowed to reorient

13      q_rot=q_rot*bohr
!      do i=1,9
!         write(30+myid,*) q_rot(i)
!      enddo
      q_temp=q
!      write(30+myid,*)
!      do i=1,9
! !        write(30+myid,*) q_temp(i)
!      enddo
      call rm_cmass2(q_temp,mass,natom)
      call rm_cmass2(q_rot,mass,natom)
      call vec_to_mat(q_temp,cart_mat,natom)
      call vec_to_mat(q_rot,cart_mat2,natom)
      call vec_to_mat(dq,grad_mat,natom)
      call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
      call rotmol(natom,grad_mat,grad_mat,U_rot)
      call mat_to_vec(grad_mat,dq,natom)


!!!!!!!!!!




10 continue
 return
end subroutine rd_molpro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wr_molpro9(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num)

implicit none
real*8 :: q(natom*3)
character(len=3) :: nuclei(natom)
character :: line*80,str*5,filename*18,filename2*16,filename5*13,st1*2,checkpoint*24,checkpoint4*19,checkpoint5*19,checkpoint6*23,str2*4,str4*4,str3*3,st2*4,st3*4,inp*10,out*9,str5*3,st5*4
character :: punch*18,filename10*35,filename2b*26
integer :: f2,natom,i,k,id,myid,chk_ind,rename,inp_num,i3,i5,i2,i4


!!!!!!!input keywords etc
write(unit=str5, fmt='(I3)') inp_num
if(inp_num.gt.9.and.inp_num.le.99) then
   i5=2
elseif(inp_num.gt.99) then
   i5=1 
else
   i5=3
endif
st5=".abi"
filename5="molpro" // str5(i5:3) // st5

!!!processor dep. script
write(unit=str, fmt='(I5)') myid+10000
st1=".x"
filename="./Molscript" // str(1:5) // st1
filename10="./Molscript" // str(1:5) // st1
filename2="Molscript" // str(1:5) // st1
filename2b="chmod 777 "// filename2
inp="fort" // str(1:5)

punch="punch,pun"//str(1:5)//",new"


open(unit=9,file=filename)
write(9,'(a11)') '#!/bin/bash'
!write(9,*) 'ldd ./IMLS.x'
!write(9,'(a31)') '. /etc/profile.d/env-modules.sh'
!write(9,'(a37)') '. /etc/profile.d/zenv-modules-load.sh'
!write(9,'(a18)') 'module load molpro'
!write(9,'(a42)') 'rm /state/partition1/richard/molpro/*'

!write(9,*) '/nethome/users/dawesr/bin/molpro_binary -n 1 ',inp,' -s '
write(9,*) 'molpros ',inp,'-s'
write(9,*) 'rm ', filename
close(9)

 open(unit=7,file=filename5)
 open(unit=myid+15,file=inp)
 70   read(7,80)line
 id=INDEX(line,'xyz')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 70
 endif
    write(myid+15,80)line
    if(chk_ind>0)then
       write(myid+15,*) checkpoint
    endif
    write(myid+15,*) punch
!    write(myid+15,*)'geometry={nosym;'
write(myid+15,*)'noorient'
    write(myid+15,*)'geometry={'
!    write(myid+15,*)'nosym'
!    write(myid+15,*)'noorient'
    write(myid+15,*)natom
    write(myid+15,*)'MOLPRO is running for IMLS'
 do i=1,natom
    write(myid+15,107) nuclei(i),(q(k),k=i*3-2,i*3)
 enddo
    write(myid+15,*)'}'


 90   read(7,80)line
 id=INDEX(line,'---')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 90
 endif
if(f2.eq.2) write(myid+15,*)'forces'
if(f2.eq.2) write(myid+15,*)'VARSAV'
! if(f2.eq.2) write(myid+15,*)'force,numerical,central,varsav,variable=energd,proc=myproc'
!if(f2.eq.2) write(myid+15,*)'force,numerical,central,varsav,variable=energd'
 if(f2.eq.2) write(myid+15,*) 'show,gradx'
 if(f2.eq.2) write(myid+15,*) 'show,grady'
 if(f2.eq.2) write(myid+15,*) 'show,gradz'
 
 100   read(7,80)line
 id=INDEX(line,'----')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 100
 endif


   write(myid+15,*) '---'
 close(7)
 close(myid+15)
 80   format(a80)
 107  format(a2,1x,3(f24.14,1x))
 call system(filename2b)
! call system('sleep10')
 call system(filename10)
 return


end subroutine wr_molpro9
subroutine wr_molpro(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num)

implicit none
real*8 :: q(natom*3)
character(len=3) :: nuclei(natom)
character :: line*80,str*3,filename*16,filename2*14,filename5*13,st1*2,checkpoint*24,checkpoint4*19,checkpoint5*19,checkpoint6*23,str2*4,str4*4,str3*3,st2*4,st3*4,inp*8,out*9,str5*3,st5*4
character :: punch*16,filename10*33
integer :: f2,natom,i,k,id,myid,chk_ind,rename,inp_num,i3,i5,i2,i4


!!!!!!!input keywords etc
write(unit=str5, fmt='(I3)') inp_num
if(inp_num.gt.9.and.inp_num.le.99) then
   i5=2
elseif(inp_num.gt.99) then
   i5=1 
else
   i5=3
endif
st5=".abi"
filename5="molpro" // str5(i5:3) // st5

!!!processor dep. script
write(unit=str, fmt='(I3)') myid+15
if(myid+15.gt.9.and.myid+15.le.99) then
   i=2
elseif(myid+15.gt.99) then
   i=1 
else
   i=3
endif
st1=".x"
filename="./Molscript" // str(i:3) // st1
filename10="./Molscript" // str(i:3) // st1
filename2="Molscript" // str(i:3) // st1
inp="fort." // str(i:3)

write(unit=str3, fmt='(I3)') myid+115
if(myid+115.gt.9.and.myid+115.le.99) then
   i3=2
elseif(myid+115.gt.99) then
   i3=1 
else
   i3=3
endif
out="fort." // str3(i3:3)
punch="punch,pun"//str3(i3:3)//",new"

write(unit=str2, fmt='(I4)') chk_ind
if(chk_ind.gt.9.and.chk_ind.le.99) then
   i2=3
elseif(chk_ind.gt.99.and.chk_ind.le.999) then
   i2=2
elseif(chk_ind.gt.999) then
   i2=1 
else
   i2=4
endif
st2=".wfu"
st3=".bak"

write(unit=str4, fmt='(I4)') rename
if(rename.gt.9.and.rename.le.99) then
   i4=3
elseif(rename.gt.99.and.rename.le.999) then
   i4=2
elseif(rename.gt.999) then
   i4=1
else
   i4=4
endif

checkpoint="file,2,ckekk" // str2(i2:4) // "_" // str(i:3) // st2


checkpoint4="$path/ckekk" // str4(i4:4) // st3
checkpoint5="$path/ckekk" // str2(i2:4) // st3
checkpoint6="$path/ckekk" // str2(i2:4) // "_" // str(i:3) // st2



open(unit=9,file=filename)
write(9,'(a11)') '#!/bin/bash'
write(9,'(a42)') 'rm /state/partition1/richard/molpro/*'
if(chk_ind>0)then
   if(chk_ind.ne.rename)then
      write(9,*) 'cp '//checkpoint5//' '//checkpoint6
   endif
endif

write(9,*) 'molpro  < ',inp,' > ',out
if(chk_ind>0)then   
   if(rename>0) then   
      write(9,*) 'mv '//checkpoint6//' '//checkpoint4
   else
      write(9,*) 'rm ', checkpoint6
   endif
endif
close(9)





 open(unit=7,file=filename5)
 open(unit=myid+15)
 70   read(7,80)line
 id=INDEX(line,'xyz')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 70
 endif
    write(myid+15,80)line
    if(chk_ind>0)then
       write(myid+15,*) checkpoint
    endif
    write(myid+15,*) punch
    write(myid+15,*)'geometry={nosym;'
!    write(myid+15,*)'nosym'
    write(myid+15,*)'noorient'
    write(myid+15,*)natom
    write(myid+15,*)'MOLPRO is running for IMLS'
 do i=1,natom
    write(myid+15,107) nuclei(i),(q(k),k=i*3-2,i*3)
 enddo
    write(myid+15,*)'}'


 90   read(7,80)line
 id=INDEX(line,'---')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 90
 endif
if(f2.eq.2) write(myid+15,*)'force'
! if(f2.eq.2) write(myid+15,*)'force,numerical,central,varsav,variable=energd,proc=myproc'
!if(f2.eq.2) write(myid+15,*)'force,numerical,central,varsav,variable=energd'
! if(f2.eq.2) write(myid+15,*) 'show,gradx'
! if(f2.eq.2) write(myid+15,*) 'show,grady'
! if(f2.eq.2) write(myid+15,*) 'show,gradz'
    write(myid+15,80)line	
 close(7)
 close(myid+15)
 80   format(a80)
 107  format(a2,1x,3(f24.14,1x))

 call system(filename10)
 return


end subroutine wr_molpro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wr_aces(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num)

implicit none
real*8 :: q(natom*3)
character(len=3) :: nuclei(natom)
integer :: f2,natom,k,i,myid,chk_ind,rename,inp_num,i5,i4,i3,i2,i6
character :: line*80,str*3,filename*17,st1*2,checkpoint*9,str2*3,str3*3,str4*3,st2*4,st3*4,inp*8,out*11,checkpoint2*25,checkpoint3*25,checkpoint4*26,checkpoint5*25,filename2*15
character :: checkpoint6*17,str5*3,filename5*12,st5*4,str6*3,checkpoint7*13,checkpoint8*8,checkpoint9*27,checkpoint10*10,checkpoint11*10

write(unit=str5, fmt='(I3)') inp_num
if(inp_num.gt.9.and.inp_num.le.99) then
   i5=2
elseif(inp_num.gt.99) then
   i5=1 
else
   i5=3
endif
st5=".abi"
filename5="aces" // str5(i5:3) // st5


write(unit=str6, fmt='(I3)') myid+1
if(myid+1.gt.9.and.myid+1.le.99) then
   i6=2
elseif(myid+1.gt.99) then
   i6=1 
else
   i6=3
endif


write(unit=str, fmt='(I3)') myid+15
if(myid+15.gt.9.and.myid+15.le.99) then
   i=2
elseif(myid+15.gt.99) then
   i=1 
else
   i=3
endif
st1=".x"
filename="./Acesscript" // str(i:3) // st1
filename2="Acesscript" // str(i:3) // st1
inp="fort." // str(i:3)


write(unit=str3, fmt='(I3)') myid+115
if(myid+115.gt.9.and.myid+115.le.99) then
   i3=2
elseif(myid+115.gt.99) then
   i3=1 
else
   i3=3
endif
out="../fort." // str3(i3:3)

write(unit=str2, fmt='(I3)') chk_ind
if(chk_ind.gt.9.and.chk_ind.le.99) then
   i2=2
elseif(chk_ind.gt.99) then
   i2=1 
else
   i2=3
endif
st2=".chk"
st3=".bak"

write(unit=str4, fmt='(I3)') rename
if(rename.gt.9.and.rename.le.99) then
   i4=2
elseif(rename.gt.99) then
   i4=1 
else
   i4=3
endif
checkpoint="./run" // str6(i6:3) //"/"
!checkpoint2="./stored_scf/aceck" // str4(i4:3) // st2
!checkpoint3="./stored_scf/aceck" // str2(i2:3) // st2
checkpoint4="../stored_scf/aceck" // str4(i4:3) // st3
checkpoint5="./stored_scf/aceck" // str2(i2:3) // st3
checkpoint6="./run" // str6(i6:3) //"/OLDAOMOS"
checkpoint7="./run" // str6(i6:3) //"/ZMAT"
checkpoint8="./GENBAS"
checkpoint9="./run" // str6(i6:3) //"/SAVEDIR"
checkpoint10="./run" // str6(i6:3) //"/*"
checkpoint11="./AOBASMOS"
open(unit=9,file=filename)
write(9,'(a11)') '#!/bin/bash'
write(9,*) 'rm -r ', checkpoint9
write(9,*) 'rm ', checkpoint10

if(chk_ind>0)then
   if(chk_ind.ne.rename)then
      write(9,*) 'cp ', checkpoint5, '   ',checkpoint6
   endif
endif
write(9,*) 'mv ',inp, checkpoint7
write(9,*) 'cp ',checkpoint8, '   ', checkpoint
write(9,*) 'cd ', checkpoint
write(9,*) 'xaces2 > ',out
if(chk_ind>0)then   
   if(rename>0) then
      !      if(chk_ind.ne.rename)then
      write(9,*) 'mv ', checkpoint11, '  ',checkpoint4
      !      endif
      !      write(9,*) 'cp ', checkpoint2, '  ', checkpoint4
   else
!      write(9,*) 'rm ', checkpoint6
   endif
!   if(chk_ind.ne.rename)then
!      write(9,*) 'cp ', checkpoint5, '  ',checkpoint3
!   endif
  
endif
close(9)

open(unit=7,file=filename5)
open(unit=myid+15)
write(myid+15,*)'ACES is running for IMLS'
do i=1,natom
   write(myid+15,107) nuclei(i),(q(k),k=i*3-2,i*3)
enddo
write(myid+15,*)' '
20 read(7,fmt='(a80)',end=100)line
write(myid+15,fmt='(a80)')line
goto 20
100 close(7)
if(f2==2) then
   backspace(myid+15)
   write(myid+15,'(a20)')'GRAD_CALC=ANALYTICAL'
   write(myid+15,'(a17)')' '
endif


if(chk_ind>0)then
   if(chk_ind.ne.rename) then
      backspace(myid+15)
      write(myid+15,'(a17)')'GUESS=READ_AO_MOS'
      write(myid+15,'(a17)')' '
   endif
endif

close(myid+15)

107 format(a2,1x,3(f17.14,1x))

call system(filename)
return

end subroutine wr_aces





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wr_gauss(q,nuclei,f2,f3,natom,myid,chk_ind,rename,inp_num)

implicit none
real*8 :: q(natom*3)
character(len=3) :: nuclei(natom)
integer :: f3,f2,natom,id,i,i2,i3,i4,i5,j,k,myid,chk_ind,rename,inp_num
character :: line*80,str*3,filename*17,st1*2,checkpoint*35,str2*4,str3*3,str4*4,st2*4,st3*4,inp*8,out*8,checkpoint2*26,checkpoint3*26,checkpoint4*26,checkpoint5*26,filename2*15
character :: checkpoint6*31,str5*3,filename5*12,st5*4

write(unit=str5, fmt='(I3)') inp_num
if(inp_num.gt.9.and.inp_num.le.99) then
   i5=2
elseif(inp_num.gt.99) then
   i5=1 
else
   i5=3
endif
st5=".abi"
filename5="gauss" // str5(i5:3) // st5

write(unit=str, fmt='(I3)') myid+15
if(myid+15.gt.9.and.myid+15.le.99) then
   i=2
elseif(myid+15.gt.99) then
   i=1 
else
   i=3
endif
st1=".x"
filename="./gausscript" // str(i:3) // st1
filename2="gausscript" // str(i:3) // st1
inp="fort." // str(i:3)


write(unit=str3, fmt='(I3)') myid+115
if(myid+115.gt.9.and.myid+115.le.99) then
   i3=2
elseif(myid+115.gt.99) then
   i3=1 
else
   i3=3
endif
out="fort." // str3(i3:3)

!!!!!
write(unit=str2, fmt='(I4)') chk_ind
if(chk_ind.gt.9.and.chk_ind.le.99) then
   i2=3
elseif(chk_ind.gt.99.and.chk_ind.le.999) then
   i2=2
elseif(chk_ind.gt.999) then
   i2=1 
else
   i2=4
endif

st2=".chk"
st3=".bak"

write(unit=str4, fmt='(I4)') rename
if(rename.gt.9.and.rename.le.99) then
   i4=3
elseif(rename.gt.99.and.rename.le.999) then
   i4=2
elseif(rename.gt.999) then
   i4=1
else
   i4=4
endif


checkpoint="%Chk=/bigdisk/SCR/dftck" // str2(i2:4) // st2 // "_" // str(i:3)
checkpoint2="/bigdisk/SCR/dftck" // str4(i4:4) // st2
checkpoint3="/bigdisk/SCR/dftck" // str2(i2:4) // st2
checkpoint4="/bigdisk/SCR/dftck" // str4(i4:4) // st3
checkpoint5="/bigdisk/SCR/dftck" // str2(i2:4) // st3
checkpoint6="/bigdisk/SCR/dftck" // str2(i2:4) // st2  // "_" // str(i:3)

! f3=1 scf
! f3=2 guess=read
! f3=3 qc

! f2=0 sp
! f2=1 force

open(unit=9,file=filename)
write(9,'(a11)') '#!/bin/bash'
if(chk_ind>0)then
   if(chk_ind.ne.rename)then
      write(9,*) 'cp ', checkpoint5, '  ',checkpoint6
   endif
endif

write(9,*) 'g09 < ',inp,' > ',out
if(chk_ind>0)then   
   if(rename>0) then
      !      if(chk_ind.ne.rename)then
      write(9,*) 'mv ', checkpoint6, '  ',checkpoint4
      !      endif
      !      write(9,*) 'cp ', checkpoint2, '  ', checkpoint4
   else
      write(9,*) 'rm ', checkpoint6
   endif
!   if(chk_ind.ne.rename)then
!      write(9,*) 'cp ', checkpoint5, '  ',checkpoint3
!   endif
  
endif
close(9)




 open(unit=7,file=filename5)
 open(unit=myid+15)
 if(chk_ind>0)then
    write(myid+15,'(a34)') checkpoint
 endif
 70   read(7,80)line
 id=INDEX(line,'xyz')
 if(id.eq.0) then
    write(myid+15,80)line
    goto 70
 endif

 if(f2.eq.2) then
    backspace(myid+15)
    write(myid+15,'(a42)')'force iop(7/33=1) iop(7/18=1) iop(7/18=10)'
    write(myid+15,'(a17)')' '
 endif
 if(chk_ind>0)then
    if(chk_ind.ne.rename) then
       backspace(myid+15)
       write(myid+15,'(a17)')'guess=read'
       write(myid+15,'(a17)')' '
    endif
 endif
! elseif(f3.eq.3) then
!    backspace(myid+15)
!    write(myid+15,'(a17)')'scf=qc'
!    write(myid+15,'(a17)')' '
! endif

 write(myid+15,80)line

 read(7,80)line
 write(myid+15,80)line

 read(7,80)line
 write(myid+15,80)line


 do i=1,natom
    write(myid+15,107) nuclei(i),(q(k),k=i*3-2,i*3)
 enddo
 write(myid+15,80)' '
 close(7)
 close(myid+15)
 80   format(a80)
 107  format(a2,1x,3(f24.14,1x))

 call system(filename)
 return
end subroutine wr_gauss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rd_aces(q,pot,dq,nuclei,f2,f3,natom,myid,mass)

implicit none
real*8 :: q(natom*3),dq(natom*3),q_rot(natom*3),q_temp(natom*3),mass(natom),bohr
real*8 :: pot,zz(4),z,cart_mat(3,natom),cart_mat2(3,natom),grad_mat(3,natom),quat(4),U_rot(3,3)
character(len=3) :: nuclei(natom)
character(len=80) :: line,st
character :: filename*8,str*3,str6*3,GRD*12
integer :: f1,f2,f3,natom,info,id,i,j,k,l,n,myid,i6,ierr
logical :: var
bohr=0.5291772083d0
write(unit=str6, fmt='(I3)') myid+1
if(myid+1.gt.9.and.myid+1.lt.99) then
   i6=2
elseif(myid+1.gt.99) then
   i6=1 
else
   i6=3
endif

GRD="./run" // str6(i6:3) //"/GRD"


write(unit=str, fmt='(I3)') myid+115
filename= "fort." // str

open(unit=9,file=filename)
  

call getrec('CCSD(T)        =',16,9,line,info)
if(info==0) then
   id=index(line,'=')
   read(line(id+1:id+21),'(f21.15)')pot
else
   rewind(9)
   call getrec('E(SCF)=',7,9,line,info)
   if(info==0) then
      id=index(line,'=')
      read(line(id+1:id+21),'(f21.15)')pot
   else
      pot=2d2
      dq=0d0
      do i=1,3*natom
         write(925+myid,*) q(i)
      enddo
      return
   endif
endif
close(9)
if(f2==2) then
   inquire(file=GRD,exist=var)
   if(var) then
      open(unit=10,file=GRD)
      read(10,'(a80)')line
      read(line,*) i,z
      
      
      do j=1,natom
         k=(j-1)*3
         read(10,'(a80)')line
         read(line,*) z,q_rot(k+1),q_rot(k+2),q_rot(k+3)
      enddo
      do j=1,natom
         k=(j-1)*3
         read(10,'(a80)')line
         read(line,*) z,dq(k+1),dq(k+2),dq(k+3)
      enddo
      q_rot=q_rot*bohr
      
      q_temp=q
      call rm_cmass(q_temp,mass)
      call rm_cmass(q_rot,mass)
      call vec_to_mat(q_temp,cart_mat,natom)
      call vec_to_mat(q_rot,cart_mat2,natom)
      call vec_to_mat(dq,grad_mat,natom)
      call qtrfit(natom,cart_mat2,cart_mat,mass,quat,U_rot,ierr)
      call rotmol(natom,grad_mat,grad_mat,U_rot)
      call mat_to_vec(grad_mat,dq,natom)
      close(10)
   else
      pot=2d2
      dq=0d0
      do i=1,3*natom
         write(925+myid,*) q(i)
      enddo
      return
   endif
endif
!!$ call getrec('Translational invariance is used',32,9,line,info) 
!!$ read(9,'(a80)')line
!!$ do i=1,natom
!!$ read(9,'(a80)')line
!!$ read(line(10:75),*)dq((i-1)*3+1:i*3)
!!$ enddo
return
end subroutine rd_aces

!/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\
!c  str=string to find
!c  ns= string dimension
!c  n=open
!c  get string
subroutine getrec(str,ns,n,line,info)
 implicit none
 integer ::  ns,n,info,is,id
 character(len=80) :: str,line
 info = 0
 21   read(n,80,end=100)line
 is=index(line,'!')
 if(is.eq.0) is = 80
 id=index(line(1:is),str(1:ns))
      if(id.eq.0) goto 21
 80   format(a80)
      return
 100  info=1
      return
end subroutine getrec




!c/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\

subroutine vec_to_mat(cart_perms,cart_mat,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_mat(kp,k)=cart_perms((k-1)*3+kp)
   enddo
enddo
return
end subroutine vec_to_mat

subroutine mat_to_vec(cart_mat,cart_perms,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_perms((k-1)*3+kp)=cart_mat(kp,k)
   enddo
enddo
return
end subroutine mat_to_vec

subroutine rm_cmass2(cart_perms,mass,natom)
integer :: k,kp,natom
real*8 :: mass(natom),cart_perms(3*natom),mtot,cmass1(3)
mtot=0d0
do k=1,natom
   mtot=mtot+mass(k)
enddo
cmass1=0d0
do k=1,natom
   do kp=1,3
      cmass1(kp)=cmass1(kp)+cart_perms((k-1)*3+kp)*mass(k)
   enddo
enddo
cmass1=cmass1/mtot

do k=1,natom
   do kp=1,3
      cart_perms((k-1)*3+kp)=cart_perms((k-1)*3+kp)-cmass1(kp)      
   enddo
enddo
return
end subroutine rm_cmass2




SUBROUTINE sobseq(x,init)
! When the optional integer "init" is present, internally initializes a set of MAXBIT 
! direction numbers for each of MAXDIM different Sobol’ sequences. Otherwise returns as 
! the vector x of length N the next values from N of these sequences. (N must not be 
! changed between initializations.)
USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(OUT) :: x
INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
INTEGER(I4B), PARAMETER :: MAXBIT=30,MAXDIM=6
REAL(SP), SAVE :: fac
INTEGER(I4B) :: i,im,ipp,j,k,l
INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE:: iu
INTEGER(I4B), SAVE :: in
INTEGER(I4B), DIMENSION(MAXDIM), SAVE :: ip,ix,mdeg
INTEGER(I4B), DIMENSION(MAXDIM*MAXBIT), SAVE :: iv
DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/

if (present(init)) then! Initialize, don’t return a vector.
  ix=0
  in=0
  if (iv(1) /= 1) RETURN
  fac=1.0_sp/2.0_sp**MAXBIT
  allocate(iu(MAXDIM,MAXBIT))
  iu=reshape(iv,shape(iu))! To allow both 1D and 2D addressing.
  do k=1,MAXDIM
    do j=1,mdeg(k)! Stored values require only normalization.
      iu(k,j)=iu(k,j)*2**(MAXBIT-j)
    end do
    do j=mdeg(k)+1,MAXBIT! Use the recurrence to get other values.
      ipp=ip(k)
      i=iu(k,j-mdeg(k))
      i=ieor(i,i/2**mdeg(k))
	do l=mdeg(k)-1,1,-1
        if (btest(ipp,0)) i=ieor(i,iu(k,j-l))
        ipp=ipp/2
      end do
      iu(k,j)=i
    end do
  end do
  iv=reshape(iu,shape(iv))
  deallocate(iu)
else! Calculate the next vector in the sequence.
  im=in
  do j=1,MAXBIT! Find the rightmost zero bit.
    if (.not. btest(im,0)) exit
    im=im/2
  end do
  if (j > MAXBIT) call nrerror('MAXBIT too small in sobseq')
  im=(j-1)*MAXDIM
  j=min(size(x),MAXDIM)
  ! XOR the appropriate direction number into each component of the vector 
  ! and convert to a floating number.
  ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
  x(1:j)=ix(1:j)*fac
  in=in+1! Increment the counter.
end if

END SUBROUTINE sobseq
