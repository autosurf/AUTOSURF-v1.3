!*********************************  A U T O S U R F  *******************************
!===================================================================================
!-----------------------------------------------------------------------------------
!-                                                                                 -
!-        AUTOSURF Package: A set of programs for the automated construction       -
!-              of Potential Energy Surfaces for van der Waals systems              -
!-                                                                                 -
!-----------------------------------------------------------------------------------
!===================================================================================
!***********************************************************************************
!-        Set of Fortran90 subroutines for "AUTOSURF-PES_rigid4D" PROGRAM          -
!***********************************************************************************


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      A B _ P A R A L L E L
! ----------------------------------------------------------------------------------
!***********************************************************************************
! This set of subroutines create the processor-dependent input files to 
! compute the ab initio energies using the information from the corresponding 
! sample-header file.


! ----------------------------------------------------------------------------------
!      a b _ i n i t i o
! ----------------------------------------------------------------------------------
! explain...
!
! *** Input ***
! q         <-- Cartesian coordinates for all atoms in the system
! nuclei    <-- Element labels for all the atom in the system
! mass      <-- Masses of all the atom in the system
! f1        <-- Code to be used: 1=Gaussian, 2=Molpro, 3= ACES II, 4=Molpro*
! f2        <-- Type of calculation: 1=single point energies, 2= also gradients
! natom     <-- Total number of atoms in the system
! myid      <-- process id in the MPI environment
! inp_num   <--  label defining the sample-header file: 1=molpro1.abi, 
!                                                       2=molpro2.abi, 3=molpro10.abi

! *** Output ***
! pot       <-- computed electronic energy
! dq        <-- computed gradients

subroutine  ab_initio(q,nuclei,mass,pot,dq,f1,f2,natom,myid,inp_num,midbond,R)

 implicit none
 real*8 :: q(natom*3),dq(natom*3),mass(natom)
 real*8 :: pot,R
 character(len=3) :: nuclei(natom)
 integer :: f1,f2,natom,myid,inp_num,midbond

 if(f1.eq.1) then
   call wr_gauss(q,nuclei,f2,natom,myid,inp_num) 
   call rd_gauss(q,pot,dq,nuclei,mass,f2,natom,myid)
 !elseif(f1.eq.2) then
 !  call wr_molpro(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num) 
 !  call rd_molpro(q,pot,dq,nuclei,f2,f3,natom,myid,mass)
 !elseif(f1.eq.3) then
 !  call wr_aces(q,nuclei,f2,natom,myid,chk_ind,rename,inp_num) 
 !  call rd_aces(q,pot,dq,nuclei,f2,f3,natom,myid,mass)
 elseif(f1.gt.3) then
   call wr_molpro(f1,q,nuclei,f2,natom,myid,inp_num,midbond,R) 
   call rd_molpro(q,pot,dq,nuclei,f2,natom,myid,mass)
 endif 

end subroutine  ab_initio

! ----------------------------------------------------------------------------------
!      w r _ m o l p r o
! ----------------------------------------------------------------------------------
! This subroutine creates the processor-dependent Molpro input file ("fort#") using 
! the coordinates "q" and information from the sample-header file "molpro#.abi"
! It also creates the executable processor-dependent script "Molscript#.x", which
! runs the Molpro-job and erases itself once it finishes.
!
! *** Input ***
! q         <-- Cartesian coordinates for all atoms in the system
! nuclei    <-- Element labels for all the atom in the system
! mass      <-- Masses of all the atom in the system
! f2        <-- Type of calculation: 1=single point energies, 2= also gradients
! natom     <-- Total number of atoms in the system
! myid      <-- process id in the MPI environment
! inp_num   <--  label defining the sample-header file: 1=molpro1.abi, 
!                                                       2=molpro2.abi, 3=molpro10.abi

subroutine wr_molpro(f1,q,nuclei,f2,natom,myid,inp_num,midbond,R)

 implicit none
 real*8 :: q(natom*3),R,dummyatom(3)
 character(len=3) :: nuclei(natom)
 character :: line*80,str*5,filename*18,script*16,filename5*13,inp*10,out*9,str5*3
 character :: exec_script*35,chmod777script*26,punch*18
 integer :: f1,f2,natom,i,k,id,myid,inp_num,i3,i2,i4,midbond

 80   format(a80)! to read/write a complete line
 107  format(a2,1x,3(f24.14,1x))! to write the Cartesian coordinates into the Molpro input file
 108  format(a3,1x,3(f24.14,1x))
 ! processor-dependent script
 write(str,'(I5)')myid+10000! #
 filename="./Molscript"//str(1:5)//".x"
 inp="fort"//str(1:5)! processor-dependent input file "fort#"
 open(unit=9,file=filename)! processor-dependent script "Molscript#.x"
 write(9,'(a11)')'#!/bin/bash'
 if(f1.eq.4) then
   write(9,*) 'molpros ',inp,'-s'
 elseif(f1.eq.5) then
   write(9,*) 'molpro ',inp,'-s'
 endif
 write(9,*) 'rm ', filename
 close(9)

 script = "Molscript"//str(1:5)//".x"
 chmod777script = "chmod 777 "//script
 exec_script = "./Molscript"//str(1:5)//".x"
 punch = "punch,pun"//str(1:5)//",new"

 write(str5,'(I3)')inp_num
 filename5='molpro'//trim(adjustl(str5))//'.abi'
 open(unit=7,file=filename5)! input header-sample
 open(unit=myid+7000,file=inp)! processor-dependent file "fort#" (input file for Molpro)
 70 read(7,80)line
 id=INDEX(line,'xyz')
 if(id.eq.0) then
   write(myid+7000,80)line
   goto 70
 endif
 write(myid+7000,80)line
 write(myid+7000,*) punch
 write(myid+7000,*)'noorient'
 write(myid+7000,*)'geometry={'
 if(midbond.gt.0) then
   write(myid+7000,*)natom+1
   dummyatom(:)=0d0
   dummyatom(3)=R/2d0
 else
   write(myid+7000,*)natom
 endif
 write(myid+7000,*)'MOLPRO is running for AUTOSURF'
 do i=1,natom
   write(myid+7000,107) nuclei(i),(q(k),k=i*3-2,i*3)
 enddo
 if(midbond.gt.0)write(myid+7000,108)'H99',dummyatom(1:3)
 write(myid+7000,*)'}'
 if(midbond.gt.0)write(myid+7000,*)'dummy,H99'
 90 read(7,80)line
 id=INDEX(line,'---')
 if(id.eq.0) then
    write(myid+7000,80)line
    goto 90
 endif
 if(f2.eq.2) write(myid+7000,*)'forces'
 if(f2.eq.2) write(myid+7000,*)'VARSAV'
 if(f2.eq.2) write(myid+7000,*) 'show,gradx'
 if(f2.eq.2) write(myid+7000,*) 'show,grady'
 if(f2.eq.2) write(myid+7000,*) 'show,gradz'
 100   read(7,80)line
 id=INDEX(line,'----')
 if(id.eq.0) then
    write(myid+7000,80)line
    goto 100
 endif
 write(myid+7000,*) '---'
 close(7)
 close(myid+7000)

 call system(chmod777script)
! call system('sleep10')
 call system(exec_script)
 return

end subroutine wr_molpro

! ----------------------------------------------------------------------------------
!      r d _ m o l p r o
! ----------------------------------------------------------------------------------
! This subroutine reads in the (processor-dependent) output file generated by 
! Molpro ("pun#") and retrieves the calculated energies and gradients (if computed) 
!
! *** Input ***
! q         <-- Cartesian coordinates for all atoms in the system
! nuclei    <-- Element labels for all the atom in the system
! mass      <-- Masses of all the atom in the system
! f2        <-- Type of calculation: 1=single point energies, 2= also gradients
! natom     <-- Total number of atoms in the system
! myid      <-- process id in the MPI environment
! mass      <-- Masses of all the atom in the system
!
! *** Output ***
! pot       <-- computed electronic energy
! dq        <-- computed gradients

subroutine rd_molpro(q,pot,dq,nuclei,f2,natom,myid,mass)

 implicit none
 real*8 :: q(natom*3),dq(natom*3),q_rot(natom*3),mass(natom),bohr
 real*8 :: pot,zz(4),z
 character(len=3) :: nuclei(natom),symbord(natom)
 character(len=80) :: line,st
 integer :: f1,f2,natom,info,id,i,j,k,l,n,myid,skip(natom),ierr,pass
 character :: str*5,punch*8
 logical :: var

 bohr=0.529177249d0 !value used by Molpro

 write(str,'(I5)')myid+10000
 punch="pun"//str(1:5)!output file generated by Molpro

 pass=0
 213 continue
 inquire(file=punch,exist=var)

 if(var)then
   open(unit=9,file=punch)
   call getrec('MOLPRO_ENERGY',13,9,line,info) 
   if(info==1)then
     if(pass==0) then
       call sleep(10)
       pass=1
       close(9)
       goto 213
     endif
     write(9000+myid,*) 'could not converge'
     do i=1,3*natom
       write(9000+myid,*) q(i)
     enddo
     pot=2d2
     dq(:)=0d0
     return
   endif
   id=index(line,'-')
   read(line(id:id+20),'(f20.15)')pot
   close(9)
 else
   if(pass==0) then
     call sleep(10)
     pass=1
     goto 213
   endif
   write(9000+myid,*) 'could not converge'
   do i=1,3*natom
     write(9000+myid,*) q(i)
   enddo
   pot=2d2
   dq(:)=0d0
   return
 endif

 if(f2.eq.2) then  
   open(unit=9,file=punch)
   rewind(9)
   do i=1,natom
     call getrec('GRADIENT',8,9,line,info)
     if(info==1)then
       write(9000+myid,*) 'file corrupted'
       do j=1,3*natom
         write(9000+myid,*) q(j)
       enddo
       pot=2d2
       dq(:)=0d0
       return
     endif
     id=index(line,':')
     read(line(id+1:80),*) dq(3*(i-1)+1:3*(i-1)+3)
   enddo 
 endif

 return

end subroutine rd_molpro



subroutine wr_gauss(q,nuclei,f2,natom,myid,inp_num)

 implicit none
 real*8 :: q(natom*3),R,dummyatom(3)
 character(len=3) :: nuclei(natom)
 character :: line*80,str*5,filename*18,script*16,filename5*13,inp*10,out*13,str5*3
 character :: exec_script*35,chmod777script*26,punch*18
 integer :: f2,natom,i,k,id,myid,inp_num,i3,i2,i4,midbond

 80   format(a80)! to read/write a complete line
 107  format(a2,1x,3(f24.14,1x))! to write the Cartesian coordinates into the Molpro input file
 108  format(a3,1x,3(f24.14,1x))
 ! processor-dependent script
 write(str,'(I5)')myid+10000
 filename="./Gauscript"//str(1:5)//".x"
 out="fort"//str(1:5)//".out"
 inp="fort"//str(1:5)! processor-dependent input file "fort#"
 open(unit=9,file=filename)! processor-dependent script "Gauscript#.x"
 write(9,'(a11)')'#!/bin/bash'
 write(9,*) 'g09 < ',inp,' > ',out
 write(9,*) 'rm ', filename
 close(9)

 script = "Gauscript"//str(1:5)//".x"
 chmod777script = "chmod 777 "//script
 exec_script = "./Gauscript"//str(1:5)//".x"

 write(str5,'(I3)')inp_num
 filename5='gaussi'//trim(adjustl(str5))//'.abi'
 open(unit=7,file=filename5)! input header-sample
 open(unit=myid+7000,file=inp)! processor-dependent file "fort#" (input file for Molpro)
 70 read(7,80)line
 id=INDEX(line,'xyz')
 if(id.eq.0) then
   write(myid+7000,80)line
   goto 70
 endif

! if(f2.eq.2) then
!    backspace(myid+7000)
!    write(myid+7000,'(a42)')'force iop(7/33=1) iop(7/18=1) iop(7/18=10)'
!    write(myid+7000,'(a17)')' '
! endif

 do i=1,natom
    write(myid+7000,107) nuclei(i),(q(k),k=i*3-2,i*3)
 enddo
 write(myid+7000,80)' '
 close(7)
 close(myid+7000)

 call system(chmod777script)
 call system(exec_script)
 return

end subroutine wr_gauss



subroutine rd_gauss(q,pot,dq,nuclei,mass,f2,natom,myid)

  implicit none
  real*8 :: q(natom*3),dq(natom*3),q_rot(natom*3),q_temp(natom*3),mass(natom)
  real*8 :: pot,cart_mat(3,natom),cart_mat2(3,natom),grad_mat(3,natom),grad_matt(3,natom),quat(4),U_rot(3,3)
  character(len=3) :: nuclei(natom)
  character(len=80) :: line,stdm
  character :: filename*13,str*5
  integer :: f1,f2,natom,info,id,i,j,k,l,flag,myid,ierr
  write(unit=str, fmt='(I5)') myid+10000
      filename= "fort"//str//".out"

  open(unit=9,file=filename)

  rewind(9)
  flag=0
  info=0
!12 call getrec(' ITN=',5,9,line,info)
!  if(info==1.and.flag==0) then
!     goto 15
!  endif
!  if(info==0) then
!     id=index(line,'=')
!     read(line(id+17:id+35),'(f18.10)')pot
!     flag=1
!     goto 12
!  endif
!  if(flag==1) then
!     rewind(9)
!     goto 17
!  endif



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
     write(9000+myid,*) 'could not converge'
     do i=1,3*natom
        write(9000+myid,*) q(i)
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
!!$     call rotmol(natom,grad_mat,grad_matt,U_rot)
!!$     call mat_to_vec(grad_matt,dq,natom)
  endif

40 continue
  return
end subroutine rd_gauss





!***********************************************************************************
!***********************************************************************************

!/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\-/-\
!c  str=string to find
!c  ns= string dimension
!c  n=open
!c  get string

subroutine getrec(str,ns,n,line,info)
 implicit none
 integer ::  ns,n,info,is,id
 character(len=80) :: line
 character(len=ns) :: str
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
