
module fact1
real*8  fact(0:301)
integer nfac
end module


!================================================================================
!Get basis between block and new added site to get the new block
!================================================================================
subroutine Get_basis(block,site,basis)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: block,site
	type(Total_Basis),intent(inout) :: basis

	integer :: i,j,ij,i1, x,y,num
	integer,allocatable :: pre_qn(:,:),res_qn(:,:)
	integer,external :: Get_Unique_QN
        integer  s1, s2, ss1

	!<1>: Get general information
	basis%len=block%len+site%len
	num=block%num*5  !!! at most enlarged by 3, site%num
        if(hubmodel==1)num=block%num*6

	allocate(pre_qn(2,num),res_qn(2,num))
	pre_qn=0
	res_qn=0

        ij=0

	do i=1,block%num
        do j=1, site%num

        s1=abs(block%sub(i)%num_down-site%sub(j)%num_down)
        s2=block%sub(i)%num_down+site%sub(j)%num_down
        do ss1=s1, s2, su
        ij=ij+1
	pre_qn(1,ij)=block%sub(i)%num_up+site%sub(j)%num_up
	pre_qn(2,ij)=ss1

        enddo
        enddo
        enddo

        num=ij

	!<2>: Sort pre_qn in descending order
	call QuickSort_Integer(pre_qn,num)
	basis%num=get_unique_qn(pre_qn,res_qn,num)

	allocate(basis%sub(basis%num))
	do i=1,basis%num
		basis%sub(i)%new_num_up=res_qn(1,i)
		basis%sub(i)%new_num_down=res_qn(2,i)
	end do
	deallocate(pre_qn,res_qn)

	!<3>: Get basis%sub(:)%num
	do i=1,basis%num
		basis%sub(i)%num=0
		do x=1,block%num
		do y=1,site%num
			if((block%sub(x)%num_up+site%sub(y)%num_up)==basis%sub(i)%new_num_up) then
			s1=abs(block%sub(x)%num_down-site%sub(y)%num_down)
			s2=block%sub(x)%num_down+site%sub(y)%num_down
			if(basis%sub(i)%new_num_down.ge.s1.and.basis%sub(i)%new_num_down.le.s2) then
				basis%sub(i)%num=basis%sub(i)%num+1
			endif
			endif
		end do
		101 continue
		end do
	end do

	!<4>: Get all information
	do i=1,basis%num
		allocate(basis%sub(i)%sub(basis%sub(i)%num))
		basis%sub(i)%idx=0
		basis%sub(i)%dim=0
		do x=1,block%num
		do y=1,site%num
			if((block%sub(x)%num_up+site%sub(y)%num_up)==basis%sub(i)%new_num_up) then
			s1=abs(block%sub(x)%num_down-site%sub(y)%num_down)
			s2=block%sub(x)%num_down+site%sub(y)%num_down
			if(basis%sub(i)%new_num_down.ge.s1.and.basis%sub(i)%new_num_down.le.s2) then
				basis%sub(i)%idx=basis%sub(i)%idx+1
				basis%sub(i)%sub(basis%sub(i)%idx)%bl_num_up=block%sub(x)%num_up
				basis%sub(i)%sub(basis%sub(i)%idx)%bl_num_down=block%sub(x)%num_down
				basis%sub(i)%sub(basis%sub(i)%idx)%st_num_up=site%sub(y)%num_up
				basis%sub(i)%sub(basis%sub(i)%idx)%st_num_down=site%sub(y)%num_down
				basis%sub(i)%sub(basis%sub(i)%idx)%sdim=block%sub(x)%sdim*site%sub(y)%sdim
				basis%sub(i)%sub(basis%sub(i)%idx)%spos=basis%sub(i)%dim
                                        ! current dim == spos
				basis%sub(i)%dim=basis%sub(i)%dim+block%sub(x)%sdim*site%sub(y)%sdim
				!!!goto 102
			endif
			endif
		end do
		102 continue
		end do
	end do

	!<5>: Get basis%dim
	basis%dim=0
	do i=1,basis%num
		basis%dim=basis%dim+basis%sub(i)%dim
	end do

end subroutine Get_basis


!==============================================================================
!Get super_basis from sys_bs and env_bs, according totoal quantum number
!==============================================================================
subroutine Get_super_basis(sys_bs,env_bs,super,total_up,total_down)
	use pubdata
	implicit none

	integer,intent(in) :: total_up,total_down
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Super_Basis),intent(inout) :: super

	integer :: i,x,y, s1, s2

	super%sys_len=sys_bs%len
	super%env_len=env_bs%len

	super%num=0
	do x=1,sys_bs%num
		do y=1,env_bs%num
		if(env_bs%sub(y)%new_num_up==(total_up-sys_bs%sub(x)%new_num_up)) then
			s1=abs(env_bs%sub(y)%new_num_down-sys_bs%sub(x)%new_num_down)
			s2=env_bs%sub(y)%new_num_down+sys_bs%sub(x)%new_num_down
			if(total_down.ge.s1.and.total_down.le.s2) then
				super%num=super%num+1

				!!!goto 101
			endif
			endif
		end do
		101 continue
	end do

	allocate(super%sub(super%num))
	super%dim=0
	super%idx=0
	do x=1,sys_bs%num
		do y=1,env_bs%num
			if(env_bs%sub(y)%new_num_up==(total_up-sys_bs%sub(x)%new_num_up)) then

			s1=abs(env_bs%sub(y)%new_num_down-sys_bs%sub(x)%new_num_down)
			s2=env_bs%sub(y)%new_num_down+sys_bs%sub(x)%new_num_down
			if(total_down.ge.s1.and.total_down.le.s2) then
				super%idx=super%idx+1
				super%sub(super%idx)%sys_num_up=sys_bs%sub(x)%new_num_up
				super%sub(super%idx)%sys_num_down=sys_bs%sub(x)%new_num_down
				super%sub(super%idx)%env_num_up=env_bs%sub(y)%new_num_up
				super%sub(super%idx)%env_num_down=env_bs%sub(y)%new_num_down
				super%sub(super%idx)%sys_dim=sys_bs%sub(x)%dim
				super%sub(super%idx)%env_dim=env_bs%sub(y)%dim
				super%dim=super%dim+sys_bs%sub(x)%dim*env_bs%sub(y)%dim
			endif
			endif
		end do
		102 continue
	end do

end subroutine Get_super_basis


!Get unique QNs from in_qn and save to out_qn
integer function Get_Unique_QN(in_qn,out_qn,num)
	implicit none

	integer,intent(in) ::num,in_qn(2,num)
	integer,intent(inout) :: out_qn(2,num)

	logical :: flag
	integer :: i,j,start,sublen

	Get_Unique_QN=0
	flag=.false.
	start=0
	sublen=0
	out_qn=0


	do i=1,num
		do j=start+1,num
			if( (in_qn(1,j)==in_qn(1,i)).and.(in_qn(2,j)==in_qn(2,i)) ) then
				flag=.true.
				sublen=sublen+1
			endif
		end do
		if(flag) then
			start=start+sublen
			Get_Unique_QN=Get_Unique_QN+1
			out_qn(1:2,Get_Unique_QN)=in_qn(1:2,i)

			sublen=0
			flag=.false.
		endif
	end do

	return

end function Get_Unique_QN


!Get index from idst with total length "tlen" with dist
recursive integer function Get_mod_dist(dist,idst,tlen)
	implicit none

	integer,intent(in) :: dist,idst,tlen

	if(dist<=0) then
		get_mod_dist=idst
		return
	else if(dist==1) then
		get_mod_dist=mod(idst,tlen)+1
		return
	else if(dist==2) then
		get_mod_dist=mod(idst+1,tlen)+1
		return
	else
		get_mod_dist=mod(get_mod_dist(dist-1,idst,tlen),tlen)+1
		return
	endif

end function Get_mod_dist


!======================================================
!Sort interger array in descending order
!======================================================
!For integer data
subroutine QuickSort_Integer(arr,num)
	implicit none

	integer,intent(in) ::  num
	integer,intent(inout) :: arr(2,num)

	integer :: i,x,index,min,large(2),temp(2)
	integer :: tau_qn,spin_qn,tmp_tau_qn,tmp_spin_qn

	!Sort the second number in descending order
	do i=1,num
		do x=i+1,num
			if(arr(2,i)<arr(2,x)) then
				temp(1:2)=arr(1:2,i)
				arr(1:2,i)=arr(1:2,x)
				arr(1:2,x)=temp(1:2)
			endif
		end do
	end do

	!Sort the first number in descending order
	do i=1,num
		do x=i+1,num
			if(arr(2,i)==arr(2,x)) then
				if(arr(1,i)<arr(1,x)) then
					temp(1:2)=arr(1:2,i)
					arr(1:2,i)=arr(1:2,x)
					arr(1:2,x)=temp(1:2)
				endif
			endif
		end do
	end do

end subroutine QuickSort_Integer


!For real data in descending order
recursive subroutine QuickSort_Real(arr,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	real(8),intent(inout) :: arr(len)

	integer :: i,j,index
	real(8) :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	!first sequence
	do while(i<=j)
		do while(arr(i)>min)
			i=i+1
		end do
		do while(arr(j)<min)
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			i=i+1
			j=j-1
		endif
	end do

	!recursion
	if(lower<j) then
		call QuickSort_Real(arr,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Real(arr,i,upper,len)
	endif

end subroutine QuickSort_Real


!For real(8) data
recursive subroutine QuickSort_Double(arr,arridx,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	integer,intent(inout) :: arridx(len)
	real(8),intent(inout) :: arr(len)

	integer :: i,j,index,tmpidx
	real(8) :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	!first sequence
	do while(i<=j)
		do while(arr(i)>min)
			i=i+1
		end do
		do while(arr(j)<min)
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			tmpidx=arridx(i)
			arridx(i)=arridx(j)
			arridx(j)=tmpidx

			i=i+1
			j=j-1
		endif
	end do

	!recursion
	if(lower<j) then
		call QuickSort_Double(arr,arridx,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Double(arr,arridx,i,upper,len)
	endif

end subroutine QuickSort_Double


!For complex data
recursive subroutine QuickSort_Complex(arr,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	complex*16,intent(inout) :: arr(len)

	integer :: i,j,index
	complex*16 :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	!first sequence
	do while(i<=j)
		do while(real(arr(i))>real(min))
			i=i+1
		end do
		do while(real(arr(j))<real(min))
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			i=i+1
			j=j-1
		endif
	end do

	!recursion
	if(lower<j) then
		call QuickSort_Complex(arr,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Complex(arr,i,upper,len)
	endif

end subroutine QuickSort_Complex


!=========================================================
!Update block for diagonal operator without truncation
!=========================================================
subroutine update_block_dia(block,basis,new_block)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,k,x,spos,sdim
        real*8 coef1

	!<1>: Get general information of new_block from basis
	new_block%len=basis%len
	new_block%num=basis%num
	new_block%dim=basis%dim
	new_block%up_dif=0
	new_block%down_dif=0
		
	!<2>: Get basis info for new_block
	allocate(new_block%sub(new_block%num))
	do i=1,new_block%num
		new_block%sub(i)%num_up=basis%sub(i)%new_num_up
		new_block%sub(i)%num_down=basis%sub(i)%new_num_down
		new_block%sub(i)%down_dif=0
		new_block%sub(i)%row_dim=basis%sub(i)%dim
		new_block%sub(i)%sdim=basis%sub(i)%dim

		sdim=new_block%sub(i)%sdim
		allocate(new_block%sub(i)%mat(sdim,sdim))
		new_block%sub(i)%mat=0.0d0
	end do

	!<3>: Get new_block%sub(:)%mat
	do k=1,new_block%num
		do x=1,basis%sub(k)%num
			do i=1,block%num
				if(block%sub(i)%num_up==basis%sub(k)%sub(x)%bl_num_up) then
				if(block%sub(i)%num_down==basis%sub(k)%sub(x)%bl_num_down) then
					spos=basis%sub(k)%sub(x)%spos
					sdim=basis%sub(k)%sub(x)%sdim

                coef1=dsqrt((1.0d0+new_block%sub(k)%num_down)/(1.0d0+block%sub(i)%num_down))
					new_block%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&=new_block%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&+block%sub(i)%mat(1:sdim,1:sdim)*coef1
					goto 101
				endif
				endif
			end do
			101 continue
		end do
	end do

end subroutine update_block_dia


!==========================================================
!Update site for diagonal operator without truncation
!==========================================================
subroutine update_site_dia(site,basis,new_site)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	integer :: i,k,x,y,spos,sdim
	double precision,allocatable :: uni_eye(:,:)
        real*8 coef1

	!<1>: Get general information of new_site from basis
	new_site%len=basis%len
	new_site%num=basis%num
	new_site%dim=basis%dim
	new_site%up_dif=0
	new_site%down_dif=0
		
	!<2>: Get basis info for new_site
	allocate(new_site%sub(new_site%num))
	do i=1,new_site%num
		new_site%sub(i)%num_up=basis%sub(i)%new_num_up
		new_site%sub(i)%num_down=basis%sub(i)%new_num_down
		new_site%sub(i)%row_dim=basis%sub(i)%dim
		new_site%sub(i)%sdim=basis%sub(i)%dim

		sdim=new_site%sub(i)%sdim
		allocate(new_site%sub(i)%mat(sdim,sdim))
		new_site%sub(i)%mat=0.0d0
	end do

	!<2>: Get new_site%sub(:)%mat
	do k=1,new_site%num
		do x=1,basis%sub(k)%num
			do i=1,site%num
				if(site%sub(i)%num_up==basis%sub(k)%sub(x)%st_num_up) then
				if(site%sub(i)%num_down==basis%sub(k)%sub(x)%st_num_down) then
					spos=basis%sub(k)%sub(x)%spos
					sdim=basis%sub(k)%sub(x)%sdim

					allocate(uni_eye(sdim,sdim))
					uni_eye=0.0d0
					do y=1,sdim
						uni_eye(y,y)=1.0d0
					end do

                coef1=dsqrt((1.0d0+new_site%sub(k)%num_down)/(1.0d0+site%sub(i)%num_down))
		new_site%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
			&=new_site%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
				&+uni_eye(1:sdim,1:sdim)*site%sub(i)%mat(1,1)*coef1

					deallocate(uni_eye)
					goto 101
				endif
				endif
			end do
			101 continue
		end do
	end do

end subroutine update_site_dia


!================================================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_block_ndia_sys(block,basis,new_block,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,x,y,new_id,block_id
        integer j1,j2,j3,j4,j5,j6
	integer :: lhs_id,lhs_num_up,lhs_num_down0, lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down

	integer :: lhs_dim,rhs_dim,up_dif,down_dif,down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag
        integer st_num_down
        real*8 coef
  real(8),external :: w6js
        

	!<1>: Get basis info for new_block
	new_block%len=basis%len
	new_block%up_dif=block%up_dif
	new_block%down_dif=block%down_dif
	up_dif=block%up_dif
	down_dif1=block%down_dif

	new_block%num=0
	do i=1,basis%num

        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule new sys
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
        if(lhs_num_down.lt.0)go to 101

!! we do option one where we do not revise up_dif, and down_dif,  but matching
!use triangle rule

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_block%num=new_block%num+1
			goto 101
			endif
			endif
		end do
		101 continue
        enddo
	end do

	!<2>: Get new_block%sub(:)%qn and new_block%sub(:)%dim
	allocate(new_block%sub(new_block%num))
	new_id=0
	new_block%dim=0
	do x=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		rhs_dim=basis%sub(x)%dim
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_block%dim=new_block%dim+rhs_dim
				
				new_block%sub(new_id)%num_up=rhs_num_up
				new_block%sub(new_id)%num_down=rhs_num_down
				new_block%sub(new_id)%down_dif=down_dif
				new_block%sub(new_id)%row_dim=lhs_dim
				new_block%sub(new_id)%sdim=rhs_dim
				allocate(new_block%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_block%sub(new_id)%mat=0.0d0
				goto 102
			endif
			end if
		end do
		102 continue
	end do
	end do

	!<3>: Get new_block%sub(:)%mat
	do i=1,new_block%num
		rhs_num_up=new_block%sub(i)%num_up
		rhs_num_down=new_block%sub(i)%num_down
		lhs_dim=new_block%sub(i)%row_dim
		rhs_dim=new_block%sub(i)%sdim
		
		rhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==rhs_num_up) then
			if(basis%sub(x)%new_num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue

		down_dif=new_block%sub(i)%down_dif
		
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		if(rhs_flag.and.lhs_flag) then
			do x=1,basis%sub(rhs_id)%num
                         !!!st_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
                         st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
			rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
			rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
			rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
			rhs_pos=basis%sub(rhs_id)%sub(x)%spos

                do sub_down_dif=-down_dif1, down_dif1,su !! diff for new block
				block_flag=.false.
				do y=1,block%num
					if(block%sub(y)%num_up==rhs_sub_num_up) then
					if(block%sub(y)%num_down==rhs_sub_num_down) then
					if(block%sub(y)%down_dif==sub_down_dif) then
						block_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						block_flag=.true.
						goto 106
					endif
                                        endif
                                        endif
				end do
				106 continue

				lhs_sub_flag=.false.
                        if(block_flag)then
				do y=1,basis%sub(lhs_id)%num
				if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
				if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
				if(basis%sub(lhs_id)%sub(y)%st_num_down==st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.block_flag) then


        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}

        j1=lhs_num_down !(J')
        j2=down_dif1     ! tensor rank 
        j3=rhs_num_down  !J
        j4=rhs_sub_num_down  ! j
        j5=st_num_down     ! s
        j6=lhs_sub_num_down  ! j'
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*dsqrt((1.0d0+j1)*(1.0d0+j3))
        coef=coef*(-1)**((j6+j5+j3+j2)/2)

        if(coef.ne.0.0)then
	new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&=new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&+coef*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)
                        endif
				endif
			end do
                enddo
		endif
        enddo

end subroutine update_block_ndia_sys

!================================================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_block_ndia_env(block,basis,new_block,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,x,y,new_id,block_id
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down
        integer j1,j2,j3,j4,j5,j6

	real(8) :: coef,signs
  real(8),external :: w6js
	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1,sub_down_dif
	integer :: st_num_up,st_num_down,st_num
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag

	!<1>: Get basis info for new_block
	new_block%len=basis%len
	new_block%up_dif=block%up_dif
	new_block%down_dif=block%down_dif
	up_dif=block%up_dif
	down_dif1=block%down_dif

	new_block%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_block%num=new_block%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_block%sub(:)%qn and new_block%sub(:)%dim
	allocate(new_block%sub(new_block%num))
	new_id=0
	new_block%dim=0
	do x=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		rhs_dim=basis%sub(x)%dim
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_block%dim=new_block%dim+rhs_dim
				
				new_block%sub(new_id)%num_up=rhs_num_up
				new_block%sub(new_id)%num_down=rhs_num_down
				new_block%sub(new_id)%down_dif=down_dif
				new_block%sub(new_id)%row_dim=lhs_dim
				new_block%sub(new_id)%sdim=rhs_dim
				allocate(new_block%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_block%sub(new_id)%mat=0.0d0
				goto 102
			endif
			end if
		end do
		102 continue
	end do
	end do

	!<3>: Get new_block%sub(:)%mat
	do i=1,new_block%num
		rhs_num_up=new_block%sub(i)%num_up
		rhs_num_down=new_block%sub(i)%num_down
		lhs_dim=new_block%sub(i)%row_dim
		rhs_dim=new_block%sub(i)%sdim
		
		rhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==rhs_num_up) then
			if(basis%sub(x)%new_num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue
		
         down_dif=new_block%sub(i)%down_dif  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		if(rhs_flag.and.lhs_flag) then
			do x=1,basis%sub(rhs_id)%num
				st_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down

				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su

				block_flag=.false.
				do y=1,block%num
					if(block%sub(y)%num_up==rhs_sub_num_up) then
					if(block%sub(y)%num_down==rhs_sub_num_down) then
					if(block%sub(y)%down_dif==sub_down_dif) then
						block_flag=.true.
						block_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
		lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
					endif
                                        endif
				end do
				106 continue

                        
			
        	lhs_sub_flag=.false.
                if(block_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue

                        endif



				if(lhs_sub_flag.and.block_flag) then

					Signs=0.0d0
					if(FlagSign=='B') then !For Boson
						Signs=1.0d0
					else if(FlagSign=='F') then !For Fermion
						st_num=st_num_up !!+st_num_down
						if(mod(st_num,2)==0) then
							Signs=1.0d0
						else
							Signs=-1.0d0
						endif
					endif

        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        ! j1=left tot J', j2=tensor rank, J3=righ tot
        ! j4=right block, j5=left site, j6=left block

        j1=lhs_num_down !(J')
        j2=down_dif1     ! tensor rank 
        j3=rhs_num_down  !J
        j4=rhs_sub_num_down  ! j
        j5=st_num_down     ! s
        j6=lhs_sub_num_down  ! j'
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*(-1)**((j6+j5+j3+j2)/2)*dsqrt((1.0d0+j1)*(1.0d0+j3))

                if(coef.ne.0.0)Then
					coef=coef*Signs
	new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				       &=new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)
                        endif
				endif
			end do
			end do
		endif
	end do

end subroutine update_block_ndia_env




!================================================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_site_ndia_sys(site,basis,new_site,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	real(8) :: coef,signs
  real(8),external :: w6js
        integer j1,j2,j3,j4,j5,j6
	integer :: i,x,y,new_id,site_id
	integer :: bl_num_up,bl_num_down,bl_num
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down

	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,site_flag
	real(8),allocatable :: uni_eye(:,:)

	!<1>: Get basis information for new_site
	new_site%len=basis%len
	new_site%up_dif=site%up_dif
	new_site%down_dif=site%down_dif
	up_dif=site%up_dif
	down_dif1=site%down_dif

	new_site%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_site%num=new_site%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_site%sub(:)%new_qn,%dim
	allocate(new_site%sub(new_site%num))
	new_id=0
	new_site%dim=0
	do x=1,basis%num
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		rhs_dim=basis%sub(x)%dim

        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_site%dim=new_site%dim+rhs_dim

				new_site%sub(new_id)%num_up=rhs_num_up
				new_site%sub(new_id)%num_down=rhs_num_down
				new_site%sub(new_id)%down_dif=down_dif
				new_site%sub(new_id)%row_dim=lhs_dim
				new_site%sub(new_id)%sdim=rhs_dim
				allocate(new_site%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_site%sub(new_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get new_site%sub(:)%mat
	do i=1,new_site%num
		rhs_num_up=new_site%sub(i)%num_up
		rhs_num_down=new_site%sub(i)%num_down
		rhs_dim=new_site%sub(i)%sdim

		rhs_flag=.false.

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==rhs_num_up) then
			if(basis%sub(x)%new_num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue
		
         down_dif=new_site%sub(i)%down_dif   !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		if(rhs_flag.and.lhs_flag) then
			do x=1,basis%sub(rhs_id)%num
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
				bl_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				bl_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su
				site_flag=.false.
				do y=1,site%num
					if(site%sub(y)%num_up==rhs_sub_num_up) then
					if(site%sub(y)%num_down==rhs_sub_num_down) then
					if(site%sub(y)%down_dif==sub_down_dif) then
						site_flag=.true.
						site_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
                                        endif
					endif
				end do
				106 continue

				lhs_sub_flag=.false.
                if(site_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%st_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==bl_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.site_flag) then
					allocate(uni_eye(rhs_sub_dim,rhs_sub_dim))
					uni_eye=0.0d0
					do y=1,rhs_sub_dim
						uni_eye(y,y)=1.0d0
					end do

					Signs=0.0d0
					if(FlagSign=='B') then !For Boson
						Signs=1.0d0
					else if(FlagSign=='F') then !For Fermion
						bl_num=bl_num_up !!+bl_num_down
						if(mod(bl_num,2)==0) then
							Signs=1.0d0
						else
							Signs=-1.0d0
						endif
					endif

        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        ! j1=left tot J', j2=tensor rank, J3=righ tot
        ! j4=right block, j5=left site, j6=left block

        j1=lhs_num_down !(J')
        j2=down_dif1     ! tensor rank 
        j3=rhs_num_down  !J
        j4=rhs_sub_num_down  ! s
        j5=bl_num_down     ! j
        j6=lhs_sub_num_down  ! s'
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*(-1)**((j5+j4+j1+j2)/2)*dsqrt((1.0d0+j1)*(1.0d0+j3))

                if(coef.ne.0.0)then


					coef=coef*Signs
					new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&=new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*uni_eye(1:rhs_sub_dim,1:rhs_sub_dim)*site%sub(site_id)%mat(1,1)
                endif

					deallocate(uni_eye)
				endif
			end do
			end do
		endif
	end do

       !! call block_to_disk1(new_site, 109)

end subroutine update_site_ndia_sys


!================================================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_site_ndia_env(site,basis,new_site,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	integer :: i,x,y,new_id,site_id
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down
        integer bl_num_down
        integer j1,j2,j3,j4,j5,j6
  real(8),external :: w6js

	real(8) :: coef,signs
	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1, sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,site_flag
	real(8),allocatable :: uni_eye(:,:)

	!<1>: Get basis information for new_site
	new_site%len=basis%len
	new_site%up_dif=site%up_dif
	new_site%down_dif=site%down_dif
	up_dif=site%up_dif
	down_dif1=site%down_dif

	new_site%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_site%num=new_site%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_site%sub(:)%new_qn,%dim
	allocate(new_site%sub(new_site%num))
	new_id=0
	new_site%dim=0
	do x=1,basis%num
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		rhs_dim=basis%sub(x)%dim
		
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_site%dim=new_site%dim+rhs_dim

				new_site%sub(new_id)%num_up=rhs_num_up
				new_site%sub(new_id)%num_down=rhs_num_down
				new_site%sub(new_id)%down_dif=down_dif
				new_site%sub(new_id)%row_dim=lhs_dim
				new_site%sub(new_id)%sdim=rhs_dim
				allocate(new_site%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_site%sub(new_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get new_site%sub(:)%mat
	do i=1,new_site%num
		rhs_num_up=new_site%sub(i)%num_up
		rhs_num_down=new_site%sub(i)%num_down
		rhs_dim=new_site%sub(i)%sdim

		rhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==rhs_num_up) then
			if(basis%sub(x)%new_num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue
		
         down_dif=new_site%sub(i)%down_dif !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		if(rhs_flag.and.lhs_flag) then
			do x=1,basis%sub(rhs_id)%num
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
				bl_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                        do sub_down_dif=-down_dif1, down_dif1, su
				site_flag=.false.
				do y=1,site%num
					if(site%sub(y)%num_up==rhs_sub_num_up) then
					if(site%sub(y)%num_down==rhs_sub_num_down) then
					if(site%sub(y)%down_dif==sub_down_dif) then
						site_flag=.true.
						site_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
					endif
					endif
				end do
				106 continue
				lhs_sub_flag=.false.
                        if(site_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%st_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==bl_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.site_flag) then
					allocate(uni_eye(rhs_sub_dim,rhs_sub_dim))
					uni_eye=0.0d0
					do y=1,rhs_sub_dim
						uni_eye(y,y)=1.0d0
					end do
        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        ! j1=left tot J', j2=tensor rank, J3=righ tot
        ! j4=right block, j5=left site, j6=left block

        j1=lhs_num_down !(J')
        j2=down_dif1     ! tensor rank 
        j3=rhs_num_down  !J
        j4=rhs_sub_num_down  ! s
        j5=bl_num_down     ! j
        j6=lhs_sub_num_down  ! s'
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*dsqrt((1.0d0+j1)*(1.0d0+j3))*(-1)**((j5+j4+j1+j2)/2)
                if(coef.ne.0.0)then
					new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&=new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*uni_eye(1:rhs_sub_dim,1:rhs_sub_dim)*site%sub(site_id)%mat(1,1)

                        endif
					deallocate(uni_eye)
				endif
			end do
			end do
		endif
	end do
end subroutine update_site_ndia_env


!================================================
!Update with truncation for number operator
!================================================
subroutine update_trun_dia(ori,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: ori,trun
	type(Total_Block),intent(inout) :: eff

	integer :: i,x,ori_dim,eff_dim
	integer :: eff_num_up,eff_num_down,trun_id
	double complex,allocatable :: mid(:,:)

	!<1>: Get eff%len,%num and eff%sub(:)%qn,%dim
	eff%len=ori%len
	eff%num=trun%num
	eff%dim=trun%dim
	eff%up_dif=ori%up_dif
	eff%down_dif=ori%down_dif
	allocate(eff%sub(eff%num))
	do i=1,eff%num
		eff%sub(i)%num_up=trun%sub(i)%num_up
		eff%sub(i)%num_down=trun%sub(i)%num_down
		eff%sub(i)%down_dif=0
		eff%sub(i)%row_dim=trun%sub(i)%sdim
		eff%sub(i)%sdim=trun%sub(i)%sdim
		
		eff_dim=eff%sub(i)%sdim
		allocate(eff%sub(i)%mat(eff_dim,eff_dim))
		eff%sub(i)%mat=0.0d0
	end do

	!<2>: Get eff%sub(:)%mat
	do i=1,eff%num
		eff_num_up=eff%sub(i)%num_up
		eff_num_down=eff%sub(i)%num_down
		eff_dim=eff%sub(i)%sdim
		do x=1,ori%num
			if(ori%sub(x)%num_up==eff_num_up) then
			if(ori%sub(x)%num_down==eff_num_down) then
				ori_dim=ori%sub(x)%row_dim
				allocate(mid(ori_dim,eff_dim))
 if(realcode)then
                                call DGEMM('N','N',ori_dim,eff_dim,ori_dim,1.0d0,ori%sub(x)%mat&
                                                &,ori_dim,trun%sub(i)%mat,ori_dim,0.0d0,mid,ori_dim)
                                call DGEMM('T','N',eff_dim,eff_dim,ori_dim,1.0d0,trun%sub(i)%mat&
                                                &,ori_dim,mid,ori_dim,0.0d0,eff%sub(i)%mat,eff_dim)
                                else
   call ZGEMM('N','N',ori_dim,eff_dim,ori_dim,cone,ori%sub(x)%mat&
                                                &,ori_dim,trun%sub(i)%mat,ori_dim,czero,mid,ori_dim)
                                call ZGEMM('C','N',eff_dim,eff_dim,ori_dim,cone,trun%sub(i)%mat&
                                                &,ori_dim,mid,ori_dim,czero,eff%sub(i)%mat,eff_dim)
                                endif



				deallocate(mid)
				goto 101
			endif
			endif
		end do
		101 continue
	end do

end subroutine update_trun_dia


!==========================================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!==========================================================
subroutine update_trun_ndia(ori,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: ori,trun
	type(Total_Block),intent(inout) :: eff

	integer :: i,j,k,x,y
	logical :: ori_flag,lhs_flag,rhs_flag
	integer :: rhs_dim,lhs_dim,up_dif,down_dif, down_dif1
	integer :: rhs_id,rhs_num_up,rhs_num_down
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: ori_id,eff_id,ori_row,ori_col
	double complex,allocatable :: mid(:,:)

        
	!<1>: Get eff%len and eff%num
	eff%len=ori%len
	eff%up_dif=ori%up_dif
	eff%down_dif=ori%down_dif
	up_dif=ori%up_dif
	down_dif=ori%down_dif
	down_dif1=ori%down_dif

	eff%num=0
	do x=1,trun%num
		rhs_num_up=trun%sub(x)%num_up
		rhs_num_down=trun%sub(x)%num_down
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do y=1,trun%num
			if(trun%sub(y)%num_up==lhs_num_up) then
			if(trun%sub(y)%num_down==lhs_num_down) then
				eff%num=eff%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get eff%sub(:)%qn,%dim
	allocate(eff%sub(eff%num))
	eff_id=0
	eff%dim=0
	do x=1,trun%num
		rhs_num_up=trun%sub(x)%num_up
		rhs_num_down=trun%sub(x)%num_down
		rhs_dim=trun%sub(x)%sdim

        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,trun%num
			if(trun%sub(y)%num_up==lhs_num_up) then
			if(trun%sub(y)%num_down==lhs_num_down) then
				eff_id=eff_id+1
				lhs_dim=trun%sub(y)%sdim
				eff%dim=eff%dim+rhs_dim

				eff%sub(eff_id)%num_up=rhs_num_up
				eff%sub(eff_id)%num_down=rhs_num_down
				eff%sub(eff_id)%down_dif=down_dif
				eff%sub(eff_id)%row_dim=lhs_dim
				eff%sub(eff_id)%sdim=rhs_dim
				allocate(eff%sub(eff_id)%mat(lhs_dim,rhs_dim))
				eff%sub(eff_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get eff%sub(:)%mat
	do i=1,eff%num
		rhs_num_up=eff%sub(i)%num_up
		rhs_num_down=eff%sub(i)%num_down

		rhs_flag=.false.
		do x=1,trun%num
			if(trun%sub(x)%num_up==rhs_num_up) then
			if(trun%sub(x)%num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue
         down_dif=eff%sub(i)%down_dif

		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,trun%num
			if(trun%sub(x)%num_up==lhs_num_up) then
			if(trun%sub(x)%num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		ori_flag=.false.
		do x=1,ori%num
			if(ori%sub(x)%num_up==rhs_num_up) then
			if(ori%sub(x)%num_down==rhs_num_down) then
			if(ori%sub(x)%down_dif==down_dif) then
				ori_flag=.true.
				ori_id=x
				goto 105
			endif
			endif
			endif
		end do
		105 continue

		if(rhs_flag.and.lhs_flag.and.ori_flag) then
			ori_row=ori%sub(ori_id)%row_dim
			ori_col=ori%sub(ori_id)%sdim
			rhs_dim=trun%sub(rhs_id)%sdim
			lhs_dim=trun%sub(lhs_id)%sdim
			allocate(mid(ori_row,rhs_dim))
           if(realcode)then
                        call DGEMM('N','N',ori_row,rhs_dim,ori_col,1.0d0,ori%sub(ori_id)%mat&
                                        &,ori_row,trun%sub(rhs_id)%mat,ori_col,0.0d0,mid,ori_row)
                        call DGEMM('T','N',lhs_dim,rhs_dim,ori_row,1.0d0,trun%sub(lhs_id)%mat&
                                        &,ori_row,mid,ori_row,0.0d0,eff%sub(i)%mat,lhs_dim)
                                else

 call ZGEMM('N','N',ori_row,rhs_dim,ori_col,cone,ori%sub(ori_id)%mat&
                                        &,ori_row,trun%sub(rhs_id)%mat,ori_col,czero,mid,ori_row)
                        call ZGEMM('C','N',lhs_dim,rhs_dim,ori_row,cone,trun%sub(lhs_id)%mat&
                                        &,ori_row,mid,ori_row,czero,eff%sub(i)%mat,lhs_dim)


			deallocate(mid)
		endif
                endif
	end do

end subroutine update_trun_ndia


!main

SUBROUTINE FACTRL
use fact1
NFAC=160
FACT(0)=1.0
DO 10 I=1,NFAC
10 FACT(I)=FACT(I-1)*FLOAT(I)
RETURN
END





real*8 FUNCTION W3JS(J1,J2,J3,M1,M2,M3)
use fact1
INTEGER Z,ZMIN,ZMAX
real*8 denom, cc, cc1, cc2
W3JS=0.0d0
!w3js=1
!go to 1000
IF(M1+M2+M3.NE.0) GOTO 1000
IA=J1+J2
IF(J3.GT.IA) GOTO 1000
IB=J1-J2
IF(J3.LT.IABS(IB)) GOTO 1000
JSUM=J3+IA
IC=J1-M1
ID=J2-M2
IF(MOD(JSUM,2).NE.0) GOTO 1000
IF(MOD(IC,2).NE.0) GOTO 1000
IF(MOD(ID,2).NE.0) GOTO 1000
IF(IABS(M1).GT.J1) GOTO 1000
IF(IABS(M2).GT.J2) GOTO 1000
IF(IABS(M3).GT.J3) GOTO 1000
IE=J3-J2+M1
IF=J3-J1-M2
ZMIN=MAX0(0,-IE,-IF)
IG=IA-J3
IH=J2+M2
ZMAX=MIN0(IG,IH,IC)
CC=0.0
DO 200 Z=ZMIN,ZMAX,2
DENOM=FACT(Z/2)*FACT((IG-Z)/2)*FACT((IC-Z)/2)*FACT((IH-Z)/2)*FACT((IE+Z)/2)*FACT((IF+Z)/2)
IF(MOD(Z,4).NE.0) DENOM=-DENOM
CC=CC+1.0/DENOM
200 CONTINUE
CC1=FACT(IG/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2)/FACT((JSUM+2)/2)
CC2=FACT((J1+M1)/2)*FACT(IC/2)*FACT(IH/2)*FACT(ID/2)*FACT((J3-M3)/2)*FACT((J3+M3)/2)
CC=CC*SQRT(CC1*CC2)
IF(MOD(IB-M3,4).NE.0) CC=-CC
W3JS=CC
1000 RETURN
END


real*8 FUNCTION W6JS(J1,J2,J3,L1,L2,L3)
use fact1
INTEGER W,WMIN,WMAX
INTEGER SUM1,SUM2,SUM3,SUM4
real*8 denom, omega, theta, theta1, theta2, theta3, theta4
W6JS=0.0
!w6js=1
!go to 1000

IA=J1+J2
IF(IA.LT.J3) GOTO 1000
IB=J1-J2
IF(IABS(IB).GT.J3) GOTO 1000
IC=J1+L2
IF(IC.LT.L3) GOTO 1000
ID=J1-L2
IF(IABS(ID).GT.L3) GOTO 1000
IE=L1+J2
IF(IE.LT.L3) GOTO 1000
IF=L1-J2
IF(IABS(IF).GT.L3) GOTO 1000
IG=L1+L2
IF(IG.LT.J3) GOTO 1000
IH=L1-L2
IF(IABS(IH).GT.J3) GOTO 1000
SUM1=IA+J3
SUM2=IC+L3
SUM3=IE+L3
SUM4=IG+J3
IF(MOD(SUM1,2).NE.0) GOTO 1000
IF(MOD(SUM2,2).NE.0) GOTO 1000
IF(MOD(SUM3,2).NE.0) GOTO 1000
IF(MOD(SUM4,2).NE.0) GOTO 1000
WMIN=MAX0(SUM1,SUM2,SUM3,SUM4)
II=IA+IG
IJ=J2+J3+L2+L3
IK=J3+J1+L3+L1
WMAX=MIN0(II,IJ,IK)
OMEGA=0.0
DO 200 W=WMIN,WMAX,2
DENOM=FACT((W-SUM1)/2)*FACT((W-SUM2)/2)*FACT((W-SUM3)/2)*FACT((W-SUM4)/2)*FACT((II-W)/2)*FACT((IJ-W)/2)*FACT((IK-W)/2)
IF(MOD(W,4).NE.0) DENOM=-DENOM
OMEGA=OMEGA+FACT(W/2+1)/DENOM
200 CONTINUE
THETA1=FACT((IA-J3)/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2)/FACT(SUM1/2+1)
THETA2=FACT((IC-L3)/2)*FACT((L3+ID)/2)*FACT((L3-ID)/2)/FACT(SUM2/2+1)
THETA3=FACT((IE-L3)/2)*FACT((L3+IF)/2)*FACT((L3-IF)/2)/FACT(SUM3/2+1)
THETA4=FACT((IG-J3)/2)*FACT((J3+IH)/2)*FACT((J3-IH)/2)/FACT(SUM4/2+1)
THETA=THETA1*THETA2*THETA3*THETA4
W6JS=OMEGA*SQRT(THETA)
1000 RETURN
END
real*8 FUNCTION W9JS(J1,J2,J3,J4,J5,J6,J7,J8,J9)
 real(8),external :: w6js
real*8 x, x1, x2, x3, s
!w9js=1
!go to 11
KMIN=ABS(J1-J9)
KMAX=J1+J9
I=ABS(J4-J8)
IF(I.GT.KMIN) KMIN=I
I=J4+J8
IF(I.LT.KMAX) KMAX=I
I=ABS(J2-J6)
IF(I.GT.KMIN) KMIN=I
I=J2+J6
IF(I.LT.KMAX) KMAX=I
X=0.
DO 1 K=KMIN,KMAX,2
S=1.
IF(MOD(K,2).NE.0) S=-1.
X1=W6JS(J1,J9,K,J8,J4,J7)
X2=W6JS(J2,J6,K,J4,J8,J5)
X3=W6JS(J1,J9,K,J6,J2,J3)
X=X+S*X1*X2*X3*FLOAT(K+1)
1 CONTINUE
W9JS=X
11      continue
RETURN
END








