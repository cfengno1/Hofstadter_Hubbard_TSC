subroutine block_transfer_basis(basis,outmat)
        use pubdata
        implicit none

        type(Total_Basis),intent(in) :: basis
        type(Total_Block),intent(inout) :: outmat

        integer :: i,row_dim,col_dim

        outmat%len=basis%len
        outmat%num=basis%num
        outmat%dim=basis%dim
        outmat%up_dif=0!%up_dif
        outmat%down_dif=0!inmat%down_dif

        if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
        do i=1,outmat%num
                outmat%sub(i)%num_up=basis%sub(i)%new_num_up
                outmat%sub(i)%num_down=basis%sub(i)%new_num_down
                outmat%sub(i)%row_dim=basis%sub(i)%dim
                outmat%sub(i)%sdim=basis%sub(i)%dim

        if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
                outmat%sub(i)%mat=0.d0
        end do

end subroutine block_transfer_basis


subroutine block_transfer_trans(inmat,outmat)
        use pubdata
        implicit none

        type(Total_Block),intent(in) :: inmat
        type(Total_Block),intent(inout) :: outmat

        integer :: i,row_dim,col_dim, up_dif, down_dif

        outmat%len=inmat%len
        outmat%num=inmat%num
        outmat%dim=0!!!inmat%dim
        outmat%up_dif=-inmat%up_dif
        outmat%down_dif=inmat%down_dif
        up_dif=inmat%up_dif
        down_dif=inmat%down_dif

        if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
        do i=1,outmat%num
                outmat%sub(i)%num_up=inmat%sub(i)%num_up+up_dif
                outmat%sub(i)%num_down=inmat%sub(i)%num_down+inmat%sub(i)%down_dif
                outmat%sub(i)%down_dif=-inmat%sub(i)%down_dif
                outmat%sub(i)%row_dim=inmat%sub(i)%sdim
                outmat%sub(i)%sdim=inmat%sub(i)%row_dim
                outmat%dim=outmat%dim+outmat%sub(i)%sdim
        if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
                outmat%sub(i)%mat=transpose(inmat%sub(i)%mat)
        if(.not.realcode) outmat%sub(i)%mat=dconjg(outmat%sub(i)%mat)
        end do
end subroutine block_transfer_trans

subroutine block_pass_info_trans(inmat,outmat)
        use pubdata
        implicit none

        type(Total_Block),intent(in) :: inmat
        type(Total_Block),intent(inout) :: outmat

        integer :: i,row_dim,col_dim, up_dif, down_dif

        outmat%len=inmat%len
        outmat%num=inmat%num
        outmat%dim=0!!!inmat%dim
        outmat%up_dif=-inmat%up_dif
        outmat%down_dif=-inmat%down_dif
        up_dif=inmat%up_dif
        down_dif=inmat%down_dif

        if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
        do i=1,outmat%num
                outmat%sub(i)%num_up=inmat%sub(i)%num_up+up_dif
                outmat%sub(i)%num_down=inmat%sub(i)%num_down+down_dif
                outmat%sub(i)%down_dif=-inmat%sub(i)%down_dif
                outmat%sub(i)%row_dim=inmat%sub(i)%sdim
                outmat%sub(i)%sdim=inmat%sub(i)%row_dim
                outmat%dim=outmat%dim+outmat%sub(i)%sdim
        if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
                outmat%sub(i)%mat=0.0d0  !!transpose(inmat%sub(i)%mat)
        end do
end subroutine block_pass_info_trans



!Transfer the block by allocating new space
subroutine block_transfer(inmat,outmat)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: inmat
	type(Total_Block),intent(inout) :: outmat

	integer :: i,row_dim,col_dim

	outmat%len=inmat%len
	outmat%num=inmat%num
	outmat%dim=inmat%dim
	outmat%up_dif=inmat%up_dif
	outmat%down_dif=inmat%down_dif

	if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
	do i=1,outmat%num
		outmat%sub(i)%num_up=inmat%sub(i)%num_up
		outmat%sub(i)%num_down=inmat%sub(i)%num_down
		outmat%sub(i)%down_dif=inmat%sub(i)%down_dif
		outmat%sub(i)%row_dim=inmat%sub(i)%row_dim
		outmat%sub(i)%sdim=inmat%sub(i)%sdim

	if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
		outmat%sub(i)%mat=inmat%sub(i)%mat
	end do

end subroutine block_transfer
!Transfer general information from old_block to new_block
subroutine block_pass_info(old_block,new_block)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: old_block
	type(Total_Block),intent(inout) :: new_block

	integer :: i,row_dim,col_dim

	new_block%len=old_block%len
	new_block%num=old_block%num
	new_block%dim=old_block%dim
	new_block%up_dif=old_block%up_dif
	new_block%down_dif=old_block%down_dif
	allocate(new_block%sub(new_block%num))

	do i=1,new_block%num
		new_block%sub(i)%num_up=old_block%sub(i)%num_up
		new_block%sub(i)%num_down=old_block%sub(i)%num_down
		new_block%sub(i)%down_dif=old_block%sub(i)%down_dif
		new_block%sub(i)%row_dim=old_block%sub(i)%row_dim
		new_block%sub(i)%sdim=old_block%sub(i)%sdim

		allocate(new_block%sub(i)%mat(new_block%sub(i)%row_dim,new_block%sub(i)%sdim))
		new_block%sub(i)%mat=0.0d0
	end do

end subroutine block_pass_info


!Set block data to zero
subroutine block_zero(block)
	use pubdata
	implicit none

	type(Total_Block),intent(inout) :: block

	integer :: i

	do i=1,block%num
		block%sub(i)%mat=0.0d0
	end do

end subroutine block_zero

subroutine block_add_block_two(block1,coef1,new_block)
        use pubdata
        implicit none

        double complex,intent(in) :: coef1!,coef2,coef
        type(Total_block),intent(inout) :: block1 !!!!,block2
        type(Total_block),intent(inout) :: new_block

        logical :: flag1,flag2
        integer :: i,x,num_up,num_down,idx1,idx2, x1,down_dif

!        if(not.allocated)

                x1=1
        do i=1,new_block%num
                num_up=new_block%sub(i)%num_up
                num_down=new_block%sub(i)%num_down
                down_dif=new_block%sub(i)%down_dif
                flag1=.false.
                do x=x1,block1%num
                        if(block1%sub(x)%num_up==num_up) then
                        if(block1%sub(x)%num_down==num_down) then
                        if(block1%sub(x)%down_dif==down_dif) then
                                flag1=.true.
                                idx1=x
                                x1=x+1
                                goto 101
                        endif
                        endif
                        endif
                end do
                101 continue

                if(flag1) then
                                new_block%sub(i)%mat=new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat
                endif
        end do

end subroutine block_add_block_two

!Assume block1 and block2 have the save configuration
!new_block=(block1+block2)
subroutine block_add_block(block1,coef1,block2,coef2,new_block,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef1,coef2,coef
	type(Total_block),intent(inout) :: block1,block2
	type(Total_block),intent(inout) :: new_block

	logical :: flag1,flag2
	integer :: i,x,num_up,num_down,down_dif,idx1,idx2

	do i=1,new_block%num

		num_up=new_block%sub(i)%num_up
		num_down=new_block%sub(i)%num_down
		down_dif=new_block%sub(i)%down_dif

		flag1=.false.
		do x=1,block1%num
			if(block1%sub(x)%num_up==num_up) then
			if(block1%sub(x)%num_down==num_down) then
			if(block1%sub(x)%down_dif==down_dif) then
				flag1=.true.
				idx1=x
				goto 101
			endif
			endif
			endif
		end do
		101 continue


		flag2=.false.
		do x=1,block2%num
			if(block2%sub(x)%num_up==num_up) then
			if(block2%sub(x)%num_down==num_down) then
			if(block2%sub(x)%down_dif==down_dif) then
				flag2=.true.
				idx2=x
				goto 102
			endif
			endif
			endif
		end do
		102 continue

		if(flag1) then
			if(flag2) then
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat+coef2*block2%sub(idx2)%mat
			else
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat
			endif
		else
			if(flag2) then
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef2*block2%sub(idx2)%mat
			else
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat
			endif
		endif
	end do

end subroutine block_add_block


!Transpose the block
subroutine block_transpose(block,new_block)
	use pubdata
	implicit none

	type(Total_block),intent(in) :: block
	type(Total_Block),intent(inout) :: new_block

	integer :: i

	do i=1,block%num
		new_block%sub(i)%mat=transpose(block%sub(i)%mat)
	end do

end subroutine block_transpose


!new_block = coef*old_block
subroutine block_coef(old_block,coef,new_block)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Block),intent(in) :: old_block
	type(Total_Block),intent(inout) :: new_block

	integer :: i

	do i=1,new_block%num
		new_block%sub(i)%mat=coef*old_block%sub(i)%mat
	end do

end subroutine block_coef




!====================================================================
!According to basis information, qn's are sort in descending order
!Both oper1 and oper2 have the same configuration, and are all
!diagonal operators, when acting on the states, it does not change
!the states, as well as the quantum number used to label the state
!====================================================================
!Flag='T' means transopose; Flag='N' means no transpose
subroutine block_mul_block_dia(oper1,Flag1,oper2,Flag2,coef,new_oper)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: oper1,oper2
	type(Total_Block),intent(inout) :: new_oper
	integer :: i,sdim,num_up,num_down

	new_oper%len=oper1%len
	new_oper%num=oper1%num
	new_oper%dim=oper1%dim
	new_oper%up_dif=0
	new_oper%down_dif=0

	allocate(new_oper%sub(new_oper%num))
	do i=1,new_oper%num
		new_oper%sub(i)%num_up=oper1%sub(i)%num_up
		new_oper%sub(i)%num_down=oper1%sub(i)%num_down
		new_oper%sub(i)%down_dif=oper1%sub(i)%down_dif
		new_oper%sub(i)%row_dim=oper1%sub(i)%row_dim
		new_oper%sub(i)%sdim=oper1%sub(i)%sdim

		sdim=new_oper%sub(i)%sdim
		allocate(new_oper%sub(i)%mat(sdim,sdim))

         if(realcode)then
                call DGEMM(Flag1,Flag2,sdim,sdim,sdim,coef,oper1%sub(i)%mat,sdim&
                                 &,oper2%sub(i)%mat,sdim,0.0d0,new_oper%sub(i)%mat,sdim)
                                else
           call ZGEMM(Flag1,Flag2,sdim,sdim,sdim,coef,oper1%sub(i)%mat,sdim&
                                 &,oper2%sub(i)%mat,sdim,czero,new_oper%sub(i)%mat,sdim)
                                endif


	end do

end subroutine block_mul_block_dia
!===========================================================================
!According to basis information, qn's are sort in descending order
!new_oper = (C^+_spin. (C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!Note: new_oper=oper2(2nd,on the left)*oper1(1st,on the right)
!===========================================================================
!new_oper=coef*(oper1^Flag1* oper2^Flag2)
subroutine block_mul_block_ndia(oper1,Flag1,oper2,Flag2,coef,new_oper,basis)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: oper1,oper2
	type(Total_Block),intent(inout) :: new_oper
	type(Total_Basis),intent(in) :: basis

	logical :: Flags,bl_flag1,bl_flag2,bs_flag
	integer :: up_dif1,down_dif1,up_dif2,down_dif2,up_dif,down_dif
	integer :: i,j,x,y,new_dim,mid_dim,old_dim,idx1,idx2,bs_idx,new_id,LDA,LDB,LDC
	integer :: num_up,num_down,mid_num_up,mid_num_down,new_num_up,new_num_down
	integer :: mid_num_up1,mid_num_down1, down_dif11, down_dif22
        double complex :: coef1
        double precision ::  coef11, coef12
        integer j1, j2,j3,j4,j5,j6,j11,j22,j33,j44,j55,j66, kq
        integer ji1, ji2,ji3,ji4,ji5,ji6
         real(8), external :: w3js


	!<1>: Get general information
	!<1-1>: For oper1
	up_dif1=0
	down_dif1=0
	if(Flag1=='N') then
		up_dif1=oper1%up_dif
		down_dif1=oper1%down_dif
	else if(Flag1=='T') then
		up_dif1=-oper1%up_dif
		down_dif1=oper1%down_dif
	else
		write(*,*) "Flag1 is wrong!"
		return
	endif

	!<1-2>: For oper2
	up_dif2=0
	down_dif2=0
	if(Flag2=='N') then
		up_dif2=oper2%up_dif
		down_dif2=oper2%down_dif
	else if(Flag2=='T') then
		up_dif2=-oper2%up_dif
		down_dif2=oper2%down_dif
	else
		write(*,*) "Fiag2 is wrong!"
		return
	endif

	!<1-3>: For new_oper (up_dif and down_dif)
	new_oper%len=basis%len
	new_oper%up_dif=up_dif1+up_dif2
	!!!new_oper%down_dif=down_dif1+down_dif2   we define it before calling

	up_dif=new_oper%up_dif
	down_dif=new_oper%down_dif

!! from irreducible to regular matrix <J'J'|(Op1 Op2)_q^K| JJ>, q==J'-J

	!<2>: Get general information of new_oper
	new_oper%num=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down

		do new_num_down=abs(num_down-down_dif), num_down+down_dif, su
		new_num_up=num_up+up_dif

		!For new_num_up and new_num_down
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_oper%num=new_oper%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do


	!<3>: Get general information of new_oper
	allocate(new_oper%sub(new_oper%num))
	new_id=0
	new_oper%dim=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down
		new_num_up=num_up+up_dif

		do new_num_down=abs(num_down-down_dif), num_down+down_dif, su
		!For new_num_up and new_num_down
		bs_flag=.false.
		bs_idx=0
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_id=new_id+1
				bs_flag=.true.
				bs_idx=x
				goto 102
			endif
			endif
		end do
		102 continue

		if(bs_flag) then
			!<c>: =========================== For new_oper ==================================
			new_oper%sub(new_id)%num_up=num_up
			new_oper%sub(new_id)%num_down=num_down
			new_oper%sub(new_id)%down_dif=new_num_down-num_down
			new_oper%sub(new_id)%row_dim=basis%sub(bs_idx)%dim
			new_oper%sub(new_id)%sdim=basis%sub(i)%dim

			old_dim=basis%sub(i)%dim
			new_dim=basis%sub(bs_idx)%dim
                
			allocate(new_oper%sub(new_id)%mat(new_dim,old_dim))
			new_oper%sub(new_id)%mat=0.0d0
			new_oper%dim=new_oper%dim+new_oper%sub(new_id)%sdim



        do mid_num_down= abs(num_down-down_dif2), num_down+down_dif2, su

			!<a>: ============ For oper2 ===================
			if(Flag2=='N') then
				mid_num_up1=num_up
				mid_num_down1=num_down
                                down_dif22= mid_num_down-num_down

			else if(Flag2=='T') then
				mid_num_up1=num_up+up_dif2
				mid_num_down1=mid_num_down
                                down_dif22= num_down-mid_num_down
			endif

			idx2=0
			bl_flag2=.false.
			do x=1,oper2%num
				if(oper2%sub(x)%num_up==mid_num_up1) then
				if(oper2%sub(x)%num_down==mid_num_down1) then
				if(oper2%sub(x)%down_dif==down_dif22) then
					bl_flag2=.true.
					idx2=x
					goto 103
				endif
                                endif
				endif
			end do
			103 continue

			!<b>: ============ For oper1 ===================
			if(Flag1=='N') then
				mid_num_up1=num_up+up_dif2
				mid_num_down1=mid_num_down
                                down_dif11= new_num_down-mid_num_down
			else if(Flag1=='T') then
				mid_num_up1=new_num_up
				mid_num_down1=new_num_down
                                down_dif11= mid_num_down-new_num_down
			endif

			idx1=0
			bl_flag1=.false.
			do x=1,oper1%num
				if(oper1%sub(x)%num_up==mid_num_up1) then
				if(oper1%sub(x)%num_down==mid_num_down1) then
				if(oper1%sub(x)%down_dif==down_dif11) then
					bl_flag1=.true.
					idx1=x
					goto 104
				endif
				endif
				endif
			end do
			104 continue

			!Multpilcation
			if(bl_flag1.and.bl_flag2) then
				!<a>: For oper1
				if(Flag1=='N') then
					mid_dim=oper1%sub(idx1)%sdim
					LDA=new_dim
				else if(Flag1=='T') then
					mid_dim=oper1%sub(idx1)%row_dim
					LDA=mid_dim
				endif

				!<b>: For oper2
				if(Flag2=='N') then
					LDB=mid_dim
				else if(Flag2=='T') then
					LDB=old_dim
				endif

        j1=new_num_down    !! from irreducible to regular
        j4=-j1
        j3=num_down
        j6=j3
        j2=new_oper%down_dif
        j5=-j4-j6  !! only choice !! j4=-(j5+j6)
        kq=j5    !! component for tensor
        coef11=w3js(j1,j2,j3,j4, j5, j6)
        if(coef11.eq.0.0)write(*,*)'w3js wrong', j1, j2, j3,j4,j5,j6
        if(coef11.eq.0.0)stop

        coef1=0.d0
        do j6=-mid_num_down, mid_num_down, su !! mid jm component 
        j1=new_num_down    !! first op
        j4=-j1

        j2=oper1%down_dif
        j3=mid_num_down
        j5=-j4-j6
        coef12=w3js(j1,j2,j3,j4,j5,j6)

        j11=mid_num_down    !! second op
        j44=-j6
        j22=oper2%down_dif
        j33=num_down
        j66=j33
        j55=-j44-j66
        coef12=coef12*w3js(j11,j22,j33,j44,j55,j66)*(-1)**((j11+j44)/2)
!! anotehr factor of tensor expansion !(-1)^(-k1+k2+q)sqrt(2k+1)T_q1^k1U_q2^k2w3js(k1,k2,k, q1,q2, -q)

        ji1=oper1%down_dif    !! 
        ji2=oper2%down_dif    !!
        ji3=new_oper%down_dif    !! 
        ji4=j5  !!  q1
        ji5=j55
        ji6=-kq
        coef12=coef12*w3js(ji1,ji2,ji3,ji4, ji5, ji6)*(-1)**((-ji1+ji2+kq)/2)*dsqrt(ji3+1.0d0)
        coef1=coef1+coef12
        enddo

        if(coef1.ne.0.0)then
        coef1=coef*coef1/coef11

           if(realcode)then
              call DGEMM(Flag1,Flag2,new_dim,old_dim,mid_dim,coef1,oper1%sub(idx1)%mat,LDA,&
                                                  &oper2%sub(idx2)%mat,LDB,1.0d0,new_oper%sub(new_id)%mat,new_dim)
                                else
        call ZGEMM(Flag1,Flag2,new_dim,old_dim,mid_dim,coef1,oper1%sub(idx1)%mat,LDA,&
                                                  &oper2%sub(idx2)%mat,LDB,cone,new_oper%sub(new_id)%mat,new_dim)
                                endif

                endif
		endif
	end do
        endif

	end do
	end do

end subroutine block_mul_block_ndia
!new_oper=coef*(oper1^Flag1* oper2^Flag2)
subroutine block_mul_block_ndia1(oper1,Flag1,oper2,Flag2,coef,new_oper,basis)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: oper1,oper2
	type(Total_Block),intent(inout) :: new_oper
	type(Total_Basis),intent(in) :: basis

	logical :: Flags,bl_flag1,bl_flag2,bs_flag
	integer :: up_dif1,down_dif1,up_dif2,down_dif2,up_dif,down_dif
	integer :: i,j,x,y,new_dim,mid_dim,old_dim,idx1,idx2,bs_idx,new_id,LDA,LDB,LDC
	integer :: num_up,num_down,mid_num_up,mid_num_down,new_num_up,new_num_down


	!<1>: Get general information
	!<1-1>: For oper1
	up_dif1=0
	down_dif1=0
	if(Flag1=='N') then
		up_dif1=oper1%up_dif
		down_dif1=oper1%down_dif
	else if(Flag1=='T') then
		up_dif1=-oper1%up_dif
		down_dif1=oper1%down_dif
	else
		write(*,*) "Flag1 is wrong!"
		return
	endif

	!<1-2>: For oper2
	up_dif2=0
	down_dif2=0
	if(Flag2=='N') then
		up_dif2=oper2%up_dif
		down_dif2=oper2%down_dif
	else if(Flag2=='T') then
		up_dif2=-oper2%up_dif
		down_dif2=oper2%down_dif
	else
		write(*,*) "Fiag2 is wrong!"
		return
	endif

	!<1-3>: For new_oper (up_dif and down_dif)
	new_oper%len=basis%len
	new_oper%up_dif=up_dif1+up_dif2
	new_oper%down_dif=down_dif1+down_dif2

	up_dif=new_oper%up_dif
	down_dif=new_oper%down_dif


	!<2>: Get general information of new_oper
	new_oper%num=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down
		new_num_up=num_up+up_dif
		new_num_down=num_down+down_dif

		!For new_num_up and new_num_down
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_oper%num=new_oper%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do


	!<3>: Get general information of new_oper
	allocate(new_oper%sub(new_oper%num))
	new_id=0
	new_oper%dim=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down
		new_num_up=num_up+up_dif
		new_num_down=num_down+down_dif

		!For new_num_up and new_num_down
		bs_flag=.false.
		bs_idx=0
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_id=new_id+1
				bs_flag=.true.
				bs_idx=x
				goto 102
			endif
			endif
		end do
		102 continue

		if(bs_flag) then
			!<a>: ============ For oper2 ===================
			if(Flag2=='N') then
				mid_num_up=num_up
				mid_num_down=num_down
			else if(Flag2=='T') then
				mid_num_up=num_up+up_dif2
				mid_num_down=num_down+down_dif2
			endif

			idx2=0
			bl_flag2=.false.
			do x=1,oper2%num
				if(oper2%sub(x)%num_up==mid_num_up) then
				if(oper2%sub(x)%num_down==mid_num_down) then
					bl_flag2=.true.
					idx2=x
					goto 103
				endif
				endif
			end do
			103 continue

			!<b>: ============ For oper1 ===================
			if(Flag1=='N') then
				mid_num_up=num_up+up_dif2
				mid_num_down=num_down+down_dif2
			else if(Flag1=='T') then
				mid_num_up=new_num_up
				mid_num_down=new_num_down
			endif

			idx1=0
			bl_flag1=.false.
			do x=1,oper1%num
				if(oper1%sub(x)%num_up==mid_num_up) then
				if(oper1%sub(x)%num_down==mid_num_down) then
					bl_flag1=.true.
					idx1=x
					goto 104
				endif
				endif
			end do
			104 continue

			!<c>: =========================== For new_oper ==================================
			new_oper%sub(new_id)%num_up=num_up
			new_oper%sub(new_id)%num_down=num_down
			new_oper%sub(new_id)%row_dim=basis%sub(bs_idx)%dim
			new_oper%sub(new_id)%sdim=basis%sub(i)%dim

			old_dim=basis%sub(i)%dim
			new_dim=basis%sub(bs_idx)%dim
			allocate(new_oper%sub(new_id)%mat(new_dim,old_dim))
			new_oper%sub(new_id)%mat=0.0d0
			new_oper%dim=new_oper%dim+new_oper%sub(new_id)%sdim

			!Multpilcation
			if(bl_flag1.and.bl_flag2) then
				!<a>: For oper1
				if(Flag1=='N') then
					mid_dim=oper1%sub(idx1)%sdim
					LDA=new_dim
				else if(Flag1=='T') then
					mid_dim=oper1%sub(idx1)%row_dim
					LDA=mid_dim
				endif

				!<b>: For oper2
				if(Flag2=='N') then
					LDB=mid_dim
				else if(Flag2=='T') then
					LDB=old_dim
				endif
				call DGEMM(Flag1,Flag2,new_dim,old_dim,mid_dim,coef,oper1%sub(idx1)%mat,LDA,&
						  &oper2%sub(idx2)%mat,LDB,1.0d0,new_oper%sub(new_id)%mat,new_dim)
			endif
		endif
	end do

end subroutine block_mul_block_ndia1


subroutine deallocate_block(block)
	use pubdata
	implicit none

	type(Total_block),intent(inout) :: block

	integer :: i

	do i=1,block%num
		deallocate(block%sub(i)%mat)
	end do
	deallocate(block%sub)
end subroutine deallocate_block

!Transfer model data by allocating new space
subroutine model_transfer(inmod,outmod)
	use pubdata
	implicit none

	type(Total_Model),intent(in) :: inmod
	type(Total_Model),intent(inout) :: outmod

	integer :: ii

	!Transfer index and Hamiltonian
	outmod%len=inmod%len
	call block_transfer(inmod%ham,outmod%ham)

	!For operators in the middle of system
	do ii=1,nleg11(inmod%len) !!!Ny

        if(tjmodel.ne.11)then
		call block_transfer(inmod%sub(ii)%elec_up,outmod%sub(ii)%elec_up)
		call block_transfer(inmod%sub(ii)%elec_down,outmod%sub(ii)%elec_down)
                        endif
        
                if(tjmodel==1)then
		call block_transfer(inmod%sub(ii)%spin_sd,outmod%sub(ii)%spin_sd)
                        endif
		if(v123==1)call block_transfer(inmod%sub(ii)%num_sn,outmod%sub(ii)%num_sn)
	end do


           if(lring.eq.1)then
                do ii=1, nleg2(inmod%len, ms+1)
                if(allocated(inmod%spin1(ii)%spin_sd%sub))then
        call block_transfer(inmod%spin1(ii)%spin_sd,outmod%spin1(ii)%spin_sd)
                endif
                enddo
        endif


end subroutine model_transfer


!deallocate model data
subroutine deallocate_model(model)
	use pubdata
	implicit none

	type(Total_Model),intent(inout) :: model
	integer :: ii

	!Deallocate Hamiltonian
   !Deallocate Hamiltonian
        if(allocated(model%ham%sub))then
        call deallocate_block(model%ham)
                endif


	!For operators in the middle of system
	do ii=1,nleg11(model%len)

                if(allocated(model%sub(ii)%elec_up%sub))then
		call deallocate_block(model%sub(ii)%elec_up)
		call deallocate_block(model%sub(ii)%elec_down)
		if(v123==1)call deallocate_block(model%sub(ii)%num_sn)
                        endif
                if(allocated(model%sub(ii)%spin_sd%sub))then
		call deallocate_block(model%sub(ii)%spin_sd)
                        endif
	end do

           if(lring.eq.1)then
                do ii=1, nleg2(model%len, ms+1)
                if(allocated(model%spin1(ii)%spin_sd%sub))then
            call deallocate_block(model%spin1(ii)%spin_sd)
                endif
                enddo
        endif
end subroutine deallocate_model


!=========================================================================
!Allocate new wavefunction according to the superbasis information
!=========================================================================
subroutine allocate_wave(super,wave)
	use pubdata
	implicit none

	type(Super_Basis),intent(in) :: super
	type(Wavefunction),intent(inout) :: wave

	integer :: i,sys_dim,env_dim

	wave%sys_len=super%sys_len
	wave%env_len=super%env_len
	wave%num=super%num
	wave%dim=super%dim

	allocate(wave%sub(wave%num))
	do i=1,wave%num
		wave%sub(i)%sys_num_up=super%sub(i)%sys_num_up
		wave%sub(i)%sys_num_down=super%sub(i)%sys_num_down
		wave%sub(i)%env_num_up=super%sub(i)%env_num_up
		wave%sub(i)%env_num_down=super%sub(i)%env_num_down
		wave%sub(i)%sys_dim=super%sub(i)%sys_dim
		wave%sub(i)%env_dim=super%sub(i)%env_dim

		allocate(wave%sub(i)%vec(wave%sub(i)%sys_dim,wave%sub(i)%env_dim))
		wave%sub(i)%vec=0.0d0
	end do

end subroutine allocate_wave


!Transfer wavefunction by allocating new space
subroutine wave_transfer(inwave,outwave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: inwave
	type(Wavefunction),intent(inout) :: outwave

	integer :: i,sys_dim,env_dim

	outwave%sys_len=inwave%sys_len
	outwave%env_len=inwave%env_len
	outwave%num=inwave%num
	outwave%dim=inwave%dim

	allocate(outwave%sub(outwave%num))
	do i=1,outwave%num
		outwave%sub(i)%sys_num_up=inwave%sub(i)%sys_num_up
		outwave%sub(i)%sys_num_down=inwave%sub(i)%sys_num_down
		outwave%sub(i)%env_num_up=inwave%sub(i)%env_num_up
		outwave%sub(i)%env_num_down=inwave%sub(i)%env_num_down
		outwave%sub(i)%sys_dim=inwave%sub(i)%sys_dim
		outwave%sub(i)%env_dim=inwave%sub(i)%env_dim

		allocate(outwave%sub(i)%vec(outwave%sub(i)%sys_dim,outwave%sub(i)%env_dim))
		outwave%sub(i)%vec=inwave%sub(i)%vec
	end do

end subroutine wave_transfer


!Notice that inwave and outwave have the same configuration
subroutine wave_pass(inwave,outwave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: inwave
	type(Wavefunction),intent(inout) :: outwave

	integer :: i

	outwave%sys_len=inwave%sys_len
	outwave%env_len=inwave%env_len
	outwave%num=inwave%num
	outwave%dim=inwave%dim

	do i=1,outwave%num
		outwave%sub(i)%sys_num_up=inwave%sub(i)%sys_num_up
		outwave%sub(i)%sys_num_down=inwave%sub(i)%sys_num_down
		outwave%sub(i)%env_num_up=inwave%sub(i)%env_num_up
		outwave%sub(i)%env_num_down=inwave%sub(i)%env_num_down
		outwave%sub(i)%sys_dim=inwave%sub(i)%sys_dim
		outwave%sub(i)%env_dim=inwave%sub(i)%env_dim
		outwave%sub(i)%vec=inwave%sub(i)%vec
	end do

end subroutine wave_pass

!Deallocate wavefunction
subroutine deallocate_wave(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(inout) :: wave

	integer :: i
	
	do i=1,wave%num
		deallocate(wave%sub(i)%vec)
	end do
	deallocate(wave%sub)

end subroutine deallocate_wave


!Note: left and right have the same configuration
real(8) function wave_product(left,right)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: left,right

	integer :: i,x

	wave_product=0.0d0
	do i=1,left%num
		do x=1,left%sub(i)%env_dim
			wave_product=wave_product+dot_product(left%sub(i)%vec(:,x),right%sub(i)%vec(:,x))
		end do
	end do
	return

end function wave_product


!Get the norm of the wavefunction
real(8) function wave_norm(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave

	integer :: i,x

	wave_norm=0.0d0
	do i=1,wave%num
		do x=1,wave%sub(i)%env_dim
			wave_norm=wave_norm+dot_product(wave%sub(i)%vec(:,x),wave%sub(i)%vec(:,x))
		end do
	end do

	wave_norm=dsqrt(wave_norm)
	return

end function wave_norm


!Nomalize the wavefunction
subroutine wave_normalize(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(inout) :: wave

	integer :: i,j
	real(8) :: vec_norm
	real(8),external :: wave_norm

	vec_norm=wave_norm(wave)
	do i=1,wave%num
		wave%sub(i)%vec=wave%sub(i)%vec/vec_norm
	end do

end subroutine wave_normalize


!Initialize the wavefunction
subroutine wave_initialize(wave)
        use pubdata
        implicit none

        type(Wavefunction),intent(inout) :: wave
!!!! final random WF
  double precision,allocatable :: mid_r(:,:),mid_r1(:,:)
        integer :: i, sys_dim, env_dim
      do i=1,wave%num
            sys_dim=wave%sub(i)%sys_dim
                env_dim=wave%sub(i)%env_dim

                allocate(mid_r(sys_dim,env_dim))
                allocate(mid_r1(sys_dim,env_dim))
                call random_number(mid_r)
                call random_number(mid_r1)
                if(.not.realcode)then
     wave%sub(i)%vec=dcmplx(mid_r,mid_r1)  !!!!mid_r) !,(thetay-thetax)*mid_r)
                        else
     wave%sub(i)%vec=mid_r  !!!!mid_r) !,(thetay-thetax)*mid_r)
                endif
                deallocate(mid_r, mid_r1)
        end do

        call wave_normalize(wave)

end subroutine wave_initialize


subroutine wave_initialize_zero(wave)
        use pubdata
        implicit none

       double precision,allocatable :: mid_r(:,:),mid_r1(:,:)
        type(Wavefunction),intent(inout) :: wave
        integer :: i, sys_dim, env_dim
      do i=1,wave%num
     wave%sub(i)%vec=0.0d0  !!!!mid_r) !,(thetay-thetax)*mid_r)
        end do
end subroutine wave_initialize_zero

double complex function wave_overlap(left, right)
        use pubdata
        implicit none

        type(Wavefunction),intent(in) :: left, right
        integer :: i,j, x
        wave_overlap=0.0d0
        do i=1,left%num
                do j=1, right%num
                   if(left%sub(i)%sys_num_up==right%sub(j)%sys_num_up)then
                   if(left%sub(i)%sys_num_down==right%sub(j)%sys_num_down)then
                do x=1,   left%sub(i)%env_dim
                if(left%sub(i)%env_dim.ne.right%sub(j)%env_dim)write(*,*)'overlap wrong'
        wave_overlap=wave_overlap+dot_product(left%sub(i)%vec(:,x),right%sub(i)%vec(:,x))
                end do
                endif
                endif
        end do
        end do
        return
end function wave_overlap



!Get the Gram¨CSchmidt projection operator
!proj_u(v) = <u,v>*u/<u,u> = wave_new
subroutine Gram_Schmidt_WF(wave_u,wave_v,wave_new)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave_u,wave_v
	type(Wavefunction) :: wave_new

	integer :: i,j,k
	real(8) :: prod_uv,prod_uu
	real(8),external :: wave_product

	!<1>: Get <u,v> and <u,u>
	prod_uv=wave_product(wave_u,wave_v)
	prod_uu=wave_product(wave_u,wave_u)

	!<2>: Get proj_u(v)=<u,v>*u/<u,u> = wave_new
	wave_new%sys_len=wave_u%sys_len
	wave_new%env_len=wave_u%env_len
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_u%sub(i)%vec*(prod_uv/prod_uu)
	end do

end subroutine Gram_Schmidt_WF


!Get (level+1)'s wave_new orthogonal to all level low-lying states
!level=1 for the ground state, level=2 for the first excited state
!level=3 for the second excited state (level=k+1)
!Input: wave_new = v_k with (level+1 = k)
!Output: wave_new = u_k = v_k -\sum_{j=1}^{k-1} proj_{u_j}(v_k)
subroutine Gram_Schmidt_process(wave_prev,levels,wave_new)
	use pubdata
	implicit none

	integer,intent(in) :: levels
	type(Wavefunction) :: wave_prev(levels)
	type(Wavefunction) :: wave_new

	integer :: i,j,x,y
	real(8) :: prod_uv(levels),prod_uu(levels)
	real(8),external :: wave_product
	type(Wavefunction) :: tmp_vec

	!<1>: Get <u_j,v> (coef_uv) and <u_j,u_j> (coef_uu)
	do j=1,levels
		prod_uv(j) = wave_product(wave_prev(j),wave_new)
		prod_uu(j) = wave_product(wave_prev(j),wave_prev(j))
	end do

	!<2>: Get (level+1)'s excited state: wave_new (u_k)
	call wave_transfer(wave_new,tmp_vec)
	do i=1,tmp_vec%num
		tmp_vec%sub(i)%vec=0.0d0
	end do

	do j=1,levels
	do i=1,tmp_vec%num
		tmp_vec%sub(i)%vec = tmp_vec%sub(i)%vec + (prod_uv(j)/prod_uu(j))*wave_prev(j)%sub(i)%vec
	end do
	end do

	!<3>: Normalize wave_new
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_new%sub(i)%vec - tmp_vec%sub(i)%vec
	end do
	call wave_normalize(wave_new)
	call deallocate_wave(tmp_vec)

end subroutine Gram_Schmidt_process


!Get (level+1)'s wave_new orthogonal to all level low-lying states
!level=1 for the ground state, level=2 for the first excited state
!level=3 for the second excited state (level=k+1)
!Input: wave_new = v_k with (level+1 = k)
!Output: wave_new = u_k = v_k -\sum_{j=1}^{k-1} proj_{u_j}(v_k)
subroutine Gram_Schmidt_Lanczos(wave_prev,levels,wave_lanc,lanc_num,wave_new)
	use pubdata
	implicit none

	integer,intent(in) :: levels,lanc_num
	type(Wavefunction) :: wave_prev(levels),wave_lanc(lanc_num)
	type(Wavefunction) :: wave_new

	integer :: i,j,x,y
	real(8) :: prod_uv(levels+lanc_num),prod_uu(levels+lanc_num)
	real(8),external :: wave_product
	type(Wavefunction) :: prev_vec,lanc_vec

	!<1>: Get <u_j,v> (coef_uv) and <u_j,u_j> (coef_uu)
	do j=1,levels
		prod_uv(j) = wave_product(wave_prev(j),wave_new)
		prod_uu(j) = wave_product(wave_prev(j),wave_prev(j))
	end do

	do j=1,lanc_num
		prod_uv(levels+j) = wave_product(wave_lanc(j),wave_new)
		prod_uu(levels+j) = wave_product(wave_lanc(j),wave_lanc(j))
	end do

	!<2>: Get (lanc_num+1)'s excited state: wave_new (u_k)
	!<2-1>: For wave_prev
	call wave_transfer(wave_new,prev_vec)
	do i=1,prev_vec%num
		prev_vec%sub(i)%vec=0.0d0
	end do

	do j=1,levels
	do i=1,prev_vec%num
		prev_vec%sub(i)%vec = prev_vec%sub(i)%vec + wave_prev(j)%sub(i)%vec*(prod_uv(j)/prod_uu(j))
	end do
	end do

	!<2-2>: For wave_lanc
	call wave_transfer(wave_new,lanc_vec)
	do i=1,lanc_vec%num
		lanc_vec%sub(i)%vec=0.0d0
	end do

	do j=1,lanc_num
	do i=1,lanc_vec%num
		lanc_vec%sub(i)%vec = lanc_vec%sub(i)%vec + wave_lanc(j)%sub(i)%vec*(prod_uv(levels+j)/prod_uu(levels+j))
	end do
	end do

	!<3>: Normalize wave_new
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_new%sub(i)%vec - prev_vec%sub(i)%vec - lanc_vec%sub(i)%vec
	end do
	call wave_normalize(wave_new)
	call deallocate_wave(prev_vec)
	call deallocate_wave(lanc_vec)

end subroutine Gram_Schmidt_Lanczos


!==================================================================
!Transfer basis information from the old one to the new one
!==================================================================
subroutine basis_transfer(old_basis,new_basis)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: old_basis
	type(Total_Basis),intent(inout) :: new_basis

	integer :: i,x

         if(allocated(new_basis%sub))call deallocate_basis(new_basis)

	new_basis%len=old_basis%len
	new_basis%num=old_basis%num
	new_basis%dim=old_basis%dim

	allocate(new_basis%sub(new_basis%num))
	do i=1,new_basis%num
		new_basis%sub(i)%new_num_up=old_basis%sub(i)%new_num_up
		new_basis%sub(i)%new_num_down=old_basis%sub(i)%new_num_down
		new_basis%sub(i)%idx=old_basis%sub(i)%idx
		new_basis%sub(i)%num=old_basis%sub(i)%num
		new_basis%sub(i)%dim=old_basis%sub(i)%dim
		allocate(new_basis%sub(i)%sub(new_basis%sub(i)%num))

		do x=1,new_basis%sub(i)%num
			new_basis%sub(i)%sub(x)%spos=old_basis%sub(i)%sub(x)%spos
			new_basis%sub(i)%sub(x)%sdim=old_basis%sub(i)%sub(x)%sdim
			new_basis%sub(i)%sub(x)%bl_num_up=old_basis%sub(i)%sub(x)%bl_num_up
			new_basis%sub(i)%sub(x)%bl_num_down=old_basis%sub(i)%sub(x)%bl_num_down
			new_basis%sub(i)%sub(x)%st_num_up=old_basis%sub(i)%sub(x)%st_num_up
			new_basis%sub(i)%sub(x)%st_num_down=old_basis%sub(i)%sub(x)%st_num_down
		end do
	end do

end subroutine basis_transfer

subroutine deallocate_basis(basis)
	use pubdata
	implicit none

	type(Total_basis),intent(inout) :: basis

	integer :: i

	do i=1,basis%num
		deallocate(basis%sub(i)%sub)
	end do
	deallocate(basis%sub)

end subroutine deallocate_basis


!Deallocate the superbasis information
subroutine deallocate_super_basis(super)
	use pubdata
	implicit none

	type(Super_basis),intent(inout) :: super

	deallocate(super%sub)

end subroutine deallocate_super_basis


!==================================
!For filename operation
!==================================
!Translate number into string
subroutine num_to_str(num,num_str)
	implicit none

	integer,intent(in) :: num
	character(len=6) :: num_str

	integer :: i,idx,ra,res
	character(len=6) :: str1

	!<1>: Initiate
	num_str(1:6)=" "

	!<2>: Get num_str
	idx=1
	res=num
	do while(res>0)
		if(abs(res)<10) then
			str1(idx:idx)=char(res+48)
			goto 111
		endif
		ra=res-res/10*10
		str1(idx:idx)=char(ra+48)
		res=res/10
		idx=idx+1
	end do
	111 continue

	do i=1,idx
		num_str(i:i)=str1(idx-i+1:idx-i+1)
	end do

end subroutine num_to_str


!Link string with string
subroutine str_link_str(str_name,str_len,num)
	implicit none

	integer,intent(in) :: str_len,num
	character(len=20) :: str_name

	integer :: num_len
	character(len=6) :: num_str

	call num_to_str(num,num_str)

	num_len=len_trim(num_str)
	str_name(str_len+1:str_len+num_len)=num_str(1:num_len)

end subroutine str_link_str


!Save wavefunction to disk
subroutine wave_to_disk(wave,fileid,filename,name_len)
	use pubdata
	implicit none

	integer,intent(in) :: fileid,name_len
	character(len=name_len) :: filename
	type(Wavefunction),intent(in) :: wave

	integer :: i,x,y
	
	open(fileid,file=filename,form='unformatted')

	write(fileid) wave%sys_len
	write(fileid) wave%env_len
	write(fileid) wave%num
	write(fileid) wave%dim

	do i=1,wave%num
		write(fileid)
		write(fileid) wave%sub(i)%sys_num_up,wave%sub(i)%sys_num_down
		write(fileid) wave%sub(i)%env_num_up,wave%sub(i)%env_num_down
		write(fileid) wave%sub(i)%sys_dim,wave%sub(i)%env_dim

		!write(fileid) wave%sub(i)%vec
				write(fileid) wave%sub(i)%vec  !!!(x,y)
	end do

	close(fileid)

end subroutine wave_to_disk


!Read wavefunction from disk
subroutine wave_from_disk(wave,fileid,filename,name_len)
	use pubdata
	implicit none

	integer,intent(in) :: fileid,name_len
	character(len=name_len) :: filename
	type(Wavefunction),intent(inout) :: wave

	integer :: i,x,y

	open(fileid,file=filename,form='unformatted')

	read(fileid) wave%sys_len
	read(fileid) wave%env_len
	read(fileid) wave%num
	read(fileid) wave%dim

	allocate(wave%sub(wave%num))
	do i=1,wave%num
		read(fileid)
		read(fileid) wave%sub(i)%sys_num_up,wave%sub(i)%sys_num_down
		read(fileid) wave%sub(i)%env_num_up,wave%sub(i)%env_num_down
		read(fileid) wave%sub(i)%sys_dim,wave%sub(i)%env_dim

		allocate(wave%sub(i)%vec(wave%sub(i)%sys_dim,wave%sub(i)%env_dim))
				read(fileid) wave%sub(i)%vec  !!!!(x,y)
	end do

	close(fileid)
end subroutine wave_from_disk

subroutine block_to_disk(block,fileid)
	use pubdata
	implicit none

	integer,intent(in) :: fileid
	type(Total_Block),intent(in) :: block

	integer :: i,x,y,row_dim,col_dim

	write(fileid) block%len
	write(fileid) block%num
	write(fileid) block%dim
	write(fileid) block%up_dif
	write(fileid) block%down_dif

	do i=1,block%num
		write(fileid)
		write(fileid) block%sub(i)%num_up,block%sub(i)%num_down, block%sub(i)%down_dif
		write(fileid) block%sub(i)%row_dim,block%sub(i)%sdim
				write(fileid) block%sub(i)%mat  !!!(x,y)
	end do

end subroutine block_to_disk


!Read block from disk
subroutine block_from_disk(block,fileid)
	use pubdata
	implicit none

	integer,intent(in) :: fileid
	type(Total_Block),intent(inout) :: block

	integer :: i,x,y

	read(fileid) block%len
	read(fileid) block%num
	read(fileid) block%dim
	read(fileid) block%up_dif
	read(fileid) block%down_dif

	allocate(block%sub(block%num))
	do i=1,block%num
		read(fileid)
		read(fileid) block%sub(i)%num_up,block%sub(i)%num_down, block%sub(i)%down_dif
		read(fileid) block%sub(i)%row_dim,block%sub(i)%sdim

		allocate(block%sub(i)%mat(block%sub(i)%row_dim,block%sub(i)%sdim))
				read(fileid) block%sub(i)%mat  !!(x,y)
	end do

end subroutine block_from_disk

subroutine model_to_disk(model,fileid,mod_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,mod_idx
	type(Total_Model),intent(in) :: model
	
	integer :: ii

	!<1>: Save index
	if(flag) then
		open(fileid,file=sysname(mod_idx,1),form='unformatted')
	else
		open(fileid,file=envname(mod_idx,1),form='unformatted')
	endif

	write(fileid) model%len

	!<2>: Save general data
	call block_to_disk(model%ham,fileid)

	!For operators in the middle of the system
	do ii=1,nleg11(model%len)
		call block_to_disk(model%sub(ii)%elec_up,fileid)
		call block_to_disk(model%sub(ii)%elec_down,fileid)
                if(tjmodel==1)then
		call block_to_disk(model%sub(ii)%spin_sd,fileid)
                endif
		if(v123==1)call block_to_disk(model%sub(ii)%num_sn,fileid)
	end do

           if(lring.eq.1)then
             do ii=1, nleg2(model%len, ms+1)
             if(allocated(model%spin1(ii)%spin_sd%sub))then
	call block_to_disk(model%spin1(ii)%spin_sd,fileid)
                endif
                enddo
        endif

	close(fileid)

end subroutine model_to_disk

!Read model from disk
!Note: flag=.true. for system, flag=.false. for environment
subroutine model_from_disk(model,fileid,mod_idx,Flag)
	use pubdata
	implicit none

	logical,intent(in) :: Flag
	integer,intent(in) :: fileid,mod_idx
	type(Total_Model),intent(inout) :: model
	
	integer :: ii

	!<1>: Read index
	if(flag) then
		open(fileid,file=sysname(mod_idx,1),form='unformatted')
	else
		open(fileid,file=envname(mod_idx,1),form='unformatted')
	endif

	read(fileid) model%len

	!<2>: Save general data
	call block_from_disk(model%ham,fileid)

	!For operators in the middle of system
	do ii=1, nleg11(model%len)
		call block_from_disk(model%sub(ii)%elec_up,fileid)
		call block_from_disk(model%sub(ii)%elec_down,fileid)
                if(tjmodel==1)then
		call block_from_disk(model%sub(ii)%spin_sd,fileid)
                endif
		if(v123==1)call block_from_disk(model%sub(ii)%num_sn,fileid)
	end do

           if(lring.eq.1)then
             do ii=1, nleg2(model%len, ms+1)
	call block_from_disk(model%spin1(ii)%spin_sd,fileid)
                enddo
        endif


	!==================== For space saving ======================

        if((infi_del.ne.1).or.(mod_idx.gt.num_site0/2+4))then
	if(Flag) then
		open(fileid,file=sysname(mod_idx+1,1),form='unformatted')
		close(fileid,status='delete')
	else
		open(fileid,file=envname(mod_idx+1,1),form='unformatted')
		close(fileid,status='delete')
	endif
        endif
end subroutine model_from_disk

subroutine truns_to_disk(truns,fileid,tr_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,tr_idx
	type(Total_Block),intent(in) :: truns

	!Open file
	if(flag) then
		open(fileid,file=sysname(tr_idx,2),form='unformatted')
	else
		open(fileid,file=envname(tr_idx,2),form='unformatted')
	endif
	call block_to_disk(truns,fileid)
	close(fileid)
end subroutine truns_to_disk

!Read truncation matrix to disk:
!Note: flag=.true. for system
!		 flag=.false. for environment
subroutine truns_from_disk(truns,fileid,tr_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,tr_idx
	type(Total_Block),intent(inout) :: truns

	!Open file
	if(flag) then
		open(fileid,file=sysname(tr_idx,2),form='unformatted')
	else
		open(fileid,file=envname(tr_idx,2),form='unformatted')
	endif

	call block_from_disk(truns,fileid)

	close(fileid)

end subroutine truns_from_disk


!Save basis to disk
subroutine basis_to_disk(basis,fileid,bs_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,bs_idx
	type(Total_Basis),intent(in) :: basis

	integer :: i,x,y

	!Open file
	if(flag) then
		open(fileid,file=sysname(bs_idx,3))
	else
		open(fileid,file=envname(bs_idx,3))
	endif

	write(fileid,101) "len=",basis%len
	write(fileid,101) "num=",basis%num
	write(fileid,101) "dim=",basis%dim

	do i=1,basis%num
		write(fileid,*)
		write(fileid,101) "idx=",basis%sub(i)%idx
		write(fileid,101) "num=",basis%sub(i)%num
		write(fileid,101) "dim=",basis%sub(i)%dim
		write(fileid,102) "new_num_up=",basis%sub(i)%new_num_up
		write(fileid,103) "new_num_down=",basis%sub(i)%new_num_down

		do x=1,basis%sub(i)%num
			write(fileid,101) "pos=",basis%sub(i)%sub(x)%spos
			write(fileid,101) "dim=",basis%sub(i)%sub(x)%sdim
			write(fileid,104) "bl_num_up=",basis%sub(i)%sub(x)%bl_num_up,"bl_num_down=",basis%sub(i)%sub(x)%bl_num_down
			write(fileid,104) "st_num_up=",basis%sub(i)%sub(x)%st_num_up,"st_num_down=",basis%sub(i)%sub(x)%st_num_down
		end do
	end do
	close(fileid)

101 format(A4,I6)
102 format(A11,I4)
103 format(A13,I4)
104 format(A10,I4,2X,A12,I4)

end subroutine basis_to_disk


!Read basis from disk
subroutine basis_from_disk(basis,fileid,bs_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,bs_idx
	type(Total_Basis),intent(inout) :: basis

	integer :: i,x,y
	character(len=4) :: str4
	character(len=10) :: str10
	character(len=11) :: str11
	character(len=12) :: str12
	character(len=13) :: str13

	!Open file
	if(flag) then
		open(fileid,file=sysname(bs_idx,3))
	else
		open(fileid,file=envname(bs_idx,3))
	endif

	read(fileid,101) str4,basis%len
	read(fileid,101) str4,basis%num
	read(fileid,101) str4,basis%dim

	allocate(basis%sub(basis%num))
	do i=1,basis%num
		read(fileid,*)
		read(fileid,101) str4,basis%sub(i)%idx
		read(fileid,101) str4,basis%sub(i)%num
		read(fileid,101) str4,basis%sub(i)%dim
		read(fileid,102) str11,basis%sub(i)%new_num_up
		read(fileid,103) str13,basis%sub(i)%new_num_down

		allocate(basis%sub(i)%sub(basis%sub(i)%num))
		do x=1,basis%sub(i)%num
			read(fileid,101) str4,basis%sub(i)%sub(x)%spos
			read(fileid,101) str4,basis%sub(i)%sub(x)%sdim
			read(fileid,104) str10,basis%sub(i)%sub(x)%bl_num_up,str12,basis%sub(i)%sub(x)%bl_num_down
			read(fileid,104) str10,basis%sub(i)%sub(x)%st_num_up,str12,basis%sub(i)%sub(x)%st_num_down
		end do
	end do
	close(fileid)

101 format(A4,I6)
102 format(A11,I4)
103 format(A13,I4)
104 format(A10,I4,2X,A12,I4)

end subroutine basis_from_disk
