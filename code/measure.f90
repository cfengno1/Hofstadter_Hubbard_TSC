!======================================================================
!Get all <n_i> in system and environment blocks
subroutine measure_operator_dia(ave_dia,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	double precision,intent(inout) :: ave_dia(Num_site)
        double complex  ave1(num_site)
        double precision a11

	integer :: i,sys_len,env_len,dx,dy
	type(Total_Block) :: new_oper

	!<1>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	ave_dia=0.0d0
        ave1=dcmplx(0.0, 0.0)

	!<2>: Get <n_i> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(ave1(i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do

	!Get <n_i> in the middle of the system
	call measure_sys_site(ave1(sys_len),st_oper,sys_bsm(sys_len),wave)

	!<2>: Get <n_i> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(ave1(Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do

	!Get <n_i> in the middle of the system
	call measure_env_site(ave1(Num_site-env_len+1),st_oper,env_bsm(env_len),wave)

	!<3>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max

        ave_dia=ave1
        a11=0.0d0
	do i=1,Num_site
		dx=Lattice(1,i)
		dy=Lattice(2,i)
        a11=a11+dreal(ave1(i))
		write(10,101) (i-1)/ny+1,mod(i-1, ny)+1,dreal(ave1(i)), a11
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I2,1X,2F16.12)
101     format(2i6, 3f14.8)
end subroutine measure_operator_dia

!Get all diagonal operator correlatioins using memory directly
!============================================================================
subroutine Get_oper_cor_dia(oper_dia,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	complex*16 :: value
	integer :: i,j,sys_len,env_len,ti,tj
	type(Total_Basis) :: sys_bs,env_bs
	type(Total_Block) :: st_oper2,new_oper

	!<1>: For information saving to the disk
	open(10,file="oper_dia_tmp.dat",position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max,Filename
	close(10)


	!<2>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)
	call block_mul_block_dia(st_oper,'N',st_oper,'N',cone,st_oper2)
	oper_dia=0.0d0

	!On-site density correlation
	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(oper_dia(i,i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(oper_dia(sys_len,sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(oper_dia(Num_site-i+1,Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(oper_dia(Num_site-env_len+1,Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)
	
	!Save to the disk
	do i=1,Num_site
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,i,oper_dia(i,i)
		close(10)
	end do


	!<2>: Get sysopers and envopers
	do i=1,sys_len-1
		call block_transfer(st_oper,sysopers(i))
	end do

	do i=1,env_len-1
		call block_transfer(st_oper,envopers(i))
	end do


	!<3>: Get sys_bl_bl and env_bl_bl operator correlation
	call Get_sys_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	call Get_env_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)


	!<4>: Get sys_bl_st and sys_bl_env_st
	do i=1,sys_len-1
		call sys_block_site_cor_dia(oper_dia(i,sys_len),sysopers(i),st_oper,sys_bs,wave)
		oper_dia(sys_len,i)=dconjg(oper_dia(i,sys_len))

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,sys_len,oper_dia(i,sys_len)
		write(10,*) sys_len,i,oper_dia(sys_len,i)
		close(10)

		call sys_block_env_site_cor_dia(oper_dia(i,Num_site-env_len+1),sysopers(i),st_oper,sys_bs,env_bs,wave)
		oper_dia(Num_site-env_len+1,i)=oper_dia(i,Num_site-env_len+1)

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,Num_site-env_len+1,oper_dia(i,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,i,oper_dia(Num_site-env_len+1,i)
		close(10)
	end do

	!<5>: Get env_bl_st and sys_st_env_bl
	do i=1,env_len-1
		call env_block_site_cor_dia(oper_dia(Num_site-i+1,Num_site-env_len+1),envopers(i),st_oper,env_bs,wave)
		oper_dia(Num_site-env_len+1,Num_site-i+1)=oper_dia(Num_site-i+1,Num_site-env_len+1)

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,Num_site-env_len+1,real(oper_dia(Num_site-i+1,Num_site-env_len+1))
		write(10,*) Num_site-env_len+1,Num_site-i+1,real(oper_dia(Num_site-env_len+1,Num_site-i+1))
		close(10)

		call sys_site_env_block_cor_dia(oper_dia(sys_len,Num_site-i+1),st_oper,envopers(i),sys_bs,env_bs,wave)
		oper_dia(Num_site-i+1,sys_len)=dconjg(oper_dia(sys_len,Num_site-i+1))

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,sys_len,real(oper_dia(Num_site-i+1,sys_len))
		write(10,*) sys_len,Num_site-i+1,real(oper_dia(sys_len,Num_site-i+1))
		close(10)
	end do

	!<6>: Get sys_bl_env_bl
	do i=1,sys_len-1
	do j=1,env_len-1
		call sys_block_env_block_cor_dia(oper_dia(i,Num_site-j+1),sysopers(i),envopers(j),sys_bs,env_bs,wave)
		oper_dia(Num_site-j+1,i)=dconjg(oper_dia(i,Num_site-j+1))
			
		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,Num_site-j+1,real(oper_dia(i,Num_site-j+1))
		write(10,*) Num_site-j+1,i,real(oper_dia(Num_site-j+1,i))
		close(10)
	end do
	end do

	call sys_site_env_site_cor_dia(oper_dia(sys_len,Num_site-env_len+1),st_oper,st_oper,sys_bs,env_bs,wave)
	oper_dia(Num_site-env_len+1,sys_len)=dconjg(oper_dia(sys_len,Num_site-env_len+1))

	!Save to disk
	open(10,file="oper_dia_tmp.dat",position='append')
	write(10,*) sys_len,Num_site-env_len+1,real(oper_dia(sys_len,Num_site-env_len+1))
	write(10,*) Num_site-env_len+1,sys_len,real(oper_dia(Num_site-env_len+1,sys_len))
	write(10,*)
	close(10)


	!Free space
	do i=1,sys_len-1
		call deallocate_block(sysopers(i))
	end do

	do i=1,env_len-1
		call deallocate_block(envopers(i))
	end do
	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)


	!<5>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		do j=1,Num_site
			write(10,111) i,j,real(oper_dia(i,j)),imag(oper_dia(I,j))
		end do
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I4,1X,2F16.12)

end subroutine Get_oper_cor_dia
!==============================================================
!Get all density correlcations in system block using memory
!==============================================================
subroutine Get_sys_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
        complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,sys_oper


	!<1>: For general information
	sys_len=wave%sys_len

	!Get g(i,j) for (i<j)
	do j=2,sys_len-1
		call update_site_dia(st_oper,sys_bsm(j),new_oper2)

		!Update sysopers(sys_len-1)
		if(j==sys_len-1) then
			call deallocate_block(sysopers(j))
			if(j<=trun_idx) then
				call update_site_dia(st_oper,sys_bsm(j),sysopers(j))
			else
				call update_site_dia(st_oper,sys_bsm(j),mid_oper)
				call update_trun_dia(mid_oper,sysopers(j),systruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_dia(st_oper,sys_bsm(i),mid_oper)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_dia(mid_oper,tmp_oper,systruns(i))
					endif
					call update_block_dia(tmp_oper,sys_bsm(j),new_oper1)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_dia(sysopers(i),sys_bsm(j),new_oper1)
				endif
			else
				call update_block_dia(sysopers(i),sys_bsm(j),new_oper1)
			endif

			!Update sysopers(i)
			call deallocate_block(sysopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,sysopers(i))
			else
				call update_trun_dia(new_oper1,sysopers(i),systruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
			call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,tmp_oper)
			call deallocate_block(new_oper1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<sys_len-1)
			if(j<sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,systruns(j))
				endif
				call change_sys_oper_dia(mid_oper,j,sys_oper,sys_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==sys_len-1)
			else if(j==sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,sys_oper)
				else
					call update_trun_dia(tmp_oper,sys_oper,systruns(sys_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_sys_block(oper_dia(i,j),sys_oper,sys_bsm(sys_len),wave)
			oper_dia(j,i)=dconjg(oper_dia(i,j))

			!Save to disk
			open(10,file="oper_dia_tmp.dat",position='append')
			write(10,111) i,j,real(oper_dia(i,j))
			write(10,111) j,i,real(oper_dia(j,i))
			close(10)
			
			call deallocate_block(sys_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)
end subroutine Get_sys_oper_cor_dia

!====================================================================
!Get all density correlcations in environment block using memory
!====================================================================
subroutine Get_env_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,env_oper


	!<1>: For general information
	env_len=wave%env_len

	!Get g(i,j) for (i<j)
	do j=2,env_len-1
		call update_site_dia(st_oper,env_bsm(j),new_oper2)

		!Update envopers(env_len-1)
		if(j==env_len-1) then
			call deallocate_block(envopers(j))
			if(j<=trun_idx) then
				call update_site_dia(st_oper,env_bsm(j),envopers(j))
			else
				call update_site_dia(st_oper,env_bsm(j),mid_oper)
				call update_trun_dia(mid_oper,envopers(j),envtruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_dia(st_oper,env_bsm(i),mid_oper)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_dia(mid_oper,tmp_oper,envtruns(i))
					endif
					call update_block_dia(tmp_oper,env_bsm(j),new_oper1)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_dia(envopers(i),env_bsm(j),new_oper1)
				endif
			else
				call update_block_dia(envopers(i),env_bsm(j),new_oper1)
			endif

			!Update envopers(i)
			call deallocate_block(envopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,envopers(i))
			else
				call update_trun_dia(new_oper1,envopers(i),envtruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
			call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,tmp_oper)
			call deallocate_block(new_oper1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<env_len-1)
			if(j<env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,envtruns(j))
				endif
				call change_env_oper_dia(mid_oper,j,env_oper,env_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==env_len-1)
			else if(j==env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,env_oper)
				else
					call update_trun_dia(tmp_oper,env_oper,envtruns(env_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_env_block(oper_dia(Num_site-i+1,Num_site-j+1),env_oper,env_bsm(env_len),wave)
			oper_dia(Num_site-j+1,Num_site-i+1)=dconjg(oper_dia(Num_site-i+1,Num_site-j+1))

			!Save to disk
			open(10,file="oper_dia_tmp.dat",position='append')
			write(10,111) Num_site-i+1,Num_site-j+1,real(oper_dia(Num_site-i+1,Num_site-j+1))
			write(10,111) Num_site-j+1,Num_site-i+1,real(oper_dia(Num_site-j+1,Num_site-i+1))
			close(10)
			
			call deallocate_block(env_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_env_oper_cor_dia
!==========================================================================
!Get the correlation functions between arbitrary two sites: (idx1<=idx2)
!==========================================================================
subroutine Get_oper_cor_bond_dia(value,oper1,idx1,oper2,idx2,trun_idx,wave)
	use pubdata
	implicit none

	integer,intent(in) :: idx1,idx2,trun_idx
	type(Total_Block),intent(in) :: oper1,oper2
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,j
	complex*16 :: value_sd,value_sz
	integer :: idx_one,idx_two,one,two,sys_len,env_len
	type(Total_Block) :: st_oper1,st_oper2,tmp_oper
	type(Total_Block) :: mid_oper,oper,new_oper
	type(Total_Block) :: sys_sz,sys_sd,env_sz,env_sd
	logical :: Sys_Flag,Env_Flag


	!<1>: For general information
	sys_len=wave%sys_len
	env_len=wave%env_len

	!Get proper operator for specific site
	if(idx1<=idx2) then
		idx_one=idx1
		idx_two=idx2
		call block_transfer(oper1,st_oper1)
		call block_transfer(oper2,st_oper2)
	else !For (idx1>idx2)
		idx_one=idx2
		idx_two=idx1
		call block_transfer(oper1,st_oper2)
		call block_transfer(oper2,st_oper1)
	endif

	value=0.0d0
	!Get truncation flag
	Sys_Flag=.true.
	if(trun_idx==(sys_len-1)) then
		Sys_Flag=.false.
	endif
	Env_Flag=.true.
	if(trun_idx==(env_len-1)) then
		Env_Flag=.false.
	endif


	!<2>: For (idx_one<sys_len)
	if(idx_one<sys_len) then
		!<2-1>: For (idx_two<sys_len)
		if(idx_two<sys_len) then
			call Get_sys_oper_dia(st_oper1,idx_one,tmp_oper,idx_two,trun_idx,.false.)
			call Get_sys_oper_dia(st_oper2,idx_two,mid_oper,idx_two,trun_idx,.false.)
			call block_mul_block_dia(tmp_oper,'N',mid_oper,'N',cone,new_oper)
			call deallocate_block(tmp_oper)
			call deallocate_block(mid_oper)

			!Change to (sys_len-1)'s configuration
			if(idx_two<=trun_idx) then
				call block_transfer(new_oper,tmp_oper)
			else
				call update_trun_dia(new_oper,tmp_oper,systruns(idx_two))
			endif
			call deallocate_block(new_oper)
			if(idx_two<sys_len-1) then
				call Change_sys_oper_dia(tmp_oper,idx_two,new_oper,sys_len-1,trun_idx,Sys_Flag)
				call deallocate_block(tmp_oper)
				call block_transfer(new_oper,tmp_oper)
				call deallocate_block(new_oper)
			endif

			call measure_sys_block(value,tmp_oper,sys_bsm(sys_len),wave)
			call deallocate_block(tmp_oper)

		!<2-2>: For (idx_two>=sys_len)
		else if(idx_two>=sys_len) then
			call Get_sys_oper_dia(st_oper1,idx_one,new_oper,sys_len-1,trun_idx,Sys_Flag)

			!For (idx_two=sys_len)
			if(idx_two==sys_len) then
				call sys_block_site_cor_dia(value,new_oper,st_oper2,sys_bsm(sys_len),wave)

			!For (idx_two=Num_site-env_len+1)
			else if(idx_two==Num_site-env_len+1) then
				call sys_block_env_site_cor_dia(value,new_oper,st_oper2,sys_bsm(sys_len),env_bsm(env_len),wave)

			!For (idx_two>Num_site-env_len+1)
			else if(idx_two>Num_site-env_len+1) then
				call Get_env_oper_dia(st_oper2,Num_site-idx_two+1,tmp_oper,env_len-1,trun_idx,Env_Flag)
				call sys_block_env_block_cor_dia(value,new_oper,tmp_oper,sys_bsm(sys_len),env_bsm(env_len),wave)
				call deallocate_block(tmp_oper)
			endif
			call deallocate_block(new_oper)
		endif


	!<3>: For (idx_one=sys_len)
	else if(idx_one==sys_len) then
		!For (idx_two=Num_site-env_len+1)
		if(idx_two==Num_site-env_len+1) then
			call sys_site_env_site_cor_dia(value,st_oper1,st_oper2,sys_bsm(sys_len),env_bsm(env_len),wave)

		!For (idx_two>Num_site-env_len+1)
		else if(idx_two>Num_site-env_len+1) then
			call Get_env_oper_dia(st_oper2,Num_site-idx_two+1,new_oper,env_len-1,trun_idx,Env_Flag)
			call sys_site_env_block_cor_dia(value,st_oper1,new_oper,sys_bsm(sys_len),env_bsm(env_len),wave)
			call deallocate_block(new_oper)
		endif


	!<4>: For (idx_one>sys_len)
	else if(idx_one>sys_len) then
		one=Num_site-idx_two+1
		two=Num_site-idx_one+1

		call block_transfer(st_oper1,tmp_oper)
		call deallocate_block(st_oper1)
		call block_transfer(st_oper2,st_oper1)
		call deallocate_block(st_oper2)
		call block_transfer(tmp_oper,st_oper2)
		call deallocate_block(tmp_oper)

		!For (two<env_len)
		if(two<env_len) then
			call Get_env_oper_dia(st_oper1,one,tmp_oper,two,trun_idx,.false.)
			call Get_env_oper_dia(st_oper2,two,mid_oper,two,trun_idx,.false.)
			call block_mul_block_dia(tmp_oper,'N',mid_oper,'N',cone,new_oper)
			call deallocate_block(tmp_oper)
			call deallocate_block(mid_oper)

			!Change to (env_len-1)'s configuration
			if(two<=trun_idx) then
				call block_transfer(new_oper,tmp_oper)
			else
				call update_trun_dia(new_oper,tmp_oper,envtruns(two))
			endif
			call deallocate_block(new_oper)

			if(two<env_len-1) then
				call Change_env_oper_dia(tmp_oper,two,new_oper,env_len-1,trun_idx,Env_Flag)
				call deallocate_block(tmp_oper)
				call block_transfer(new_oper,tmp_oper)
				call deallocate_block(new_oper)
			endif

			call measure_env_block(value,tmp_oper,env_bsm(env_len),wave)
			call deallocate_block(tmp_oper)

		!For (two=env_len)
		else if(two==env_len) then
			call Get_env_oper_dia(st_oper1,one,new_oper,env_len-1,trun_idx,Env_Flag)
			call env_block_site_cor_dia(value,new_oper,st_oper2,env_bsm(env_len),wave)
			call deallocate_block(new_oper)
		endif
	endif

	call deallocate_block(st_oper1)
	call deallocate_block(st_oper2)

end subroutine Get_oper_cor_bond_dia



!==========================================================
!Measure <n_i> and <n_i^2> in system block
!==========================================================
!Get <n_i> in the system block
subroutine measure_sys_block(value,sys_bl,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_num_up,bl_num_down,spos,sdim
	complex*16,allocatable :: mid(:,:)
        double complex coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					bl_num_up=sys_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=sys_bs%sub(k)%sub(x)%bl_num_down
					spos=sys_bs%sub(k)%sub(x)%spos
					sdim=sys_bs%sub(k)%sub(x)%sdim
                 coef1=1.0d0/dsqrt(1.0d0+bl_num_down)


					do y=1,sys_bl%num
						if(sys_bl%sub(y)%num_up==bl_num_up) then
						if(sys_bl%sub(y)%num_down==bl_num_down) then
							allocate(mid(sdim,env_dim))

                                                                if(realcode)then
                                   call DGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(y)%mat,sdim&
                                                                        &,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,0.0d0,mid,sdim)


                                                                else
							call ZGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(y)%mat,sdim&
									&,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,czero,mid,sdim)
                                                                endif
							do z=1,env_dim
								value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,z),mid(1:sdim,z))
							end do

							deallocate(mid)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_sys_block
!Get <n_i> in the system block
subroutine measure_sys_site(value,sys_st,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: st_num_up,st_num_down,spos,sdim
	complex*16 :: yita
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					st_num_up=sys_bs%sub(k)%sub(x)%st_num_up
					st_num_down=sys_bs%sub(k)%sub(x)%st_num_down
					spos=sys_bs%sub(k)%sub(x)%spos
					sdim=sys_bs%sub(k)%sub(x)%sdim

					do y=1,sys_st%num
						if(sys_st%sub(y)%num_up==st_num_up) then
						if(sys_st%sub(y)%num_down==st_num_down) then

                 coef1=1.0d0/dsqrt(1.0d0+st_num_down)
							yita=0.0d0
							do z=1,env_dim
								yita=yita+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,z)&
													&,wave%sub(i)%vec(spos+1:spos+sdim,z))
							end do

							value=value+yita*sys_st%sub(y)%mat(1,1)*coef1
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do
end subroutine measure_sys_site


!==========================================================
!Measure <n_i> and <n_i^2> in environment block
!==========================================================

!==========================================================
!Measure <n_i> and <n_i^2> in environment block
!==========================================================
!Get <n_i> in the environment block
subroutine measure_env_block(value,env_bl,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: env_num_up,env_num_down,sys_dim
	integer :: bl_num_up,bl_num_down,spos,sdim
	complex*16,allocatable :: mid(:,:)
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					bl_num_up=env_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=env_bs%sub(k)%sub(x)%bl_num_down
					spos=env_bs%sub(k)%sub(x)%spos
					sdim=env_bs%sub(k)%sub(x)%sdim

					do y=1,env_bl%num
						if(env_bl%sub(y)%num_up==bl_num_up) then
						if(env_bl%sub(y)%num_down==bl_num_down) then
							allocate(mid(sys_dim,sdim))
         coef1=1.0d0/dsqrt(1.0d0+bl_num_down)



                if(realcode)then
							call DGEMM('N','N',sys_dim,sdim,sdim,coef1&
									&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
									&,env_bl%sub(y)%mat,sdim,czero,mid,sys_dim)

                else
							call ZGEMM('N','N',sys_dim,sdim,sdim,coef1&
									&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
									&,env_bl%sub(y)%mat,sdim,czero,mid,sys_dim)
                        endif
							do z=1,sdim
								value=value+dot_product(mid(1:sys_dim,z),wave%sub(i)%vec(1:sys_dim,spos+z))
							end do

							deallocate(mid)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_env_block
!Get <n_i> in the environment block
subroutine measure_env_site(value,env_st,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: env_num_up,env_num_down,sys_dim
	integer :: st_num_up,st_num_down,spos,sdim
	complex*16 :: yita
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					st_num_up=env_bs%sub(k)%sub(x)%st_num_up
					st_num_down=env_bs%sub(k)%sub(x)%st_num_down
					spos=env_bs%sub(k)%sub(x)%spos
					sdim=env_bs%sub(k)%sub(x)%sdim

					do y=1,env_st%num
						if(env_st%sub(y)%num_up==st_num_up) then
						if(env_st%sub(y)%num_down==st_num_down) then

        
                                                                coef1=1.0d0/dsqrt(1.0d0+st_num_down)

							yita=0.0d0
							do z=1,sdim
								yita=yita+dot_product(wave%sub(i)%vec(1:sys_dim,spos+z)&
													&,wave%sub(i)%vec(1:sys_dim,spos+z))
							end do

							value=value+yita*env_st%sub(y)%mat(1,1)*coef1
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_env_site




!Note: value=<psi|n_1^+.n_2|psi> in system block
subroutine sys_block_site_cor_dia(value,sys_bl,sys_st,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Total_Basis),intent(in) :: sys_bs
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: value

	integer :: i,k,x,y,z
	logical :: bl_flag,st_flag
	integer :: bl_id,st_id,spos,sdim
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down
	double complex,allocatable :: mid(:,:)
	double complex :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim
					
					bl_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==bl_num_up) then
						if(sys_bl%sub(k)%num_down==bl_num_down) then
							bl_flag=.true.
							bl_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef=sys_st%sub(st_id)%mat(1,1)
                           coef=coef*dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))

						allocate(mid(sdim,env_dim))
                        if(realcode)then
						call DGEMM('N','N',sdim,env_dim,sdim,coef,sys_bl%sub(bl_id)%mat&
								&,sdim,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
								&,sdim,0.0d0,mid,sdim)
                                else



						call zGEMM('N','N',sdim,env_dim,sdim,coef,sys_bl%sub(bl_id)%mat&
								&,sdim,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
								&,sdim,czero,mid,sdim)


                                endif

						do k=1,env_dim
							value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,k),mid(1:sdim,k))
						end do
						deallocate(mid)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_cor_dia



!Note: value=<psi|n_1^+.n_2|psi> in environment block
subroutine env_block_site_cor_dia(value,env_bl,env_st,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: env_bl,env_st
	type(Total_Basis),intent(in) :: env_bs
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	logical :: bl_flag,st_flag
	integer :: bl_id,st_id,spos,sdim
	integer :: env_num_up,env_num_down,sys_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down
	complex*16,allocatable :: mid(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					bl_num_up=env_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=env_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=env_bs%sub(x)%sub(y)%st_num_up
					st_num_down=env_bs%sub(x)%sub(y)%st_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim
					
					bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==bl_num_up) then
						if(env_bl%sub(k)%num_down==bl_num_down) then
							bl_flag=.true.
							bl_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef=env_st%sub(st_id)%mat(1,1)
                           coef=coef*dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))

						allocate(mid(sys_dim,sdim))

                                                        if(realcode)then
						call DGEMM('N','N',sys_dim,sdim,sdim,coef&
								&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
								&,env_bl%sub(bl_id)%mat,sdim,0.0d0,mid,sys_dim)
                                                                else
						call ZGEMM('N','N',sys_dim,sdim,sdim,coef&
								&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
								&,env_bl%sub(bl_id)%mat,sdim,czero,mid,sys_dim)
                                                        endif
						do k=1,sdim
							value=value+dot_product(mid(1:sys_dim,k),wave%sub(i)%vec(1:sys_dim,spos+k))
						end do
						deallocate(mid)
					endif
				end do
			endif
			endif
		end do
	end do
end subroutine env_block_site_cor_dia




!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_env_block_cor_dia(value,sys_bl,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16, allocatable :: mat1(:,:),mat2(:,:)
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

         coef1=1.0d0/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_bl_num_down))
                                                        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,cone&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                                        else 
									call ZGEMM('N','N',sdim,edim,edim,cone&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                
                                                        endif

									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))*coef1
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_block_cor_dia



!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_env_site_cor_dia(value,sys_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									allocate(mat(sdim,edim))
									coef=env_st%sub(env_id)%mat(1,1)

                                  coef=coef/dsqrt((1.0d0+sys_bl_num_down)*(1.0d0+env_st_num_down))
                                                                                if(realcode)then
									call DGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat,sdim)


                                                                                        else

									call ZGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat,sdim)

                                                                                endif

									do z=1,edim
										value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),mat(1:sdim,z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_site_cor_dia
!Note: value=<psi|n_s^+.n_e|psi>


subroutine sys_site_env_block_cor_dia(value,sys_st,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)
									allocate(mat(sdim,edim))
                                                                        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat,sdim)


                        else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat,sdim)

                        endif

									do z=1,edim
										value=value+dot_product(mat(1:sdim,z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_block_cor_dia

!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_site_env_site_cor_dia(value,sys_bl,sys_st,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_st
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,st_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	double complex,allocatable :: mat(:,:)
	double complex :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat(sdim,edim))
									coef=sys_st%sub(st_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)
                                !!  coef=coef/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_st_num_down))

                        if(realcode)then

									call DGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat,sdim)
                        else

			call ZGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat,sdim)


                        endif

									do z=1,edim
										value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),mat(1:sdim,z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_site_cor_dia


!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_site_env_block_site_cor_dia(value,sys_st,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										st_flag=.true.
										st_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.env_flag.and.st_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)*env_st%sub(st_id)%mat(1,1)
                                coef=sys_st%sub(st_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)

									allocate(mat(sdim,edim))

        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat,sdim)

                                else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat,sdim)

                                                                        endif
									do z=1,edim
										value=value+dot_product(mat(1:sdim,z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do
end subroutine sys_site_env_block_site_cor_dia

!Note: value=<psi|n_s.n_e|psi>



!Note: value=<psi|n_s.n_e|psi>
subroutine sys_site_env_site_cor_dia(value,sys_st,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16 :: coef,tmp_gl

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)

									tmp_gl=0.0d0
									do z=1,edim
										tmp_gl=tmp_gl+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									value=value+tmp_gl*coef
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_site_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_site_env_block_cor_dia(value,sys_bl,sys_st,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=sys_st%sub(st_id)%mat(1,1)
                                                if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                                        else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)

                                                endif
									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_block_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_env_block_site_cor_dia(value,sys_bl,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										st_flag=.true.
										st_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=env_st%sub(st_id)%mat(1,1)

                                if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                        else
        
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                        endif

									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_block_site_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_site_env_block_site_cor_dia(value,sys_bl,sys_st,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,sys_st_flag,env_st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,sys_st_id
	integer :: env_st_num_up,env_st_num_down,env_st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					sys_st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_st_flag=.true.
							sys_st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_st_flag=.true.
										env_st_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 104
									endif
									endif
								end do
								104 continue

								if(sys_flag.and.sys_st_flag.and.env_st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
                                                        if(realcode)then 
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)

                                                                else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                                                endif
									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_block_site_cor_dia


!===========================================================
!Get structure factor for honeycomb Lattice
!Primary vector: e_1=(1,0), e_2=(-1/2,sqrt(3)/2)
!with e_y=(e_1-e_2)/2, e_x=(2*e_1+e_2)/3, e_z=(e_1+2*e_2)/3
!===========================================================
subroutine Get_structure_factor(oper_cor,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	real(8),intent(in) :: oper_cor(Num_site,Num_site)

	integer :: i,j,x,y,dx1,dy1,dx2,dy2
	integer :: x1p,x2p,kxn,kyn,xidx,yidx
	real(8) :: mkx(0:Nx),mky(0:Ny),Stri,Strinew,Spi0
	complex(8) :: nk(0:Nx,0:Ny),nk_mid

	kxn=Nx
	kyn=Ny

	mkx=0.0d0
	do i=0,kxn
		mkx(i)=2.0d0*pi*i/kxn
	end do

	mky=0.0d0
	do i=0,kyn
		mky(i)=2.0d0*pi*i/kyn
	end do

	!Get structure factor
	nk=0.0d0
	do i=0,kxn
	do j=0,kyn
		do x=1,Num_site
			nk(i,j)=nk(i,j)+dcmplx(oper_cor(x,x),0.0d0)
		end do

		nk_mid=0.0d0
		do x=1,Num_site
		do y=1,Num_site
			if(x/=y) then
				dx1=Lattice(1,x)
				dy1=Lattice(2,x)
				dx2=Lattice(1,y)
				dy2=Lattice(2,y)
				nk_mid=nk_mid+cdexp(dcmplx(0.0d0,mkx(i)*(dx1-dx2)+mky(j)*(dy1-dy2)))*dcmplx(oper_cor(x,y),0.0d0)
			endif
		end do
		end do

		nk(i,j)=nk(i,j)+nk_mid
	end do
	end do
	nk=nk/(Num_site)

	!Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do j=0,kyn
		do i=0,kxn
			write(10,111) mkx(i),mky(j),dreal(nk(i,j))
		end do
	end do
	write(10,*)
	close(10)
111 format(F16.12,1X,F16.12,1X,F16.12)

end subroutine Get_structure_factor

subroutine Get_structure_factor_cut(oper_cor,Filename,NxCut)
	use pubdata
	implicit none

	character(len=30) :: Filename
	real(8),intent(in) :: oper_cor(Num_site,Num_site)
    integer,intent(in) :: NxCut

	integer :: i,j,x,y,dx1,dy1,dx2,dy2
	integer :: x1p,x2p,kxn,kyn,xidx,yidx
	real(8) :: mkx(0:(Nx-NxCut*2)),mky(0:Ny),Stri,Strinew,Spi0
	complex(8) :: nk(0:(Nx-NxCut*2),0:Ny),nk_mid

	kxn=Nx-NxCut*2
	kyn=Ny

	mkx=0.0d0
	do i=0,kxn
		mkx(i)=2.0d0*pi*i/kxn
	end do

	mky=0.0d0
	do i=0,kyn
		mky(i)=2.0d0*pi*i/kyn
	end do

	!Get structure factor
	nk=0.0d0
	do i=0,kxn
	do j=0,kyn
		do x=NxCut*Ny+1,Num_site-NxCut*Ny
			nk(i,j)=nk(i,j)+dcmplx(oper_cor(x,x),0.0d0)
		end do

		nk_mid=0.0d0
		do x=NxCut*Ny+1,Num_site-NxCut*Ny
		do y=NxCut*Ny+1,Num_site-NxCut*Ny
			if(x/=y) then
				dx1=Lattice(1,x)
				dy1=Lattice(2,x)
				dx2=Lattice(1,y)
				dy2=Lattice(2,y)
				nk_mid=nk_mid+cdexp(dcmplx(0.0d0,mkx(i)*(dx1-dx2)+mky(j)*(dy1-dy2)))*dcmplx(oper_cor(x,y),0.0d0)
			endif
		end do
		end do

		nk(i,j)=nk(i,j)+nk_mid
	end do
	end do
	nk=nk/(Num_site-(NxCut*2)*Ny)

	!Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do j=0,kyn
		do i=0,kxn
			write(10,111) mkx(i)/PI,mky(j)/PI,dreal(nk(i,j))
		end do
	end do
	write(10,*)
	close(10)
    open(11,file="-PI_to_PI",position='append')
    	    do j=0,kyn
		    do i=0,kxn
                  if(mkx(i)<=1.0001)then
			      write(11,111) mkx(i)/PI,mky(j)/PI,dreal(nk(i,j))
                  endif
                  if(mkx(i)>=0.99999)then
			      write(11,111) (mkx(i)/PI)-2,mky(j)/PI,dreal(nk(i,j))
                  endif
		    end do
	       end do
    close(11)
111 format(F16.12,1X,F16.12,1X,F16.12)

end subroutine Get_structure_factor_cut
!==================================================================================
!<1>: Get bond-bond correlation function (idx1<idx2): (n_i.n_j)
!==================================================================================
subroutine Get_Operator_Bond_Sys(st_oper,idx1,idx2,new_oper,new_idx,trun_idx,truns)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2,sys_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_Operator_Bond_Sys"
		return
	endif

	!<1>: For bond operator: (n_i.n_j)
	call Get_sys_oper_dia(st_oper,idx1,new_oper1,idx2,trun_idx,.false.)
	call Get_sys_oper_dia(st_oper,idx2,new_oper2,idx2,trun_idx,.false.)
	call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,new_oper)
	call deallocate_block(new_oper1)
	call deallocate_block(new_oper2)

	!<2>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,sys_oper)
	else
		call update_trun_dia(new_oper,sys_oper,systruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_sys_oper_dia(sys_oper,idx2,new_oper,new_idx,trun_idx,Truns)
	call deallocate_block(sys_oper)

end subroutine Get_Operator_Bond_Sys


!==================================================================================
!<1>: Get bond-bond correlation function (idx1<idx2): (n_i.n_j)
!==================================================================================
subroutine Get_Operator_Bond_Env(st_oper,idx1,idx2,new_oper,new_idx,trun_idx,truns)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2,env_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_Operator_Bond_Env"
		return
	endif

	!<1>: For bond operator: (n_i.n_j)
	call Get_env_oper_dia(st_oper,idx1,new_oper1,idx2,trun_idx,.false.)
	call Get_env_oper_dia(st_oper,idx2,new_oper2,idx2,trun_idx,.false.)
	call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,new_oper)
	call deallocate_block(new_oper1)
	call deallocate_block(new_oper2)

	!<2>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,env_oper)
	else
		call update_trun_dia(new_oper,env_oper,envtruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_env_oper_dia(env_oper,idx2,new_oper,new_idx,trun_idx,Truns)
	call deallocate_block(env_oper)

end subroutine Get_Operator_Bond_Env
