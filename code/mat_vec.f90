
subroutine new_wave_overlap(sys_bs,  sys_bs_new,  env_bs, env_bs_new,sys_ov,env_ov,ground,ground_new, label, total_up, total_down,aa)
        use pubdata

	implicit none
	type(Total_Basis),intent(in) :: sys_bs, sys_bs_new !! for invec and outvec, resp.
	type(Total_Basis),intent(in) :: env_bs, env_bs_new !! for invec and outvec, resp.
	type(Total_Block),intent(in) :: sys_ov, env_ov
	type(super_basis) :: super,  super1
	type(Wavefunction) :: ground, ground_new, rground,  rground_new   !!! ground1 is new
	integer::i,k,x,y,sys_qn,env_dim,bl_qn,spos,sdim,sdim1
	integer::spos1, x1, y1, j, label, label1,total_up, total_down
	double precision  aa
	double precision, external :: wave_overlap

        label1=num_site-label-2
         call get_super_basis(sys_bs,  env_bs_new, super1, total_up, total_down)
           call allocate_wave(super1,rground_new)! new-wave changed sys
           call allocate_wave(super1,rground)! old-wave changed env
       call mat_diag_vec1(sys_ov,sys_bs_new,sys_bs,ground_new,rground_new,cone) !! new_wave to old sys
       call vec_mat_diag1(env_ov,env_bs,env_bs_new,ground,rground, cone)  !! old WF to new env
        aa=wave_overlap(rground_new,rground)
        call deallocate_wave(rground)
        call deallocate_wave(rground_new)
        call deallocate_super_basis(super1)
end subroutine new_wave_overlap

!Notice : outvec=alpha*sys_bl*invec,  outvec may use different basis than invec
subroutine mat_diag_vec1(sys_bl,sys_bs,sys_bs1,invec,outvec,alpha)
	use pubdata
	implicit none
	double complex,intent(in) :: alpha
	type(Total_Basis),intent(in) :: sys_bs, sys_bs1 !! for invec and outvec, resp.
	type(Total_Block),intent(in) :: sys_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec
	integer::i,k,x,y,sys_qn(2),env_dim,bl_qn(2),spos,sdim,sdim1
	integer::spos1, x1, y1, j

	do i=1,invec%num
		sys_qn(1)=invec%sub(i)%sys_num_up
		sys_qn(2)=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim
        do j=1,outvec%num
                if(sys_qn(1)==outvec%sub(j)%sys_num_up)then ! match the qn
                if(sys_qn(2)==outvec%sub(j)%sys_num_down)then ! match the qn
		if(env_dim.ne.outvec%sub(j)%env_dim)write(*,*)'wave r wrong'
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_qn(1)) then
			if(sys_bs%sub(x)%new_num_down==sys_qn(2)) then
				do y=1,sys_bs%sub(x)%num
				bl_qn(1)=sys_bs%sub(x)%sub(y)%bl_num_up
				bl_qn(2)=sys_bs%sub(x)%sub(y)%bl_num_down
				spos=sys_bs%sub(x)%sub(y)%spos
				sdim=sys_bs%sub(x)%sub(y)%sdim
                do x1=1, sys_bs1%num
                        if(sys_bs1%sub(x1)%new_num_up==sys_qn(1))then
                        if(sys_bs1%sub(x1)%new_num_down==sys_qn(2))then
!!! in each basis sub,  spos is related to site qn. they may be matched differently 
                do y1=1,sys_bs1%sub(x1)%num
                if(bl_qn(1)==sys_bs1%sub(x1)%sub(y1)%bl_num_up)then 
                if(bl_qn(2)==sys_bs1%sub(x1)%sub(y1)%bl_num_down)then 
				spos1=sys_bs1%sub(x1)%sub(y1)%spos

				do k=1,sys_bl%num
				if(sys_bl%sub(k)%num_up==bl_qn(1)) then
				if(sys_bl%sub(k)%num_down==bl_qn(2)) then
        sdim1=size(sys_bl%sub(k)%mat,1)
if(sys_bs1%sub(x1)%sub(y1)%sdim.ne.sdim1)write(*,*)'sys_bs1 wrong',sys_bs1%sub(x1)%sub(y1)%sdim, sdim1
 call zgemm('n','n',sdim1,env_dim,sdim,alpha,sys_bl%sub(k)%mat,sdim1&
        &,invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,cone&
        &,outvec%sub(j)%vec(spos1+1:spos1+sdim1,1:env_dim),sdim1)

							goto 101
						endif
						endif
					end do
					101 continue
                endif
                endif
                enddo  !y1
		endif
		endif
				end do
                        enddo ! do x1
			endif
			endif
		end do !! do x
        endif
	endif
	end do
	enddo

end subroutine mat_diag_vec1


!Notice : outvec=alpha*env_bl*invec
subroutine vec_mat_diag1(env_bl,env_bs,env_bs1,invec,outvec,alpha)
	use pubdata
	implicit none
	
	double complex,intent(in) :: alpha
	type(Total_Basis),intent(in) :: env_bs,env_bs1 ! for invec and outvec
	type(Total_Block),intent(in) :: env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec
	integer :: i,k,x,y,env_qn(2),sys_dim,bl_qn(2),spos,sdim, sdim1
	integer ::j, x1, y1, spos1
	do i=1,invec%num
		env_qn(1)=invec%sub(i)%env_num_up
		env_qn(2)=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim
	do j=1, outvec%num
	if(env_qn(1)==outvec%sub(j)%env_num_up)then
	if(env_qn(2)==outvec%sub(j)%env_num_down)then
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_qn(1)) then
			if(env_bs%sub(x)%new_num_down==env_qn(2)) then
				do y=1,env_bs%sub(x)%num
					bl_qn(1)=env_bs%sub(x)%sub(y)%bl_num_up
					bl_qn(2)=env_bs%sub(x)%sub(y)%bl_num_down
					spos=env_bs%sub(x)%sub(y)%spos
				!	sdim=env_bs%sub(x)%sub(y)%sdim
                do x1=1,env_bs1%num
			if(env_bs1%sub(x1)%new_num_up==env_qn(1)) then
			if(env_bs1%sub(x1)%new_num_down==env_qn(2)) then
				do y1=1,env_bs1%sub(x1)%num
	if(bl_qn(1)==env_bs1%sub(x1)%sub(y1)%bl_num_up)then
	if(bl_qn(2)==env_bs1%sub(x1)%sub(y1)%bl_num_down)then
					spos1=env_bs1%sub(x1)%sub(y1)%spos
					sdim1=env_bs1%sub(x1)%sub(y1)%sdim

			do k=1,env_bl%num
				if(env_bl%sub(k)%num_up==bl_qn(1)) then
				if(env_bl%sub(k)%num_down==bl_qn(2)) then
        sdim=size(env_bl%sub(k)%mat,1)!older qn of env, matching to old vec
        if(sdim.ne.env_bs%sub(x)%sub(y)%sdim)write(*,*)'env_bs1 wrong'
        call zgemm('n','n',sys_dim,sdim1,sdim,alpha&
                                &,invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
                                                &,sys_dim,env_bl%sub(k)%mat,sdim,cone&
                                &,outvec%sub(j)%vec(1:sys_dim,spos1+1:spos1+sdim1),sys_dim)
							goto 101
						endif
						endif
					end do
                                endif
                                endif
                                        enddo
                                endif   
                                endif   
                                        enddo
					101 continue
				end do
			endif
			endif
		end do
	endif
	endif
	end do
	enddo
end subroutine vec_mat_diag1

!==========For the diagonal part (diag¡ÁIn) or (In¡Ádiag)=========
!Notice : outvec=coef*sys_bl*invec
subroutine block_vec_dia(sys_bl,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none
	
	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,sys_num_up,sys_num_down
	integer :: env_dim,bl_num_up,bl_num_down,spos,sdim
        double complex :: coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

        coef1=coef/dsqrt(1.0d0+bl_num_down)

					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==bl_num_up) then
						if(sys_bl%sub(k)%num_down==bl_num_down) then
	
        call ZGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(k)%mat,sdim&
                                                                        &,invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,cone&
                                                                        &,outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim)

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

end subroutine block_vec_dia


!Notice : outvec=coef*env_bl*invec
subroutine vec_block_dia(env_bl,env_bs,invec,outvec,coef)
	use pubdata
	implicit none
	
	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,env_num_up,env_num_down
	integer :: sys_dim,bl_num_up,bl_num_down,spos,sdim
        double complex :: coef1

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					bl_num_up=env_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=env_bs%sub(x)%sub(y)%bl_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim

        coef1=coef/dsqrt(1.0d0+bl_num_down)
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==bl_num_up) then
						if(env_bl%sub(k)%num_down==bl_num_down) then
	
   call ZGEMM('N','C',sys_dim,sdim,sdim,coef1&
                                                                        &,invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
                                                                        &,sys_dim,env_bl%sub(k)%mat,sdim,cone&
                                                                        &,outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim)

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

end subroutine vec_block_dia


!Notice : outvec=alpha*sys_st*invec
subroutine site_vec_dia(sys_st,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: st_flag
	integer :: i,k,x,y,st_id,spos,sdim,env_dim
	integer :: sys_num_up,sys_num_down,st_num_up,st_num_down
	double complex :: coef_st, coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(st_flag) then
						coef_st=coef*sys_st%sub(st_id)%mat(1,1)
                                                coef_st=coef_st/dsqrt(1.d0+st_num_down)
						outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
							&= outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
							&+ coef_st*invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine site_vec_dia


!Notice : outvec=alpha*env_st*invec
subroutine vec_site_dia(env_st,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: st_flag
	integer :: i,k,x,y,st_id,spos,sdim,sys_dim
	integer :: env_num_up,env_num_down,st_num_up,st_num_down
	double complex  :: coef_st

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					st_num_up=env_bs%sub(x)%sub(y)%st_num_up
					st_num_down=env_bs%sub(x)%sub(y)%st_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					if(st_flag) then
						coef_st=coef*env_st%sub(st_id)%mat(1,1)
                                                coef_st=coef_st/dsqrt(1.d0+st_num_down)
						outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
							&= outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
							&+ coef_st*invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine vec_site_dia


!==========For the diagonal part (diag¡Ádiag)========================
!Notice : outvec=coef*sys_bl*sys_st*invec
subroutine block_site_vec_dia(sys_bl,sys_st,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bl_flag,st_flag
	integer :: i,k,x,y,spos,sdim
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_id,bl_num_up,bl_num_down
	integer :: st_id,st_num_up,st_num_down
	double complex :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

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
							bl_id=k
							bl_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef_tmp=coef*sys_st%sub(st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))
              call ZGEMM('N','N',sdim,env_dim,sdim,coef_tmp,sys_bl%sub(bl_id)%mat&
                                                                &,sdim,invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim&
                                                                &,cone,outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim)

					endif
				end do
			endif
			endif
		end do
	end do

end subroutine block_site_vec_dia


!Notice : outvec=coef*env_bl*env_st*invec
subroutine vec_block_site_dia(env_bl,env_st,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bl_flag,st_flag
	integer :: i,k,x,y,spos,sdim,sys_dim
	integer :: env_num_up,env_num_down
	integer :: bl_id,bl_num_up,bl_num_down
	integer :: st_id,st_num_up,st_num_down
	double complex :: coef_tmp

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

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
							bl_id=k
							bl_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef_tmp=coef*env_st%sub(st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))
                          call ZGEMM('N','C',sys_dim,sdim,sdim,coef_tmp&
                                                                &,invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
                                                                &,sys_dim,env_bl%sub(bl_id)%mat,sdim,cone&
                                                                &,outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim)


					endif
				end do
			endif
			endif
		end do
	end do

end subroutine vec_block_site_dia


!Notice : outvec=coef*sys_st*env_st*invec
subroutine site_vec_site_dia(sys_st,env_st,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_st_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_st_id,env_st_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_st_num_up,sys_st_num_down,env_st_num_up,env_st_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
        double complex :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_st_num_up=sys_bs%sub(sys_bs_id)%sub(x)%st_num_up
				sys_st_num_down=sys_bs%sub(sys_bs_id)%sub(x)%st_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_st_flag=.false.
				do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
						sys_st_id=k
						sys_st_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_st_num_up=env_bs%sub(env_bs_id)%sub(y)%st_num_up
					env_st_num_down=env_bs%sub(env_bs_id)%sub(y)%st_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==env_st_num_up) then
						if(env_st%sub(k)%num_down==env_st_num_down) then
							env_st_id=k
							env_st_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_st_flag.and.env_st_flag) then
						coef_tmp=coef*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+sys_st_num_down)*(1.0d0+env_st_num_down))
						outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
							&=outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
							&+coef_tmp*invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)
					endif
				end do
			end do
		endif
	end do

end subroutine site_vec_site_dia


!Notice : outvec=coef*sys_bl*env_st*invec
subroutine block_vec_site_dia(sys_bl,env_st,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_st_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_bl_id,env_st_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,env_st_num_up,env_st_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
        double complex :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_bl_num_up=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_up
				sys_bl_num_down=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_bl_flag=.false.
				do k=1,sys_bl%num
					if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
					if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
						sys_bl_id=k
						sys_bl_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_st_num_up=env_bs%sub(env_bs_id)%sub(y)%st_num_up
					env_st_num_down=env_bs%sub(env_bs_id)%sub(y)%st_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==env_st_num_up) then
						if(env_st%sub(k)%num_down==env_st_num_down) then
							env_st_id=k
							env_st_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_bl_flag.and.env_st_flag) then
						coef_tmp=coef*env_st%sub(env_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+env_st_num_down)*(1.0d0+sys_bl_num_down))

                   call ZGEMM('N','N',sys_dim,env_dim,sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat,sys_dim&
                                                                &,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim&
                                                                &,cone,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim)


					endif
				end do
			end do
		endif
	end do

end subroutine block_vec_site_dia


!Notice : outvec=coef*sys_st*env_bl*invec
subroutine site_vec_block_dia(sys_st,env_bl,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_bl_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_st_id,env_bl_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_st_num_up,sys_st_num_down,env_bl_num_up,env_bl_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	double complex :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_st_num_up=sys_bs%sub(sys_bs_id)%sub(x)%st_num_up
				sys_st_num_down=sys_bs%sub(sys_bs_id)%sub(x)%st_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_st_flag=.false.
				do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
						sys_st_id=k
						sys_st_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_bl_num_up=env_bs%sub(env_bs_id)%sub(y)%bl_num_up
					env_bl_num_down=env_bs%sub(env_bs_id)%sub(y)%bl_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==env_bl_num_up) then
						if(env_bl%sub(k)%num_down==env_bl_num_down) then
							env_bl_id=k
							env_bl_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_st_flag.and.env_bl_flag) then
						coef_tmp=coef*sys_st%sub(sys_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_st_num_down))

      call ZGEMM('N','C',sys_dim,env_dim,env_dim,coef_tmp&
                                                                &,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
                                                                &,sys_dim,env_bl%sub(env_bl_id)%mat,env_dim,cone&
                                                                &,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim)
                                        endif

				end do
			end do
		endif
	end do

end subroutine site_vec_block_dia


!Notice : outvec=coef*sys_bl*env_bl*invec
subroutine block_vec_block_dia(sys_bl,env_bl,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_bl_id,env_bl_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,env_bl_num_up,env_bl_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	double complex,allocatable :: mat(:,:)
        double complex coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_bl_num_up=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_up
				sys_bl_num_down=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_bl_flag=.false.
				do k=1,sys_bl%num
					if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
					if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
						sys_bl_id=k
						sys_bl_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_bl_num_up=env_bs%sub(env_bs_id)%sub(y)%bl_num_up
					env_bl_num_down=env_bs%sub(env_bs_id)%sub(y)%bl_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==env_bl_num_up) then
						if(env_bl%sub(k)%num_down==env_bl_num_down) then
							env_bl_id=k
							env_bl_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_bl_flag.and.env_bl_flag) then
        if(sys_dim*env_dim==0)write(*,*)'dimen wrong'
        if(sys_dim*env_dim==0)stop
						allocate(mat(sys_dim,env_dim))
						coef1=coef/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_bl_num_down))

  call ZGEMM('N','N',sys_dim,env_dim,sys_dim,cone,sys_bl%sub(sys_bl_id)%mat,sys_dim&
                                                                &,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
                                                                &,sys_dim,czero,mat,sys_dim)
                                                call ZGEMM('N','C',sys_dim,env_dim,env_dim,coef1,mat,sys_dim,env_bl%sub(env_bl_id)%mat&
                                                                &,env_dim,cone,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim&
                                                                &,env_pos+1:env_pos+env_dim),sys_dim)

						deallocate(mat)
					endif
				end do
			end do
		endif
	end do

end subroutine block_vec_block_dia


!=============For the non-diagonal part=============================================
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!===================================================================================
!Outvec=coef*(C^+_i(sys_bl)*C_j(sys_st)+h.c.)*invec Matrix-vector multiplication
subroutine block_site_vec_ndia(sys_bl,FlagBL,sys_st,FlagST,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagBL,FlagST
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bs_flag,bl_flag,st_flag
	integer :: i,k,x,y,up_dif,down_dif
	integer :: sys_num_up,sys_num_down,bl_num,st_num,tot_num
	integer :: bl_num_up,bl_num_down,new_bl_num_up,new_bl_num_down
	integer :: st_num_up,st_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,old_pos,old_dim,new_pos,new_dim,env_dim
	double complex :: coef_tmp
	real(8) :: Signs,SignBL,SignST
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=sys_st%up_dif).or.(down_dif1/=sys_st%down_dif)) then
		write(*,*) "Quantum number error in block_site_vec!"
		stop
	endif

	!Start multiplication mat_vec
	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					bl_num_up=sys_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=sys_bs%sub(k)%sub(x)%bl_num_down
					st_num_up=sys_bs%sub(k)%sub(x)%st_num_up
					st_num_down=sys_bs%sub(k)%sub(x)%st_num_down
					old_pos=sys_bs%sub(k)%sub(x)%spos
					old_dim=sys_bs%sub(k)%sub(x)%sdim
					
        do down_dif=-down_dif1, down_dif1, su  !! for sys_bl only
        do st_down_dif=-down_dif1, down_dif1, su  !! for sys_st

						bl_flag=.false.
						do y=1,sys_bl%num
							if(sys_bl%sub(y)%num_up==bl_num_up) then
							if(sys_bl%sub(y)%num_down==bl_num_down) then
							if(sys_bl%sub(y)%down_dif==down_dif) then
								bl_id=y
								bl_flag=.true.
								goto 102
							endif
							endif
							endif
						end do
						102 continue

					new_st_num_up=st_num_up-up_dif
					new_st_num_down=st_num_down-st_down_dif

						st_flag=.false.

						do y=1,sys_st%num
							if(sys_st%sub(y)%num_up==new_st_num_up) then
							if(sys_st%sub(y)%num_down==new_st_num_down) then
							if(sys_st%sub(y)%down_dif==st_down_dif) then
								st_id=y
								st_flag=.true.
								goto 103
							endif
							endif
							endif
						end do
						103 continue

					!For C^+_i(block).C_j(site) case
					new_bl_num_up=bl_num_up+up_dif
					new_bl_num_down=bl_num_down+down_dif
                
                if(bl_flag.and.st_flag)then

					bs_flag=.false.
					do y=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(sys_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
						if(sys_bs%sub(k)%sub(y)%st_num_up==new_st_num_up) then
						if(sys_bs%sub(k)%sub(y)%st_num_down==new_st_num_down) then
							bs_flag=.true.
							new_pos=sys_bs%sub(k)%sub(y)%spos
							new_dim=sys_bs%sub(k)%sub(y)%sdim
							goto 101
						endif
						endif
						endif
						endif
					end do
					101 continue

                                                if(bs_flag)then
							!<a>: For sys_block operator
							SignBL=cone
							
							!<b>: For sys_site operator
							SignST=czero
							if(FlagST=='B') then !For Boson
								SignST=cone
							else if(FlagST=='F') then !For Fermion
								bl_num=bl_num_up !!!+bl_num_down
								tot_num=bl_num
								if(mod(tot_num,2)==0) then
									SignST=cone
								else
									SignST=-cone
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down !! diag new_sys_st_num
        j6=sys_num_down
              if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef1=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef1
        endif

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)/dsqrt(1.0d0+down_dif1)

				Signs=SignBL*SignST*coef1
			coef_tmp=Signs* coef* sys_st%sub(st_id)%mat(1,1)
	call ZGEMM('N','N',new_dim,env_dim,old_dim,coef_tmp,sys_bl%sub(bl_id)%mat,new_dim&
		&,invec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,cone&
				&,outvec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim)

if(up_dif.ne.0)then
        coef_tmp=dconjg(coef_tmp)
                call ZGEMM('C','N',old_dim,env_dim,new_dim,coef_tmp,sys_bl%sub(bl_id)%mat,new_dim&
                        &,invec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim,cone&
                                &,outvec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim)

        endif



						endif
					endif

                        endif
        go to 110

					!================================================
					!For C^+_j(site).C_i(block) case
					new_bl_num_up=bl_num_up-up_dif
					new_bl_num_down=bl_num_down-down_dif
					new_st_num_up=st_num_up+up_dif
					new_st_num_down=st_num_down+down_dif

					bs_flag=.false.
					do y=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(sys_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
							bs_flag=.true.
							new_pos=sys_bs%sub(k)%sub(y)%spos
							new_dim=sys_bs%sub(k)%sub(y)%sdim
							goto 104
						endif
						endif
					end do
					104 continue

					if(bs_flag) then
						bl_flag=.false.
						do y=1,sys_bl%num
							if(sys_bl%sub(y)%num_up==new_bl_num_up) then
							if(sys_bl%sub(y)%num_down==new_bl_num_down) then
								bl_id=y
								bl_flag=.true.
								goto 105
							endif
							endif
						end do
						105 continue

						st_flag=.false.
						do y=1,sys_st%num
							if(sys_st%sub(y)%num_up==st_num_up) then
							if(sys_st%sub(y)%num_down==st_num_down) then
								st_id=y
								st_flag=.true.
								goto 106
							endif
							endif
						end do
						106 continue
						
						if(bl_flag.and.st_flag) then
							!<a>: For sys_block operator
							SignBL=cone

							!<b>: For sys_site operator
							SignST=czero
							if(FlagST=='B') then !For Boson
								SignST=cone
							else if(FlagST=='F') then !For Fermion
								bl_num=new_bl_num_up !!!+new_bl_num_down
								tot_num=bl_num
								if(mod(tot_num,2)==0) then
									SignST=cone
								else
									SignST=-cone
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down !! diag new_sys_st_num
        j6=sys_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)
							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef* sys_st%sub(st_id)%mat(1,1)
							call ZGEMM('C','N',new_dim,env_dim,old_dim,coef_tmp,sys_bl%sub(bl_id)%mat,old_dim&
									&,invec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,cone&
									&,outvec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim)
						endif
						endif
					endif
110     continue
                        
				end do
                        enddo
                        enddo
			endif
			endif
		end do
	end do


end subroutine block_site_vec_ndia


subroutine vec_block_site_ndia(env_bl,FlagBL,env_st,FlagST,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagBL,FlagST
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bs_flag,bl_flag,st_flag
	integer :: i,k,x,y,up_dif,down_dif,env_num_up,env_num_down
	integer :: sys_num_up,sys_num_down,sys_num,bl_num,st_num,tot_num
	integer :: bl_num_up,bl_num_down,new_bl_num_up,new_bl_num_down
	integer :: st_num_up,st_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,old_pos,old_dim,new_pos,new_dim,sys_dim
	double complex :: coef_tmp
        double precision ::Signs,SignBL,SignST
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=env_bl%up_dif
	down_dif1=env_bl%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in vec_block_site!"
		stop
	endif

	!Start multiplication mat_vec
	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					bl_num_up=env_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=env_bs%sub(k)%sub(x)%bl_num_down
					st_num_up=env_bs%sub(k)%sub(x)%st_num_up
					st_num_down=env_bs%sub(k)%sub(x)%st_num_down
					old_pos=env_bs%sub(k)%sub(x)%spos
					old_dim=env_bs%sub(k)%sub(x)%sdim

					!================================================
					!For C^+_i(block).C_j(site) case

                        do bl_down_dif=-down_dif1, down_dif1, su

						bl_flag=.false.
						do y=1,env_bl%num
							if(env_bl%sub(y)%num_up==bl_num_up) then
							if(env_bl%sub(y)%num_down==bl_num_down) then
							if(env_bl%sub(y)%down_dif==bl_down_dif) then
								bl_id=y
								bl_flag=.true.
								goto 102
							endif
							endif
							endif
						end do
						102 continue

                        do st_down_dif=-down_dif1, down_dif1, su
						st_flag=.false.
					new_st_num_up=st_num_up-up_dif
					new_st_num_down=st_num_down-st_down_dif

						do y=1,env_st%num
							if(env_st%sub(y)%num_up==new_st_num_up) then
							if(env_st%sub(y)%num_down==new_st_num_down) then
							if(env_st%sub(y)%down_dif==st_down_dif) then
								st_id=y
								st_flag=.true.
								goto 103
							endif
							endif
							endif
						end do
						103 continue
					new_bl_num_up=bl_num_up+up_dif
					new_bl_num_down=bl_num_down+bl_down_dif

                if(st_flag.and.bl_flag)then
					bs_flag=.false.
					do y=1,env_bs%sub(k)%num
						if(env_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(env_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
						if(env_bs%sub(k)%sub(y)%st_num_up==new_st_num_up) then
						if(env_bs%sub(k)%sub(y)%st_num_down==new_st_num_down) then
							bs_flag=.true.
							new_pos=env_bs%sub(k)%sub(y)%spos
							new_dim=env_bs%sub(k)%sub(y)%sdim
							goto 101
						endif
						endif
						endif
						endif
					end do
					101 continue

					if(bs_flag) then
							!<a>: For env_block operator
							SignBL=czero
							if(FlagBL=='B') then !For Boson
								SignBL=cone
							else if(FlagBL=='F') then !For Fermion
								sys_num=sys_num_up !!+sys_num_down
								st_num=new_st_num_up !+new_st_num_down
								tot_num=sys_num+st_num
								if(mod(tot_num,2)==0) then
									SignBL=cone
								else
									SignBL=-cone
								endif
							endif

							!<b>: For env_site operator
							SignST=czero
							if(FlagST=='B') then !For Boson
								SignST=cone
							else if(FlagST=='F') then !For Fermion
								sys_num=sys_num_up !!+sys_num_down
								tot_num=sys_num
								if(mod(tot_num,2)==0) then
									SignST=cone
								else
									SignST=-cone
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down !! diag new_sys_st_num
        j6=env_num_down

      if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef1=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef1
        endif

                if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)/dsqrt(1.0d0+down_dif1)


							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef *env_st%sub(st_id)%mat(1,1)
				call ZGEMM('N','C',sys_dim,new_dim,old_dim,coef_tmp&
								&,invec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim&
									&,env_bl%sub(bl_id)%mat,new_dim,cone&
									&,outvec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim)

   if(up_dif.ne.0)then
        coef_tmp=dconjg(coef_tmp)
                        call ZGEMM('N','N',sys_dim,old_dim,new_dim,coef_tmp&
                                     &,invec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim&
                                     &,env_bl%sub(bl_id)%mat,new_dim,cone&
                                           &,outvec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim)
                endif
        




						endif
					endif
                                        endif
                go to 110


					!For C_j^+(site).C_i(block) case
					new_bl_num_up=bl_num_up-up_dif
					new_bl_num_down=bl_num_down-down_dif
					new_st_num_up=st_num_up+up_dif
					new_st_num_down=st_num_down+down_dif

					bs_flag=.false.
					do y=1,env_bs%sub(k)%num
						if(env_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(env_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
							bs_flag=.true.
							new_pos=env_bs%sub(k)%sub(y)%spos
							new_dim=env_bs%sub(k)%sub(y)%sdim
							goto 104
						endif
						endif
					end do
					104 continue

					if(bs_flag) then
						bl_flag=.false.
						do y=1,env_bl%num
							if(env_bl%sub(y)%num_up==new_bl_num_up) then
							if(env_bl%sub(y)%num_down==new_bl_num_down) then
								bl_id=y
								bl_flag=.true.
								goto 105
							endif
							endif
						end do
						105 continue

						st_flag=.false.
						do y=1,env_st%num
							if(env_st%sub(y)%num_up==st_num_up) then
							if(env_st%sub(y)%num_down==st_num_down) then
								st_id=y
								st_flag=.true.
								goto 106
							endif
							endif
						end do
						106 continue
						
						if(bl_flag.and.st_flag) then
							!<a>: For env_block operator
							SignBL=czero
							if(FlagBL=='B') then !For Boson
								SignBL=cone
							else if(FlagBL=='F') then !For Fermion
								sys_num=sys_num_up!!+sys_num_down
								st_num=st_num_up!!+st_num_down
								tot_num=sys_num+st_num
								if(mod(tot_num,2)==0) then
									SignBL=cone
								else
									SignBL=-cone
								endif
							endif

							!<b>: For env_site operator
							SignST=czero
							if(FlagST=='B') then !For Boson
								SignST=cone
							else if(FlagST=='F') then !For Fermion
								sys_num=sys_num_up!!+sys_num_down
								tot_num=sys_num
								if(mod(tot_num,2)==0) then
									SignST=cone
								else
									SignST=-cone
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down !! diag new_sys_st_num
        j6=env_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)
							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef* env_st%sub(st_id)%mat(1,1)
							call ZGEMM('N','N',sys_dim,new_dim,old_dim,coef_tmp&
									&,invec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim&
									&,env_bl%sub(bl_id)%mat,old_dim,cone&
									&,outvec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim)
						endif
					endif
					endif
110             continue
				end do
				end do
				end do
			endif
			endif
		end do
	end do

end subroutine vec_block_site_ndia

subroutine block_vec_block_ndia(sys_bl,FlagSL,env_bl,FlagEL,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagEL
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_bl_id,sys_num
	integer :: sys_st_num_down, env_st_num_up,env_st_num_down,env_st_num,sys_bl_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
        integer :: sys_st_num_up
	integer :: old_sys_pos,old_env_pos,old_sys_dim,old_env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim,new_env_dim
	double complex :: coef_tmp
        double precision::SignSL,SignEL,Signs
	double complex,allocatable :: mat(:,:)
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j24
        integer j5, j6,j55, j66,  env_down_dif, bl_down_dif, bl_down_dif1
        integer ji1,ji2,ji3,ji4, ji5, ji6,ji7,ji8,ji9, jj, j7,j8,j9,j77,j88,j99
        integer x1, y1
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=env_bl%up_dif).or.(down_dif1/=env_bl%down_dif)) then
		write(*,*) "Quantum number error in block_vec_block!"
		stop
	endif


        !if(env_bl%len==9) call block_to_disk1(env_bl, 101)
        !if(env_bl%len==9) stop

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim

                do down_dif=-down_dif1, down_dif1, su
							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(k)%down_dif==down_dif) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                        
				new_sys_bl_num_up=sys_bl_num_up+up_dif
				new_sys_bl_num_down=sys_bl_num_down+down_dif


                do env_down_dif=-down_dif1, down_dif1, su
				new_env_bl_num_up=env_bl_num_up-up_dif
				new_env_bl_num_down=env_bl_num_down-env_down_dif

							env_bl_flag=.false.
							do k=1,env_bl%num
					if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
					if(env_bl%sub(k)%down_dif==env_down_dif) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue


					if(sys_bl_flag.and.env_bl_flag) then
        do new_sys_num_down= abs(sys_num_down-down_dif1), sys_num_down+down_dif1, su
						new_sys_num_up=sys_num_up+up_dif
						new_env_num_up=env_num_up-up_dif
					new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
					do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
						do l=1,sys_bs%sub(k)%num
					if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
					if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_down==sys_st_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_up==sys_st_num_up) then
								sys_bs_flag=.true.
					new_sys_pos=sys_bs%sub(k)%sub(l)%spos
						new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
                                        x1=k
											goto 103
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
							if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
							if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_down==env_st_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_up==env_st_num_up) then 
											env_bs_flag=.true.
									new_env_pos=env_bs%sub(k)%sub(l)%spos
									new_env_dim=env_bs%sub(k)%sub(l)%sdim
                                                        y1=k
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=cone

										!<b>: For env_block operator
										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
								sys_num=sys_num_up!+sys_num_down
								env_st_num=env_st_num_up !!+env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif
        coef1=czero
        
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down

        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down      
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down

      if(iw6j3(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j3(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6) 
        iw6j3(j1,j2,j3-j1,j4,j5,j6)=1
        w6j3(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 105
      if(iw6j3(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j3(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j3(j11,j22,j33-j11,j44,j55,j66)=1
        w6j3(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif

        if(coef11.eq.0.0)go to 105
        coef1=coef1*coef11
        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        if(mod((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2, 2)==1)coef1=-coef1

										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1

        go to 1051

										allocate(mat(new_sys_dim,old_env_dim))
						call ZGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												&,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												&,czero,mat,new_sys_dim)
						call ZGEMM('N','C',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
									&,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)

1051    continue

                    allocate(mat(new_sys_dim,old_env_dim))
                                                       call ZGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
                                                                   &,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
                                                                          &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                                               &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
                                                                                    &,czero,mat,new_sys_dim)
                                            call ZGEMM('N','N',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
                                                    &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
                                                          &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                                                 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
                                                                                deallocate(mat)
  if(up_dif.ne.0)then
        coef_tmp=dconjg(coef_tmp)
                allocate(mat(old_sys_dim,new_env_dim))
                        call ZGEMM('C','N',old_sys_dim,new_env_dim,new_sys_dim,cone&
                                        &,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
                                                &,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim&
                                                                &,czero,mat,old_sys_dim)
                                      call ZGEMM('N','C',old_sys_dim,old_env_dim,new_env_dim,coef_tmp&
                                        &,mat,old_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
                                                &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim)
                                                                                deallocate(mat)
        endif




										goto 105
									endif
									endif
                                                                endif
								end do
							endif
							105 continue
        go to 110

							
							new_sys_bl_num_up=sys_bl_num_up-up_dif
							new_sys_bl_num_down=sys_bl_num_down-down_dif
							new_env_bl_num_up=env_bl_num_up+up_dif
							new_env_bl_num_down=env_bl_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down
        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down    
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down
        coef1=dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2)
										SignSL=cone

										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
								sys_num=new_sys_num_up 
								env_st_num=env_st_num_up 
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif

										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1
										allocate(mat(new_sys_dim,old_env_dim))
										call ZGEMM('C','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												 &,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												 &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												 &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												 &,czero,mat,new_sys_dim)
										call ZGEMM('N','T',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												 &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
												 &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
                                enddo
12      continue
                endif
						end do
						end do
                                enddo
                                enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine block_vec_block_ndia



!Outvec=coef*(C^+_i(sys_bl)*C_j(env_bl)+h.c.)*invec Matrix-wavefunction multiplication
subroutine block_vec_block_ndia01(sys_bl,FlagSL,env_bl,FlagEL,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagEL
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_bl_id,sys_num
	integer :: sys_st_num_down, env_st_num_up,env_st_num_down,env_st_num,sys_bl_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
        integer :: sys_st_num_up
	integer :: old_sys_pos,old_env_pos,old_sys_dim,old_env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim,new_env_dim
	double complex :: coef_tmp
        double precision::SignSL,SignEL,Signs
	double complex,allocatable :: mat(:,:)
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, env_down_dif, bl_down_dif, bl_down_dif1
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=env_bl%up_dif).or.(down_dif1/=env_bl%down_dif)) then
		write(*,*) "Quantum number error in block_vec_block!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim

                do down_dif=-down_dif1, down_dif1, su
							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(k)%down_dif==down_dif) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                        
				new_sys_bl_num_up=sys_bl_num_up+up_dif
				new_sys_bl_num_down=sys_bl_num_down+down_dif

                do env_down_dif=-down_dif1, down_dif1, su
				new_env_bl_num_up=env_bl_num_up-up_dif
				new_env_bl_num_down=env_bl_num_down-env_down_dif

							env_bl_flag=.false.
							do k=1,env_bl%num
					if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
					if(env_bl%sub(k)%down_dif==env_down_dif) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue

                if(sys_bl_flag.and.env_bl_flag)then


        do new_sys_num_down= abs(sys_num_down-down_dif1), sys_num_down+down_dif1, su
							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
							new_env_num_down=new_sys_num_down
                !! total spin=0
							!<1>: For C_i^+(sys_block).C_j(env_block) case

							sys_bs_flag=.false.
					do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
						do l=1,sys_bs%sub(k)%num
					if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
					if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_down==sys_st_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_up==sys_st_num_up) then
								sys_bs_flag=.true.
					new_sys_pos=sys_bs%sub(k)%sub(l)%spos
						new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 103
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
							if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
							if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_down==env_st_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_up==env_st_num_up) then 
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=cone

										!<b>: For env_block operator
										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
								sys_num=sys_num_up!+sys_num_down
								env_st_num=env_st_num_up !!+env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif


        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down

        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down     
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down

      if(iw6j3(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j3(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6) 
        iw6j3(j1,j2,j3-j1,j4,j5,j6)=1
        w6j3(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j3(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j3(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j3(j11,j22,j33-j11,j44,j55,j66)=1
        w6j3(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11      

                if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        if(mod((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2,2).ne.0)coef1=-coef1


										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1

										allocate(mat(new_sys_dim,old_env_dim))
										call ZGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												&,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												&,czero,mat,new_sys_dim)
										call ZGEMM('N','N',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												&,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 105
									endif
									endif
                                                                endif
								end do
							endif
							105 continue
        go to 110

							
							!<2>: For C^+_j(env_block)).C_i(sys_block) case
							new_sys_bl_num_up=sys_bl_num_up-up_dif
							new_sys_bl_num_down=sys_bl_num_down-down_dif
							new_env_bl_num_up=env_bl_num_up+up_dif
							new_env_bl_num_down=env_bl_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							!!new_sys_num_down=sys_num_down-down_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down !! diag
        j6=new_sys_bl_num_down
        j11=new_env_num_down  !! new-env
        j22=down_dif1
        j33=env_num_down       !! new-env-
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down
        coef1=dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2)
										SignSL=cone

										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
								sys_num=new_sys_num_up 
								env_st_num=env_st_num_up 
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif

										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1
										allocate(mat(new_sys_dim,old_env_dim))
										call ZGEMM('T','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												 &,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												 &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												 &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												 &,czero,mat,new_sys_dim)
										call ZGEMM('N','T',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												 &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
												 &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
                                enddo
12      continue
                                endif
						end do
						end do
                                enddo
                                enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do


end subroutine block_vec_block_ndia01


!Outvec=coef*(C^+_i(sys_bl)*C_j(env_st)+h.c.)*invec Matrix-wavefunction multiplication
subroutine block_vec_site_ndia(sys_bl,FlagSL,env_st,FlagET,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagET
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_st_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_st_id,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	integer :: old_sys_pos,old_env_pos,old_sys_dim,env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim
        integer env_bl_num_down    !! jenv=jenv'

        integer:: sys_st_num_down
	double complex :: coef_tmp
	real(8) :: SignSL,SignET,Signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in block_vec_site!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
						sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							env_dim=env_bs%sub(y)%sub(n)%sdim
							
							!<1>: For C_i^+(sys_block).C_j(env_site) case

                        do bl_down_dif=-sys_bl%down_dif, sys_bl%down_dif, su
							sys_bl_flag=.false.
				new_sys_bl_num_up=sys_bl_num_up+up_dif
				new_sys_bl_num_down=sys_bl_num_down+bl_down_dif
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(k)%down_dif==bl_down_dif) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                do st_down_dif=-env_st%down_dif, env_st%down_dif, su
							env_st_flag=.false.
				new_env_st_num_up=env_st_num_up-up_dif
				new_env_st_num_down=env_st_num_down-st_down_dif

							do k=1,env_st%num
								if(env_st%sub(k)%num_up==new_env_st_num_up) then
								if(env_st%sub(k)%num_down==new_env_st_num_down) then
								if(env_st%sub(k)%down_dif==st_down_dif) then
									env_st_id=k
									env_st_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue

							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
							!new_sys_num_down=sys_num_down+down_dif
							!new_env_num_down=env_num_down-down_dif
  if(sys_bl_flag.and.env_st_flag) then
			do new_sys_num_down=abs(sys_num_down-sys_bl%down_dif), sys_num_down+sys_bl%down_dif, su
                           new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
						if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
						if(sys_bs%sub(k)%sub(l)%st_num_down==sys_st_num_down) then
									sys_bs_flag=.true.
								new_sys_pos=sys_bs%sub(k)%sub(l)%spos
								new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 103
										endif
                                                                                endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 104
										endif
                                                                                endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_bl_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=cone

										!<b>: For env_site operator
										SignET=czero
										if(FlagET=='B') then !For Boson
											SignET=cone
										else if(FlagET=='F') then !For Fermion
								sys_num=sys_num_up !!+sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=cone
											else
												SignET=-cone
											endif
										endif
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down !! diag
        j6=new_sys_bl_num_down
        j11=new_env_num_down  
        j22=down_dif1        
        j33=env_num_down     
        j44= env_st_num_down 
        j55=env_bl_num_down    
        j66=new_env_st_num_down  

      if(iw6j3(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j3(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j3(j1,j2,j3-j1,j4,j5,j6)=1
        w6j3(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j2(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j2(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j2(j11,j22,j33-j11,j44,j55,j66)=1
        w6j2(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11  

                if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.d0+j11)*(1.d0+j33)/(1.0d0+down_dif1))
        if(mod((j6+j5+j3+j2)/2+(j55++j44+j11+j22)/2,2).ne.0) coef1=-coef1
							Signs=SignSL*SignET*coef1
					coef_tmp=Signs* coef* env_st%sub(env_st_id)%mat(1,1)


										call ZGEMM('N','N',new_sys_dim,env_dim,old_sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat&
												&,new_sys_dim,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+env_dim),old_sys_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+env_dim),new_sys_dim)

        if(up_dif.ne.0)then
        coef_tmp=dconjg(coef_tmp)
  call ZGEMM('C','N',old_sys_dim,env_dim,new_sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat&
                                             &,new_sys_dim,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                                              &,new_env_pos+1:new_env_pos+env_dim),new_sys_dim,cone&
                                                                    &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                                            &,old_env_pos+1:old_env_pos+env_dim),old_sys_dim)

                                endif
                                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
							

                go to 110
							new_sys_bl_num_up=sys_bl_num_up-up_dif
							new_sys_bl_num_down=sys_bl_num_down-down_dif
							new_env_st_num_up=env_st_num_up+up_dif
							new_env_st_num_down=env_st_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==env_st_num_up) then
								if(env_st%sub(k)%num_down==env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_bl_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										SignSL=cone

										SignET=czero
										if(FlagET=='B') then 
											SignET=cone
										else if(FlagET=='F') then 
							sys_num=new_sys_num_up  
							tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=cone
											else
												SignET=-cone
											endif
										endif
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down
        j11=new_env_num_down  
        j22=down_dif1        
        j33=env_num_down    
        j44= env_st_num_down  
        j55=env_bl_num_down  
        j66=new_env_st_num_down  
        coef1=dsqrt((1.d0+j11)*(1.d0+j33)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j6+j5+j3+j2)/2+(j55++j44+j11+j22)/2)
										Signs=SignSL*SignET
										coef_tmp=Signs* coef* env_st%sub(env_st_id)%mat(1,1)
										call ZGEMM('T','N',new_sys_dim,env_dim,old_sys_dim,coef_tmp&
												&,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+env_dim),old_sys_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+env_dim),new_sys_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                                endif
						end do
                                        enddo
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

        

end subroutine block_vec_site_ndia


!Outvec=coef*(C^+_i(sys_st)*C_j(env_bl)+h.c.)*invec Matrix-wavefunction multiplication
subroutine site_vec_block_ndia(sys_st,FlagST,env_bl,FlagEL,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagST,FlagEL
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_bl_flag
	integer :: sys_bs_id,env_bs_id,sys_st_id,env_bl_id
	integer :: env_st_num_up,env_st_num_down,env_st_num,tot_num
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num,sys_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
	integer :: old_sys_pos,old_env_pos,old_sys_dim,sys_dim
	integer :: new_sys_pos,new_env_pos,old_env_dim,new_env_dim
	double complex :: coef_tmp
	real(8) :: SignST,SignEL,Signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif
	if((up_dif/=env_bl%up_dif).or.(down_dif1/=env_bl%down_dif)) then
		write(*,*) "Quantum number error in site_vec_block!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							sys_dim=sys_bs%sub(x)%sub(m)%sdim

							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim
							

                do st_down_dif=-sys_st%down_dif, sys_st%down_dif, su
							sys_st_flag=.false.
							new_sys_st_num_up=sys_st_num_up+up_dif
							new_sys_st_num_down=sys_st_num_down+st_down_dif
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==sys_st_num_up) then
								if(sys_st%sub(k)%num_down==sys_st_num_down) then
								if(sys_st%sub(k)%down_dif==st_down_dif) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                do bl_down_dif=-env_bl%down_dif, env_bl%down_dif, su
							env_bl_flag=.false.
							new_env_bl_num_up=env_bl_num_up-up_dif
							new_env_bl_num_down=env_bl_num_down-bl_down_dif
							do k=1,env_bl%num
				if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
						if(env_bl%sub(k)%down_dif==bl_down_dif) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 102
								endif
								endif
                                                                endif
							end do
							102 continue

				new_sys_num_up=sys_num_up+up_dif
				new_env_num_up=env_num_up-up_dif
                 if(sys_st_flag.and.env_bl_flag) then
		do new_sys_num_down=abs(sys_num_down-sys_st%down_dif),sys_num_down+sys_st%down_dif, su
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
							if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
							if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
							if(sys_bs%sub(k)%sub(l)%bl_num_down==sys_bl_num_down) then
							if(sys_bs%sub(k)%sub(l)%bl_num_up==sys_bl_num_up) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 103
										endif
										endif
										endif
                                                        endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
							if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
							if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_up==env_st_num_up) then
							if(env_bs%sub(k)%sub(l)%st_num_down==env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
		if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
		if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down     
        j44= env_bl_num_down
        j55=env_st_num_down 
        j66=new_env_bl_num_down

      if(iw6j2(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j2(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j2(j1,j2,j3-j1,j4,j5,j6)=1
        w6j2(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j3(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j3(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j3(j11,j22,j33-j11,j44,j55,j66)=1
        w6j3(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11       

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(j33+1.d0)/(1.0d0+down_dif1))
        if(mod((j66+j55+j33+j22)/2+(j5+j4+j1+j2)/2,2).ne.0)coef1=-coef1


						SignST=czero
					if(FlagST=='B') then 
								SignST=cone
							else if(FlagST=='F') then 
				sys_bl_num=sys_bl_num_up 
								tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=cone
											else
												SignST=-cone
											endif
										endif

										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
						sys_num=sys_num_up  
						env_st_num=env_st_num_up  
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif

							Signs=SignST*SignEL*coef1
						coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)

         call ZGEMM('N','N',sys_dim,new_env_dim,old_env_dim,coef_tmp&
                            &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
                                       &,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
                                               &,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
                                                    &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
                                                          &,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)


          if(up_dif.ne.0)then
        coef_tmp=dconjg(coef_tmp)
                        call ZGEMM('N','C',sys_dim,old_env_dim,new_env_dim,coef_tmp&
                              &,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
                                               &,new_env_pos+1:new_env_pos+new_env_dim),sys_dim&
                                                      &,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
                                                               &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
                                                                   &,old_env_pos+1:old_env_pos+old_env_dim),sys_dim)
                        endif

                go to 105
	call ZGEMM('N','T',sys_dim,new_env_dim,old_env_dim,coef_tmp&
			&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
				&,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
				&,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
					&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
						&,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)
                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
                go to 110 

							new_sys_st_num_up=sys_st_num_up-up_dif
							new_sys_st_num_down=sys_st_num_down-down_dif
							new_env_bl_num_up=env_bl_num_up+up_dif
							new_env_bl_num_down=env_bl_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==new_sys_st_num_up) then
								if(sys_st%sub(k)%num_down==new_sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										SignST=czero
										if(FlagST=='B') then !For Boson
											SignST=cone
										else if(FlagST=='F') then !For Fermion
							sys_bl_num=sys_bl_num_up 
								tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=cone
											else
												SignST=-cone
											endif
										endif

										SignEL=czero
										if(FlagEL=='B') then !For Boson
											SignEL=cone
										else if(FlagEL=='F') then !For Fermion
							sys_num=new_sys_num_up 
							env_st_num=env_st_num_up 
									tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=cone
											else
												SignEL=-cone
											endif
										endif

        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down       
        j44= env_bl_num_down
        j55=env_st_num_down 
        j66=new_env_bl_num_down

        coef1=dsqrt((1.0d0+j11+1)*(j33+1)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j66+j55+j33+j22)/2+(j5+j4+j1+j2)/2)

										Signs=SignST*SignEL*coef1
										coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)
										call ZGEMM('N','T',sys_dim,new_env_dim,old_env_dim,coef_tmp&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
												&,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                                        endif
						end do
						end do
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine site_vec_block_ndia


!Outvec=coef*(C^+_i(sys_st)*C_j(env_st)+h.c.)*invec Matrix-wavefunction multiplication
subroutine site_vec_site_ndia(sys_st,FlagST,env_st,FlagET,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagST,FlagET
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n
	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_st_flag
	integer :: sys_bs_id,env_bs_id,sys_st_id,env_st_id,up_dif,down_dif
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_bl_num_down,env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	integer :: old_sys_pos,old_env_pos,sys_dim,env_dim,new_sys_pos,new_env_pos
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44
        integer j5, j6,j55, j66, st_down_dif1, st_down_dif2
          real(8),external :: w6js, w3js, w9js
	double complex :: coef_tmp, coef1
	real(8) :: SignST,SignET,Signs, coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in site_vec_site!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
        if(sys_num_down.ne.env_num_down)write(51,*)'wrong1'

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							!For environment part
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							env_dim=env_bs%sub(y)%sub(n)%sdim
							

        do st_down_dif1=-sys_st%down_dif, sys_st%down_dif, su
							sys_st_flag=.false.
							do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
					if(sys_st%sub(k)%down_dif==st_down_dif1) then
									sys_st_id=k
									sys_st_flag=.true.
				new_sys_st_num_up=sys_st_num_up+up_dif
				new_sys_st_num_down=sys_st_num_down+st_down_dif1
									goto 101
								endif
                                                                endif
								endif
							end do
							101 continue

        do st_down_dif2=-env_st%down_dif, env_st%down_dif, su
							env_st_flag=.false.
					new_env_st_num_up=env_st_num_up-up_dif
					new_env_st_num_down=env_st_num_down-st_down_dif2

							do k=1,env_st%num
					if(env_st%sub(k)%num_up==new_env_st_num_up) then
					if(env_st%sub(k)%num_down==new_env_st_num_down) then
					if(env_st%sub(k)%down_dif==st_down_dif2) then
									env_st_id=k
									env_st_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue

				new_sys_num_up=sys_num_up+up_dif
				new_env_num_up=env_num_up-up_dif
				if(sys_st_flag.and.env_st_flag) then
                do new_sys_num_down=abs(sys_num_down-sys_st%down_dif), sys_num_down+sys_st%down_dif, su
							!<1>: For C_i^+(sys_site).C_j(env_site) case
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
					if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
					if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
					if(sys_bs%sub(k)%sub(l)%bl_num_down==sys_bl_num_down) then
								new_sys_pos=sys_bs%sub(k)%sub(l)%spos
									sys_bs_flag=.true.
											goto 103
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
							do l=1,env_bs%sub(k)%num
								if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
								if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
								if(env_bs%sub(k)%sub(l)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 104
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

				if(sys_st_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
						do k=1,outvec%num
			if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
			if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_site operator
										SignST=czero
										if(FlagST=='B') then !For Boson
											SignST=cone
										else if(FlagST=='F') then !For Fermion
							sys_bl_num=sys_bl_num_up  
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=cone
											else
												SignST=-cone
											endif
										endif

										!<b>: For env_block operator
										SignET=czero
										if(FlagET=='B') then !For Boson
											SignET=cone
										else if(FlagET=='F') then !For Fermion
								sys_num=sys_num_up !!+sys_num_down 
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=cone
											else
												SignET=-cone
											endif
										endif

				Signs=SignST*SignET
				coef_tmp=Signs* coef
                
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1        
        j33=env_num_down    
        j44= env_st_num_down 
        j55=env_bl_num_down    
        j66=new_env_st_num_down   


      if(iw6j2(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j2(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j2(j1,j2,j3-j1,j4,j5,j6)=1
        w6j2(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j2(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j2(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j2(j11,j22,j33-j11,j44,j55,j66)=1
        w6j2(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11    

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        if(mod((j55+j44+j11+j22)/2+(j5+j4+j1+j2)/2, 2).ne.0)coef1=-coef1
       coef1=coef_tmp*coef1*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)

		outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
		&=outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
	&+coef1*invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)

        if(up_dif.ne.0)then
        coef1=dconjg(coef1)
 outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)&
                &=outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)&
        &+coef1*invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)

        endif




                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
        go to 110
							
							new_sys_st_num_up=sys_st_num_up-up_dif
							new_sys_st_num_down=sys_st_num_down-down_dif
							new_env_st_num_up=env_st_num_up+up_dif
							new_env_st_num_down=env_st_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==new_sys_st_num_up) then
								if(sys_st%sub(k)%num_down==new_sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==env_st_num_up) then
								if(env_st%sub(k)%num_down==env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue

							if(sys_st_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_site operator
										SignST=czero
										if(FlagST=='B') then !For Boson
											SignST=cone
										else if(FlagST=='F') then !For Fermion
								sys_bl_num=sys_bl_num_up !!+sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=cone
											else
												SignST=-cone
											endif
										endif

										SignET=czero
										if(FlagET=='B') then !For Boson
											SignET=cone
										else if(FlagET=='F') then !For Fermion
						sys_num=new_sys_num_up  
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=cone
											else
												SignET=-cone
											endif
										endif

										Signs=SignST*SignET
										coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
										outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
											&=outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
											&+coef_tmp*invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                        endif
                                                enddo
						end do
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine site_vec_site_ndia




!=============================================================
!For sweep process with full system size
!=============================================================
subroutine model_wave(sys,env,sys_bs,env_bs,invec,outvec)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Model),intent(in) :: sys,env
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,idx,sys_len,env_len
        integer :: ii, is, js, jc, jc1, i1, j, j1,k1,j12, j11,j2, ij
	integer :: sys_xidx,sys_yidx,env_xidx,env_yidx
	double complex :: coefx,coefy,coefsd,coefsz,coefsn
	real(8) :: maghx,maghy,maghz,bg_time,end_time
        double complex alpha, beta
	type(Total_Block) :: sys_st,env_st, tmp_bl, tmp_bl1

	call cpu_time(bg_time)

	!<1>: Initiate outvec
	sys_len=outvec%sys_len
	env_len=outvec%env_len
	do i=1,outvec%num
		outvec%sub(i)%vec=czero
	end do

	!<2>: For block-hamiltonian term
	call block_vec_dia(sys%ham,sys_bs,invec,outvec,cone)
	call vec_block_dia(env%ham,env_bs,invec,outvec,cone)
	!<3>: For all the two-body interaction term
                        is=nposi_lat(sys%len+1)
                        js=nposi_lat(sys%len+env%len+2)
                        ii=neiba(is,js)
 if(hubbard(is).ne.0.0)then
                coefx=hubbard(is)
        call site_vec_site_dia(st_double,st_si,sys_bs,env_bs,invec,outvec,coefx)
                endif

      if(hubbard(js).ne.0.0)then
                coefx=hubbard(js)
        call site_vec_site_dia(st_si,st_double,sys_bs,env_bs,invec,outvec,coefx)
                endif


                 if(ii.ne.0)then
                                jz1=jz(ii,is)
                                jd1=jd(ii, is)
                                    coefx=jt(ii, is) !! hoping in either x or y
                           coefx=jt(ii, is)*pint(is,js)
                                      jn1=jn(ii, is)


                if(coefx.ne.0.0)then
                call block_transfer_trans(st_elec_down, tmp_bl)
call site_vec_site_ndia(st_elec_up,'F',tmp_bl,'F',sys_bs,env_bs,invec,outvec,coefx)
                endif

                if(jd1.ne.0.0)then
	call site_vec_site_ndia(st_sd,'B',st_sd,'B',sys_bs,env_bs,invec,outvec,jd1)
                endif

                if(jn1.ne.0.0)then
	call site_vec_site_dia(num_elec,num_elec,sys_bs,env_bs,invec,outvec,jn1)
                        endif
			endif

                        if(coup_bs(1,1).eq.1)then  !! coefs assorbed
                call block_transfer_trans(st_elec_down, tmp_bl)
	call block_site_vec_ndia(sys%snc(1)%elec_up,'F',tmp_bl,'F',sys_bs,invec,outvec,cone)!coefx)

                if(tjmodel==1)then
	call block_site_vec_ndia(sys%snc(1)%spin_sd,'B',st_sd,'B',sys_bs,invec,outvec,cone)!coefsd)
                endif
	if(v123==1)call block_site_vec_dia(sys%snc(1)%num_sn,num_elec,sys_bs,invec,outvec,cone)!coefsn)
	endif

                        if(coup_bs(1,2).eq.1)then                
        call block_transfer_trans(st_elec_down, tmp_bl)
  call block_vec_site_ndia(sys%snc(2)%elec_up,'F',tmp_bl,'F',sys_bs,env_bs,invec,outvec,cone)!coefy)
                if(tjmodel==1)then
  call block_vec_site_ndia(sys%snc(2)%spin_sd,'B',st_sd,'B',sys_bs,env_bs,invec,outvec,cone)!!coefsd)
                                endif
  if(v123==1)call block_vec_site_dia(sys%snc(2)%num_sn,num_elec,sys_bs,env_bs,invec,outvec,cone)!coefsn)
			endif


                        if(coup_bs(2,1).eq.1)then
        call block_transfer_trans(env%snc(1)%elec_down, tmp_bl1)
     call site_vec_block_ndia(st_elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                if(tjmodel==1)then
                call block_transfer_trans(env%snc(1)%spin_sd, tmp_bl)
          call site_vec_block_ndia(st_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)!coefsd)
                        call deallocate_block(tmp_bl)
                                endif
           if(v123==1)call site_vec_block_dia(num_elec,env%snc(1)%num_sn,sys_bs,env_bs,invec,outvec,cone)!coefsn)
                        call deallocate_block(tmp_bl1)
			endif


                                if(coup_bs(2,2).eq.1)then
        call block_transfer_trans(st_elec_down, tmp_bl)	
call vec_block_site_ndia(env%snc(2)%elec_up,'F',tmp_bl,'F',env_bs,invec,outvec,cone)!1coefx)
                if(tjmodel==1)then
	call vec_block_site_ndia(env%snc(2)%spin_sd,'B',st_sd,'B',env_bs,invec,outvec,cone)!coefsd)
                        endif
	if(v123==1)call vec_block_site_dia(env%snc(2)%num_sn,num_elec,env_bs,invec,outvec,cone)!coefsn)
			endif
                                jc=0
                        do i=1,sys%len+1+env%len
                        if(inb(i).ne.0)then
                                jc=jc+1
                        if(i.le.sys%len)then  !!
                                i1=inb1(i,sys%len)

        call block_transfer_trans(env%sub_s(jc)%elec_down, tmp_bl1)
       call block_vec_block_ndia(sys%sub(i1)%elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                if(tjmodel==1)then
                call block_transfer_trans(env%sub_s(jc)%spin_sd, tmp_bl)
	call block_vec_block_ndia(sys%sub(i1)%spin_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl)
                                endif
	if(v123==1)call block_vec_block_dia(sys%sub(i1)%num_sn,env%sub_s(jc)%num_sn,sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl1)
                                else
                i1=inb1(i-sys%len-1,env%len)
                        jc1=jc-coup_bb(1)
        call block_transfer_trans(env%sub(i1)%elec_down, tmp_bl1)
	call block_vec_block_ndia(sys%sub_s(jc1)%elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                if(tjmodel==1)then
                call block_transfer_trans(env%sub(i1)%spin_sd, tmp_bl)
	call block_vec_block_ndia(sys%sub_s(jc1)%spin_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl)
                                endif
	if(v123==1)call block_vec_block_dia(sys%sub_s(jc1)%num_sn,env%sub(i1)%num_sn,sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl1)
                               endif
                                endif
                                enddo
101     continue

111             continue

                
	call cpu_time(end_time)
end subroutine model_wave


   subroutine hmerge(sys,env,sys_bs,env_bs)
	use pubdata
	implicit none
	type(Total_Basis),intent(in)::sys_bs,env_bs
	type(Total_Model),intent(inout)::sys,env
	type(Total_block)::  tmp_bl
	integer i,k,j,mbsite,ib,ib1,i1,k1,k11,imm, jc, jc1
	integer ia1, ia, syssite, totsite, dmsite 
	real*8 h1
			totsite=(sys%len+env%len+2+nleg-1)/nleg*nleg
                        dmsite=sys%len+env%len+2
			syssite=sys%len+1
		coup_bs(1:2,1:2)=0
		do k=1,neibt
		j=neib(syssite,k) 
		if(nposi(nposi_dm(j)).eq.3)coup_bs(2,1)=1 
		if(nposi(nposi_dm(j)).eq.1)coup_bs(1,1)=1 

		j=neib(nposi_lat(dmsite),k) 
		if(nposi(nposi_dm(j)).eq.3)coup_bs(2,2)=1 
		if(nposi(nposi_dm(j)).eq.1)coup_bs(1,2)=1  
		enddo

				do ia=1,2
				if(coup_bs(1,ia).eq.1)then
        if(allocated(sys%snc(ia)%elec_up%sub))call deallocate_block(sys%snc(ia)%elec_up)
  call block_pass_info(sys%sub(1)%elec_up,sys%snc(ia)%elec_up)
                if(tjmodel==1)then
	if(allocated(sys%snc(ia)%spin_sd%sub))call deallocate_block(sys%snc(ia)%spin_sd)
  call block_pass_info(sys%sub(1)%spin_sd,sys%snc(ia)%spin_sd)
                                endif
	if(allocated(sys%snc(ia)%num_sn%sub))call deallocate_block(sys%snc(ia)%num_sn)
  if(v123==1)call block_pass_info(sys%sub(1)%num_sn,  sys%snc(ia)%num_sn)
					endif

				if(coup_bs(2,ia).eq.1)then
        if(allocated(env%snc(ia)%elec_up%sub))call deallocate_block(env%snc(ia)%elec_up)
        if(allocated(env%snc(ia)%elec_down%sub))call deallocate_block(env%snc(ia)%elec_down)
  call block_pass_info(env%sub(1)%elec_up,env%snc(ia)%elec_up)
  call block_pass_info(env%sub(1)%elec_down,env%snc(ia)%elec_down)
                if(tjmodel==1)then
	if(allocated(env%snc(ia)%spin_sd%sub))call deallocate_block(env%snc(ia)%spin_sd)
  call block_pass_info(env%sub(1)%spin_sd,env%snc(ia)%spin_sd)
					endif
	if(allocated(env%snc(ia)%num_sn%sub))call deallocate_block(env%snc(ia)%num_sn)
  if(v123==1)call block_pass_info(env%sub(1)%num_sn,  env%snc(ia)%num_sn)
					endif
                                        enddo
		
		do ia=1,2 
				j=sys%len+1+(env%len+1)*(ia-1) !! site in DMRG
			do k1=1,neibt
			Ib=neib(nposi_lat(j),k1) !! latt posi.
			if(nposi(nposi_dm(ib)).eq.3)then	
			Ib1=nposi_dm(Ib)-sys%len-1  
			Ib1=inb1(Ib1,env%len)
		if(coup_bs(2,ia).ne.0)then
                jz1=jz(k1, nposi_lat(j))
                jd1=jd(k1, nposi_lat(j))   
                jt1=jt(k1,nposi_lat(j))*pint(nposi_lat(j),ib)
                jn1=jn(k1,nposi_lat(j))

                if(jt1.ne.0.0)then                        
  call block_add_block_two(env%sub(ib1)%elec_up,jt1,env%snc(ia)%elec_up)
  call block_add_block_two(env%sub(ib1)%elec_down,dconjg(jt1),env%snc(ia)%elec_down)
                        endif
                if(tjmodel==1)then
                        if(jd1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%spin_sd,jd1,env%snc(ia)%spin_sd)
                                endif
                                endif
                                if(jn1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%num_sn,jn1,env%snc(ia)%num_sn)
                                        endif
				endif
				endif
			if(nposi(nposi_dm(Ib)).eq.1)then	
                        ib1=nposi_dm(ib)
			Ib1=inb1(Ib1,sys%len)
	if(coup_bs(1,ia).eq.1)then
                jz1=jz(k1, nposi_lat(j))
                jd1=jd(k1, nposi_lat(j))   
                jt1=jt(k1,nposi_lat(j))*pint(nposi_lat(j),ib)
                jn1=jn(k1, nposi_lat(j))
                if(jt1.ne.0.0)then

                jt1=dconjg(jt1)

	call block_add_block_two(sys%sub(ib1)%elec_up,jt1,sys%snc(ia)%elec_up)
                        endif
                if(tjmodel==1)then
                        if(jd1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%spin_sd,jd1,sys%snc(ia)%spin_sd)
                                endif
				endif
                        if(jn1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%num_sn,jn1,sys%snc(ia)%num_sn)
                                endif
				endif
				endif
				enddo
				enddo

		coup_bb(1:2)=0
		k11=0
           do I=1, sys%len   
		k1=0
		inb(I)=0
		do k=1,neibt
		J=neib(I,k)
		if(nposi(nposi_dm(J)).eq.3)then 
		k1=k1+1
		endif
		enddo
		inb(I)=k1
		if(nleg11(env%len).lt.nleg11(sys%len))inb(i)=0
		if(inb(I).ne.0)k11=k11+1
		enddo
        
			coup_bb(1)=k11!
				inb(sys%len+1)=0
				inb(sys%len+1+env%len+1)=0
      			do I=1, env%len   !
			k1=0
			inb(i+sys%len+1)=0
			do k=1,neibt
			J=neib(nposi_lat(I+sys%len+1),k) ! in first block
			if(nposi(nposi_dm(J)).eq.1)then

			if(inb(J).eq.0)then
				k1=k1+1!!
			endif
				endif
			enddo
			inb(I+sys%len+1)=k1
		if(inb(I+sys%len+1).ne.0)k11=k11+1
		enddo
		coup_bb(2)=k11-coup_bb(1)

			if(coup_bb(1).ne.0)then
			do i=1,coup_bb(1)
       if(allocated(env%sub_s(i)%elec_up%sub)) call deallocate_block(env%sub_s(i)%elec_up)
       if(allocated(env%sub_s(i)%elec_down%sub))call deallocate_block(env%sub_s(i)%elec_down)
       if(allocated(env%sub_s(i)%spin_sd%sub))call deallocate_block(env%sub_s(i)%spin_sd)
       if(allocated(env%sub_s(i)%num_sn%sub))call deallocate_block(env%sub_s(i)%num_sn)
	call block_pass_info(env%sub(1)%elec_up,env%sub_s(i)%elec_up)
	call block_pass_info(env%sub(1)%elec_down,env%sub_s(i)%elec_down)
                if(tjmodel==1)then
	call block_pass_info(env%sub(1)%spin_sd,env%sub_s(i)%spin_sd)
                        endif
	if(v123==1)call block_pass_info(env%sub(1)%num_sn,env%sub_s(i)%num_sn)
			enddo
			endif
			if(coup_bb(2).ne.0)then
			do i=1,coup_bb(2)
	if(allocated(sys%sub_s(i)%elec_up%sub))then
                call deallocate_block(sys%sub_s(i)%elec_up)
                if(tjmodel==1)then
                call deallocate_block(sys%sub_s(i)%spin_sd)
                        endif
                if(v123==1)call deallocate_block(sys%sub_s(i)%num_sn)
                        endif
	call block_pass_info(sys%sub(1)%elec_up,sys%sub_s(i)%elec_up)
                if(tjmodel==1)then
	call block_pass_info(sys%sub(1)%spin_sd,sys%sub_s(i)%spin_sd)
                        endif
	if(v123==1)call block_pass_info(sys%sub(1)%num_sn,sys%sub_s(i)%num_sn)
				enddo
				endif

			jc=0
		do I=1,sys%len+env%len+1
		if(inb(I).ne.0)then
		Ia=nposi_lat(i)
		if(I.le.sys%len)I1=inb1(I,sys%len)
		if(I.gt.sys%len)I1=inb1(I-sys%len-1,env%len)
		jc=jc+1
			do k1=1,neibt
			Ib=neib(Ia,k1)
		if(nposi(nposi_dm(Ib)).eq.3.and.I.le.sys%len)then
			Ib1=inb1(nposi_dm(Ib)-sys%len-1,env%len)

                jt1=jt(k1, ia)*pint(ia,ib)
                jd1=jd(k1, ia)  
                jz1=jz(k1, ia)
                jn1=jn(k1, ia)  


        if(jt1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%elec_up,jt1,env%sub_s(jc)%elec_up)
	call block_add_block_two(env%sub(ib1)%elec_down,dconjg(jt1),env%sub_s(jc)%elec_down)
                endif
                if(tjmodel==1)then
        if(jd1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%spin_sd,jd1,env%sub_s(jc)%spin_sd)
                endif
                        endif
        if(jn1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%num_sn,jn1,env%sub_s(jc)%num_sn)
                endif
				endif


			if(nposi(nposi_dm(Ib)).eq.1.and.I.gt.sys%len+1)then	
				if(inb(nposi_dm(ib)).eq.0)then
			Ib1=inb1(Ib,sys%len)
				jc1=jc-coup_bb(1)
                jt1=jt(k1, ia)*pint(ia,ib)
                jd1=jd(k1, ia)   
                jz1=jz(k1, ia)
                jn1=jn(k1, ia)  
                if(jt1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%elec_up,dconjg(jt1),sys%sub_s(jc1)%elec_up)
                        endif

                if(tjmodel==1)then
                if(jd1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%spin_sd,jd1,sys%sub_s(jc1)%spin_sd)
                        endif
			endif
                                if(jn1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%num_sn,jn1,sys%sub_s(jc1)%num_sn)
                                endif
                        endif
			endif
			enddo

	endif
	enddo
   end subroutine hmerge


subroutine block_site_vec_block_site_ndia(sys_bl,FlagSL,sys_st,FlagST,env_bl,FlagEL,env_st,FlagET,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	double complex,intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagST,FlagEL,FlagET
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n
	integer :: sl_up_dif,sl_down_dif,st_up_dif,st_down_dif
	integer :: el_up_dif,el_down_dif,et_up_dif,et_down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag,sys_st_flag,env_st_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_bl_id,sys_st_id,env_st_id
	integer :: sys_st_num_up,sys_st_num_down,sys_st_num,sys_bl_num,sys_num
	integer :: new_sys_st_num_up,new_sys_st_num_down,new_env_st_num_up,new_env_st_num_down
	integer :: env_st_num_up,env_st_num_down,env_st_num,env_bl_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
	integer :: old_sys_pos,old_env_pos,old_sys_dim,old_env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim,new_env_dim
	double complex :: coef_tmp,SignSL,SignST,SignEL,SignET,Signs
	double complex,allocatable :: mat(:,:)

	!Get general information
	sl_up_dif=sys_bl%up_dif
	sl_down_dif=sys_bl%down_dif
	st_up_dif=sys_st%up_dif
	st_down_dif=sys_st%down_dif
	el_up_dif=env_bl%up_dif
	el_down_dif=env_bl%down_dif
	et_up_dif=env_st%up_dif
	et_down_dif=env_st%down_dif
	if(((sl_up_dif-st_up_dif)/=(el_up_dif-et_up_dif)).or.((sl_down_dif-st_down_dif)/=(el_down_dif-et_down_dif))) then
		write(*,*) "Quantum number error in block_site_vec_block_site!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim
							
							!<1>: For C_1^+(sys_bl).C_2(sys_st).C_3^+(env_st).C_4(env_bl) case
							new_sys_bl_num_up=sys_bl_num_up+sl_up_dif
							new_sys_bl_num_down=sys_bl_num_down+sl_down_dif
							new_sys_st_num_up=sys_st_num_up-st_up_dif
							new_sys_st_num_down=sys_st_num_down-st_down_dif

							new_env_bl_num_up=env_bl_num_up-el_up_dif
							new_env_bl_num_down=env_bl_num_down-el_down_dif
							new_env_st_num_up=env_st_num_up+et_up_dif
							new_env_st_num_down=env_st_num_down+et_down_dif

							new_sys_num_up=sys_num_up+sl_up_dif-st_up_dif
							new_sys_num_down=sys_num_down+sl_down_dif-st_down_dif
							new_env_num_up=env_num_up-el_up_dif+et_up_dif
							new_env_num_down=env_num_down-el_down_dif+et_down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
							end do
							101 continue

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==new_sys_st_num_up) then
								if(sys_st%sub(k)%num_down==new_sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 102
								endif
								endif
							end do
							102 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
								if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 103
								endif
								endif
							end do
							103 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==env_st_num_up) then
								if(env_st%sub(k)%num_down==env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 104
								endif
								endif
							end do
							104 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 105
										endif
										endif
									end do
								endif
								endif
							end do
							105 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 106
										endif
										endif
									end do
								endif
								endif
							end do
							106 continue

							if(sys_bl_flag.and.sys_st_flag.and.env_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=1.0d0

										!<b>: For sys_site operator
										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
											sys_bl_num=sys_bl_num_up+sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<c>: For env_site operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
											sys_num=sys_num_up+sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif

										!<d>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
											sys_num=sys_num_up+sys_num_down
											env_st_num=env_st_num_up+env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif

										Signs=SignSL*SignST*SignET*SignEL
						coef_tmp=Signs* coef*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
									allocate(mat(new_sys_dim,old_env_dim))
								call ZGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												&,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												&,czero,mat,new_sys_dim)
										call ZGEMM('N','N',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												&,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,cone&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 107
									endif
									endif
								end do
							endif
							107 continue

							
							!<2>: For C_4^+(env_bl).C_3(env_st).C_2^+(sys_st).C_1(sys_bl) case
							new_sys_bl_num_up=sys_bl_num_up-sl_up_dif
							new_sys_bl_num_down=sys_bl_num_down-sl_down_dif
							new_sys_st_num_up=sys_st_num_up+st_up_dif
							new_sys_st_num_down=sys_st_num_down+st_down_dif

							new_env_bl_num_up=env_bl_num_up+el_up_dif
							new_env_bl_num_down=env_bl_num_down+el_down_dif
							new_env_st_num_up=env_st_num_up-et_up_dif
							new_env_st_num_down=env_st_num_down-et_down_dif

							new_sys_num_up=sys_num_up-sl_up_dif+st_up_dif
							new_sys_num_down=sys_num_down-sl_down_dif+st_down_dif
							new_env_num_up=env_num_up+el_up_dif-et_up_dif
							new_env_num_down=env_num_down+el_down_dif-et_down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 108
								endif
								endif
							end do
							108 continue

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==sys_st_num_up) then
								if(sys_st%sub(k)%num_down==sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 109
								endif
								endif
							end do
							109 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 110
								endif
								endif
							end do
							110 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==new_env_st_num_up) then
								if(env_st%sub(k)%num_down==new_env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 111
								endif
								endif
							end do
							111 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 112
										endif
										endif
									end do
								endif
								endif
							end do
							112 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 113
										endif
										endif
									end do
								endif
								endif
							end do
							113 continue


							if(sys_bl_flag.and.sys_st_flag.and.env_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=1.0d0

										!<b>: For sys_site operator
										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
											sys_bl_num=new_sys_bl_num_up+new_sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<c>: For env_site operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
											sys_num=new_sys_num_up+new_sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif

										!<d>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
											sys_num=new_sys_num_up+new_sys_num_down
											env_st_num=new_env_st_num_up+new_env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif

										Signs=SignSL*SignST*SignET*SignEL

                if(realcode)then
		coef_tmp=Signs* dconjg(coef)*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
                        else
		coef_tmp=Signs* dconjg(coef*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1))
                        endif
									allocate(mat(new_sys_dim,old_env_dim))
		call ZGEMM('C','N',new_sys_dim,old_env_dim,old_sys_dim,cone&
												 &,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												 &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												 &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												 &,czero,mat,new_sys_dim)
										call ZGEMM('N','C',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												 &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,cone&
												 &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 114
									endif
									endif
								end do
							endif
							114 continue
						end do
						end do
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine block_site_vec_block_site_ndia
