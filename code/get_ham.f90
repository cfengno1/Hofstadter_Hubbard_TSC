!==========================================================
!Get system hamiltonian from block and one site
!==========================================================
subroutine Get_Hamiltonian_Sys(pre,next,basis, trun)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: basis
	type(Total_Model),intent(inout) :: pre
	type(Total_Model),intent(inout) :: next
	type(Total_Block),intent(in) :: trun

	integer :: i,x,y,old_xidx,old_yidx,new_xidx,new_yidx
        integer :: i1, j, j1,ij, k, nleg1,ii

        integer k11, k12, k13
	double complex :: coefx,coefy,coefsd,coefsz,coefsn
	type(Total_Block) :: tmp_oper, tmp_oper1
        type(Block_set) tmp_sys !!!!_up,tmp_down,tmp_sd, tmp_sz,tmp_sn


	!<1>: Update the Hamiltonian
	next%len=pre%len+1
	old_xidx=Lattice(1,pre%len)
	old_yidx=Lattice(2,pre%len)
	new_xidx=Lattice(1,next%len)
	new_yidx=Lattice(2,next%len)

	call update_block_dia(pre%ham,basis,next%ham)


	!<2>: Update all the operators

        kc=0

        if(tjmodel==1)k11=5
        if(tjmodel==1)k12=2
        if(hubmodel==1)k11=6
        if(hubmodel==1)k12=5


        do ij=1, 6  !!!k11, k12
        if(ij==2)go to 22
        if(ij==4)go to 22
        if(tjmodel==1)then
        if(ij==6)go to 22
        endif
        if(hubmodel==1)then
        if(ij==3)go to 22
        if(ij==5.and.v123.eq.0)go to 22
                !! ij=1 and 6 left for Hubbard model

        endif

                j=pre%len+1
                        do k=1, neibt
                I=neib(J,k)
          if(I.le.pre%len.and.I.gt.0)then!!! coupling of this site with system
                        i1=inb1(i, pre%len) !!index for the oper
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, j)
         coefx = dconjg(jt(k,j)*pint(j,i))  !!!*pind(j,i)
	coefsd= jd(k, j)!*pind(j,i)
	coefsz= jz(k, j)!*pinz(j,i)
	coefsn=  jn(k, j)
                if(ij==1)then
                if(coefx.ne.0.0)then
                call block_transfer_trans(st_elec_down, tmp_oper)
	call Get_hamiltonian_ndia(pre%sub(i1)%elec_up,tmp_oper, basis,coefx,next%ham,'F','S')
                        call deallocate_block(tmp_oper)
                endif
                endif

                if(ij==2)then
                if(coefx.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%elec_down,st_elec_down,basis,coefx,next%ham,'F','S')
                endif
                endif
                if(ij==3)then
                if(coefsd.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%spin_sd,st_sd, basis,coefsd,next%ham, 'B','S')
                endif
                endif

	endif
	endif
        enddo  !! enddo k=neibt

	do y=1,pre%len   !!!!Ny
		!Update center-site operators
                I1=INB1(y, pre%len)
                if(I1.le.nleg11(pre%len))then !! useful operator related to y site

                if(ij==1)then
    call update_block_ndia_sys(pre%sub(i1)%elec_up,basis,next%sub(i1)%elec_up,'F')
    call update_block_ndia_sys(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')

                       ! call block_to_disk1(next%sub(i1)%elec_up, 301)
                       ! call block_to_disk1(next%sub(i1)%elec_down, 302)
                        call deallocate_block(pre%sub(i1)%elec_up)
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==2)then
   call update_block_ndia_sys(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
           call deallocate_block(pre%sub(i1)%elec_down)
                    endif
                 if(ij==3)then
       call update_block_ndia_sys(pre%sub(i1)%spin_sd,basis,next%sub(i1)%spin_sd,'B')
                  call deallocate_block(pre%sub(i1)%spin_sd)
                endif
               if(ij==4)then
             call update_block_dia(pre%sub(i1)%spin_sz,basis,next%sub(i1)%spin_sz)
               call deallocate_block(pre%sub(i1)%spin_sz)
                  endif
                  if(ij==5)then
             call update_block_dia(pre%sub(i1)%num_sn,basis,next%sub(i1)%num_sn)
              call deallocate_block(pre%sub(i1)%num_sn)
                        endif
           !       if(ij==6)then
            ! call update_block_dia(pre%sub(i1)%double,basis,next%sub(i1)%double)
             ! call deallocate_block(pre%sub(i1)%double)
              !          endif

		        endif
                        enddo

        !!if(ij==3.and.next%len.eq.4)stop

!!! adding the diagonal to the ham
        go to 11
                if(esite(pre%len+1).ne.0.0)then
        call block_transfer(next%ham, tmp_oper)
call block_add_block(tmp_oper,0.d0,tmp_sys%num_sn,esite(pre%len+1),next%ham,1.d0)
                call deallocate_block(tmp_oper)
                endif
11      continue


	!Get Hamiltonian from all the bond interaction terms
	!<3>: Get Hamiltonian: x-bond interaction
		!For hopping interaction


                if(ij==1)then
          call update_site_ndia_sys(st_elec_up,basis,tmp_oper,'F')
          call update_site_ndia_sys(st_elec_down,basis,tmp_oper1,'F')
                endif
                if(ij==2)then
          call update_site_ndia_sys(st_elec_down,basis,tmp_oper,'F')
                endif
        if(ij==3)then
          call update_site_ndia_sys(st_sd,basis,tmp_oper,'B')
                       ! call block_to_disk1(tmp_oper, 105)
                endif
                if(ij==4)then
          call update_site_dia(st_sz,basis,tmp_oper)
                endif
                if(ij==5)then
          call update_site_dia(num_elec,basis,tmp_oper)
                endif
                if(ij==6)then
          call update_site_dia(st_double,basis,tmp_oper)
                endif


        if(ij==3)then
   if(lring.eq.1)then
                ms=0
        call block_transfer(tmp_oper, tmp_sys%spin_sd)
        call get_operator(pre,next,tmp_sys,basis,trun)
                endif
                endif




                j=pre%len+1
                        do k=1, neibt
                I=neib(J,k)
          if(I.le.pre%len.and.I.gt.0)then!!! coupling of this site with system
                        i1=inb1(i, pre%len) !!index for the oper
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, j)
        
	coefsd= jd(k, j)!*pind(j,i)
	coefsz= jz(k, j)!*pinz(j,i)
	coefsn=  jn(k, j)
                if(ij==4)then
		call Get_hamiltonian_dia(next%sub(i1)%spin_sz,'N', tmp_oper,'N',coefsz,next%ham,.false.)
                endif
                if(ij==5)then
                if(coefsn.ne.0.0)then
		call Get_hamiltonian_dia(next%sub(i1)%num_sn,'N', tmp_oper,'N',coefsn,next%ham,.false.)
                endif
                endif

	endif
	endif
        enddo  !! do k


                if(ij==6)then
                if(hubbard(pre%len+1).ne.0.0)then
                coefx=hubbard(pre%len+1)
call block_add_block_two(tmp_oper,coefx,next%ham)
                endif
                endif

!!! reorder the operators

                  nleg1=nleg11(pre%len+1)
         do I=1,pre%len   !!! this the middle block has more than
                I1=INB1(I,pre%len) !! this is the I1 before merge
                J = INB1(I,pre%len+1)
                              
                        if(J.le.nleg1.and.j.ne.i1)then

                if(ij==1)then
        call block_transfer(next%sub(i1)%elec_up,next%sub(j)%elec_up)
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==3)then
        call block_transfer(next%sub(i1)%spin_sd,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(next%sub(i1)%spin_sz,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(next%sub(i1)%num_sn,next%sub(j)%num_sn)
                endif
       !         if(ij==6)then
       ! call block_transfer(next%sub(i1)%double,next%sub(j)%double)
        !        endif
                                endif
                                enddo

                    j=inb1(pre%len+1,pre%len+1)
                if(ij==1)then
        call block_transfer(tmp_oper,next%sub(j)%elec_up)
        call block_transfer(tmp_oper1,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(tmp_oper,next%sub(j)%elec_down)
                        endif
                if(ij==3)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(tmp_oper,next%sub(j)%num_sn)
                endif
        !        if(ij==6)then
        !call block_transfer(tmp_oper,next%sub(j)%double)
        !        endif

                        if(nleg1.lt.nleg11(pre%len))then
                        do j=nleg1+1,nleg11(pre%len)
                        if(ij==1)then
                call deallocate_block(next%sub(j)%elec_up)
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==2)then
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==3)then
                call deallocate_block(next%sub(j)%spin_sd)
                        endif
                        if(ij==4)then
                call deallocate_block(next%sub(j)%spin_sz)
                        endif
                        if(ij==5)then
                call deallocate_block(next%sub(j)%num_sn)
                        endif
        !                if(ij==6)then
         !       call deallocate_block(next%sub(j)%double)
          !              endif
                                enddo
                                endif
                
		call deallocate_block(tmp_oper)
	!For Hamiltonian


	!For operators in the middle of system
	do ii=1, nleg11(next%len)

                        if(ij==1)then
		call update_trun_ndia(next%sub(ii)%elec_up,pre%sub(ii)%elec_up,trun) 
		call update_trun_ndia(next%sub(ii)%elec_down,pre%sub(ii)%elec_down,trun) 
                        call deallocate_block(next%sub(ii)%elec_up)
                        call deallocate_block(next%sub(ii)%elec_down)
                        endif

                        if(ij==2)then
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                                endif

                        if(ij==3)then
		call update_trun_ndia(next%sub(ii)%spin_sd,pre%sub(ii)%spin_sd,trun)
                        call deallocate_block(next%sub(ii)%spin_sd)
                        endif

                                if(ij==4)then
		call update_trun_dia(next%sub(ii)%spin_sz,pre%sub(ii)%spin_sz,trun)
                        call deallocate_block(next%sub(ii)%spin_sz)
                                endif

                        if(ij==5)then
		call update_trun_dia(next%sub(ii)%num_sn,pre%sub(ii)%num_sn,trun)
                        call deallocate_block(next%sub(ii)%num_sn)
                                endif
          !              if(ij==6)then
	!	call update_trun_dia(next%sub(ii)%double,pre%sub(ii)%double,trun)
         !               call deallocate_block(next%sub(ii)%double)
	end do

22     continue
                        enddo
	call update_trun_dia(next%ham,tmp_oper,trun)
                        call block_transfer(tmp_oper, pre%ham)
                        call deallocate_block(tmp_oper)
                        call deallocate_block(next%ham)
                pre%len=next%len

end subroutine Get_Hamiltonian_Sys


!==========================================================
!Get system hamiltonian from block and one site
!==========================================================
subroutine Get_Hamiltonian_Env(pre,next,basis, trun)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: basis
	type(Total_Model),intent(inout) :: pre
	type(Total_Model),intent(inout) :: next
	type(Total_Block),intent(in) :: trun

	integer :: i,x,y,old_xidx,old_yidx,new_xidx,new_yidx
        integer k11, k12, k13
        integer :: i1, j, j1,ij, k, nleg1,ii
	double complex :: coefx,coefy,coefsd,coefsz,coefsn
	type(Total_Block) :: tmp_oper, tmp_oper1
        type(Block_set) tmp_sys !!!!_up,tmp_down,tmp_sd, tmp_sz,tmp_sn


	!<1>: Update the Hamiltonian
	next%len=pre%len+1
        old_xidx=Lattice(1,nxc*nleg-pre%len+1)
        old_yidx=Lattice(2,nxc*nleg-pre%len+1)
        new_xidx=Lattice(1,nxc*nleg-pre%len+2)
        new_yidx=Lattice(2,nxc*nleg-pre%len+2)

	call update_block_dia(pre%ham,basis,next%ham)

	!<2>: Update all the operators
        if(tjmodel==1)k11=5
        if(tjmodel==1)k12=2
        if(hubmodel==1)k11=6
        if(hubmodel==1)k12=5
        do ij=1, 6
        if(ij==2)go to 22
        if(ij==4)go to 22
        if(tjmodel==1)then
        if(ij==6)go to 22
        endif
        if(hubmodel==1)then
        if(ij==3)go to 22
        if(ij==5.and.v123.eq.0)go to 22
        endif
                        do k=1, neibt
                j=pre%len+1+mshift
                I=neib(nposi_lat(J),k)
         II=nposi_dm(I)-mshift  
          if(II.le.pre%len.and.II.gt.0)then!!! coupling of this site with system
                        i1=inb1(ii, pre%len) !!index for the oper
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, nposi_lat(j))
        coefx = jt(k, nposi_lat(j))*pint(nposi_lat(j),i)
	coefsd= jd(k,nposi_lat(j))!*pind(j,i)
	coefsz= jz(k,nposi_lat(j))!*pinz(j,i)
	coefsn= jn(k,nposi_lat(j))
                if(ij==1)then
                if(coefx.ne.0.0)then
                call block_transfer_trans(st_elec_up, tmp_oper)
        coefx=dconjg(coefx)
	call Get_hamiltonian_ndia(pre%sub(i1)%elec_down,tmp_oper,basis,coefx,next%ham,'F','E')
!! makechange
                        call deallocate_block(tmp_oper)
                endif
                endif
                if(ij==2)then
                if(coefx.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%elec_down,st_elec_down,basis, coefx,next%ham,'F','E')
                endif
                endif
                if(ij==3)then
                if(coefsd.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%spin_sd,st_sd, basis,coefsd,next%ham,'B','E')
                endif
                endif
	endif
	endif
        enddo

!!! reorder the operators
	do y=1,pre%len   !!!!Ny
		!Update center-site operators
                I1=INB1(y, pre%len)
                if(I1.le.nleg11(pre%len))then !! useful operator related to y site

                if(ij==1)then
                call update_block_ndia_env(pre%sub(i1)%elec_up,basis,next%sub(i1)%elec_up,'F')
                call update_block_ndia_env(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
                        call deallocate_block(pre%sub(i1)%elec_up)
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==2)then
                        call update_block_ndia_env(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==3)then
                        call update_block_ndia_env(pre%sub(i1)%spin_sd,basis,next%sub(i1)%spin_sd,'B')
                        !call block_to_disk1(next%sub(i1)%spin_sd, 104)
                        call deallocate_block(pre%sub(i1)%spin_sd)
                        endif
                        if(ij==4)then
                        call update_block_dia(pre%sub(i1)%spin_sz,basis,next%sub(i1)%spin_sz)
                        call deallocate_block(pre%sub(i1)%spin_sz)
                        endif
                        if(ij==5)then
                        call update_block_dia(pre%sub(i1)%num_sn,basis,next%sub(i1)%num_sn)
                        call deallocate_block(pre%sub(i1)%num_sn)
                        endif
                     !   if(ij==6)then
                      !  call update_block_dia(pre%sub(i1)%double,basis,next%sub(i1)%double)
                       ! call deallocate_block(pre%sub(i1)%double)
                       ! endif

		        endif
                        enddo

!!! adding the diagonal to the ham

        go to 11
                if(esite(nposi_lat(j)).ne.0.0)then
        call block_transfer(next%ham, tmp_oper)
call block_add_block(tmp_oper,0.d0,tmp_sys%num_sn,esite(nposi_lat(j)),next%ham,1.d0)
                call deallocate_block(tmp_oper)
                endif
11      continue


	!Get Hamiltonian from all the bond interaction terms
	!<3>: Get Hamiltonian: x-bond interaction
		!For hopping interaction


                if(ij==1)then
          call update_site_ndia_env(st_elec_up,basis,tmp_oper,'F')
          call update_site_ndia_env(st_elec_down,basis,tmp_oper1,'F')
                endif
                if(ij==2)then
          call update_site_ndia_env(st_elec_down,basis,tmp_oper,'F')
                endif
        if(ij==3)then
          call update_site_ndia_env(st_sd,basis,tmp_oper,'B')
                endif
                if(ij==4)then
          call update_site_dia(st_sz,basis,tmp_oper)
                endif
                if(ij==5)then
          call update_site_dia(num_elec,basis,tmp_oper)
                endif
                if(ij==6)then
          call update_site_dia(st_double,basis,tmp_oper)
                endif


        if(ij==3)then
   if(lring.eq.1)then
                ms=1
                call block_transfer(tmp_oper, tmp_sys%spin_sd)
        call get_operator(pre,next,tmp_sys,basis,trun)
                endif
                endif

                        do k=1, neibt
                I=neib(nposi_lat(J),k)
         II=nposi_dm(I)-mshift  
          if(II.le.pre%len.and.II.gt.0)then!!! coupling of this site with system
                        i1=inb1(ii, pre%len) !!index for the oper
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, nposi_lat(j))
        coefx = jt(k, nposi_lat(j))*pint(nposi_lat(j),i)
	coefsd= jd(k,nposi_lat(j))!*pind(j,i)
	coefsz= jz(k,nposi_lat(j))!*pinz(j,i)
	coefsn= jn(k,nposi_lat(j))
                if(ij==4)then
		call Get_hamiltonian_dia(next%sub(i1)%spin_sz,'N', tmp_oper,'N',coefsz,next%ham,.false.)
                endif

                if(ij==5)then
                if(coefsn.ne.0.0)then
		call Get_hamiltonian_dia(next%sub(i1)%num_sn,'N', tmp_oper,'N',coefsn,next%ham,.false.)
                endif
                endif

                endif











	endif
        enddo

                if(ij==6)then
                if(hubbard(nposi_lat(j)).ne.0.0)then
                coefx=hubbard(nposi_lat(j))
call block_add_block_two(tmp_oper,coefx,next%ham)
                endif
                endif
!!! reorder the operators

                  nleg1=nleg11(pre%len+1)
         do I=1,pre%len   !!! this the middle block has more than
                I1=INB1(I,pre%len) !! this is the I1 before merge
                J = INB1(I,pre%len+1)
                        if(J.le.nleg1.and.j.ne.i1)then

                if(ij==1)then
        call block_transfer(next%sub(i1)%elec_up,next%sub(j)%elec_up)
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==3)then
        call block_transfer(next%sub(i1)%spin_sd,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(next%sub(i1)%spin_sz,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(next%sub(i1)%num_sn,next%sub(j)%num_sn)
                endif
      !          if(ij==6)then
       ! call block_transfer(next%sub(i1)%double,next%sub(j)%double)
        !        endif
                                endif
                                enddo

                    j=inb1(pre%len+1,pre%len+1)
                if(ij==1)then
        call block_transfer(tmp_oper,next%sub(j)%elec_up)
        call block_transfer(tmp_oper1,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(tmp_oper,next%sub(j)%elec_down)
                        endif
                if(ij==3)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(tmp_oper,next%sub(j)%num_sn)
                endif
       !         if(ij==6)then
       ! call block_transfer(tmp_oper,next%sub(j)%double)
        !        endif

                        if(nleg1.lt.nleg11(pre%len))then
                        do j=nleg1+1,nleg11(pre%len)
                        if(ij==1)then
                call deallocate_block(next%sub(j)%elec_up)
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==2)then
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==3)then
                call deallocate_block(next%sub(j)%spin_sd)
                        endif
                        if(ij==4)then
                call deallocate_block(next%sub(j)%spin_sz)
                        endif
                        if(ij==5)then
                call deallocate_block(next%sub(j)%num_sn)
                        endif
        !                if(ij==6)then
         !       call deallocate_block(next%sub(j)%double)
          !              endif
                                enddo
                                endif
                
		call deallocate_block(tmp_oper)
	!For Hamiltonian


	!For operators in the middle of system
	do ii=1, nleg11(next%len)

                        if(ij==1)then
		call update_trun_ndia(next%sub(ii)%elec_up,pre%sub(ii)%elec_up,trun) 
                        call deallocate_block(next%sub(ii)%elec_up)
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                        endif

                        if(ij==2)then
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                                endif

                        if(ij==3)then
		call update_trun_ndia(next%sub(ii)%spin_sd,pre%sub(ii)%spin_sd,trun)
                        call deallocate_block(next%sub(ii)%spin_sd)
                        endif

                                if(ij==4)then
		call update_trun_dia(next%sub(ii)%spin_sz,pre%sub(ii)%spin_sz,trun)
                        call deallocate_block(next%sub(ii)%spin_sz)
                                endif

                        if(ij==5)then
		call update_trun_dia(next%sub(ii)%num_sn,pre%sub(ii)%num_sn,trun)
                        call deallocate_block(next%sub(ii)%num_sn)
                                endif
         !               if(ij==6)then
	!	call update_trun_dia(next%sub(ii)%double,pre%sub(ii)%double,trun)
         !               call deallocate_block(next%sub(ii)%double)
          !                      endif
	end do

22     continue
                        enddo
	call update_trun_dia(next%ham,tmp_oper,trun)
                        call block_transfer(tmp_oper, pre%ham)
                        call deallocate_block(tmp_oper)
                        call deallocate_block(next%ham)
                pre%len=next%len

end subroutine Get_Hamiltonian_Env




subroutine get_operator(pre,next,sys_set, basis,trun)
	use pubdata
	implicit integer (i-n)
	type(Total_Basis),intent(in) :: basis
	type(Total_Model),intent(inout) :: pre
	type(Total_Model),intent(inout) :: next
	type(Total_Block) :: tmp_oper
	type(Total_Block) :: sys_sn,sys_sd,tmp_bl, tmp_bl1, tmp_bl2
	type(Block_set) :: sys_set
	type(Total_Block),intent(in) :: trun

	double complex :: alpha, beta
        double complex ::  cc1, cc2, ph1(3) 
        integer  qn1, qn_dif

	num20=nleg2(pre%len, ms+1) !! previous number of operators

		if(num20.gt.0)then
         do I=1, num20   !!! updating them existing spin coupling terms
        
        if(allocated(tmp_bl%sub))call deallocate_block(tmp_bl)
        if(ms.eq.0)then
  call update_block_ndia_sys(pre%spin1(i)%spin_sd,basis,tmp_bl,'B')
                else
  call update_block_ndia_env(pre%spin1(i)%spin_sd,basis,tmp_bl,'B')
                endif
        	call deallocate_block(pre%spin1(i)%spin_sd)
         call block_transfer(tmp_bl, pre%spin1(i)%spin_sd)
        call deallocate_block(tmp_bl)
			enddo  !!! end of do i
			endif !!! num20.ne.0

        if(pre%len.ge.2)then
                        do j21=1, num20
                        j=tri3b(j21, pre%len, ms+1)  
                if(j.eq.pre%len+1)then
                        alpha=1.0d0
                          alpha=jring(1)*tri3s(j21, pre%len, ms+1)*ci


  if(cdabs(alpha).gt.0.0)then
                        if(ms.eq.0)then
                tmp_bl%down_dif=0
        call block_mul_block_ndia(pre%spin1(j21)%spin_sd,'N',sys_set%spin_sd,'N',cone,tmp_bl,basis)
         call block_add_block_two(tmp_bl,alpha,next%ham)
                        else
                tmp_bl%down_dif=0
        call block_mul_block_ndia(pre%spin1(j21)%spin_sd,'N',sys_set%spin_sd,'N',cone,tmp_bl,basis)
         call block_add_block_two(tmp_bl,alpha,next%ham)
                        endif
                 call deallocate_block(tmp_bl)
                        endif


                call deallocate_block(pre%spin1(j21)%spin_sd)
                        endif
                        enddo
                        endif


                num21=num20  !!! making new spin2 operators from here
!!!! new spin2 called spin1%spin0(i, j)
                !I=neib(nposi_lat(J),k)             !! we use fixed coordinate
                !II=nposi_dm(I)-mshift  !!! dmrg position - system-length

                do ij=1, kk1 !! kk1 number of different rings, 
                 j = pre%len+1 

                        do k=1,neibt1  !!!! NN 
                    j1=neib(nposi_lat(J+mshift),k)
                I=nposi_dm(j1)-mshift  !!! dmrg position - system-length

                       if(I.le.pre%len.and.I.gt.0)then 

                                k1=0
                        if(ms.eq.0)then !! from system
                        j2=tri3(j,j1,ij)
                       if(j2.gt.pre%len+1.and.j2.le.num_site)k1=1
                                else
                        j2=tri3(nposi_lat(j+mshift),j1,ij)
                       if(j2.lt.nposi_lat(j+mshift).and.j2.gt.0)k1=1 !! not in env
                                endif

                        if(k1.ne.0)then !! adding a term
                        i1=inb1(i, pre%len) !! j1
                        ji=num_spin2(j, i, ij)  !! the same for sys end env
                num21=num21+1
        ind2e(ji,pre%len, ms+1)=num21!!!! temperarily storing spin2 operators

                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
   !  call block_pass_info_trans(next%sub(i1)%spin_sd, pre%spin1(num21)%spin_sd)
        k11=0
                if(allocated(tmp_bl%sub))call deallocate_block(tmp_bl)
                call block_transfer(next%sub(i1)%spin_sd, tmp_bl) 

                if(allocated(tmp_bl1%sub))call deallocate_block(tmp_bl1)
                call block_transfer(sys_set%spin_sd, tmp_bl1) 

                beta=1.0d0
                if(allocated(tmp_bl2%sub))call deallocate_block(tmp_bl2)
                tmp_bl2%down_dif=su
      call block_mul_block_ndia(tmp_bl,'N', tmp_bl1,'N',cone, tmp_bl2, basis)
             call block_transfer(tmp_bl2,pre%spin1(num21)%spin_sd)
         !!call block_add_block_two(tmp_bl2,beta,pre%spin1(num21)%spin_sd)

                call deallocate_block(tmp_bl2)
                call deallocate_block(tmp_bl)
                call deallocate_block(tmp_bl1)
                endif  !if(i1
                endif  !if(k1
                endif  !if(I
                enddo
                enddo


!! now we transfer all operators to new model

        do ij=1,kk1
   do j=1,pre%len
                do j1=j+1, pre%len+1
                i=num_spin2(j,j1, ij)
                ki=ind2e(i,pre%len, ms+1) !!added new f pre%len

                if(ki.gt.0)then
                kk2=ind2e(i,pre%len+1, ms+1)!!!!=ddkk2
          if(kk2.le.nleg2(pre%len+1, ms+1).and.kk2.gt.0)then!new ones needed
        call block_transfer(pre%spin1(ki)%spin_sd,next%spin1(kk2)%spin_sd)
        if(allocated(pre%spin1(ki)%spin_sd%sub))then
        call deallocate_block(pre%spin1(ki)%spin_sd)
                endif
   if(ki.gt.nleg2(pre%len, ms+1))ind2e(i,pre%len, ms+1)=0
                endif
                endif
                enddo
                enddo
                enddo


        do j=1,nleg2(next%len, ms+1)               !!!!(pre%tlen, ms+1)
                if(allocated(next%spin1(j)%spin_sd%sub))then
                call update_trun_ndia(next%spin1(j)%spin_sd,pre%spin1(j)%spin_sd,trun)
                call deallocate_block(next%spin1(j)%spin_sd)
                        endif
        end do

end subroutine get_operator



!Note: trun=get_trun_oper(den,trun,keep,pm)
subroutine Get_Hamiltonian_Trun(pre,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: trun
	type(Total_Model),intent(in) :: pre
	type(Total_Model),intent(inout) :: eff

	integer :: ii

	!For Hamiltonian
	call update_trun_dia(pre%ham,eff%ham,trun)
	eff%len=pre%len

	!For operators in the middle of system
	do ii=1, nleg11(pre%len)

		call update_trun_ndia(pre%sub(ii)%elec_up,eff%sub(ii)%elec_up,trun)
                        call deallocate_block(pre%sub(ii)%elec_up)

                                if(tjmodel==1)then
		call update_trun_ndia(pre%sub(ii)%spin_sd,eff%sub(ii)%spin_sd,trun)
                        call deallocate_block(pre%sub(ii)%spin_sd)
                                endif
                if(v123==1)then
		call update_trun_dia(pre%sub(ii)%num_sn,eff%sub(ii)%num_sn,trun)
                        call deallocate_block(pre%sub(ii)%num_sn)
                                endif
	end do

	!For operators in the boundary of system

end subroutine Get_Hamiltonian_Trun


!=========================================================================
!Get hamiltonian from diagonal interaction term n_i.n_i+1
!=========================================================================
subroutine Get_hamiltonian_dia(Oper1,Flag1,Oper2,Flag2,Coef,Ham,FlagConjg)
	use pubdata
	implicit none

	double complex,intent(in) :: Coef
	logical,intent(in) :: FlagConjg
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: Oper1,Oper2
	type(Total_Block),intent(inout) :: Ham

	integer :: i,x,y
	logical :: oper1_flag,oper2_flag
	integer :: num_up,num_down,row_dim,sdim,oper1_id,oper2_id
	double complex,allocatable :: mat(:,:)
        double complex coef1

	do i=1,ham%num
		num_up=ham%sub(i)%num_up
		num_down=ham%sub(i)%num_down
		row_dim=ham%sub(i)%row_dim
		sdim=ham%sub(i)%sdim

		oper1_flag=.false.
		do x=1,oper1%num
			if(oper1%sub(x)%num_up==num_up) then
			if(oper1%sub(x)%num_down==num_down) then
				oper1_id=x
				oper1_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		oper2_flag=.false.
		do x=1,oper2%num
			if(oper2%sub(x)%num_up==num_up) then
			if(oper2%sub(x)%num_down==num_down) then
				oper2_id=x
				oper2_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(oper1_flag.and.oper2_flag) then
			allocate(mat(sdim,sdim))

        coef1=coef/dsqrt(1.0d0+num_down)
	!		call DGEMM(Flag1,Flag2,sdim,sdim,sdim,coef1,oper1%sub(oper1_id)%mat&
	!				 &,sdim,oper2%sub(oper2_id)%mat,sdim,0.0d0,mat,sdim)

			!For the case with conjugate part
	!		if(FlagConjg) then
	!			ham%sub(i)%mat=ham%sub(i)%mat+(mat+transpose(mat))
	!		else
	!			ham%sub(i)%mat=ham%sub(i)%mat+mat
	!		endif


     if(realcode)then
                        call DGEMM(Flag1,Flag2,sdim,sdim,sdim,coef1,oper1%sub(oper1_id)%mat&
                                         &,sdim,oper2%sub(oper2_id)%mat,sdim,0.0d0,mat,sdim)
                        else
     call ZGEMM(Flag1,Flag2,sdim,sdim,sdim,coef1,oper1%sub(oper1_id)%mat&
                                         &,sdim,oper2%sub(oper2_id)%mat,sdim,czero,mat,sdim)
                        endif

                        !For the case with conjugate part
                        if(FlagConjg) then
         if(realcode)then
                                ham%sub(i)%mat=ham%sub(i)%mat+(mat+transpose(mat))
                                        else
                                ham%sub(i)%mat=ham%sub(i)%mat+(mat+dconjg(transpose(mat)))
                                        endif

                        else
                                ham%sub(i)%mat=ham%sub(i)%mat+mat
                        endif





			deallocate(mat)
		endif
	end do

end subroutine Get_hamiltonian_dia


!================================================================
!Get hamiltonian from non-diagonal interaction term (C^+_i.C_i+1)
!new_oper = (C^+_spin. C_spin) with spin=(up, down)
!(1) if spin=(up spin), then Up_dif=Up_bias,Down_dif=0
!(2) if spin=(down spin), then Up_dif=0, Down_dif=Down_bias
!================================================================
!Get Ham= Coef* [Oper1*(Oper2^T)+ h.c.) (i.e. C^+_1.C_2 + h.c.)
!================================================================
!subroutine Get_hamiltonian_ndia1(Oper1,Oper2,Coef,Ham)
!	call Get_hamiltonian_ndia(pre%sub(i1)%elec_up,basis,coefx,next%ham,'F')
subroutine Get_hamiltonian_ndia(block,site,basis,coef,Ham, FlagSign, Flags)
	use pubdata
        use fact1
	implicit none

	character(len=1),intent(in) :: FlagSign, Flags
	type(Total_Block),intent(in) :: block, site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: Ham
        double complex, intent (in) :: coef

	integer :: i,x,y,new_id,block_id
        integer j1,j2,j3,j4,j5,j6
	integer :: lhs_id,lhs_num_up,lhs_num_down0, lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down, st_id

	integer :: lhs_dim,rhs_dim,up_dif,down_dif,down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag, st_flag
        integer st_num_down, st_num_up, new_st_num_up, new_st_num_down, st_down_dif
        real*8  coef11, Sign
        double complex coef1
	double complex,allocatable :: mat(:,:)
  real(8),external :: w6js

	up_dif=block%up_dif
	down_dif1=block%down_dif

	!<3>: Get new_block%sub(:)%mat
	do i=1,ham%num
		rhs_num_up=ham%sub(i)%num_up
		rhs_num_down=ham%sub(i)%num_down
		lhs_dim=ham%sub(i)%row_dim
		rhs_dim=ham%sub(i)%sdim
		
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

		down_dif=ham%sub(i)%down_dif
		lhs_num_up=rhs_num_up
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
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
                                st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
                                st_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su
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

                        if(block_flag)then

        do new_st_num_down=abs(st_num_down-site%down_dif), st_num_down+site%down_dif,su

                  st_flag=.false.
                        st_down_dif=st_num_down-new_st_num_down
                        new_st_num_up=st_num_up-site%up_dif
                                             do y=1,site%num
                                       if(site%sub(y)%num_up==new_st_num_up) then
                                              if(site%sub(y)%num_down==new_st_num_down) then
                                                        if(site%sub(y)%down_dif==st_down_dif) then
                                                                st_id=y
                                                                st_flag=.true.
                                                                goto 1031
                                                        endif
                                                        endif
                                                        endif
                                                end do
                                                1031 continue


				lhs_sub_flag=.false.

				do y=1,basis%sub(lhs_id)%num
				if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
				if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
				if(basis%sub(lhs_id)%sub(y)%st_num_down==new_st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
                                        endif
                                        endif
				end do
				105 continue

				if(lhs_sub_flag.and.block_flag.and.st_flag) then
        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        ! j1=left tot J', j2=tensor rank, J3=righ tot
        ! j4=right block, j5=left site, j6=left block


        j1=lhs_sub_num_down !(j')
        j2=rhs_sub_num_down  ! j
        j3=down_dif1     ! tensor rank 
        j4=st_num_down  !J
        j5=new_st_num_down  ! j
        j6=rhs_num_down     !  J

        if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef11=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef11=w6js(j1,j2,j3,j4, j5,j6)  !! eq51
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef11
        endif

        
        coef11=coef11*dsqrt((1.0d0+j6)/(1.d0+j3))*(-1)**((j2+j5+j3+j6)/2)*site%sub(st_id)%mat(1,1)

                                                        !<a>: For sys_block
                                                        !operator
                if(Flags=='S')then !!! sys,  add st sign
                                              Sign=1.0d0
                                                        if(FlagSign=='F') then
                                                         if(mod(rhs_sub_num_up, 2)==1)Sign=-1.0d0
                                                                endif

                        else !! env and st
                        Sign=1.0d0
                                    if(FlagSign=='F') then
                                            if(mod(st_num_up, 2)==0)Sign=-1.0d0
                                                          endif
                                                                endif



        if(coef11.ne.0.0)then
                coef1=coef11*coef*Sign
	ham%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&=ham%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&+coef1*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)


        if(site%up_dif.ne.0)then
                allocate(mat(rhs_sub_dim, lhs_sub_dim))
                mat= transpose(block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim))
                if(.not.realcode)mat=dconjg(mat)
	ham%sub(i)%mat(rhs_pos+1:rhs_pos+rhs_sub_dim,lhs_pos+1:lhs_pos+lhs_sub_dim)&
		&=ham%sub(i)%mat(rhs_pos+1:rhs_pos+rhs_sub_dim,lhs_pos+1:lhs_pos+lhs_sub_dim)&
				&+dconjg(coef1)*mat

			deallocate(mat)
                        endif
                                endif 
				endif
			end do
        endif
                enddo
                enddo
		endif
        enddo

       ! call block_to_disk1(ham,103)
       ! call block_to_disk1(block,103)

end subroutine Get_hamiltonian_ndia


subroutine Get_hamiltonian_ndia1(Oper1,Oper2,Coef,Ham)
	use pubdata
	implicit none

	real(8),intent(in) :: Coef
	type(Total_Block),intent(in) :: Oper1,Oper2
	type(Total_Block),intent(inout) :: Ham

	logical :: oper1_flag,oper2_flag
	integer :: i,x,y,up_dif,down_dif
	integer :: num_up,num_down,mid_num_up,mid_num_down
	integer :: row_dim,col_dim,oper1_id,oper2_id
        integer down_dif1, q1,m1
        real*8 coef1
        real*8, external ::w3js

	!Get next%ham from oper1 and oper2
	up_dif=Oper1%up_dif
	down_dif1=Oper1%down_dif
	if((up_dif/=Oper2%up_dif).or.(down_dif1/=Oper2%down_dif)) then
		write(*,*) "up_dif or down_dif error in Get_hamiltonian_ndia!"
		stop
	endif

	do i=1,ham%num
		num_up=ham%sub(i)%num_up
		num_down=ham%sub(i)%num_down

        do down_dif=-down_dif1, down_dif1, su !! j change
		mid_num_up=num_up-up_dif
		mid_num_down=num_down-down_dif

		oper1_flag=.false.
		do x=1,oper1%num
			if(oper1%sub(x)%num_up==mid_num_up) then
			if(oper1%sub(x)%num_down==mid_num_down) then
			if(oper1%sub(x)%down_dif==down_dif) then
				oper1_id=x
				oper1_flag=.true.
				goto 101
			endif
			endif
			end if
		end do
		101 continue

		oper2_flag=.false.
		do x=1,oper2%num
			if(oper2%sub(x)%num_up==mid_num_up) then
			if(oper2%sub(x)%num_down==mid_num_down) then
			if(oper2%sub(x)%down_dif==down_dif) then
				oper2_id=x
				oper2_flag=.true.
				goto 102
			endif
			endif
			endif
		end do
		102 continue
		
		if(oper1_flag.and.oper2_flag) then

			row_dim=ham%sub(i)%row_dim
			col_dim=oper1%sub(oper1_id)%sdim

                coef1=0.d0
        do q1=-down_dif1, down_dif1, su       !! q1 for tensor operator Jz
         m1=num_down-q1
                coef1=coef1+(-1)**(1-q1/2)*w3js(num_down,down_dif1,mid_num_down,-num_down,q1,num_down-q1)**2
                enddo
                coef1=coef1*coef*dsqrt((1.d0+num_down)/(1.d0+down_dif1))
			!<1>: For the case with (C_1^+.C_2)


        if(coef1.ne.0.0)then
			call DGEMM('N','T',row_dim,row_dim,col_dim,coef1,oper1%sub(oper1_id)%mat&
					&,row_dim,oper2%sub(oper2_id)%mat,row_dim,1.0d0,ham%sub(i)%mat,row_dim)

                endif
		endif
        enddo
	end do

end subroutine Get_hamiltonian_ndia1

!===========================================================================
!Get operator in the new representation
!===========================================================================
!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_sys_oper_dia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,sys_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_dia(st_oper,sys_bsm(idx_one),sys_oper)
			else
				call update_site_dia(st_oper,sys_bsm(idx_one),mid_oper)
				call update_trun_dia(mid_oper,sys_oper,systruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(sys_oper,sys_bsm(idx+1),mid_oper)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(sys_oper,sys_bsm(idx+1),new_oper)
		call deallocate_block(sys_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_dia(st_oper,sys_bsm(idx_two),new_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

end subroutine Get_sys_oper_dia


!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_env_oper_dia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,env_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_dia(st_oper,env_bsm(idx_one),env_oper)
			else
				call update_site_dia(st_oper,env_bsm(idx_one),mid_oper)
				call update_trun_dia(mid_oper,env_oper,envtruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(env_oper,env_bsm(idx+1),mid_oper)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(env_oper,env_bsm(idx+1),new_oper)
		call deallocate_block(env_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_dia(st_oper,env_bsm(idx_two),new_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

end subroutine Get_env_oper_dia


!Change operator from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Change_sys_oper_dia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,sys_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(sys_oper,sys_bsm(idx+1),mid_oper)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(sys_oper,sys_bsm(idx+1),new_oper)
		call deallocate_block(sys_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_sys_oper_dia


!Change operator from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Change_env_oper_dia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,env_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(env_oper,env_bsm(idx+1),mid_oper)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(env_oper,env_bsm(idx+1),new_oper)
		call deallocate_block(env_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif
	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif

end subroutine Change_env_oper_dia


!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_sys_oper_dia2(st_oper1,st_oper2,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper1,st_oper2
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: sys_oper,tmp_oper,mid_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,sys_oper)
		else
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,tmp_oper)
			call update_site_dia(tmp_oper,sys_bsm(idx_one),mid_oper)
			call deallocate_block(tmp_oper)

			if(idx_one<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx_one))
			endif
			call deallocate_block(mid_oper)
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(sys_oper,sys_bsm(idx+1),mid_oper)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(sys_oper,sys_bsm(idx+1),new_oper)
		call deallocate_block(sys_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,new_oper)
		else
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,mid_oper)
			call update_site_dia(mid_oper,sys_bsm(idx_two),new_oper)
			call deallocate_block(mid_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

end subroutine Get_sys_oper_dia2


!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_env_oper_dia2(st_oper1,st_oper2,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper1,st_oper2
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: env_oper,tmp_oper,mid_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,env_oper)
		else
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,tmp_oper)
			call update_site_dia(tmp_oper,env_bsm(idx_one),mid_oper)
			call deallocate_block(tmp_oper)

			if(idx_one<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx_one))
			endif
			call deallocate_block(mid_oper)
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(env_oper,env_bsm(idx+1),mid_oper)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(env_oper,env_bsm(idx+1),new_oper)
		call deallocate_block(env_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,new_oper)
		else
			call block_mul_block_dia(st_oper1,'N',st_oper2,'N',1.0d0,mid_oper)
			call update_site_dia(mid_oper,env_bsm(idx_two),new_oper)
			call deallocate_block(mid_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

end subroutine Get_env_oper_dia2


!==================================================================================
!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_sys_oper_ndia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,sys_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_ndia_sys(st_oper,sys_bsm(idx_one),sys_oper,Types)
			else
				call update_site_ndia_sys(st_oper,sys_bsm(idx_one),mid_oper,Types)
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),mid_oper,Types)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),new_oper,Types)
		call deallocate_block(sys_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_ndia_sys(st_oper,sys_bsm(idx_two),new_oper,Types)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

end subroutine Get_sys_oper_ndia


!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_env_oper_ndia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,env_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_ndia_env(st_oper,env_bsm(idx_one),env_oper,Types)
			else
				call update_site_ndia_env(st_oper,env_bsm(idx_one),mid_oper,Types)
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_env(env_oper,env_bsm(idx+1),mid_oper,Types)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_env(env_oper,env_bsm(idx+1),new_oper,Types)
		call deallocate_block(env_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_ndia_env(st_oper,env_bsm(idx_two),new_oper,Types)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

end subroutine Get_env_oper_ndia


!Change operator from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Change_sys_oper_ndia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,sys_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),mid_oper,Types)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),new_oper,Types)
		call deallocate_block(sys_oper)
		

	!If flag=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif
	endif

	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_sys_oper_ndia


!Change operator from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Change_env_oper_ndia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,env_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_env(env_oper,env_bsm(idx+1),mid_oper,Types)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_env(env_oper,env_bsm(idx+1),new_oper,Types)
		call deallocate_block(env_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_env_oper_ndia
