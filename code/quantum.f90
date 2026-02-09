!=================================================
!For the t-J model on the square Lattice
!=================================================
program Quantum
!include "cxml_include.f90"
	use pubdata
	implicit none

	real(8),parameter :: trun_err=1.0E-18
	integer :: i,x,y,tindex,trun_idx,levels,upindex
	real(8) :: bg_time,end_time
        !!!!integer :: sys_len_start

  call factrl

        tjmodel=0
        hubmodel=1

        realcode=.false.
        if(hubmodel==1)N_maX=3
        if(tjmodel==11)N_max =1
        if(tjmodel==1)N_maX=2
        ss=2  !!! spin=1
        ss=1

        density1=.true.
        measbl=.true.
        measli=.true.
        superc=.true.
        


	!<1>: Set Model parameter
        call dmrg_setup(num_site/2-1, num_site/2-1)
        call dmrg_setup0()
	call Get_Lattice
	call Get_site_operator_su2
	call Get_single_site_basis
	call Get_model_name(sysname,envname)
	call Get_wave_name
!	call Read_parameter

        open(1212,file='restart.dat')
        !!!!read(1212,*)total, restart2
        read(1212,*)restart1,sys_len_start
        read(1212,*)restart2
        close(1212)        

        iw6j1=0
        w6j1=0.d0
        iw6j2=0
        w6j2=0.d0
        iw6j3=0
        w6j3=0.d0


        infi_del=1
        cone=1.0d0
        czero=0.0d0
        openy=0
        !!!   2nd layer (even column has positive bias esite=Delta)
        esite=0
        do x=2,nx,2
        do y=1,ny
        i=y+(x-1)*ny
        esite(i)=6
        enddo
        enddo
        esite=0

        open(103, file='block.dat')

	!Start program
        pinz=1
        pind=1
                jd=0
                jz=0
                jn=0
                jt=0

        ci=dcmplx(0.0d0, 1.0d0)

        nrr=3
        nrr=9
        nrr=num_site-ny
        nrr=8+num_site !! pick up one triangle
        nrr= 1+3*num_site   !! all tri
        nrr=0

        pindr=0.0d0

      !  if(hubmodel==1)then
        jz=0.d0
        jt=0.d0
        jd=0.d0
        v123=0
        if(v123==0)jn=0.0d0

                jz(1:6, 1:num_site)=1.0d0
                jz(7:12, 1:num_site)=0.1d00
        jd(1:18,1:num_site)=jz(1:18,1:num_site)*(-dsqrt(3.0d0))
        jn=0

        t2=0.0d0

        if(tjmodel==1)then !!! real t-J
        hubbard=0.0d0
        jt(1:6,1:num_site )=-3.0d0*dsqrt(2.0d0)
        jt(7:12, 1:num_site)=-3.0d0*dsqrt(cdabs(jz(7,1)))*dsqrt(2.0d0)
        jn(1:12, 1:num_site)=-0.250d0*jz(1:12, 1:num_site)
        v123=1
        endif

       !!! jt=jt*0.001
        if(hubmodel==1)then
        jz=0.0d0
        jd=0.0d0
        jn=0.0d0
        v123=0
        hubbard=14.0d0
        jt(1:6,1:num_site )=-1.0d0*dsqrt(2.0d0)
        endif

        lring=0
        jring(1:num_site)=0.50d0*0.0d0!!!0.30d0
        jring(1:num_site)=jring(1:num_site)*dsqrt(6.0d0)
        if(jring(1).ne.0.0)lring=1

        
        call dmrg_setup(num_site/2-1, num_site/2-1)
        call dmrg_setup0()
        call dmrg_setup(num_site/2-1, num_site/2-1)

        theta=pi/2.0d0
                gauge=1
                msym=0
         call get_phase()
                if(msym==1)stop

		!Get the Wavefunction
		kept_max=2400
		sweep_num=4
                  iter1=0
                if(restart1.ne.0)sweep_num=3  !!! this is adjustable 
                if(restart1.ne.0)iter1=1  !!! this is adjustable  betweenn 0 
		trun_idx=Int(DLOG(kept_min*1.0d0+0.001d0)/DLOG(N_max*1.0d0))
        trun_idx=1
                 kept_min=min(3*kept_max/4, 3600)
                kept_maxf=kept_max

        if(restart1==0)then
                kept_min=min(1400, kept_max*3/4)
                kept_max=min(2400, kept_max)
                        endif


		!<1>: For the ground state
		tot_num_up=num_site/2-Num_site/16    !upindex*num_up_int
                tot_num_up=2*tot_num_up
		tot_num_down=0    !!2*6!!Num_site/4+num_site/4/4+11 !!/2-1 !Num_site/2-tot_num_up
		!Get ground state
        if(restart1.ne.4)then
		levels=1
		call Get_Ground_State(trun_idx,trun_err,levels)
		call Save_result_to_disk(Energys(1:levels),levels)
                endif

		!<2>: Measurement
		call Measurement(trun_idx)

end program Quantum


!Get ground_state wavefunction



!Get ground_state wavefunction
subroutine Get_Ground_State(trun_idx,trun_err,levels)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	integer,intent(in) :: trun_idx,levels

	logical :: newflag
	integer :: i,iter,tmp_keep,point, idx1, jk, jk1
	real(8) :: bg_time,end_time,tmp_err
	real(8),parameter :: lanerr=20.0e-6

	!<1>: For the Warm_up process
	call cpu_time(bg_time)
	open(7,file="Energy_All.dat",position='append')
	write(7,*) "Warm_up:","point=",trun_idx,kept_min, kept_max
	write(7,*) "jz, jt",jz(1:8,1), jt(1:8,1),tot_num_up, num_site
	write(7,*) "kept",kept_min, kept_max, hubbard(1)
	close(7)

	spectrum=.false.
	lanzero=lanerr*(10**3)
	tmp_err=trun_err*(10**2)

        theta=pi/2.0d0
                gauge=1
         call get_phase()

        if(restart1.eq.0)then
	call warm_up_point(trun_idx,Num_site/2,tmp_err,levels)
                endif

                 if(restart1.le.-1)then  !! call again
        call Get_Hamiltonian_Again(levels,trun_idx)
        restart1=-restart1
        endif


	newflag=.false.
	lanzero=lanerr*(15)
	tmp_err=trun_err*(15)

        if(restart1.ne.0)then
	lanzero=lanerr*2 !!*(10**2)
	tmp_err=trun_err*2 !!*(10**2)
                endif

        if(restart2.ne.0)then
	lanzero=lanerr !!*(10**2)
	tmp_err=trun_err !!*(10**2)
                endif

    if(restart1.eq.0)kept_max=kept_max*2/3+kept_maxf/3
                kept_min=min(3*kept_max/4, 3600)

        idx1=trun_idx
        if(restart1.le.1)then
        if(restart1.eq.0)sys_len_start=num_site/2-1
        restart1=restart2
call left_to_right(sys_len_start,Num_site-idx1-1,tmp_err,levels,.true.)
        restart1=1
        restart2=1
                endif

       if(restart1.eq.0)kept_max=kept_max*2/3+kept_maxf/3
                kept_min=min(3*kept_max/4, 3600)

                if(restart1.le.2)then
        if(restart1.le.1)sys_len_start=num_site-idx1
        restart1=restart2
	call right_to_left(sys_len_start-1,idx1,tmp_err,levels,.true.)
        restart1=2
                        endif

        restart1=1

	!<2>: Sweep
	iter=1


                kept_max=kept_maxf
                kept_min=min(3*kept_max/4, 3600)
        do while(iter<=sweep_num)

        open(7,file="Energy_All.dat",position='append')
                write(7,*) "Sweep=",iter, ' t2=', t2
        write(7,*)"jz, jt",jz(1,1),jz(5,1),jt(1,1),jt(5,1)
        write(7,*) "kept",kept_min, kept_max, iter
        close(7)

        lanzero=lanerr
        tmp_err=trun_err

        call left_to_right(idx1,Num_site-idx1-1,tmp_err,levels,newflag)
        call right_to_left(Num_site-idx1-1,idx1,tmp_err,levels,newflag)

                call cpu_time(end_time)
                open(7,file="Energy_All.dat",position='append')
                write(7,"(A6,F12.4,2X,A3)") "Time =",end_time-bg_time,"Sec"
                write(7,*)
                close(7)

                iter=iter+1
        end do
        spectrum=.true.
        lanzero=lanerr !*(0.1d0)
        tmp_err=trun_err !*(0.1d0)
        call left_to_right(idx1,Num_site-ny*2-1,tmp_err,levels,newflag)

end subroutine Get_Ground_State




        subroutine get_phase()
        use pubdata
                real*8 sc1(6+6), ph1, ph2
                integer jm(4)
        character ca

        pint=1
                pint0=0
        phasey=dcmplx(cos(theta),  sin(theta))
        phasey1=-1.0d0
               
        do i1=1, ny
        do i2=1, nx
        i=i1+(i2-1)*ny
        do k=1, 6,2
        j=neib(i,k)
        if(j.ne.0)then
        if(k.eq.1)then
        pint(i,j)=phasey1*(-1)**(i2) !(-1)**(i2+1) j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
             if(mod(i2,2)==0)pint0(i,k)=pi
             if(mod(i2,2)==1)pint0(i,k)=0.0d0
        endif
        if(k.eq.3)then
        pint(i,j)=phasey   !! j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
             pint0(i,k)=-pi/2.0d0 !! phasey was used, but -pi/2 is required to
                                  !! have the same flux as symmetric gauge used
        endif
        if(k.eq.5)then
        pint(i,j)=phasey1*(-1)**(i2+1)   !! j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
             if(mod(i2,2)==1)pint0(i,k)=pi
             if(mod(i2,2)==0)pint0(i,k)=0.0d0
        endif
        endif
        enddo
        enddo
        enddo

        if(gauge==2)pint=1
        phasey=dcmplx(cos(theta),  sin(theta))
        phasey1=phasey
               
        do i1=1, ny
        do i2=1, nx
        i=i1+(i2-1)*ny
        do k=1, 6,2
        j=neib(i,k)
        if(j.ne.0)then
        if(k.eq.1)then
                        if(gauge==2)then
        pint(i,j)=phasey*(-1)**(i1+i2)   !! j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
                        endif
                    if(mod(i1+i2,2)==0)then
          pint0(i,k)=pi/2.0d0-pint0(i,k)
                else
          pint0(i,k)=-pi/2.0d0-pint0(i,k)
        endif
        endif
        if(k.eq.3)then
                        if(gauge==2)then
        pint(i,j)=phasey*(-1)**(i2)   !! j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
                endif
                    if(mod(i2,2)==0)then
          pint0(i,k)=pi/2.0d0-pint0(i,k)
                else
          pint0(i,k)=-pi/2.0d0-pint0(i,k)
        endif
        endif
        if(k.eq.5)then
                        if(gaUGE==2)THEN
        pint(i,j)=phasey*(-1)**(i1)   !! j-site crosses boundary
        pint(j,i)=dconjg(pint(i,j))
                endif
                    if(mod(i1,2)==0)then
          pint0(i,k)=pi/2.0d0-pint0(i,k)
                else
          pint0(i,k)=-pi/2.0d0-pint0(i,k)
        endif
        endif
        pint0(j,k+1)=-pint0(i,k)        
                endif

        enddo
        enddo
        enddo

                open(24, file='uphase.dat')
                open(26, file='Sc_phase_sym.dat')
        uphase(1)=0.0d0
        do i2=1, nx
        do i1=1, ny
        i=i1+(i2-1)*ny
        if(i1==1.and.i2.ne.1)then
                Uphase(i)=Uphase(i-ny)+pint0(i-ny,3)
        else
                        if(i.ne.1)then
                Uphase(i)=Uphase(i-1)+pint0(i-1,1)
                        endif
        endif
        if(uphase(i).ge.2*pi)uphase(i)=uphase(i)-2*pi 
        if(uphase(i).ge.2*pi)uphase(i)=uphase(i)-2*pi 
        if(uphase(i).le.-2*pi)uphase(i)=uphase(i)+2*pi 
        if(uphase(i).le.-2*pi)uphase(i)=uphase(i)+2*pi 
          write(24,*)i1, i2, i, uphase(i)
          write(25,14)i1, i2, i, pint0(i,1), pint0(i,3)
        enddo
        enddo
        14      format(3i6, 2f16.8)


        if(msym.ne.1)return        
                read(20,*)ca
                read(20,*)ca
                read(20,*)ca
                read(20,*)ca
                read(20,*)ca
                do i=1,220
                read(20,*)jm, sc1
                j1=neib(jm(1), 1)
                        sc1(7:12)=-sc1(7:12) !! opposite flux of original model
                do k=1,6
                        if(k.le.3)then
                j2=neib(jm(2), 2*k-1)
                                else
                j2=neib(jm(2), 2*(k-3))
                                       endif
                sc1(6+k)=sc1(6+k)+uphase(jm(1))
                sc1(6+k)=sc1(6+k)+uphase(j1)
                sc1(6+k)=sc1(6+k)-uphase(jm(2))
                sc1(6+k)=sc1(6+k)-uphase(j2)
        if(sc1(6+k).le.-pi+1e-5)sc1(6+k)=sc1(6+k)+2*pi 
        if(sc1(6+k).le.-pi+1e-5)sc1(6+k)=sc1(6+k)+2*pi 
        if(sc1(6+k).gt.pi+1e-5)sc1(6+k)=sc1(6+k)-2*pi 
        if(sc1(6+k).gt.pi+1e-5)sc1(6+k)=sc1(6+k)-2*pi 
                        enddo
                write(26,12)jm, sc1(7:12)
12              format(4i6, 6f20.12)
        enddo

        stop
        return 
        end subroutine get_phase        







!Perform measurement
subroutine Measurement(trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	real(8) :: bg_time,end_time,value
	integer :: i,j,i1, j1,label,sys_len,env_len,x,y

	integer :: name_len,knum
	character(len=30) :: Filename
	type(Total_Model) :: sys,env
	type(Total_Basis) :: sys_bs,env_bs
	type(Wavefunction) :: eigvec
        character ca

	double precision,dimension(Num_site) :: ave_spin,ave_spin2
	double precision,dimension(Num_site) :: ave_num_elec_up,ave_num_elec_down,ave_spin_sd2,ave_spin_sd21

	double precision, dimension(Num_site) :: ave_num_elec,ave_num_hole,ave_num_elec2,ave_num_elec_up2
	double precision,dimension(Num_site) :: ave_spin_sz,ave_num_double,ave_spin_sz2,ave_num_elec_down2

	double complex ,dimension(Num_site,Num_site) :: num_elec_cor,spin_sz_cor,spin_sd_cor,num_hole_cor
	double complex,dimension(Num_site,Num_site) :: elec_cor,elec_down_cor,spin_cor
	double precision,dimension(Num_site,Num_site) :: elec_cor1,spin_cor1
        integer re_meas
        real*8 d11, d12
        re_meas=1


	!<1>: ============ Read data from disk =================
	!<1-1>: Read eigen_state from disk (Ground state only)
	call wave_from_disk(eigvec,1001,wave_name(1),7)
	sys_len=eigvec%sys_len
	env_len=eigvec%env_len

	!<1-2>: Read truncation operators from disk
	do i=trun_idx+1,sys_len
		call truns_from_disk(systruns(i),1001,i,.true.)
	end do

	do i=trun_idx+1,env_len
		call truns_from_disk(envtruns(i),1001,i,.false.)
	end do

	!<1-3>: Read basis from disk
	do i=2,sys_len
		call basis_from_disk(sys_bsm(i),1001,i,.true.)
	end do

	do i=2,env_len
		call basis_from_disk(env_bsm(i),1001,i,.false.)
	end do

      xi1=nx/6*ny-ny+1
      xj1=xi1
        xii=num_site
        xff=0
        if(superc)then
        call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor1.dat", xi1, xi1,0, 1)
        endif

	Filename="num_elec.dat"
	call measure_operator_dia(ave_num_elec,num_elec,eigvec,trun_idx,Filename)


                xj1=xi1
        if(measli)then
        spin_cor=0.0d0
        elec_cor=0.0d0
        num_elec_cor=0.0d0
filename='spin_cor11.dat'
        call Get_density_cor(spin_cor, st_sd,st_sd,eigvec,trun_idx,'B',Filename,xi1, xj1)
        open(251,file='spin_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(spin_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(spin_cor(i, j)), imag(spin_cor(i,j))
        endif
        enddo
        enddo
28     format(3i8, 3f21.14)

filename='elec_cor11.dat'
        call Get_density_cor(elec_cor, st_elec_up,st_elec_down,eigvec,trun_idx,'F',Filename, xi1, xj1)

        open(251,file='elec_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(elec_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(elec_cor(i, j)), imag(elec_cor(i,j))
        endif
        enddo
        enddo
filename='density_cor11.dat'
        call Get_density_cor(num_elec_cor, num_elec,num_elec,eigvec,trun_idx,'B',Filename, xi1, xj1)

        open(251,file='density_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(num_elec_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(num_elec_cor(i, j)), real(ave_num_elec(i)),real(ave_num_elec(j))
        endif
        enddo
        enddo
        endif


        if(measbl)then
        xi1=1
        xj1=1
        Filename="spin_cor1.dat"
        call Get_oper_cor_ndia1(spin_cor,st_sd,st_sd,eigvec,trun_idx,'B',Filename)
        spin_cor1=spin_cor
        Filename="strfac_spin.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)

        Filename="density_cor1.dat"
        call Get_oper_cor_ndia1(num_elec_cor,num_elec, num_elec,eigvec,trun_idx,'B',Filename)
        spin_cor1=num_elec_cor
        Filename="strfac_spin.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)

        Filename="elec_cor1.dat"
        call Get_oper_cor_ndia1(elec_cor,st_elec_up,st_elec_down,eigvec,trun_idx,'F',Filename)
        spin_cor1=elec_cor
        Filename="strfac_elec.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)


        endif


	!<6>: ========== Free space ===========
	do i=trun_idx+1,sys_len
		call deallocate_block(systruns(i))
	end do

	do i=trun_idx+1,env_len
		call deallocate_block(envtruns(i))
	end do

	do i=2,sys_len
		call deallocate_basis(sys_bsm(i))
	end do

	do i=2,env_len
		call deallocate_basis(env_bsm(i))
	end do
	call deallocate_wave(eigvec)

end subroutine Measurement

subroutine Get_Hamiltonian_Again(levels,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	integer :: i,label,sys_len,env_len,slen,elen,levels
	type(Total_Model) :: sys,env,new_sys,new_env
	type(Wavefunction) :: eigvec


	!<1>: ============ Read data from disk =================
	!<1-1>: Read eigen_state from disk (Ground state only)
	call wave_from_disk(eigvec,1001,wave_name(1),7)
	sys_len=eigvec%sys_len
	env_len=eigvec%env_len
        slen=max(num_site0/2,sys_len)
        elen=env_len

	!<1-2>: Read truncation operators from disk
	do i=trun_idx+1,slen
		call truns_from_disk(systruns(i),1001,i,.true.)
	end do

	do i=trun_idx+1,elen
		call truns_from_disk(envtruns(i),1001,i,.false.)
	end do

	!<1-3>: Read basis from disk
	do i=2,slen
		call basis_from_disk(sys_bsm(i),1001,i,.true.)
	end do

	do i=2,elen
		call basis_from_disk(env_bsm(i),1001,i,.false.)
	end do
	
	!<2>: Get Hamiltonian for system part
	!<2-1>: Get Hamiltonian without truncation
	label=174

                if(label==1)then
	call Get_hamiltonian_site(sys)
        write(*,*)sys%ham%sub(1)%mat, sys%ham%len,'bl1'
          do i=1,sys%ham%num
      sys%ham%sub(i)%mat=esite(1)*num_elec%sub(i)%mat
      sys%ham%sub(i)%mat=sys%ham%sub(i)%mat+hubbard(num_site)*st_double%sub(i)%mat
        end do
	call model_to_disk(sys,1001,label,.true.)
                        endif
                if(label.ne.1)then
	call model_from_disk(sys,1001,label,.true.)
                        endif

	do while((label+1)<=trun_idx)
                          call dmrg_setup(label,num_site-2-label)
                                mshift=0
		call Get_hamiltonian_sys(sys,new_sys,sys_bsm(label+1),systruns(label+1))
		call deallocate_model(new_sys)
		call model_to_disk(sys,1001,label+1,.true.)
		label=label+1
                end do
	
	!<2-2>: Get Hamiltonian with truncation
	do while((label+1)<=slen)
                          call dmrg_setup(label,num_site-2-label)
                                mshift=0
		call Get_hamiltonian_sys(sys,new_sys,sys_bsm(label+1),systruns(label+1))
                                write(*,*)sys%len, sys%ham%sub(1)%mat,'bl2'
!!		call deallocate_model(sys)
!!		call Get_hamiltonian_trun(new_sys,sys,systruns(label+1))
		call deallocate_model(new_sys)
		call model_to_disk(sys,1001,label+1,.true.)
		label=label+1
	end do
	call deallocate_model(sys)
        stop



	
	!<3>: Get Hamiltonian for environment part
	!<3-1>: Get Hamiltonian without truncation
22      continue
	label=83

                if(label.eq.1)then
	call Get_hamiltonian_site(env)
        if(esite(num_site).ne.0.0)then
          do i=1,env%ham%num
                env%ham%sub(i)%mat=esite(num_site)*num_elec%sub(i)%mat
!!!old_block%sub(i)%mat
        end do
        endif
          do i=1,env%ham%num
                env%ham%sub(i)%mat=env%ham%sub(i)%mat+hubbard(num_site)*st_double%sub(i)%mat
        end do
                else
	call Get_hamiltonian_site(env)
                call deallocate_model(env)
	call model_from_disk(env,1001,label,.false.)
                endif

                        write(*,*)'here', label
	do while((label+1)<=trun_idx)
                          call dmrg_setup(num_site-2-label, label)
                                mshift=num_site-2-label+1
	call Get_hamiltonian_env(env,new_env,env_bsm(label+1), envtruns(label+1))
!!		call deallocate_model(env)
!!		call model_transfer(new_env,env)
		call deallocate_model(new_env)
		call model_to_disk(env,1001,label+1,.false.)
		label=label+1
	end do
	
	!<3-2>: Get Hamiltonian with truncation
	do while((label+1)<=108)  !!!env_len)
                          call dmrg_setup(num_site-2-label, label)
                                mshift=num_site-2-label+1
	call Get_hamiltonian_env(env,new_env,env_bsm(label+1), envtruns(label+1))
!		call deallocate_model(env)
!		call Get_hamiltonian_trun(new_env,env,envtruns(label+1))
		call deallocate_model(new_env)
		call model_to_disk(env,1001,label+1,.false.)
		label=label+1
	end do
	call deallocate_model(env)
	

	!<4>: ========== Free space ===========
        call deallocate_wave(eigvec)
	do i=trun_idx+1,sys_len
		call deallocate_block(systruns(i))
	end do

	do i=trun_idx+1,env_len
		call deallocate_block(envtruns(i))
	end do

	do i=2,sys_len
		call deallocate_basis(sys_bsm(i))
	end do

	do i=2,env_len
		call deallocate_basis(env_bsm(i))
	end do

end subroutine Get_Hamiltonian_Again


subroutine Get_Hamiltonian_Again1(levels,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	integer :: i,label,sys_len,env_len, levels
	type(Total_Model) :: sys,env,new_sys,new_env
	type(Wavefunction) :: eigvec


	!<1>: ============ Read data from disk =================
	!<1-1>: Read eigen_state from disk (Ground state only)

	!<1-2>: Read truncation operators from disk
	do i=trun_idx+1,num_site/2+2
		call truns_from_disk(envtruns(i),1001,i,.true.)
	end do
	do i=2,num_site/2+2
		call basis_from_disk(env_bsm(i),1001,i,.true.)
	end do

	label=1
                do label=1, num_site/2+2
	call model_from_disk(env,1001,label,.true.)
	call model_to_disk(env,1001,label,.false.)
	call deallocate_model(env)
                i=label
                        if(i.ge.2)then
		call truns_to_disk(envtruns(i),1001,i,.false.)
		call basis_to_disk(env_bsm(i),1001,i,.false.)
                        endif
	end do
	
	

	do i=trun_idx+1,num_site/2+2
		call deallocate_block(envtruns(i))
		call deallocate_basis(env_bsm(i))
	end do
end subroutine Get_Hamiltonian_Again1


!Save results to the disk
subroutine Save_result_to_disk(Eigvals,levels)
	use pubdata
	implicit none

	integer,intent(in) :: levels
	real(8),intent(in) :: Eigvals(levels)

	integer :: i

	open(10,file="Spectrum.dat",position='append')
		write(10,*) "t=",jt(1,1),"N_up=",tot_num_up,"N_down=",tot_num_down&
					&,"M1=",kept_min,"M2=",kept_max,"E0=",(Eigvals(i),i=1,levels),"S=",Entropys
	close(10)

	open(19,file="Spectrum_copy.dat",position='append')
            write(19,*) Eigvals(1)
    close(19)

end subroutine Save_result_to_disk
