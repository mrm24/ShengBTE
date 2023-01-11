!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2023 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2023 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2023 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2023 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!  Copyright (C) 2021-2022 Fanchen Meng <fanchem@g.clemson.edu>
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Main program.
program ShengBTE
  use data
  use misc
  use config
  use input
  use processes
  use wedgetc
  use iterations
  use conductivity
  use scaling
  use phonon_routines
  use dos_routines
  use gruneisen
  use integrals
  use mpi_f08
  !$ use omp_lib
  !$ use omp_lib_kinds
  implicit none

  real(kind=dp) :: kappa_sg(3,3),kappa_old(3,3),relchange
  integer :: i,j,ii,jj,kk,ll,mm,nn
  integer :: Tcounter
  real(kind=dp),allocatable :: energy(:,:),q0(:,:),q0_reduced(:,:),velocity(:,:,:)
  complex(kind=dp),allocatable :: eigenvect(:,:,:)

  real(kind=dp),allocatable :: grun(:,:)

  real(kind=dp),allocatable :: rate_scatt(:,:),rate_scatt_plus(:,:),rate_scatt_minus(:,:)
  real(kind=dp),allocatable :: tau_zero(:,:),tau(:,:),tau_b(:,:),tau2(:,:)
  real(kind=dp),allocatable :: dos(:),pdos(:,:),rate_scatt_isotope(:,:)
  real(kind=dp),allocatable :: F_n(:,:,:),F_n_0(:,:,:),F_n_aux(:,:)
  real(kind=dp),allocatable :: ThConductivity(:,:,:)
  real(kind=dp),allocatable :: ThConductivityMode(:,:,:,:)
  real(kind=dp),allocatable :: ticks(:),cumulative_kappa(:,:,:,:)

  integer :: Ntri
  ! IJK are vectors of product(ngrid) points in the lattice coordinate
  ! system, with each component ranging from 0 to Ngrid(:)-1.  Index_N
  ! is a mapping of 3 indices for an individual q-point into one index.
  integer,allocatable :: Index_i(:),Index_j(:),Index_k(:),IJK(:,:),Index_N(:,:,:)
  real(kind=dp),allocatable :: Phi(:,:,:,:),R_j(:,:),R_k(:,:)

  integer :: nlist,Ntotal_plus,Ntotal_minus
  integer,allocatable :: nequi(:),list(:)
  integer,allocatable :: AllEquiList(:,:),TypeofSymmetry(:,:),eqclasses(:)
  integer,allocatable :: N_plus(:),N_minus(:)
  integer :: Naccum_plus, Naccum_minus
  real(kind=dp) :: radnw,kappa_or_old
  real(kind=dp),allocatable :: Pspace_plus_total(:,:)
  real(kind=dp),allocatable :: Pspace_minus_total(:,:)
  real(kind=dp),allocatable :: ffunc(:,:),radnw_range(:),v_or(:,:),F_or(:,:)
  real(kind=dp),allocatable :: kappa_or(:),kappa_wires(:,:)

  integer :: iorient,ierr
  character(len=4) :: aux,aux2
  character(len=1024) :: path
  character(len=128) :: sorientation
  integer :: max_nthreads

  real(kind=dp) :: dnrm2

  call MPI_INIT_thread(MPI_THREAD_FUNNELED,i,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  ! Read CONTROL and initialize variables such as the
  ! direct/reciprocal lattice vectors, the symmetry operations, the
  ! Pearson deviation coefficients related to isotopic disorder and so
  ! on and so forth.
  call read_config()

  !$ max_nthreads = omp_get_max_threads()
  !$ if (max_nthreads < nthreads) then
  !$   if (myid.eq.0) then
  !$      write(*,'(2(A,I0),A)') " Info: Number of threads requested (", nthreads, &
  !$                 &") larger than system value of OMP_NUM_THREADS (",max_nthreads,")"
  !$      write(*,*) "Info: OMP_NUM_THREADS has overriden requested number of threads"
  !$   end if
  !$   nthreads = max_nthreads
  !$ end if
  !$ ! if negative value passed for nthreads default to OMP_NUM_THREADS
  !$ if (nthreads > 0) then
  !$    call omp_set_num_threads(nthreads)
  !$ else
  !$    nthreads = max_nthreads
  !$ end if

  if((myid == 0) .and. (nthreads > 1)) then
     write(*,'(A,I0,A)') " Info: Each MPI process will use ", nthreads," thread(s)"
  end if

  allocate(energy(nbands,nptk))
  allocate(eigenvect(nbands,nbands,nptk))
  allocate(grun(nbands,nptk))
  allocate(q0(nptk,3),q0_reduced(nptk,3))
  allocate(velocity(3,nbands,nptk))
  allocate(IJK(3,nptk))
  allocate(Index_N(0:(ngrid(1)-1),0:(ngrid(2)-1),0:(ngrid(3)-1)))
  allocate(F_n(3,nbands,nptk))
  allocate(F_n_0(3,nbands,nptk),F_n_aux(nbands,nptk))
  allocate(ThConductivity(3,3,nbands))
  allocate(ThConductivityMode(3,3,nbands,nptk))
  allocate(kappa_wires(nbands,nwires))
  allocate(Nequi(nptk))
  allocate(List(nptk))
  allocate(AllEquiList(nsymm_rot,nptk))
  allocate(TypeOfSymmetry(nsymm_rot,nptk))
  allocate(eqclasses(nptk))
  if (nanowires) then
    allocate(ffunc(nbands,nptk))
    allocate(v_or(nbands,nptk))
    allocate(F_or(nbands,nptk))
    allocate(kappa_or(nbands))
  end if

  ! Obtain the q-point equivalence classes defined by symmetry
  ! operations.
  call Id2Ind(IJK)
  call wedge(Nlist,Nequi,List,ALLEquiList,TypeofSymmetry)
  do ll=1,Nlist
     do kk=1,Nequi(ll)
        eqclasses(ALLEquiList(kk,ll))=ll
     end do
  end do


  do ii=1,Ngrid(1)
     do jj=1,Ngrid(2)
        do kk=1,Ngrid(3)
           ll=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
           q0(ll,:)=rlattvec(:,1)*(ii-1)/real(ngrid(1),kind=dp)+&
                rlattvec(:,2)*(jj-1)/real(ngrid(2),kind=dp)+&
                rlattvec(:,3)*(kk-1)/real(ngrid(3),kind=dp)
           q0_reduced(ll,1)=real(ii-1,kind=dp)/real(ngrid(1),kind=dp)
           q0_reduced(ll,2)=real(jj-1,kind=dp)/real(ngrid(2),kind=dp)
           q0_reduced(ll,3)=real(kk-1,kind=dp)/real(ngrid(3),kind=dp)
        end do
     end do
  end do

  if(myid.eq.0) then
     open(1,file="BTE.ReciprocalLatticeVectors",status="replace")
     write(1,"(3F20.10,A12)") rlattvec(:,1), "# nm-1,b1" 
     write(1,"(3F20.10,A12)") rlattvec(:,2), "# nm-1,b2" 
     write(1,"(3F20.10,A12)") rlattvec(:,3), "# nm-1,b3" 
     close(1)
     open(1,file="BTE.qpoints",status="replace")
     do ll=1,Nlist
        write(1,"(3(I9,1x),3(E20.10,1x))") ll, List(ll),Nequi(ll),q0_reduced(List(ll),:)
     end do
     close(1)
     open(1,file="BTE.qpoints_full",status="replace")
     do ll=1,nptk
        write(1,"(2(I9,1x),3(E20.10,1x))") ll,eqclasses(ll),q0_reduced(ll,:)
     end do
     close(1)
  end if

  ! Obtain the phonon spectrum from the 2nd-order force constants (in
  ! either of the two formats suported) and the dielectric parameters
  ! in CONTROL.
  if(myid.eq.0) then
     write(*,*) "Info: about to obtain the spectrum"
     if(espresso) then
        write(*,*) "Info: expecting Quantum Espresso 2nd-order format"
     else
        write(*,*) "Info: expecting Phonopy 2nd-order format"
     end if
  end if
  call eigenDM(energy,eigenvect,velocity)
  if(myid.eq.0)write(*,*) "Info: spectrum calculation finished"

  open(101,file="BTE.cvVsT")
  open(102,file="BTE.KappaTensorVsT_sg")
  if (myid.eq.0) print*, "Info: start calculating specific heat and kappa in the small-grain limit "
  do Tcounter=1,ceiling((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.0_dp)))  exit 
     if (T.gt.(T_max+1.0_dp)) T=T_max 
     if (myid.eq.0) then
        print *, "Info: Temperature=", T
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call create_directory(trim(adjustl(path))//C_NULL_CHAR)
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
        ! Compute the harmonic integrals: lattice specific heat and
        ! small-grain-limit reduced thermal conductivity. Write out this
        ! information, as well as the spectrum itself.
        open(1,file="BTE.cv",status="replace")
        write(1,*) cv(energy)
        close(1)
        call kappasg(energy,velocity,kappa_sg)
        open(1,file="BTE.kappa_sg",status="replace")
        write(1,"(9E14.5)") kappa_sg
        close(1)
        write(101,"(F7.1,E14.5)") T,cv(energy)
        flush(101)
        write(102,"(F7.1,9E14.5)") T,kappa_sg
        flush(102)
        call change_directory(".."//C_NULL_CHAR)
     end if
  enddo
  close(101)
  close(102)

  write(aux,"(I0)") 3*Nbands
  if(myid.eq.0) then
     open(1,file="BTE.v",status="replace")
     do ii=1,Nbands
        do ll=1,Nlist
           write(1,"(3E20.10)") velocity(:,ii,list(ll))
        end do
     end do
     close(1)
     open(1,file="BTE.v_full",status="replace")
     do ii=1,Nbands
        do ll=1,nptk
           write(1,"(3E20.10)") velocity(:,ii,ll)
        end do
     end do
  end if
  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.omega",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") energy(:,list(ll))
     end do
     close(1)
  end if

  ! Compute the normalized boundary scattering rates.
  allocate(tau_b(Nbands,Nlist))
  !$OMP PARALLEL default(none) private(ll,ii) &
  !$OMP & shared(nlist,nbands,nptk,tau_b,velocity,list)

  !$OMP DO collapse(2) schedule(static)
  do ll=1,Nlist
     do ii=1,Nbands
        tau_b(ii,ll)=1.0_dp/dnrm2(3,velocity(:,ii,List(ll)),1)
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  if (myid.eq.0) then
     write(aux,"(I0)") Nbands
     open(2,file="BTE.w_boundary",status="replace")
     do i=1,Nbands
        do ll = 1,Nlist
           write(2,"(2E20.10)") energy(i,list(ll)),  1.0_dp/tau_b(i,ll)
        enddo
     end do
     close(2)
  endif
  deallocate(tau_b)

  ! Locally adaptive estimates of the total and projected densities of states.
  allocate(ticks(nticks))
  allocate(dos(nticks))
  allocate(pdos(nticks,natoms))
  call calc_dos(energy,velocity,eigenvect,ticks,dos,pdos)
  if(myid.eq.0) then
     open(1,file="BTE.dos",status="replace")
     do mm=1,nticks
        write(1,"(2E14.5)") ticks(mm),dos(mm)
     end do
     close(1)
     write(aux,"(I0)") Natoms
     open(1,file="BTE.pdos",status="replace")
     do mm=1,nticks
        write(1,"(E14.5,"//trim(adjustl(aux))//"E14.5)")&
             ticks(mm),pdos(mm,:)
     end do
     close(1)
  end if
  deallocate(ticks,dos,pdos)


  allocate(rate_scatt_isotope(Nbands,Nlist))
  if (isotopes) then
     call calc_isotopescatt(energy,velocity,eigenvect,nlist,list,&
          rate_scatt_isotope)
  else
     rate_scatt_isotope=0.0_dp
  end if

  if(myid.eq.0) then
     write(aux,"(I0)") Nbands
     open(1,file="BTE.w_isotopic",status="replace")
     do i=1,Nbands
        do ll=1,Nlist
           write(1,"(2E20.10)") energy(i,list(ll)),rate_scatt_isotope(i,ll)
        enddo
     end do
     close(1)
  end if

  nstates=ceiling(real(nlist*nbands,kind=dp)/real(numprocs,kind=dp))
  allocate(rate_scatt(Nbands,Nlist))
  allocate(rate_scatt_plus(Nbands,Nlist))
  allocate(rate_scatt_minus(Nbands,Nlist))
  allocate(tau_zero(Nbands,Nlist))
  allocate(tau(Nbands,Nlist))
  allocate(tau2(Nbands,nptk))
  allocate(radnw_range(nwires))
  do ii=1,nwires
     radnw_range(ii)=rmin+(ii-1.0_dp)*dr
  end do

  ! N_plus: number of allowed absorption processes
  ! N_minus: number of allowed emission processes.
  allocate(N_plus(Nlist*Nbands))
  allocate(N_minus(Nlist*Nbands))
  ! Phase space volume per mode and their sum.
  allocate(Pspace_plus_total(Nbands,Nlist))
  allocate(Pspace_minus_total(Nbands,Nlist))
  
  call NP_driver(energy,velocity,Nlist,List,IJK,&
       N_plus,Pspace_plus_total,N_minus,Pspace_minus_total)

  Ntotal_plus=sum(N_plus)
  Ntotal_minus=sum(N_minus)

  if(myid.eq.0)write(*,*) "Info: Ntotal_plus =",Ntotal_plus
  if(myid.eq.0)write(*,*) "Info: Ntotal_minus =",Ntotal_minus


  allocate(Naccum_plus_array(nstates))
  allocate(Naccum_minus_array(nstates))
  Naccum_plus=0
  Naccum_minus=0
  Naccum_plus_array = 0
  Naccum_minus_array = 0
  do nn=1,nstates
      mm=myid*nstates+nn
      if (mm.gt.nlist*nbands) cycle
      ! keep track of how many processes are required in the for the states 1->(nn-1)
      ! to allow accurate array accessing across MPI processes and threads
      Naccum_plus_array(nn) = Naccum_plus
      Naccum_minus_array(nn) = Naccum_minus

      Naccum_plus=Naccum_plus+N_plus(mm)
      Naccum_minus=Naccum_minus+N_minus(mm)
  enddo


  if(myid.eq.0) then
     open(1,file="BTE.P3_plus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_plus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3_minus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_minus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") 2.0_dp*(Pspace_plus_total(:,ll)+&
             Pspace_minus_total(:,ll)*0.5_dp)/3.0_dp
     end do
     close(1)
     do ii=1,Nlist
        Pspace_plus_total(:,ii)=Pspace_plus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_plus_total",status="replace")
     write(1,*) sum(Pspace_plus_total)
     close(1)
     do ii=1,Nlist
        Pspace_minus_total(:,ii)=Pspace_minus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_minus_total",status="replace")
     write(1,*) sum(Pspace_minus_total)
     close(1)
     open(1,file="BTE.P3_total",status="replace")
     write(1,*) 2.0_dp*(sum(Pspace_plus_total)+sum(Pspace_minus_total)*0.5_dp)/3.0_dp
     close(1)
  end if

  if(onlyharmonic) then
! weighted phase space (WP3) is calculated here if onlyharmonic=.true., otherwise WP3 will be calculated later together with BTE.w_anharmonic 
  if (myid.eq.0) print*, "Info: start calculating weighted phase space"
  do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.0_dp)))  exit 
     if (T.gt.(T_max+1.0_dp)) T=T_max 
     if (myid.eq.0) then
        print *, "Info: Temperature=", T
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
     endif
     call RTA_driver(energy,velocity,Nlist,List,IJK,rate_scatt,&
          rate_scatt_plus,rate_scatt_minus,Pspace_plus_total,Pspace_minus_total)
     if(myid.eq.0) then
        open(1,file="BTE.WP3_plus",status="replace")
        do ll=1,Nlist
        do i=1,Nbands
              write(1,"(2E14.5)") energy(i,list(ll)),Pspace_plus_total(i,ll)
           enddo
        end do
        close(1)
        open(1,file="BTE.WP3_minus",status="replace")
        do ll=1,Nlist
        do i=1,Nbands
              write(1,"(2E14.5)") energy(i,list(ll)),Pspace_minus_total(i,ll)
           enddo
        end do
        close(1)
        open(1,file="BTE.WP3",status="replace")
        do ll=1,Nlist
        do i=1,Nbands
              write(1,"(2E14.5)") energy(i,list(ll)),&
                   Pspace_plus_total(i,ll)+Pspace_minus_total(i,ll)
           enddo
        end do
        close(1)
     end if
  enddo
  endif  ! onlyharmonic
  deallocate(Pspace_plus_total)
  deallocate(Pspace_minus_total)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Up to this point, no anharmonic information is used.
  if(onlyharmonic) then
     write(*,*) "Info: onlyharmonic=.true., stopping here"
     call MPI_FINALIZE(ierr)
     stop
  end if

  ! Load the anharmonic IFCs from FORCE_CONSTANTS_3RD.
  call read3fc(Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)

  call mode_grun(energy,eigenvect,Ntri,Phi,R_j,R_k,&
       Index_i,Index_j,Index_k,grun)

  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.gruneisen",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") grun(:,list(ll))
     end do
     close(1)
  end if
  open(101,file="BTE.gruneisenVsT_total")
  do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.0_dp)))  exit 
     if (T.gt.(T_max+1.0_dp)) T=T_max 
     if (myid.eq.0) then
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
        open(1,file="BTE.gruneisen_total",status="replace")
        write(1,*) total_grun(energy,grun)
        close(1)
        write(101,"(F7.1,E14.5)") T,total_grun(energy,grun)
        flush(101)
        call change_directory(".."//C_NULL_CHAR)
     endif
  enddo
  close(101)
  deallocate(grun)

  if (myid.eq.0) write(*,*) "Info: max(N_plus), max(N_minus)", MAXVAL(N_plus), MAXVAL(N_minus)
  if (myid.eq.0) write(*,*) "Info: calculating Vp_plus and Vp_minus"
  call calculate_Vp(energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,Naccum_plus,Naccum_minus)

  ! This is the most expensive part of the calculation: obtaining the
  ! three-phonon scattering amplitudes for all allowed processes.
  ! When the iterative solution to the full linearized BTE is not
  ! requested (i.e., when the relaxation-time approximation is
  ! enough) we use optimized routines with a much smaller memory footprint.
  ! weighted phase space volume per mode .

  allocate(Pspace_plus_total(Nbands,Nlist))
  allocate(Pspace_minus_total(Nbands,Nlist))
  open(303,file="BTE.KappaTensorVsT_RTA")
  if(convergence) then
     open(403,file="BTE.KappaTensorVsT_CONV")
     allocate(Indof2ndPhonon_plus(Naccum_plus))
     allocate(Indof3rdPhonon_plus(Naccum_plus))
     allocate(Indof2ndPhonon_minus(Naccum_minus))
     allocate(Indof3rdPhonon_minus(Naccum_minus))
     allocate(Gamma_plus(Naccum_plus))
     allocate(Gamma_minus(Naccum_minus))
     Indof2ndPhonon_plus=0 
     Indof3rdPhonon_plus=0
     Indof2ndPhonon_minus=0
     Indof3rdPhonon_minus=0
     Gamma_plus=0.0_dp
     Gamma_minus=0.0_dp
  endif
  if (myid.eq.0) print*, "Info: start calculating kappa"
  do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.0_dp)))  exit 
     if (T.gt.(T_max+1.0_dp)) T=T_max 
     if (myid.eq.0) then
        print *, "Info: Temperature=", T
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
     endif

     rate_scatt=0.0_dp
     rate_scatt_plus=0.0_dp
     rate_scatt_minus=0.0_dp
     Pspace_plus_total=0.0_dp
     Pspace_minus_total=0.0_dp

     if(convergence) then
        !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nstates,myid,Nbands,nlist) &
        !$OMP & shared(energy,velocity,IJK,list,omega_max) &
        !$OMP & shared(N_plus,Naccum_plus_array,Gamma_plus,rate_scatt_plus,Pspace_plus_total) &
        !$OMP & shared(N_minus,Naccum_minus_array,Gamma_minus,rate_scatt_minus,Pspace_minus_total) &
        !$OMP & private(i,ll,mm,nn,Naccum_plus,Naccum_minus)
        do nn=1,nstates
            mm=myid*nstates+nn
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1

            if (mm.gt.nlist*nbands) cycle

            Naccum_plus=Naccum_plus_array(nn)
            Naccum_minus=Naccum_minus_array(nn)

            if (energy(i,list(ll)) /= 0.0_dp .and. energy(i,list(ll)) .le. omega_max) then
              call Ind_driver(nn,energy,velocity,Nlist,List,IJK,&
                 N_plus,N_minus,Naccum_plus,Naccum_minus, &
                 rate_scatt_plus(i,ll),rate_scatt_minus(i,ll),&
                 Pspace_plus_total(i,ll),Pspace_minus_total(i,ll))
            end if
        enddo 
        !$OMP END PARALLEL DO
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (myid == 0) then
          call MPI_REDUCE(MPI_IN_PLACE,rate_scatt_plus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(MPI_IN_PLACE,rate_scatt_minus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(MPI_IN_PLACE,Pspace_plus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(MPI_IN_PLACE,Pspace_minus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
        else
          call MPI_REDUCE(rate_scatt_plus,rate_scatt_plus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(rate_scatt_minus,rate_scatt_minus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(Pspace_plus_total,Pspace_plus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
          call MPI_REDUCE(Pspace_minus_total,Pspace_minus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
               MPI_SUM,0,MPI_COMM_WORLD,ll)
        end if
        rate_scatt=rate_scatt_plus+rate_scatt_minus
     else
        call RTA_driver(energy,velocity,Nlist,List,IJK,rate_scatt,rate_scatt_plus,&
             rate_scatt_minus,Pspace_plus_total,Pspace_minus_total)
     end if

     ! Passing in rate_scatt rather than rate_scatt(1,1) in some MPI implementations
     ! can result in rate_scatt being shifted from 
     !    rate_scatt(1:nbands,1:nptk) -> rate_scatt(0:nbands-1,0:nptk-1)
     ! this results in all sorts memory shenanigans which is best avoided.
     call MPI_BCAST(rate_scatt(1,1),Nbands*Nlist,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ll)

     if(myid.eq.0) then
        open(1,file="BTE.WP3_plus",status="replace")
        do i=1,Nbands
           do ll=1,Nlist
              write(1,"(2E14.5)") energy(i,list(ll)),Pspace_plus_total(i,ll)
           enddo
        end do
        close(1)
        open(1,file="BTE.WP3_minus",status="replace")
        do i=1,Nbands
           do ll=1,Nlist
              write(1,"(2E14.5)") energy(i,list(ll)),Pspace_minus_total(i,ll)
           enddo
        end do
        close(1)
        open(1,file="BTE.WP3",status="replace")
        do i=1,Nbands
           do ll=1,Nlist
              write(1,"(2E14.5)") energy(i,list(ll)),Pspace_plus_total(i,ll)+Pspace_minus_total(i,ll)
           enddo
        end do
        close(1)
     end if

     write(aux,"(I0)") Nbands
     if(myid.eq.0) then
        open(1,file="BTE.w_anharmonic",status="replace")
        open(2,file="BTE.w_anharmonic_plus",status="replace")
        open(3,file="BTE.w_anharmonic_minus",status="replace")
        do i=1,Nbands
           do ll=1,Nlist
              write(1,"(2E20.10)") energy(i,list(ll)),rate_scatt(i,ll)
              write(2,"(2E20.10)") energy(i,list(ll)),rate_scatt_plus(i,ll)
              write(3,"(2E20.10)") energy(i,list(ll)),rate_scatt_minus(i,ll)
           enddo
        end do
        close(1)
        close(2)
        close(3)
     end if

     ! Obtain the total scattering rates in the relaxation time approximation.
     rate_scatt=rate_scatt+rate_scatt_isotope
     deallocate(rate_scatt_isotope)
     if(myid.eq.0) then
        open(1,file="BTE.w",status="replace")
        do i=1,Nbands
           do ll = 1,Nlist
              write(1,"(2E20.10)") energy(i,list(ll)),rate_scatt(i,ll)
           enddo
        end do
        close(1)
     end if

     tau_zero=0.0_dp
     !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
     !$OMP & shared(nlist,nbands,rate_scatt,tau_zero) private(ll,i)
     do ll = 1,Nlist
        do i=1,Nbands
           if(rate_scatt(i,ll)>1.0e-8_dp) then
              tau_zero(i,ll)=1.0_dp/rate_scatt(i,ll)
           end if
        end do
     end do
     !$OMP END PARALLEL DO


     ! Set up everything to start the iterative process.
     call iteration0(Nlist,Nequi,ALLEquiList,energy,velocity,tau_zero,F_n)
     F_n_0=F_n
     if(myid.eq.0) then
        ! Open all output files.
        open(2001,file="BTE.kappa",status="replace")
        open(2002,file="BTE.kappa_tensor",status="replace")
        open(2003,file="BTE.kappa_scalar",status="replace")
        call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
        !$OMP PARALLEL DO default(none) schedule(static) private(ll) shared(nbands,ThConductivity)
        do ll=1,nbands
           call symmetrize_tensor(ThConductivity(:,:,ll))
        end do
        !$OMP END PARALLEL DO
        write(aux,"(I0)") 9*Nbands
        write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") 0,ThConductivity
        write(2002,"(I9,9E20.10,E20.10)") 0,sum(ThConductivity,dim=3)
        write(2003,"(I9,E20.10,E20.10)") 0,&
             sum(sum(ThConductivity,dim=3),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.0_dp
        write(303,"(F7.1,9E14.5)") T,sum(ThConductivity,dim=3)
        flush(303)
     endif
     ! parallel iteration
     ! Iterate to convergence if desired.
     if(convergence) then
        do ii=1,maxiter
           call iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
              & energy,velocity,tau_zero,F_n)
           if(myid.eq.0) then
               kappa_old=sum(ThConductivity,dim=3)
                ! Correct F_n to prevent it drifting away from the symmetry of the system.
                do ll=1,nptk
                   F_n(:,:,ll)=matmul(symmetrizers(:,:,ll),F_n(:,:,ll))
                end do
                call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
                !$OMP PARALLEL DO default(none) schedule(static) private(ll) shared(nbands,ThConductivity)
                do ll=1,nbands
                   call symmetrize_tensor(ThConductivity(:,:,ll))
                end do
                !$OMP END PARALLEL DO
                write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") ii,ThConductivity
                flush(2001)
                write(2002,"(I9,9E20.10)") ii,sum(ThConductivity,dim=3)
                flush(2002)
                write(2003,"(I9,E20.10)") ii,&
                     sum(sum(ThConductivity,dim=3),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.0_dp
                flush(2003)
                relchange=twonorm3x3(sum(ThConductivity,dim=3)-kappa_old)/&
                     twonorm3x3(kappa_old)
                write(*,*) "Info: Iteration",ii
                write(*,*) "Info:","Relative change","=",relchange
           endif
           call MPI_BCAST(relchange,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           if(relchange.lt.eps)exit
        end do
        if (myid.eq.0) then
        write(403,"(F7.1,9E14.5,I6)") T,sum(ThConductivity,dim=3),ii
        flush(403)
        endif
     end if
     if (myid.eq.0) then
        close(2001)
        close(2002)
        close(2003)

        ! Write out the converged scattering rates.
        !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
        !$OMP & shared(nlist,nbands,F_n,tau,velocity,energy,List) private(ll,ii)
        do ll=1,Nlist
           do ii=1,Nbands
              tau(ii,ll)=dot_product(F_n(:,ii,List(ll)),velocity(:,ii,List(ll)))/&
                   (dot_product(velocity(:,ii,List(ll)),velocity(:,ii,List(ll)))*energy(ii,List(ll)))
           end do
        end do
        !$OMP END PARALLEL DO
        write(aux,"(I0)") Nbands
        open(1,file="BTE.w_final",status="replace")
        do i=1,Nbands
           do ll = 1,Nlist
              write(1,"(2E20.10)") energy(i,list(ll)),  1.0_dp/tau(i,ll)
           enddo
        end do
        close(1)
        !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
        !$OMP & shared(nptk,nbands,F_n,tau2,velocity,energy) private(ll,ii)
        do ll=1,nptk
           do ii=1,Nbands
              tau2(ii,ll)=dot_product(F_n(:,ii,ll),velocity(:,ii,ll))/&
                    (dot_product(velocity(:,ii,ll),velocity(:,ii,ll))*energy(ii,ll))
           end do
        end do
        !$OMP END PARALLEL DO
        write(aux,"(I0)") Nbands

        ! If results for nanowires have been requested, obtain a lower bound
        ! for the thermal conductivity along each crystallographic orientation
        ! by using the bulk RTA results.
        if(nanowires) then
           do iorient=1,norientations
              write(sorientation,"(I128)") iorient
              open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation))//"_lower",status="replace")
              !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
              !$OMP & private(ii,jj) shared(nbands,nptk,velocity,v_or,F_or,F_n_0,uorientations,iorient)
              do ii=1,nptk
                 do jj=1,Nbands
                    v_or(jj,ii)=dot_product(velocity(:,jj,ii),uorientations(:,iorient))
                    F_or(jj,ii)=dot_product(F_n_0(:,jj,ii),uorientations(:,iorient))
                 end do
              end do
              !$OMP END PARALLEL DO
              do mm=1,Nwires
                 radnw=radnw_range(mm)
                 call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
                 !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
                 !$OMP & private(ii,jj) shared(nbands,nptk,F_n_aux,F_or,ffunc)
                 do ii=1,nptk
                    do jj=1,Nbands
                       F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(jj,ii)
                    end do
                 end do
                 !$OMP END PARALLEL DO
                 call TConductScalar(energy,v_or,F_n_aux,kappa_or)
                 write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.0_dp*radnw,&
                      kappa_or,sum(kappa_or)
              end do
              close(3001)
           end do
        end if
        allocate(ticks(nticks),cumulative_kappa(3,3,nbands,nticks))
        ! Cumulative thermal conductivity.
        call CumulativeTConduct(energy,velocity,F_n,ticks,cumulative_kappa)
        !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
        !$OMP & private(ii,ll) shared(nbands,nticks,cumulative_kappa)
        do ii=1,nticks
           do ll=1,nbands
              call symmetrize_tensor(cumulative_kappa(:,:,ll,ii))
           end do
        end do
        !$OMP END PARALLEL DO
        write(aux,"(I0)") 9*nbands+1
        open(2002,file="BTE.cumulative_kappa_tensor",status="replace")
        open(2003,file="BTE.cumulative_kappa_scalar",status="replace")
        do ii=1,nticks
           write(2002,"(10E20.10)") ticks(ii),&
                sum(cumulative_kappa(:,:,:,ii),dim=3)
           write(2003,"(2E20.10)") ticks(ii),&
                sum(sum(cumulative_kappa(:,:,:,ii),dim=3),&
                reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.0_dp
        end do
        close(2002)
        close(2003)
        deallocate(ticks,cumulative_kappa)

        ! Cumulative thermal conductivity vs angular frequency.
        allocate(ticks(nticks),cumulative_kappa(3,3,nbands,nticks))
        call CumulativeTConductOmega(energy,velocity,F_n,ticks,cumulative_kappa)
        !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
        !$OMP & private(ii,ll) shared(nbands,nticks,cumulative_kappa)
        do ii=1,nticks
           do ll=1,nbands
              call symmetrize_tensor(cumulative_kappa(:,:,ll,ii))
           end do
        end do
        !$OMP END PARALLEL DO
        write(aux,"(I0)") 9*nbands+1
        open(2002,file="BTE.cumulative_kappaVsOmega_tensor",status="replace")
        do ii=1,nticks
           write(2002,"(10E20.10)") ticks(ii),&
                sum(cumulative_kappa(:,:,:,ii),dim=3)
        end do
        close(2002)
        deallocate(ticks,cumulative_kappa)
     end if

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     ! If results for nanowires have been requested, repeat the iterative process
     ! for each desired orientation, introducing the appropriate scaling of
     ! relaxation times.
     if(nanowires) then
        kappa_wires=0.0_dp
        kk=ceiling(real(nwires,kind=dp)/real(numprocs,kind=dp))
        do iorient=1,norientations
           if(myid.eq.0) then
              write(*,"(A,I0,A,3(1x,I0))") "Info: nanowires with orientation ",&
                   iorient,":",orientations(:,iorient)
           end if
           write(sorientation,"(I128)") iorient
           !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
           !$OMP & private(ii,jj) shared(nbands,nptk,velocity,v_or,F_or,F_n_0,uorientations,iorient)
           do ii=1,nptk
              do jj=1,Nbands
                 v_or(jj,ii)=dot_product(velocity(:,jj,ii),uorientations(:,iorient))
                 F_or(jj,ii)=dot_product(F_n_0(:,jj,ii),uorientations(:,iorient))
              end do
           end do
           !$OMP END PARALLEL DO
           do mm=1,nwires
              radnw=radnw_range(mm)
              call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
              !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
              !$OMP & private(ii,jj) shared(nbands,nptk,F_n_aux,F_or,ffunc)
              do ii=1,nptk
                 do jj=1,Nbands
                    F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(jj,ii)
                 end do
              end do
              !$OMP END PARALLEL DO
              call TConductScalar(energy,v_or,F_n_aux,kappa_or)
              if(convergence) then
                 do ii=1,maxiter
                    call iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                         & Ntotal_plus,Ntotal_minus,energy,v_or,Gamma_plus,Gamma_minus,&
                         tau_zero,F_n_aux)
                    if (myid .eq. 0) then
                        kappa_or_old=sum(kappa_or)
                        !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
                        !$OMP & private(ll,jj) shared(nbands,nptk,F_n_aux,ffunc)
                        do ll=1,nptk
                           do jj=1,nbands
                              F_n_aux(jj,ll)=F_n_aux(jj,ll)*ffunc(jj,ll)
                           end do
                        end do
                        !$OMP END PARALLEL DO
                        call TConductScalar(energy,v_or,F_n_aux,kappa_or)
                        relchange=abs((sum(kappa_or)-kappa_or_old)/kappa_or_old)
                    endif
                    call MPI_BCAST(relchange,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
                    if(relchange.lt.eps)exit
                 end do
              end if
              kappa_wires(:,mm)=kappa_or
           end do
           if(myid.eq.0) then
              write(aux,"(I0)") 3*nbands
              open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation)),status="replace")
              do ii=1,Nwires
                 radnw=radnw_range(ii)
                 write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.0_dp*radnw,&
                      kappa_wires(:,ii),sum(kappa_wires(:,ii))
              end do
              close(3001)
           end if
        end do
     end if
     if (myid.eq.0) call change_directory(".."//C_NULL_CHAR)
  end do ! Tcounter


  if(myid.eq.0)write(*,*) "Info: normal exit"

  call free_config()
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
end program ShengBTE
