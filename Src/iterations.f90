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

! Code used to initialize F_n for the zeroth-order iteration of the
! BTE and to perform successive iterations.
module iterations
  use data
  use config
  use wedgetc
  use mpi_f08
  implicit none

contains

  ! Fill F_n with its initial values for the iteration, computed from the
  ! RTA approximation to tau.
  subroutine iteration0(Nlist,Nequi,ALLEquiList,omega,velocity,tau_zero,F_n)
    implicit none

    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),tau_zero(Nbands,Nlist)
    real(kind=dp),intent(out) :: F_n(3,Nbands,nptk)

    integer :: ii,kk,ll

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity) &
    !$OMP & shared(F_n,Nlist,nbands,Nequi,ALLEquiList) private(ll,ii,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do ii=1,Nbands
             F_n(:,ii,ALLEquiList(kk,ll))=tau_zero(ii,ll)*velocity(:,ii,ALLEquiList(kk,ll))*&
                  omega(ii,ALLEquiList(kk,ll))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine iteration0

  ! Advance the algorithm one iteration. F_n is used both as the input
  ! and as the output.
  subroutine iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                      & omega,velocity,tau_zero,F_n)
    implicit none
    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer,intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer,intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands)
    real(kind=dp),intent(in) :: omega(nbands,nptk),velocity(3,nbands,nptk)
    real(kind=dp),intent(in) :: tau_zero(nbands,nlist)
    real(kind=dp),intent(inout) :: F_n(3,Nbands,nptk)

    integer :: ID_equi(Nsymm_rot,nptk),Naccum_plus,Naccum_minus
    integer :: i,j,k,jj,kk,ll,mm,nn,mmm,nnn
    real(kind=dp) :: DeltaF(3,Nbands,nptk)

    call symmetry_map(ID_equi)
    DeltaF=0.0_dp
    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nstates,myid,nbands,nlist) &
    !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
    !$OMP & shared(N_plus,Naccum_plus_array,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus) &
    !$OMP & shared(N_minus,Naccum_minus_array,Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus) &
    !$OMP & private(i,j,k,kk,ll,mm,nn,mmm,nnn,Naccum_plus,Naccum_minus)
    do nnn=1,nstates
        mmm=myid*nstates+nnn
        i=modulo(mmm-1,Nbands)+1
        ll=int((mmm-1)/Nbands)+1

        if (mmm.gt.nlist*nbands) cycle

        Naccum_plus=Naccum_plus_array(nnn)
        Naccum_minus=Naccum_minus_array(nnn)

        do kk=1,Nequi(ll)
           if ((N_plus(mmm).ne.0)) then
              do jj=1,N_plus(mmm) 
                 j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                 mm=(Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                 k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                 nn=(Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                 DeltaF(:,i,ALLEquiList(kk,ll))=DeltaF(:,i,ALLEquiList(kk,ll))+&
                      Gamma_plus(Naccum_plus+jj)*(F_n(:,k,ID_equi(TypeofSymmetry(kk,ll),nn))-&
                      F_n(:,j,ID_equi(TypeofSymmetry(kk,ll),mm)))
              end do !jj
           end if
           if ((N_minus(mmm).ne.0)) then
              do jj=1,N_minus(mmm)
                 j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                 mm=(Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                 k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                 nn=(Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                 DeltaF(:,i,ALLEquiList(kk,ll))=DeltaF(:,i,ALLEquiList(kk,ll))+&
                      Gamma_minus(Naccum_minus+jj)*(F_n(:,k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                      F_n(:,j,ID_equi(TypeofSymmetry(kk,ll),mm)))*0.5_dp
              end do !jj
           end if
        end do !kk
    end do !nnn
    !$OMP END PARALLEL DO

    call MPI_BARRIER(MPI_COMM_WORLD,mmm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,DeltaF,nptk*nbands*3,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,mm)

         

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity,ALLEquiList) &
    !$OMP & shared(F_n,DeltaF,Nlist,nbands,Nequi) private(ll,i,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do i=1,Nbands
             F_n(:,i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(:,i,ALLEquiList(kk,ll))*&
                   omega(i,ALLEquiList(kk,ll))+tau_zero(i,ll)*DeltaF(:,i,ALLEquiList(kk,ll))
          enddo
       enddo
    enddo  
    !$OMP END PARALLEL DO       

  end subroutine iteration

  ! Restricted variation of the above subroutine, limited to cases where kappa
  ! is an scalar. Used in the nanowire calculation.
  subroutine iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
      & Ntotal_plus,Ntotal_minus,&
      & omega,velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer,intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer,intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    real(kind=dp),intent(in) :: omega(nbands,nptk),velocity(nbands,nptk)
    real(kind=dp),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=dp),intent(inout) :: F_n(Nbands,nptk)

    integer :: ID_equi(Nsymm_Rot,nptk),Naccum_plus,Naccum_minus
    integer :: i,j,k,jj,kk,ll,mm,nn,mmm,nnn
    real(kind=dp) :: DeltaF(Nbands,nptk)

    DeltaF=0.0_dp
    call symmetry_map(ID_equi)
    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nstates,myid,nbands,nlist) &
    !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
    !$OMP & shared(N_plus,Naccum_plus_array,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus) &
    !$OMP & shared(N_minus,Naccum_minus_array,Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus) &
    !$OMP & private(i,j,k,kk,ll,mm,nn,mmm,nnn,Naccum_plus,Naccum_minus)
    do nnn=1,nstates
        mmm=myid*nstates+nnn
        i=modulo(mmm-1,Nbands)+1
        ll=int((mmm-1)/Nbands)+1

        if (mmm.gt.nlist*nbands) cycle
        Naccum_plus=Naccum_plus_array(nnn)
        Naccum_minus=Naccum_minus_array(nnn)

        do kk=1,Nequi(ll)
           if ((N_plus(mmm).ne.0)) then
              do jj=1,N_plus(mmm) 
                  j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  mm=(Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                  k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  nn=(Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                  DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))
              end do !jj
           end if
           if ((N_minus(mmm).ne.0)) then
               do jj=1,N_minus(mmm)
                  j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  mm=(Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                  k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  nn=(Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                  DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))*0.5_dp
               end do !jj
           end if
        end do !kk
    end do ! nnn
    !$OMP END PARALLEL DO

    call MPI_BARRIER(MPI_COMM_WORLD,mmm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,DeltaF,nptk*nbands,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,mm)

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity,ALLEquiList) &
    !$OMP & shared(F_n,DeltaF,Nlist,nbands,Nequi) private(ll,i,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do i=1,Nbands
             F_n(i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(i,ALLEquiList(kk,ll))*&
                  omega(i,ALLEquiList(kk,ll))+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll))
          enddo
       enddo
    enddo   
    !$OMP END PARALLEL DO 
  end subroutine iteration_scalar
end module iterations
