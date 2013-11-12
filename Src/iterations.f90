!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2013 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2013 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2013 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2013 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
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
module iterations
  use data
  use config
  use wedgetc
  implicit none

contains

  subroutine iteration0(Nlist,Nequi,ALLEquiList,omega,velocity,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm,nptk)
    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),tau_zero(Nbands,Nlist)
    real(kind=8),intent(out) :: F_n(Nbands,nptk,3)

    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll

    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    do ll=1,Nlist
       do ii=1,Nbands
          do kk=1,Nequi(ll)
             F_n(ii,ALLEquiList(kk,ll),:)=tau_zero(ii,ll)*velocity(ALLEquiList(kk,ll),ii,:)*&
                  omega(ALLEquiList(kk,ll),ii)
          end do
       end do
    end do
  end subroutine iteration0

  subroutine iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,Ntotal_plus,Ntotal_minus,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Indof2ndPhonon_minus,Indof3rdPhonon_minus,omega,&
       velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm,nptk)
    integer(kind=4),intent(in) :: TypeofSymmetry(Nsymm,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    integer(kind=4),intent(in) :: Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus)
    integer(kind=4),intent(in) :: Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus)
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands,3)
    real(kind=8),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=8),intent(inout) :: F_n(Nbands,nptk,3)

    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ID_equi(Nsymm,nptk),Naccum_plus,Naccum_minus
    integer(kind=4) :: i,j,k,ii,jj,kk,ll,mm,nn
    real(kind=8) :: DeltaF(Nbands,nptk,3)

    DeltaF=0.d0
    call symmetry_map(ID_equi)
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    DeltaF=0.d0
    do ll=1,Nlist
       do i=1,Nbands
          if (((ll-1)*Nbands+i).eq.1) then
             Naccum_plus=0
             Naccum_minus=0
          else
             Naccum_plus=Naccum_plus+N_plus((ll-1)*Nbands+i-1)
             Naccum_minus=Naccum_minus+N_minus((ll-1)*Nbands+i-1)
          end if
          do kk=1,Nequi(ll)
             if ((N_plus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_plus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))
                end do !jj
             end if
             if ((N_minus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_minus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))*5.D-1
                end do !jj
             end if
             F_n(i,ALLEquiList(kk,ll),:)=tau_zero(i,ll)*velocity(ALLEquiList(kk,ll),i,:)*&
                  omega(ALLEquiList(kk,ll),i)+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll),:)
          end do !kk
       end do
    end do
  end subroutine iteration

  subroutine iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,Ntotal_plus,&
       Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Indof2ndPhonon_minus,&
       Indof3rdPhonon_minus,omega,velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm,nptk)
    integer(kind=4),intent(in) :: TypeofSymmetry(Nsymm,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    integer(kind=4),intent(in) :: Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus)
    integer(kind=4),intent(in) :: Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus)
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands)
    real(kind=8),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=8),intent(inout) :: F_n(Nbands,nptk)

    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ID_equi(Nsymm,nptk),Naccum_plus,Naccum_minus
    integer(kind=4) :: i,j,k,ii,jj,kk,ll,mm,nn
    real(kind=8) :: DeltaF(Nbands,nptk)

    DeltaF=0.d0
    call symmetry_map(ID_equi)
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    DeltaF=0.d0
    do ll=1,Nlist
       do i=1,Nbands
          if (((ll-1)*Nbands+i).eq.1) then
             Naccum_plus=0
             Naccum_minus=0
          else
             Naccum_plus=Naccum_plus+N_plus((ll-1)*Nbands+i-1)
             Naccum_minus=Naccum_minus+N_minus((ll-1)*Nbands+i-1)
          end if
          do kk=1,Nequi(ll)
             if ((N_plus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_plus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))
                end do !jj
             end if
             if ((N_minus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_minus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))*5.D-1
                end do !jj
             end if
             F_n(i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(ALLEquiList(kk,ll),i)*&
                  omega(ALLEquiList(kk,ll),i)+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll))
          end do !kk
       end do
    end do
  end subroutine iteration_scalar
end module iterations
