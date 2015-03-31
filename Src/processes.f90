!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2015 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2015 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2015 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2015 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
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

! Compute the number of allowed three-phonon processes, their
! scattering amplitudes and their phase-space volume.

module processes
  use iso_fortran_env
  use misc
  use data
  use config
  implicit none

  real(kind=8),parameter :: hbarp=hbar*1e22

contains

  ! Number of absorption and emission processes.
  subroutine Nprocesses(mm,N_plus,N_minus,energy,velocity,Nlist,List,IJK)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    integer(kind=4),intent(out) :: N_plus,N_minus
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus=0
    N_minus=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,List(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    ! omega, omegap and omegadp are the frequencies of the three
    ! phonons involved in each process. A locally adaptive broadening
    ! (sigma) is computed for each process and used to determine
    ! whether it is allowed by conservation of energy. The regularized
    ! version of the Dirac delta chosen for this is a gaussian.
    if (omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus=N_plus+1
                   endif
                end if
             end do ! k
             !--------END absorption process-------------
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      N_minus=N_minus+1
                   endif
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
    end if
  end subroutine Nprocesses

  ! Wrapper around Nprocesses that splits the work among processors.
  subroutine Nprocesses_driver(energy,velocity,Nlist,List,IJK,&
       N_plus,N_minus)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    integer(kind=4),intent(in) :: Nlist
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(out) :: N_plus(Nlist*Nbands)
    integer(kind=4),intent(out) :: N_minus(Nlist*Nbands)

    integer(kind=4) :: mm
    integer(kind=4) :: N_plus_reduce(Nlist*Nbands)
    integer(kind=4) :: N_minus_reduce(Nlist*Nbands)

    N_plus=0
    N_minus=0
    N_plus_reduce=0
    N_minus_reduce=0

    do mm=myid+1,Nbands*Nlist,numprocs
       call Nprocesses(mm,N_plus_reduce(mm),N_minus_reduce(mm),&
            energy,velocity,Nlist,List(1:Nlist),IJK)
    end do

    call MPI_ALLREDUCE(N_plus_reduce,N_plus,Nbands*Nlist,&
         MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(N_minus_reduce,N_minus,Nbands*Nlist,&
         MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mm)
  end subroutine Nprocesses_driver

  ! Scattering amplitudes of absorption processes.
  subroutine Ind_plus(mm,N_plus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_plus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_plus(N_plus),Indof3rdPhonon_plus(N_plus)
    real(kind=8),intent(out) :: Gamma_plus(N_plus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_plus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,rr,ss,tt
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    complex(kind=8) :: Vp,Vp0,prefactor

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    ! The following loop is very similar to the central part of
    ! Nprocesses(), with the exception that when an allowed process is
    ! detected its scattering amplitude is computed. This amplitude
    ! involves the anharmonic force constants, and its precise
    ! expression can be found in the manuscript.
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(N_plus_count)=(index_N(qprime(1),qprime(2),qprime(3))-1)*Nbands+j
                      Indof3rdPhonon_plus(N_plus_count)=(index_N(qdprime(1),qdprime(2),qdprime(3))-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      !--------BEGIN calculation of Vp-----------
                      Vp=0.
                      do ll=1,Ntri
                         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
                              masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
                              phexp(dot_product(realqprime/ngrid,R_j(:,ll)))*&
                              phexp(-dot_product(realqdprime/ngrid,R_k(:,ll)))
                         Vp0=0.
                         do rr=1,3
                            do ss=1,3
                               do tt=1,3
                                  Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                                       eigenvect(index_N(q(1),q(2),q(3)),i,tt+3*(Index_i(ll)-1))*&
                                       eigenvect(index_N(qprime(1),qprime(2),qprime(3)),j,ss+3*(Index_j(ll)-1))*&
                                       conjg(eigenvect(index_N(qdprime(1),qdprime(2),qdprime(3)),k,rr+3*(Index_k(ll)-1)))
                               end do
                            end do
                         end do
                         Vp=Vp+prefactor*Vp0
                      end do
                      !--------END calculation of Vp-------------
                      Gamma_plus(N_plus_count)=hbarp*pi/4.d0*(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)*abs(Vp)**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(N_plus_count)=Gamma_plus(N_plus_count)*5.60626442*1.d8/nptk ! THz
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
    if(N_plus_count.ne.N_plus) write(error_unit,*) "Error: in Ind_plus, N_plus_count!=N_plus"
  end subroutine Ind_plus

  ! Scattering amplitudes of emission processes. See Ind_plus() for details.
  subroutine Ind_minus(mm,N_minus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_minus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_minus(N_minus),Indof3rdPhonon_minus(N_minus)
    real(kind=8),intent(out) :: Gamma_minus(N_minus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,rr,ss,tt
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) ::  omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    complex(kind=8) :: Vp,Vp0,prefactor

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(N_minus_count)=(index_N(qprime(1),qprime(2),qprime(3))-1)*Nbands+j
                      Indof3rdPhonon_minus(N_minus_count)=(index_N(qdprime(1),qdprime(2),qdprime(3))-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      !--------BEGIN calculation of Vp-----------
                      Vp=0.
                      do ll=1,Ntri
                         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
                              masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
                              phexp(-dot_product(realqprime/ngrid,R_j(:,ll)))*&
                              phexp(-dot_product(realqdprime/ngrid,R_k(:,ll)))
                         Vp0=0.
                         do rr=1,3
                            do ss=1,3
                               do tt=1,3
                                  Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                                       eigenvect(index_N(q(1),q(2),q(3)),i,tt+3*(Index_i(ll)-1))*&
                                       conjg(eigenvect(index_N(qprime(1),qprime(2),qprime(3)),j,ss+3*(Index_j(ll)-1)))*&

                                       conjg(eigenvect(index_N(qdprime(1),qdprime(2),qdprime(3)),k,rr+3*(Index_k(ll)-1)))
                               end do
                            end do
                         end do
                         Vp=Vp+prefactor*Vp0
                      end do
                      !--------END calculation of Vp-------------
                      Gamma_minus(N_minus_count)=hbarp*pi/4.d0*(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)*abs(Vp)**2
                      Gamma_minus(N_minus_count)=Gamma_minus(N_minus_count)*5.60626442*1.d8/nptk
                   end if
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
    end if
    if(N_minus_count.ne.N_minus) write(error_unit,*) "Error: in Ind_minus, N_minus_count!=N_minus"
  end subroutine Ind_minus

  ! Wrapper around Ind_plus and Ind_minus that splits the work among processors.
  subroutine Ind_driver(energy,velocity,eigenvect,Nlist,List,IJK,N_plus,N_minus,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus,&
       Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus,rate_scatt)
       
    implicit none

    include "mpif.h"
    
    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands)
    integer(kind=4),intent(in) :: N_minus(Nlist*Nbands)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)
    integer(kind=4),intent(out) :: Indof2ndPhonon_plus(:)
    integer(kind=4),intent(out) :: Indof3rdPhonon_plus(:)
    real(kind=8),intent(out) :: Gamma_plus(:)
    integer(kind=4),intent(out) :: Indof2ndPhonon_minus(:)
    integer(kind=4),intent(out) :: Indof3rdPhonon_minus(:)
    real(kind=8),intent(out) :: Gamma_minus(:)
    real(kind=8),intent(out) :: rate_scatt(Nbands,Nlist)

    integer(kind=4) :: i
    integer(kind=4) :: ll
    integer(kind=4) :: mm
    integer(kind=4) :: maxsize
    integer(kind=4) :: Ntotal_plus
    integer(kind=4) :: Ntotal_minus
    integer(kind=4) :: Naccum_plus(Nbands*Nlist)
    integer(kind=4) :: Naccum_minus(Nbands*Nlist)
    integer(kind=4),allocatable :: Indof2ndPhonon(:)
    integer(kind=4),allocatable :: Indof3rdPhonon(:)
    integer(kind=4),allocatable :: Indof2ndPhonon_plus_reduce(:)
    integer(kind=4),allocatable :: Indof3rdPhonon_plus_reduce(:)
    integer(kind=4),allocatable :: Indof2ndPhonon_minus_reduce(:)
    integer(kind=4),allocatable :: Indof3rdPhonon_minus_reduce(:)
    real(kind=8) :: rate_scatt_reduce(Nbands,Nlist)
    real(kind=8),allocatable :: Gamma0(:)
    real(kind=8),allocatable :: Gamma_plus_reduce(:)
    real(kind=8),allocatable :: Gamma_minus_reduce(:)

    maxsize=max(maxval(N_plus),maxval(N_minus))
    allocate(Indof2ndPhonon(maxsize))
    allocate(Indof3rdPhonon(maxsize))
    allocate(Gamma0(maxsize))
    
    Naccum_plus(1)=0
    Naccum_minus(1)=0
    do mm=2,Nbands*Nlist
       Naccum_plus(mm)=Naccum_plus(mm-1)+N_plus(mm-1)
       Naccum_minus(mm)=Naccum_minus(mm-1)+N_minus(mm-1)
    end do
    Ntotal_plus=sum(N_plus)
    Ntotal_minus=sum(N_minus)

    allocate(Indof2ndPhonon_plus_reduce(Ntotal_plus))
    allocate(Indof3rdPhonon_plus_reduce(Ntotal_plus))
    allocate(Indof2ndPhonon_minus_reduce(Ntotal_minus))
    allocate(Indof3rdPhonon_minus_reduce(Ntotal_minus))
    allocate(Gamma_plus_reduce(Ntotal_plus))
    allocate(Gamma_minus_reduce(Ntotal_minus))

    Indof2ndPhonon_plus_reduce=0
    Indof3rdPhonon_plus_reduce=0
    Indof2ndPhonon_minus_reduce=0
    Indof3rdPhonon_minus_reduce=0
    Gamma_plus_reduce=0.d0
    Gamma_minus_reduce=0.d0
    rate_scatt_reduce=0.d0    
    
    do mm=myid+1,Nbands*NList,numprocs
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       if(N_plus(mm).ne.0) then
          call Ind_plus(mm,N_plus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_plus(mm)),Indof3rdPhonon(1:N_plus(mm)),&
               Gamma0(1:N_plus(mm)))
          Indof2ndPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof2ndPhonon(1:N_plus(mm))
          Indof3rdPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof3rdPhonon(1:N_plus(mm))
          Gamma_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Gamma0(1:N_plus(mm))
          rate_scatt_reduce(i,ll)=rate_scatt_reduce(i,ll)+sum(Gamma0(1:N_plus(mm)))
       end if
       if(N_minus(mm).ne.0) then
          call Ind_minus(mm,N_minus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_minus(mm)),Indof3rdPhonon(1:N_minus(mm)),&
               Gamma0(1:N_minus(mm)))
          Indof2ndPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof2ndPhonon(1:N_minus(mm))
          Indof3rdPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof3rdPhonon(1:N_minus(mm))
          Gamma_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Gamma0(1:N_minus(mm))
          rate_scatt_reduce(i,ll)=rate_scatt_reduce(i,ll)+sum(Gamma0(1:N_minus(mm)))*5.D-1
       end if
    end do

    deallocate(Gamma0)
    deallocate(Indof3rdPhonon)
    deallocate(Indof2ndPhonon)

    call MPI_ALLREDUCE(Indof2ndPhonon_plus_reduce,Indof2ndPhonon_plus,&
         Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof3rdPhonon_plus_reduce,Indof3rdPhonon_plus,&
         Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof2ndPhonon_minus_reduce,Indof2ndPhonon_minus,&
         Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof3rdPhonon_minus_reduce,Indof3rdPhonon_minus,&
         Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Gamma_plus_reduce,Gamma_plus,Ntotal_plus,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Gamma_minus_reduce,Gamma_minus,Ntotal_minus,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(rate_scatt_reduce,rate_scatt,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ll)

    deallocate(Indof2ndPhonon_plus_reduce)
    deallocate(Indof3rdPhonon_plus_reduce)
    deallocate(Indof2ndPhonon_minus_reduce)
    deallocate(Indof3rdPhonon_minus_reduce)
    deallocate(Gamma_plus_reduce)
    deallocate(Gamma_minus_reduce)
  end subroutine Ind_driver

  ! Mode contribution to the phase-space volume. The algorithm is very
  ! similar to the rest of subroutines in this module. No matrix
  ! element is calculated in this case (compare with Ind_*) when an
  ! allowed process is detected. Instead, the amplitude of each
  ! gaussian is added to the final result. This subroutine computes the
  ! absorption part.
  subroutine D_plus(mm,N_plus,energy,velocity,Nlist,List,IJK,P3_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_plus
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(out) :: P3_plus(N_plus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_plus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus_count=N_plus_count+1
                      P3_plus(N_plus_count)=exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
    if(N_plus_count.ne.N_plus) write(error_unit,*) "Error: in D_plus, N_plus_count!=N_plus"
  end subroutine D_plus

  ! Same as D_plus, but for emission processes.
  subroutine D_minus(mm,N_minus,energy,velocity,Nlist,List,IJK,P3_minus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_minus
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(out) :: P3_minus(N_minus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if(abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      N_minus_count=N_minus_count+1
                      P3_minus(N_minus_count)=exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                end if
             end do ! k
             !--------END emission process-------------!
          end do ! ii
       end do  ! j
    end if
    if(N_minus_count.ne.N_minus) write(error_unit,*) "Error: in D_minus, N_minus_count!=N_minus"
  end subroutine D_minus

  ! Wrapper around D_plus and D_minus that splits the work among processors.
  subroutine P3_driver(energy,velocity,Nlist,List,IJK,N_plus,N_minus,&
       Pspace_plus_total,Pspace_plus_partial,&
       Pspace_minus_total,Pspace_minus_partial)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands)
    integer(kind=4),intent(in) :: N_minus(Nlist*Nbands)
    real(kind=8),intent(out) :: Pspace_plus_total(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_plus_partial(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_minus_total(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_minus_partial(Nbands,Nlist)

    integer(kind=4) :: ii
    integer(kind=4) :: jj
    integer(kind=4) :: mm
    real(kind=8),allocatable :: Pspace_plus_tmp(:)
    real(kind=8),allocatable :: Pspace_minus_tmp(:)

    Pspace_plus_total=0.d0
    Pspace_plus_partial=0.d0
    Pspace_minus_total=0.d0
    Pspace_minus_partial=0.d0

    allocate(Pspace_plus_tmp(maxval(N_plus)))
    allocate(Pspace_minus_tmp(maxval(N_minus)))

    do mm=myid+1,Nbands*Nlist,numprocs
       Pspace_plus_tmp=0.d0
       if(N_plus(mm).ne.0) then
          ii=modulo(mm-1,Nbands)+1
          jj=int((mm-1)/Nbands)+1
          call D_plus(mm,N_plus(mm),energy,velocity,Nlist,List,IJK,&
               Pspace_plus_tmp)
          Pspace_plus_partial(ii,jj)=Pspace_plus_partial(ii,jj)+sum(Pspace_plus_tmp)
       end if
       Pspace_minus_tmp=0.d0
       if(N_minus(mm).ne.0) then
          ii=modulo(mm-1,Nbands)+1
          jj=int((mm-1)/Nbands)+1
          call D_minus(mm,N_minus(mm),energy,velocity,Nlist,List,IJK,&
               Pspace_minus_tmp)
          Pspace_minus_partial(ii,jj)=Pspace_minus_partial(ii,jj)+sum(Pspace_minus_tmp)
       end if
    end do

    call MPI_ALLREDUCE(Pspace_plus_partial,Pspace_plus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ii)
    call MPI_ALLREDUCE(Pspace_minus_partial,Pspace_minus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ii)

    deallocate(Pspace_minus_tmp)
    deallocate(Pspace_plus_tmp)
  end subroutine P3_driver
end module processes
