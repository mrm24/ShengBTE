!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2022 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2022 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2022 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2022 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
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

! Compute the number of allowed three-phonon processes, their
! scattering amplitudes and their phase-space volume.

module processes
  use iso_fortran_env
  use misc
  use data
  use config
  use mpi_f08
  implicit none

  real(kind=dp),parameter :: hbarp=hbar*1.0e22_dp

contains
  subroutine calculate_Vp(energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK, plus_size,minus_size)

    implicit none

    integer,intent(in) :: NList,List(Nlist),IJK(3,nptk),Ntri,plus_size,minus_size
    integer,intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk)
    real(kind=dp),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=dp),intent(in) :: eigenvect(Nbands,Nbands,nptk)

    integer :: q(3),qprime(3),qdprime_p(3),qdprime_m(3),i,j,k
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))

    integer :: ii,jj,kk,ll,ss_p,ss_m,mm,nn
    integer :: N_plus_count, N_minus_count
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: omega,omegap,omegadp
    real(kind=dp) :: realqprime(3),realqdprime_p(3),realqdprime_m(3)

    complex(kind=dp), allocatable :: eigenvect_i(:),eigenvect_j(:), eigenvect_k(:)

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do

    allocate(Vp_plus_matrix(plus_size))
    allocate(Vp_minus_matrix(minus_size))
    
    Vp_plus_matrix  = 0.0_dp
    Vp_minus_matrix = 0.0_dp
    
    call MPI_BARRIER(MPI_COMM_WORLD,mm)
    
    !$OMP parallel default(none) shared(myid,numprocs,ngrid,nlist,Nbands,nstates,nptk,IJK,rlattvec) &
    !$OMP & shared(energy,velocity,scalebroad,Vp_minus_matrix,Vp_plus_matrix,Naccum_plus_array,Naccum_minus_array) &
    !$OMP & shared(Index_N,Index_i,Index_j,Index_k,Phi,R_k,R_j,Ntri,eigenvect,omega_max,list) &
    !$OMP & private(i,j,k,ii,jj,kk,ll,q,qprime,qdprime_p,qdprime_m,realqprime,realqdprime_p,realqdprime_m) &
    !$OMP & private(omega,omegap,omegadp,vel_j) &
    !$OMP & private(ss_p,ss_m,sigma,N_plus_count,N_minus_count,nn,mm,eigenvect_i,eigenvect_j,eigenvect_k)

    allocate(eigenvect_i(nbands),eigenvect_j(nbands),eigenvect_k(nbands))

    !$OMP DO schedule(dynamic,1)
    do nn=1,nstates
       N_plus_count  = Naccum_plus_array(nn)
       N_minus_count = Naccum_minus_array(nn)
       mm=myid*nstates+nn
       if (mm.gt.nlist*nbands) cycle
         i=modulo(mm-1,Nbands)+1
         ll=int((mm-1)/Nbands)+1
         q=IJK(:,list(ll))
         omega=energy(i,list(ll))
         eigenvect_i(:) = eigenvect(:,i,list(ll))
         ! Loop over all processes, detecting those that are allowed and
         ! computing their amplitudes.
         if(omega .ne. 0.0_dp .and. omega .le. omega_max) then
            do ii=1,nptk
               qprime=IJK(:,ii)
               realqprime=matmul(rlattvec,real(qprime,kind=dp)/real(ngrid,kind=dp))
               qdprime_p=q+qprime
               qdprime_p=modulo(qdprime_p,Ngrid)
               realqdprime_p=matmul(rlattvec,real(qdprime_p,kind=dp)/real(ngrid,kind=dp))
               ss_p=Index_N(qdprime_p(1),qdprime_p(2),qdprime_p(3))

               qdprime_m=q-qprime
               qdprime_m=modulo(qdprime_m,Ngrid)
               realqdprime_m=matmul(rlattvec,real(qdprime_m,kind=dp)/real(ngrid,kind=dp))
               ss_m=Index_N(qdprime_m(1),qdprime_m(2),qdprime_m(3))
                     
               do j=1,Nbands
                  omegap=energy(j,ii)
                  if (omegap .eq. 0.0_dp) cycle
                  eigenvect_j = eigenvect(:,j,ii)
                  do k=1,Nbands
                     !--------BEGIN absorption process-----------!
                     omegadp=energy(k,ss_p)
                     if (omegadp .ne. 0.0_dp) then
                        sigma=scalebroad*base_sigma(velocity(:,j,ii)-velocity(:,k,ss_p))
                        if(abs(omega+omegap-omegadp).le.(2.0_dp*sigma)) then
                           eigenvect_k =  eigenvect(:,k,ss_p)
                           N_plus_count = N_plus_count + 1
                           Vp_plus_matrix(N_plus_count)=Vp_plus(eigenvect_i,eigenvect_j,eigenvect_k,&
                                 realqprime,realqdprime_p,&
                                 Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                        end if
                     end if
                     !--------END absorption process-------------!
                     !--------BEGIN emission process-------------!
                     omegadp=energy(k,ss_m)
                     if (omegadp .ne. 0.0_dp) then
                        sigma=scalebroad*base_sigma(velocity(:,j,ii)-velocity(:,k,ss_m))
                        if (abs(omega-omegap-omegadp).le.(2.0_dp*sigma)) then
                           eigenvect_k =  eigenvect(:,k,ss_m)
                           N_minus_count = N_minus_count + 1
                           Vp_minus_matrix(N_minus_count)=Vp_minus(eigenvect_i,eigenvect_j,eigenvect_k,&
                                 realqprime,realqdprime_m,&
                                 Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                        end if
                     end if
                     !--------END emission process---------------!
                  end do ! k
               end do ! ii
            end do  ! j
         end if
    end do ! nstates
    !$OMP END DO

    deallocate(eigenvect_i,eigenvect_j,eigenvect_k)
    !$OMP END PARALLEL
    return
  end subroutine calculate_Vp

  ! Compute one of the matrix elements involved in the calculation of Ind_plus.
  function Vp_plus(eigenvect_i,eigenvect_j,eigenvect_k,realqprime,realqdprime,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none

    complex(kind=dp),intent(in) :: eigenvect_i(Nbands),eigenvect_j(Nbands),eigenvect_k(Nbands)
    real(kind=dp),intent(in) :: realqprime(3)
    real(kind=dp),intent(in) :: realqdprime(3)
    integer,intent(in) :: Ntri
    real(kind=dp),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=dp),intent(in) :: R_j(3,Ntri)
    real(kind=dp),intent(in) :: R_k(3,Ntri)
    integer,intent(in) :: Index_i(Ntri)
    integer,intent(in) :: Index_j(Ntri)
    integer,intent(in) :: Index_k(Ntri)

    real(kind=dp) :: Vp_plus

    integer :: ll
    integer :: rr
    integer :: ss
    integer :: tt
    complex(kind=dp) :: prefactor
    complex(kind=dp) :: Vp0,Vp_plus_c

    Vp_plus_c= cmplx(0.0_dp,0.0_dp,kind=dp)
    
    do ll=1,Ntri
       prefactor=1.0_dp/sqrt(masses(types(Index_i(ll)))*&
            masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
            phexp(dot_product(realqprime,R_j(:,ll)))*&
            phexp(-dot_product(realqdprime,R_k(:,ll)))
       Vp0=cmplx(0.0_dp,0.0_dp,kind=dp)
       do rr=1,3
          do ss=1,3
             do tt=1,3
                Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                     eigenvect_i(3*(Index_i(ll)-1)+tt)*&
                     eigenvect_j(3*(Index_j(ll)-1)+ss)*&
                     conjg(eigenvect_k(3*(Index_k(ll)-1)+rr))
             end do
          end do
       end do
       Vp_plus_c=Vp_plus_c+prefactor*Vp0
    end do
    Vp_plus = abs(Vp_plus_c)
  end function Vp_plus

  ! Compute one of the matrix elements involved in the calculation of Ind_minus.
  function Vp_minus(eigenvect_i,eigenvect_j,eigenvect_k,realqprime,realqdprime,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none

    complex(kind=dp),intent(in) :: eigenvect_i(Nbands),eigenvect_j(Nbands),eigenvect_k(Nbands)
    real(kind=dp),intent(in) :: realqprime(3)
    real(kind=dp),intent(in) :: realqdprime(3)    
    integer,intent(in) :: Ntri
    real(kind=dp),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=dp),intent(in) :: R_j(3,Ntri)
    real(kind=dp),intent(in) :: R_k(3,Ntri)
    integer,intent(in) :: Index_i(Ntri)
    integer,intent(in) :: Index_j(Ntri)
    integer,intent(in) :: Index_k(Ntri)

    real(kind=dp) :: Vp_minus

    integer :: ll
    integer :: rr
    integer :: ss
    integer :: tt
    complex(kind=dp) :: prefactor
    complex(kind=dp) :: Vp0, Vp_minus_c

    Vp_minus_c=cmplx(0.0_dp,0.0_dp,kind=dp)

    do ll=1,Ntri
       prefactor=1.0_dp/sqrt(masses(types(Index_i(ll)))*&
            masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
            phexp(-dot_product(realqprime,R_j(:,ll)))*&
            phexp(-dot_product(realqdprime,R_k(:,ll)))
       Vp0=cmplx(0.0_dp,0.0_dp,kind=dp)
       do rr=1,3
          do ss=1,3
             do tt=1,3
                Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                     eigenvect_i(3*(Index_i(ll)-1)+tt)*&
                     conjg(eigenvect_j(3*(Index_j(ll)-1)+ss))*&
                     conjg(eigenvect_k(3*(Index_k(ll)-1)+rr))
             end do
          end do
       end do
       Vp_minus_c=Vp_minus_c+prefactor*Vp0
    end do

    Vp_minus = abs(Vp_minus_c)
  end function Vp_minus

 
  ! Wrapper around Ind_plus and Ind_minus that splits the work among processors.
  subroutine Ind_driver(nn, energy,velocity,Nlist,List,IJK,N_plus,N_minus, Naccum_plus, Naccum_minus,&
        rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus)
    implicit none

    real(kind=dp),intent(in) :: energy(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    integer,intent(in) :: NList,nn
    integer,intent(in) :: List(Nlist)
    integer,intent(in) :: IJK(3,nptk)
    integer,intent(in) :: N_plus(Nlist*Nbands),Naccum_plus
    integer,intent(in) :: N_minus(Nlist*Nbands),Naccum_minus
    real(kind=dp),intent(out) :: rate_scatt_plus,rate_scatt_minus
    real(kind=dp),intent(out) :: WP3_plus
    real(kind=dp),intent(out) :: WP3_minus

    integer :: i, mm
    integer :: ll

    integer :: q(3),qprime(3),qdprime_p(3),qdprime_m(3),j,k,N_plus_count,N_minus_count
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ss_p, ss_m
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: fBEprime,fBEdprime
    real(kind=dp) :: omega,omegap,omegadp
    real(kind=dp) :: WP3
    real(kind=dp) :: Vp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do

    N_plus_count=0
    N_minus_count=0
    mm=myid*nstates+nn
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(i,list(ll))
    ! Loop over all processes, detecting those that are allowed and
    ! computing their amplitudes.
    do ii=1,nptk
       qprime=IJK(:,ii)
       qdprime_p=q+qprime
       qdprime_p=modulo(qdprime_p,Ngrid)
       ss_p=Index_N(qdprime_p(1),qdprime_p(2),qdprime_p(3))

       qdprime_m=q-qprime
       qdprime_m=modulo(qdprime_m,Ngrid)
       ss_m=Index_N(qdprime_m(1),qdprime_m(2),qdprime_m(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             if (omegap .eq. 0.0_dp) cycle
             fBEprime=1.0_dp/(exp(hbar*omegap/Kb/T)-1.0_dp)
             vel_j(:) = velocity(:,j,ii)
             do k=1,Nbands
                !--------BEGIN absorption process------------!
                omegadp=energy(k,ss_p)
                if (omegadp .ne. 0.0_dp) then
                   sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss_p))
                   if(abs(omega+omegap-omegadp).le.(2.0_dp*sigma)) then
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(Naccum_plus+N_plus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_plus(Naccum_plus+N_plus_count)=(ss_p-1)*Nbands+k
                      fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                      Vp=Vp_plus_matrix(Naccum_plus+N_plus_count) 
                      WP3=(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_plus=WP3_plus+WP3
                      Gamma_plus(Naccum_plus+N_plus_count)=hbarp*pi/4.0_dp*WP3*Vp**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(Naccum_plus+N_plus_count)=&
                      Gamma_plus(Naccum_plus+N_plus_count)*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp) ! THz
                      rate_scatt_plus=rate_scatt_plus+Gamma_plus(Naccum_plus+N_plus_count)        
                   end if
                end if
                !--------END absorption process-------------!
                !--------BEGIN emission process-------------!
                omegadp=energy(k,ss_m)
                if (omegadp .ne. 0.0_dp) then
                   sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss_m))
                   if (abs(omega-omegap-omegadp).le.(2.0_dp*sigma)) then
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(Naccum_minus+N_minus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_minus(Naccum_minus+N_minus_count)=(ss_m-1)*Nbands+k
                      fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                      Vp=Vp_minus_matrix(Naccum_minus+N_minus_count)
                      WP3=(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_minus=WP3_minus+WP3
                      Gamma_minus(Naccum_minus+N_minus_count)=hbarp*pi/4.0_dp*WP3*Vp**2
                      Gamma_minus(Naccum_minus+N_minus_count)=&
                      Gamma_minus(Naccum_minus+N_minus_count)*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp)
                      rate_scatt_minus = rate_scatt_minus + Gamma_minus(Naccum_minus+N_minus_count)*0.5_dp
                   endif
                endif
                !--------END emission process---------------!
             end do ! k
          end do ! ii
    end do  ! j
    WP3_plus=WP3_plus/real(nptk,kind=dp)
    WP3_minus=WP3_minus*0.5_dp/real(nptk,kind=dp)
  end subroutine Ind_driver

  ! Compute the number of allowed absorption processes and their contribution
  ! to phase space.
  subroutine NP_plus(mm,energy,velocity,Nlist,List,IJK,N_plus,P_plus)
    implicit none

    integer,intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk)
    integer,intent(out) :: N_plus
    real(kind=dp),intent(out) :: P_plus

    integer :: q(3),qprime(3),qdprime(3),i,j,k
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: omega,omegap,omegadp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus=0
    P_plus=0.0_dp
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(i,list(ll))
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          qdprime=q+qprime
          qdprime=modulo(qdprime,Ngrid)
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if ((omegap .ne. 0.0_dp).and.(omegadp .ne. 0.0_dp)) then
                   sigma=scalebroad*base_sigma(vel_j-velocity(:,k,ss))
                   if(abs(omega+omegap-omegadp).le.(2.0_dp*sigma)) then
                      N_plus=N_plus+1
                      P_plus=P_plus+&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
  end subroutine NP_plus

  ! Same as NP_plus, but for emission processes.
  subroutine NP_minus(mm,energy,velocity,Nlist,List,IJK,N_minus,P_minus)
    implicit none

    integer,intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk)
    integer,intent(out) :: N_minus
    real(kind=dp),intent(out) :: P_minus

    integer :: q(3),qprime(3),qdprime(3),i,j,k
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: omega,omegap,omegadp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus=0
    P_minus=0.0_dp
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(i,list(ll))
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          qdprime=q-qprime
          qdprime=modulo(qdprime,Ngrid)
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if ((omegap .ne. 0.0_dp).and.(omegadp .ne. 0.0_dp)) then
                   sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss))
                   if(abs(omega-omegap-omegadp).le.(2.0_dp*sigma)) then
                      N_minus=N_minus+1
                      P_minus=P_minus+&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                end if
             end do ! k
             !--------END emission process-------------!
          end do ! ii
       end do  ! j
    end if
  end subroutine NP_minus

  ! Wrapper around NP_plus and NP_minus that splits the work among processors.
  subroutine NP_driver(energy,velocity,Nlist,List,IJK,&
       N_plus,Pspace_plus_total,N_minus,Pspace_minus_total)
    implicit none

    real(kind=dp),intent(in) :: energy(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    integer,intent(in) :: NList
    integer,intent(in) :: List(Nlist)
    integer,intent(in) :: IJK(3,nptk)
    integer,intent(out) :: N_plus(Nlist*Nbands)
    integer,intent(out) :: N_minus(Nlist*Nbands)
    real(kind=dp),intent(out) :: Pspace_plus_total(Nbands,Nlist)
    real(kind=dp),intent(out) :: Pspace_minus_total(Nbands,Nlist)

    integer :: mm

    Pspace_plus_total=0.0_dp
    Pspace_minus_total=0.0_dp
    N_plus=0
    N_minus=0

    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nbands,nlist,numprocs,myid,omega_max) &
    !$OMP & shared(energy,velocity,List,IJK) &
    !$OMP & shared(N_plus,Pspace_plus_total,N_minus,Pspace_minus_total) &
    !$OMP & private(mm)
    do mm=myid+1,Nbands*Nlist,numprocs
       if (energy(modulo(mm-1,Nbands)+1,List(int((mm-1)/Nbands)+1)).le.omega_max) then
          call NP_plus(mm,energy,velocity,Nlist,List,IJK,&
               N_plus(mm),Pspace_plus_total(modulo(mm-1,Nbands)+1, int((mm-1)/Nbands)+1))
          call NP_minus(mm,energy,velocity,Nlist,List,IJK,&
               N_minus(mm),Pspace_minus_total(modulo(mm-1,Nbands)+1, int((mm-1)/Nbands)+1))
       endif
    end do
    !$OMP END PARALLEL DO

    call MPI_ALLREDUCE(MPI_IN_PLACE,N_plus,Nbands*Nlist,MPI_INTEGER,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,N_minus,Nbands*Nlist,MPI_INTEGER,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,Pspace_plus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,Pspace_minus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
  end subroutine NP_driver

  ! RTA-only version of Ind_plus.
  subroutine RTA_plus(nn,q,omega,energy,velocity,IJK,&
       Gamma_plus,WP3_plus)
    implicit none

    integer,intent(in) :: nn,IJK(3,nptk), q(3)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk), omega
    real(kind=dp),intent(out) :: Gamma_plus,WP3_plus

    integer :: qprime(3),qdprime(3),i,j,k
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss,mm
    integer :: N_plus_count
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: fBEprime,fBEdprime
    real(kind=dp) :: omegap,omegadp
    real(kind=dp) :: realqprime(3),realqdprime(3)
    real(kind=dp) :: WP3
    real(kind=dp) :: Vp

    Gamma_plus=0.0_dp
    WP3_plus=0.0_dp
    N_plus_count = Naccum_plus_array(nn)
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    mm=myid*nstates+nn
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,real(qprime,kind=dp)/real(ngrid,kind=dp))
          qdprime=q+qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,real(qdprime,kind=dp)/real(ngrid,kind=dp))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             if (omegap .eq. 0.0_dp) cycle
             fBEprime=1.0_dp/(exp(hbar*omegap/Kb/T)-1.0_dp)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if (omegadp .eq. 0.0_dp) cycle
                sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss))
                if(abs(omega+omegap-omegadp).le.(2.0_dp*sigma)) then
                   N_plus_count = N_plus_count + 1
                   fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                   WP3=(fBEprime-fBEdprime)*&
                      exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                      (omega*omegap*omegadp)
                   WP3_plus=WP3_plus+WP3
                   if (.not.onlyharmonic) then
                   Vp=Vp_plus_matrix(N_plus_count)
                   Gamma_plus=Gamma_plus+hbarp*pi/4.0_dp*WP3*Vp**2
                   endif
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
       WP3_plus=WP3_plus/real(nptk,kind=dp)
    end if
    Gamma_plus=Gamma_plus*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp) ! THz
  end subroutine RTA_plus

  ! RTA-only version of Ind_minus.
  subroutine RTA_minus(nn,q,omega,energy,velocity,&
       IJK,Gamma_minus,WP3_minus)
    implicit none

    integer,intent(in) :: nn,IJK(3,nptk),q(3)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk), omega
    real(kind=dp),intent(out) :: Gamma_minus,WP3_minus

    integer :: qprime(3),qdprime(3),i,j,k,N_minus_count
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss,mm
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: fBEprime,fBEdprime
    real(kind=dp) :: omegap,omegadp
    real(kind=dp) :: realqprime(3),realqdprime(3)
    real(kind=dp) :: WP3
    real(kind=dp) :: Vp

    Gamma_minus=0.0_dp
    WP3_minus=0.0_dp
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=Naccum_minus_array(nn)
    mm=myid*nstates+nn
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,real(qprime,kind=dp)/real(ngrid,kind=dp))
          qdprime=q-qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,real(qdprime,kind=dp)/real(ngrid,kind=dp))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             if (omegap .eq. 0.0_dp) cycle
             fBEprime=1.0_dp/(exp(hbar*omegap/Kb/T)-1.0_dp)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if (omegadp .eq. 0.0_dp) cycle
                sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss))
                if (abs(omega-omegap-omegadp).le.(2.0_dp*sigma)) then
                   N_minus_count = N_minus_count + 1
                   fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                   WP3=(fBEprime+fBEdprime+1)*&
                      exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                      (omega*omegap*omegadp)
                   WP3_minus=WP3_minus+WP3
                   if (.not.onlyharmonic) then
                   Vp=Vp_minus_matrix(N_minus_count)
                   Gamma_minus=Gamma_minus+hbarp*pi/4.0_dp*WP3*Vp**2
                   endif
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
       WP3_minus=WP3_minus*0.5_dp/real(nptk,kind=dp)
    end if
    Gamma_minus=Gamma_minus*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp)
  end subroutine RTA_minus

  ! Wrapper around RTA_plus and RTA_minus that splits the work among processors.
  subroutine RTA_driver(energy,velocity,Nlist,List,IJK,&
       rate_scatt,rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus)
    implicit none

    real(kind=dp),intent(in) :: energy(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    integer,intent(in) :: NList
    integer,intent(in) :: List(Nlist)
    integer,intent(in) :: IJK(3,nptk)
    real(kind=dp),intent(out) :: rate_scatt(Nbands,Nlist),rate_scatt_plus(Nbands,Nlist),rate_scatt_minus(Nbands,Nlist)
    real(kind=dp),intent(out) :: WP3_plus(Nbands,Nlist)
    real(kind=dp),intent(out) :: WP3_minus(Nbands,Nlist)

    integer :: i
    integer :: ll, q(3)
    integer :: mm, nn
    real(kind=dp) :: Gamma_plus,Gamma_minus

    rate_scatt=0.0_dp
    WP3_plus=0.0_dp
    WP3_minus=0.0_dp

    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(Nbands,NList,myid,nstates) &
    !$OMP & shared(energy,velocity,List,IJK,Naccum_plus_array,Naccum_minus_array,omega_max) &
    !$OMP & shared(WP3_plus,WP3_minus,rate_scatt_plus,rate_scatt_minus) &
    !$OMP & private(nn,mm,i,ll,Gamma_plus,Gamma_minus,q)
    do nn=1,nstates
       mm=myid*nstates+nn
       if (mm.gt.nlist*nbands) cycle
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       q=IJK(:,list(ll))
       if (energy(i,List(ll)).le.omega_max) then
          call RTA_plus(nn,q,energy(i,List(ll)),energy,velocity,IJK,&
               Gamma_plus,WP3_plus(i,ll))
          rate_scatt_plus(i,ll)=Gamma_plus
          call RTA_minus(nn,q,energy(i,List(ll)),energy,velocity,IJK,&
               Gamma_minus,WP3_minus(i,ll))
          rate_scatt_minus(i,ll)=Gamma_minus*0.5_dp
       endif
    end do
    !$OMP END PARALLEL DO

    call MPI_ALLREDUCE(MPI_IN_PLACE,rate_scatt_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,rate_scatt_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,WP3_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,WP3_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    rate_scatt=rate_scatt_plus+rate_scatt_minus
  end subroutine RTA_driver
end module processes
