!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2023 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2023 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2023 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2023 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!  Copyright (C) 2021-2022 Fanchen Meng <fanchem@g.clemson.edu>
!  Copyright (C) 2022-2023 Ben Durham <bd740@york.ac.uk>
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

! Routines to calculate the density of states and related quantities.
module dos_routines
  use config
  use mpi_f08
  use data, only : dp
  implicit none

contains

  ! Compute a first estimate of the broadening for each mode.
  subroutine calc_sigma0(velocity,sigma)
    implicit none

    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    real(kind=dp),intent(out) :: sigma(nbands,nptk)

    integer :: ii,jj

    do ii=1,nptk
       do jj=1,nbands
          sigma(jj,ii)=base_sigma(velocity(:,jj,ii))
       end do
    end do
  end subroutine calc_sigma0

  ! Compute the 25th and 75th percentiles of log(sigma)
  subroutine calc_percentiles(sigma,per25,per75)
    implicit none

    real(kind=dp),intent(in) :: sigma(nbands,nptk)
    real(kind=dp),intent(out) :: per25,per75

    integer :: ii,jj,pos
    real(kind=dp) :: logsigma(nptk*nbands)

    do ii=1,nptk
       do jj=1,nbands
          pos=(ii-1)*nbands+jj
          if(sigma(jj,ii)==0.0_dp) then
             logsigma(pos)=-huge(logsigma)
          else
             logsigma(pos)=log(sigma(jj,ii))
          end if
       end do
    end do

    call dlasrt("I",nptk*nbands,logsigma,ii)
    per25=logsigma((25*nptk*nbands)/100)
    per75=logsigma((75*nptk*nbands)/100)
  end subroutine calc_percentiles

  ! Refine the initial estimate of smearing to avoid huge
  ! peaks in the DOS.
  subroutine refine_sigma(sigma)
    implicit none

    real(kind=dp),intent(inout) :: sigma(nbands,nptk)

    real(kind=dp) :: per25,per75,delta,lbound

    call calc_percentiles(sigma,per25,per75)
    delta=per75-per25
    lbound=exp(per25-1.5_dp*delta)
    sigma=max(sigma,lbound)
  end subroutine refine_sigma

  ! Compute the DOS, the projected DOS and the isotopic scattering rates.
  subroutine calc_dos(omega,velocity,eigenvect,ticks,dos,pdos)
    implicit none

    real(kind=dp),intent(in) :: omega(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    complex(kind=dp),intent(in) :: eigenvect(nbands,nbands,nptk)
    
    real(kind=dp),intent(out) :: ticks(nticks)
    real(kind=dp),intent(out) :: dos(nticks)
    real(kind=dp),intent(out) :: pdos(nticks,natoms)

    integer :: ii,jj,kk,mm
    real(kind=dp) :: sigma(nbands,nptk),thisomega,thissigma,weight,prod
    REAL(kind=dp)  EMIN,EMAX

    call calc_sigma0(velocity,sigma)
    call refine_sigma(sigma)
    dos=0.0_dp
    pdos=0.0_dp

    EMIN=1.0e10_dp
    EMAX=-1.0e10_dp
    DO ii=1,NPTK
    DO jj=1,Nbands
       emin = min(emin, omega(jj,ii))
       emax = max(emax, omega(jj,ii))
    ENDDO
    ENDDO
    emax=emax*1.1_dp

    do ii=1,nticks
       ticks(ii)=emin+(emax-emin)*real(ii,kind=dp)/real(nticks,kind=dp)
    enddo
    do mm=1,nticks
          thisomega=ticks(mm)
          if(thisomega==0.0_dp) then
             cycle
          end if
          do ii=1,nptk
             do jj=1,Nbands
                thissigma=sigma(jj,ii)
                weight=exp(-(thisomega-omega(jj,ii))**2/(thissigma**2))&
                     /thissigma/sqrt(pi)
                if(abs(thisomega-omega(jj,ii)).lt.2.5_dp*thissigma) then
                   dos(mm)=dos(mm)+weight
                   do kk=1,natoms
                      prod=(abs(dot_product(&
                           eigenvect(((kk-1)*3+1):((kk-1)*3+3),jj,ii),&
                           eigenvect(((kk-1)*3+1):((kk-1)*3+3),jj,ii))))
                      pdos(mm,kk)=pdos(mm,kk)+weight*prod
                   end do
                end if
             end do
          end do
    end do
    dos=dos/real(nptk,kind=dp)
    pdos=pdos/real(nptk,kind=dp)
  end subroutine calc_dos

  ! Compute all the isotopic information, including the matrix elements for
  ! the iterative procedure in global object (TODO: remove global object, and
  ! pass it by reference)
  subroutine calc_isotopescatt(omega,velocity,eigenvect,nlist,list,&
       rate_scatt_isotope)
    implicit none

    real(kind=dp),intent(in) :: omega(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    complex(kind=dp),intent(in) :: eigenvect(nbands,nbands,nptk)
    integer,intent(in) :: nlist
    integer,intent(in) :: list(nptk)

    real(kind=dp),intent(out) :: rate_scatt_isotope(nbands,nlist)

    integer :: ii,jj,kk,mm,nn,ib,iq, ncount
    real(kind=dp) :: sigma(nbands,nptk),thisomega,thissigma,weight,prod,GammaIso

    call calc_sigma0(velocity,sigma)
    call refine_sigma(sigma)

    rate_scatt_isotope=0.0_dp

    !Driver
    ! Get sizes to store matrix elements
    do nn=1,nstates
      mm=myid*nstates+nn
      if (mm.gt.nlist*nbands) then
         cycle
      end if
      ! Get the indexes
      ib=modulo(mm-1,Nbands)+1
      iq=int((mm-1)/Nbands)+1
      thisomega=omega(ib,list(iq))
      if(thisomega==0.0_dp) then
         cycle
      end if
      do ii=1,nptk
         do jj=1,Nbands
            thissigma=sigma(jj,ii)
            if(abs(thisomega-omega(jj,ii)).lt.2.5_dp*thissigma) then
                  IsoInfo%nisotopic = IsoInfo%nisotopic + 1
            end if
         end do
      end do
    end do ! nn

    ! Allocate isotopic information
    allocate(IsoInfo%Indof1stPhononIso(IsoInfo%nisotopic),&
             IsoInfo%Indof2ndPhononIso(IsoInfo%nisotopic),&
             IsoInfo%Gamma_isotopic(IsoInfo%nisotopic))

    ! Compute isotopic scattering information
    ncount = 0
    do nn=1,nstates
      mm=myid*nstates+nn
      if (mm.gt.nlist*nbands) then
         cycle
      end if
      ! Get the indexes
      ib=modulo(mm-1,Nbands)+1
      iq=int((mm-1)/Nbands)+1
      thisomega=omega(ib,list(iq))
      if(thisomega==0.0_dp) then
         cycle
      end if
      do ii=1,nptk
         do jj=1,Nbands
            thissigma=sigma(jj,ii)
            weight=exp(-(thisomega-omega(jj,ii))**2/(thissigma**2))&
                 /thissigma/sqrt(pi)
            GammaIso = 0.0_dp
            if(abs(thisomega-omega(jj,ii)).lt.2.5_dp*thissigma) then
                  do kk=1,natoms
                     prod=(abs(dot_product(&
                          eigenvect(((kk-1)*3+1):((kk-1)*3+3),nn,list(mm)),&
                          eigenvect(((kk-1)*3+1):((kk-1)*3+3),jj,ii))))**2
                     GammaIso=GammaIso+&
                          weight*prod*gfactors(types(kk))
                  end do
                  GammaIso = GammaIso/&
                     (2.0_dp*nptk)*pi*thisomega**2
                  ncount = ncount + 1
                  IsoInfo%Indof1stPhononIso(ncount) = mm
                  IsoInfo%Indof2ndPhononIso(ncount) = ii * Nbands + jj
                  IsoInfo%Gamma_isotopic(ncount) = GammaIso
                  rate_scatt_isotope(ib,iq) = rate_scatt_isotope(ib,iq)+&
                     GammaIso
            end if
         end do
      end do
    end do ! nn

    !Now collapse over procs
    call MPI_ALLREDUCE(MPI_IN_PLACE,rate_scatt_isotope,Nbands*Nlist,&
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,nn)

  end subroutine calc_isotopescatt
end module dos_routines
