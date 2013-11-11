!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2013 Wu Li <wu.li@cea.fr>
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

module conductivity
  use config
  use data
  implicit none

contains

  subroutine TConduct(omega,velocity,F_n,ThConductivity)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ThConductivity(Nbands,3,3)

    real(kind=8) :: fBE,tmp(3,3)
    integer(kind=4) :: ii,jj,dir1,dir2

    ThConductivity=0.d0
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          ThConductivity(jj,:,:)=ThConductivity(jj,:,:)+fBE*(fBE+1)*omega(ii,jj)*tmp
       end do
    end do
    ThConductivity=1e21*hbar**2*ThConductivity/(kB*T*T*V*nptk)
  end subroutine TConduct


  subroutine TConductScalar(omega,velocity,F_n,ThConductivity)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands),F_n(Nbands,nptk)
    real(kind=8),intent(out) :: ThConductivity(Nbands)

    real(kind=8) :: fBE
    integer(kind=4) :: ii,jj

    ThConductivity=0.d0
    do jj=1,Nbands
       do ii=2,nptk
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          ThConductivity(jj)=ThConductivity(jj)+fBE*(fBE+1)*omega(ii,jj)*&
               velocity(ii,jj)*F_n(jj,ii)
       end do
    end do
    ThConductivity=1e21*hbar**2*ThConductivity/(kB*T*T*V*nptk)
  end subroutine TConductScalar


  subroutine CumulativeTConduct(omega,velocity,F_n,ticks,results)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ticks(nticks),results(Nbands,3,3,Nticks)

    real(kind=8) :: fBE,tmp(3,3),lambda
    integer(kind=4) :: ii,jj,kk,dir1,dir2

    real(kind=8) :: dnrm2

    do ii=1,nticks
       ticks(ii)=10.**(8.*(ii-1.)/(nticks-1.)-2.)
    end do

    results=0.
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          tmp=fBE*(fBE+1)*omega(ii,jj)*tmp
          lambda=dot_product(F_n(jj,ii,:),velocity(ii,jj,:))/(&
               (omega(ii,jj)+1d-12)*dnrm2(3,velocity(ii,jj,:),1))
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(jj,:,:,kk)=results(jj,:,:,kk)+tmp
             end if
          end do
       end do
    end do
    results=1e21*hbar**2*results/(kB*T*T*V*nptk)
  end subroutine CumulativeTConduct


end module conductivity
