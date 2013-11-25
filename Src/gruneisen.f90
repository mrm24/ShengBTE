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

! Subroutines and functions related to Grüneisen parameters.
module gruneisen
  use data
  use config
  use input
  implicit none

contains

  ! Subroutine to compute the mode Grüneisen parameters.
  subroutine mode_grun(omega,eigenvect,Ntri,Phi,&
       R_j,R_k,Index_i,Index_j,Index_k,grun)
    implicit none

    integer(kind=4),intent(in) :: Ntri,Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: omega(nptk,nbands)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,nbands,nbands)
    real(kind=8),intent(out) :: grun(nptk,nbands)

    real(kind=8),parameter :: unitfactor=9.6472d4 ! From nm*eV/(amu*A^3*THz^2) to 1.

    integer(kind=4) :: ik,ii,jj,kk,iband,itri,ialpha,ibeta
    real(kind=8) :: kspace(nptk,3)
    complex(kind=8) :: factor1,factor2,factor3,g(nbands)

    do ii=1,Ngrid(1)
       do jj=1,Ngrid(2)
          do kk=1,Ngrid(3)
             ik=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
             kspace(ik,:)=rlattvec(:,1)*(ii-1.0)/ngrid(1)+&
                  rlattvec(:,2)*(jj-1.0)/ngrid(2)+&
                  rlattvec(:,3)*(kk-1.0)/ngrid(3)
          end do
       end do
    end do    
    do ik=1,nptk
       g=0.0
       do iband=1,nbands
          do itri=1,Ntri
             factor1=exp(iunit*dot_product(kspace(ik,:),R_j(:,itri)))/&
                  sqrt(masses(Index_i(itri))*masses(Index_j(itri)))
             do ialpha=1,3
                factor2=factor1*conjg(eigenvect(ik,iband,&
                     3*(Index_i(itri)-1)+ialpha))
                do ibeta=1,3
                   factor3=factor2*eigenvect(ik,iband,&
                        3*(Index_j(itri)-1)+ibeta)
                   g(iband)=g(iband)+factor3*dot_product(&
                        Phi(ialpha,ibeta,:,itri),&
                        (lattvec(:,Index_k(itri))+R_k(:,itri)))
                end do
             end do
          end do
       end do
       g=-unitfactor*g/6.d00/omega(ik,:)**2
       grun(ik,:)=real(g)
    end do
  end subroutine mode_grun

  ! Obtain the total Grüneisen parameter as a weighted sum over modes.
  function total_grun(omega,grun)
    implicit none
    real(kind=8),intent(in) :: omega(nptk,nbands),grun(nptk,nbands)

    integer(kind=4) :: ii,jj
    real(kind=8) :: total_grun,weight,dBE,x

    total_grun=0.
    weight=0.
    do jj=1,nbands
       do ii=1,nptk
          x=hbar*omega(ii,jj)/(2.*kB*T)
          if(x.eq.0.) then
             dBE=1.
          else
             dBE=(x/sinh(x))**2.
          end if
          weight=weight+dBE
          total_grun=total_grun+dBE*grun(ii,jj)
       end do
    end do
    total_grun=total_grun/weight
  end function total_grun
end module gruneisen
