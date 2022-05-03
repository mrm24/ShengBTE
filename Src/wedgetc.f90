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

! Code related to the irreducible wedge of the BZ.

module wedgetc
  use symmetry
  use config
  use misc
  use mpi_f08
  implicit none

contains

  ! Find q-points in the irreducible wedge. The Brillouin zone is
  ! discretized into product(ngrid) points through Ngrid(:) divisions
  ! along each lattice vector. List(1:Nlist) gives Nlist points in the
  ! irreducible wedge. Nequi(1:Nlist) gives the number of equivalent
  ! points in the whole Brillouin zone.  These q-points are indexed
  ! along the G1 direction first, then the G2 direction and finally
  ! the G3 direction. In other words, when the single index is
  ! increased, the G1 component changes first, then the G2 component
  ! and finally the G3 component.
  subroutine wedge(Nlist,Nequi,List,ALLEquiList,TypeofSymmetry)
    implicit none

    integer(kind=4),intent(out) :: Nlist,Nequi(nptk),List(nptk)
    integer(kind=4),intent(out) :: ALLEquiList(Nsymm_rot,nptk),TypeofSymmetry(Nsymm_rot,nptk)

    integer(kind=4) :: ID_equi(Nsymm_rot,nptk),ii,ll,Ntot,ierr

    integer(kind=4) :: NAllList,AllList(nptk),EquiList(Nsymm_rot)

    call symmetry_map(ID_equi)

    Ntot=0
    NList=0
    NAllList=0
    List=0
    AllList=0
    do ii=1,nptk
       ! Check if the current q point is equivalent to some point
       ! included in AllList, which contains all q points already
       ! considered and all their images through the symmetry
       ! operations.
       if(.not.any(AllList(1:NAllList).eq.ii)) then
          Nlist=Nlist+1
          List(Nlist)=ii
          Nequi(Nlist)=0
          
          do ll=1,Nsymm_rot
             ! Check if the point generated by the ll-th symmetric
             ! operation on q is equivalent to some points included in
             ! EquiList, which contains all the equivalent points
             ! generated by the previous operations.
             if(.not.any(EquiList(1:Nequi(Nlist)).eq.ID_equi(ll,ii))) then
                Nequi(Nlist)=Nequi(Nlist)+1
                EquiList(Nequi(Nlist))=ID_equi(ll,ii)
                NAllList=NAllList+1
                AllList(NAllList)=ID_equi(ll,ii)
                ALLEquiList(Nequi(Nlist),Nlist)=ID_equi(ll,ii)
                TypeofSymmetry(Nequi(Nlist),Nlist)=ll
             end if
          end do
          Ntot=Ntot+Nequi(Nlist)
       end if
    end do

    ! Print a summary of the findings.
    if(myid.eq.0)write(*,*) "Info: Ntot =",Ntot
    if(myid.eq.0)write(*,*) "Info: Nlist =",Nlist
    if(ntot.ne.nptk) then
       if(myid.eq.0)write(error_unit,*) "Error: ntot!=nptk"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
  end subroutine wedge
end module wedgetc
