!> \file drive.F90
!! \brief main program

!> NEK5000: Spectral Element Computational Fluid Dynamics Solver
!! COPYRIGHT (c) 2008-2010 UCHICAGO ARGONNE, LLC
!!
!! The UChicago Argonne, LLC as Operator of Argonne National
!! Laboratory holds copyright in the Software. The copyright holder
!! reserves all rights except those expressly granted to licensees,
!! and U.S. Government license rights.

!> License
!!
!!    NEK5000 is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    NEK5000 is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with NEK5000.  If not, see <http://www.gnu.org/licenses/>.

program NEKTON
  implicit none

  integer :: intracomm

  call nek_init(intracomm)
  call nek_solve()
  call nek_end()

  call exitt()
END PROGRAM
