! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Module  for FF calculation with ML correction
module xtb_gfnff_ffml
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: Tffml, set_ffml

   ! holds info for ML correction of GFN-FF calculation
   type :: Tffml
     logical  :: fixMD
! contains  ! procedures...
   end type Tffml

contains

! setup type for FF calculation with ML correction
subroutine set_ffml(self)
   class(Tffml), intent(out) :: self
   ! run deterministic MD simulation (0K temperature and fixed shifts)
   self%fixMD = .true.

end subroutine set_ffml


end module xtb_gfnff_ffml
