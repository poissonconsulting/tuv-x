! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! one dimension grid type
module tuvx_grid

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: grid_t, grid_ptr

  type ::  grid_t
    !> grid handle
    type(string_t) :: handle_
    !> units
    type(string_t) :: units_
    !> number of wavelength grid cells
    integer(musica_ik) :: ncells_
    !> cell centers
    real(musica_dk), allocatable :: mid_(:)
    !> cell edges
    real(musica_dk), allocatable :: edge_(:)
    !> cell deltas
    real(musica_dk), allocatable :: delta_(:)
  contains
    !> Returns the units for the grid
    procedure :: units
  end type grid_t

  !> Pointer type for building sets of spectral wght objects
  type :: grid_ptr
    class(grid_t), pointer :: val_ => null( )
  end type grid_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for the grid
  type(string_t) function units( this )

    class(grid_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid
