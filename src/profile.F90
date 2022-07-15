! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> Profile type
module tuvx_profile

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: profile_t, profile_ptr

  type ::  profile_t
    !> grid handle
    type(string_t) :: handle_
    !> units
    type(string_t) :: units_
    !> number of wavelength grid cells
    integer(musica_ik)           :: ncells_
    !> scale heigth
    real(musica_dk)              :: hscale_
    !> cell centers
    real(musica_dk), allocatable :: mid_val_(:)
    !> cell edges
    real(musica_dk), allocatable :: edge_val_(:)
    !> cell deltas
    real(musica_dk), allocatable :: delta_val_(:)
    !> layer densities
    real(musica_dk), allocatable :: layer_dens_(:)
    !> layer densities including "exo" model layer
    real(musica_dk), allocatable :: exo_layer_dens_(:)
    !> overhead column burden
    real(musica_dk), allocatable :: burden_dens_(:)
  contains
    !> Returns the units for the profile
    procedure :: units
  end type profile_t

  !> Pointer type for building sets of spectral wght objects
  type :: profile_ptr
    class(profile_t), pointer :: val_ => null( )
  end type profile_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for the profile
  type(string_t) function units( this )

    class(profile_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile
