! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile
  ! Profile type

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: profile_t, profile_ptr

  type profile_t
    type(string_t) :: handle_ ! grid handle
    type(string_t) :: units_ ! units
    integer(musica_ik)           :: ncells_ ! number of wavelength grid cells
    real(musica_dk)              :: hscale_ ! scale heigth
    real(musica_dk), allocatable :: mid_val_(:) ! cell centers
    real(musica_dk), allocatable :: edge_val_(:) ! cell edges
    real(musica_dk), allocatable :: delta_val_(:) ! cell deltas
    real(musica_dk), allocatable :: layer_dens_(:) ! layer densities
    real(musica_dk), allocatable :: exo_layer_dens_(:) ! layer densities including "exo" model layer
    real(musica_dk), allocatable :: burden_dens_(:) ! overhead column burden
  contains
    procedure :: units
  end type profile_t

  type profile_ptr
    ! Pointer type for building sets of profile objects
    class(profile_t), pointer :: val_ => null( )
  end type profile_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function units( this )
    ! Returns the units for the profile

    class(profile_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile
