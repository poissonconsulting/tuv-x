! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_radiative_transfer_solver
! General interface for radiative transfer solvers

   use musica_constants,               only : dk => musica_dk

   implicit none

   private
   public :: radiative_transfer_solver_t, radiation_field_t

   type :: radiation_field_t
     real(dk), allocatable :: edr_(:,:) ! \todo needs description
     real(dk), allocatable :: eup_(:,:) ! \todo needs description
     real(dk), allocatable :: edn_(:,:) ! \todo needs description
     real(dk), allocatable :: fdr_(:,:) ! \todo needs description
     real(dk), allocatable :: fup_(:,:) ! \todo needs description
     real(dk), allocatable :: fdn_(:,:) ! \todo needs description
   contains
     final :: finalize
   end type radiation_field_t

   type, abstract :: radiative_transfer_solver_t
     contains
     procedure(update_radiation_field), deferred :: update_radiation_field
   end type radiative_transfer_solver_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_streams,       &
      n_layers, spherical_geometry, grid_warehouse, profile_warehouse,        &
      radiator_warehouse ) result( radiation_field )
    ! Solves for the radiation field based on given conditions

     use musica_constants,             only : dk => musica_dk
     use tuvx_grid_warehouse,          only : grid_warehouse_t
     use tuvx_profile_warehouse,       only : profile_warehouse_t
     use tuvx_radiator_warehouse,      only : radiator_warehouse_t
     use tuvx_spherical_geometry,      only : spherical_geometry_t

     import radiative_transfer_solver_t, radiation_field_t

     class(radiative_transfer_solver_t), intent(inout) :: this
     integer, intent(in)                       :: n_streams          ! Number of streams
     integer, intent(in)                       :: n_layers           ! Number of vertical layers
     real(dk), intent(in)                      :: solar_zenith_angle ! Solar zenith angle [degrees]
     type(grid_warehouse_t), intent(inout)     :: grid_warehouse     ! Available grids
     type(profile_warehouse_t), intent(inout)  :: profile_warehouse  ! Available profiles
     type(radiator_warehouse_t), intent(inout) :: radiator_warehouse ! Set of radiators
     type(spherical_geometry_t), intent(inout) :: spherical_geometry ! Spherical geometry calculator

     class(radiation_field_t), pointer         :: radiation_field

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleans up memory for a radiative transfer solver

    type(radiation_field_t), intent(inout) :: this

    if( allocated( this%edr_ ) ) deallocate( this%edr_ )
    if( allocated( this%eup_ ) ) deallocate( this%eup_ )
    if( allocated( this%edn_ ) ) deallocate( this%edn_ )
    if( allocated( this%fdr_ ) ) deallocate( this%fdr_ )
    if( allocated( this%fup_ ) ) deallocate( this%fup_ )
    if( allocated( this%fdn_ ) ) deallocate( this%fdn_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiative_transfer_solver
