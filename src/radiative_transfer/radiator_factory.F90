! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiator_factory
! Builds :f:type:`~tuvx_radiator/radiator_t` s for
! :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`. 

  use tuvx_radiator,                   only : radiator_t
  use tuvx_radiator_aerosol,           only : radiator_aerosol_t

  implicit none

  private
  public :: radiator_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function radiator_builder( config, grid_warehouse ) result( new_radiator )
    ! Builder of :f:type:`~tuvx_radiator/radiator_t` objects

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    type(config_t),         intent(inout) :: config        ! Radiator configuration data
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(radiator_t),      pointer       :: new_radiator  ! New :f:type:`~tuvx_radiator/radiator_t` object

    ! Local variables
    type(string_t) :: radiator_type
    character(len=*), parameter :: Iam = 'Radiator builder'

    new_radiator => null()
    call config%get( 'type', radiator_type, Iam )

    select case( radiator_type%to_char() )
      case( 'base' )
        new_radiator => radiator_t( config, grid_warehouse )
      case( 'aerosol' )
        new_radiator => radiator_aerosol_t( config, grid_warehouse )
      case default
        call die_msg( 460768245, "Invalid radiator type: '"//                 &
                                 radiator_type%to_char()//"'" )
    end select

  end function radiator_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_factory
