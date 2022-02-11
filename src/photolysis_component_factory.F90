! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Photolysis component factory

!> Builder of photolysis components
module photolysis_component_factory

  use photolysis_component,     only : component_t

  implicit none
  private

  public :: component_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build a photolysis component by name
  function component_builder( config, gridWareHouse, ProfileWareHouse ) result( new_component )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use radXfer_component_core,        only : radXfer_component_core_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_Profile_warehouse,        only : Profile_warehouse_t

    !> New photolysis component
    class(component_t), pointer   :: new_component
    !> Photolysis component configuration data
    type(config_t), intent(inout) :: config
    !> grid warehouse
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse
    !> Profile warehouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

    character(len=*), parameter :: my_name = "photolysis component builder"
    type(string_t) :: component_type

    new_component => null()
    call config%get( 'type', component_type, my_name )

    select case( component_type%to_char() )
      case( 'Radiative transfer' )
        new_component => radXfer_component_core_t( config, gridWareHouse, ProfileWareHouse )
      case default
        call die_msg( 935006810, "Unsupported model component type: '"//        &
                                 component_type%to_char( )//"'" )
    end select

  end function component_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photolysis_component_factory
