! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_factory
! Provides a function which creates grids for the 
! :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`.

  use tuvx_grid,                       only : grid_t
  use tuvx_grid_equal_delta,           only : equal_delta_t
  use tuvx_grid_from_csv_file,         only : from_csv_file_t
  use tuvx_grid_from_config,           only : from_config_t

  implicit none

  private
  public :: grid_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function grid_builder( config ) result( new_grid_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> Arguments
    !> Grid configuration data
    type(config_t), intent(inout) :: config

    class(grid_t), pointer :: new_grid_t ! New A :f:type:`~tuvx_grid/grid_t` object

    !> Local variables
    type(string_t) :: grid_type
    character(len=*), parameter :: Iam = 'Grid builder: '

    new_grid_t => null()
    call config%get( 'type', grid_type, Iam )

    select case( grid_type%to_char() )
      case( 'equal interval' )
        new_grid_t => equal_delta_t( config )
      case( 'from csv file' )
        new_grid_t => from_csv_file_t( config )
      case( 'from config file' )
        new_grid_t => from_config_t( config )
      case default
        call die_msg( 460768215, "Invalid grid type: '" &
          // grid_type%to_char()//"'" )
    end select

  end function grid_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_factory
