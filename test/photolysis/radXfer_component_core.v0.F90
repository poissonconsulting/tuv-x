! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The radXfer component module

!> The radXfer_component_t type and related functions
module radXfer_component_core

  use photolysis_component,            only : component_t
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk, ik => musica_ik
  use musica_string,                   only : string_t
  use micm_grid_warehouse,             only : grid_warehouse_t
  use micm_Profile_warehouse,          only : Profile_warehouse_t
  use micm_radXfer_xsect_warehouse,    only : radXfer_xsect_warehouse_t
  use micm_radiator_warehouse,         only : radiator_warehouse_t
  use abstract_radXfer,                only : abstract_radXfer_t
  use delta_eddington,                 only : delta_eddington_t
  use disord,                          only : disord_t

  implicit none
  private

  public :: radXfer_component_core_t

  !> radXfer component core
  !!
  !! Calculates the atmospheric radiation field
  type, extends(component_t) :: radXfer_component_core_t
    private
    integer(ik)                              :: nStreams_ = 0_ik
    class(abstract_radXfer_t), pointer       :: radXferSolver_ => null()
    !> Copy of the original radXfer component configuration
    type(config_t)                           :: config_
    !> Copy of original radXfer cross section warehouse
    type(radXfer_xsect_warehouse_t), pointer :: radXferXsectWareHouse_ => null()
    !> Copy of original radiator warehouse
    type(radiator_warehouse_t), pointer      :: RadiatorWareHouse_ => null()
  contains
    !> Returns the name of the component
    procedure :: name => component_name
    !> Returns a description of the component purpose
    procedure :: description
    !> Calculate the radiation field at given solar zenith angle
    procedure :: upDate
    !> Finalize the radXfer core
    final :: finalize
  end type radXfer_component_core_t

  !> Constructor
  interface radXfer_component_core_t
    module procedure constructor
  end interface radXfer_component_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> radXfer component core constructor
  !!
  !! Sets up radXfer objects for solving
  function constructor( config, gridWareHouse, ProfileWareHouse ) result( new_component )

    use musica_assert,                 only : assert, assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New radXfer component core
    type(radXfer_component_core_t), pointer :: new_component
    !> radXfer configuration data
    type(config_t), intent(inout) :: config
    !> grid warehouse
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse
    !> Profile warehouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

    character(len=*), parameter :: Iam = 'radXfer component radXfer component constructor'
    type(string_t)       :: solver

    allocate( new_component )

    ! instantiate and initialize the radXfer cross section warehouse
    new_component%radXferXsectWareHouse_ => radXfer_xsect_warehouse_t( config, gridWareHouse, ProfileWareHouse )
    ! instantiate and initialize the radiator warehouse
    new_component%radiatorWareHouse_ => radiator_warehouse_t( config )

    ! save the configuration (used for preprocessing input data only)
    new_component%config_ = config

    ! get the radiative transfer solver
    call config%get( "Solver", solver, Iam, default="Delta Eddington" )
    if( solver == 'Discrete Ordinance' ) then
      call config%get( "nStreams", new_component%nStreams_, Iam, default=4_ik )
      allocate( disord_t :: new_component%radXferSolver_ )
    else
      allocate( delta_eddington_t :: new_component%radXferSolver_ )
    endif

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Model component name
  type(string_t) function component_name( this )

    !> Arguments
    class(radXfer_component_core_t), intent(in) :: this

    component_name = "radXfer component"

  end function component_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Model component description
  type(string_t) function description( this )

    !> Arguments
    class(radXfer_component_core_t), intent(in) :: this

    description = "calculate the radiation field"

  end function description

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the radiation field
  subroutine upDate( this, zenithAngle, SphericalGeom, GridWareHouse, ProfileWareHouse )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use micm_1d_grid,                  only : abs_1d_grid_t
    use micm_radiator_warehouse,       only : warehouse_iterator_t
    use micm_abs_radiator_type,        only : abs_radiator_t
    use micm_abs_radiator_type,        only : radiator_state_t
    use spherical_geom_type,           only : spherical_geom_t

    !> Arguments
    real(dk), intent(in)                           :: zenithAngle
    !> radXfer component radXfer component
    class(radXfer_component_core_t), intent(inout) :: this
    !> grid warehouse
    type(grid_warehouse_t), intent(inout)          :: GridWareHouse
    !> Profile warehouse
    type(Profile_warehouse_t), intent(inout)       :: ProfileWareHouse
    !> Profile warehouse
    type(spherical_geom_t), intent(inout)          :: SphericalGeom

    !> Local variables
    real(dk), parameter                  :: kfloor = 1.e-36_dk
    real(dk), parameter                  :: precis = 1.e-7_dk

    integer(ik)                          :: nlyr, nlambda
    real(dk), allocatable                :: dscat(:,:)
    real(dk), allocatable                :: dscat_accum(:,:)
    real(dk), allocatable                :: dabs_accum(:,:)
    real(dk), allocatable                :: asym_accum(:,:)
    type(string_t)                       :: Handle
    type(warehouse_iterator_t), pointer  :: iter
    class(abs_radiator_t), pointer       :: aRadiator
    type(radiator_state_t)               :: aRadiatorState
    type(radiator_state_t), allocatable  :: RadiatorStates(:)
    class(abs_1d_grid_t), pointer        :: zGrid
    class(abs_1d_grid_t), pointer        :: lambdaGrid

    Handle = 'Vertical Z'
    zGrid => GridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength'
    lambdaGrid => GridWareHouse%get_grid( Handle )

    nlyr = zGrid%ncells_ ; nlambda = lambdaGrid%ncells_
    allocate( dscat(nlyr,nlambda) )
    allocate( dscat_accum, mold=dscat )
    allocate( dabs_accum, mold=dscat )
    allocate( asym_accum, mold=dscat )
    dscat_accum = 0._dk
    dabs_accum  = 0._dk
    asym_accum  = 0._dk

    allocate( RadiatorStates(0) )
    !> iterate over radiators accumulating radiative properties
    iter => this%RadiatorWareHouse_%get_iterator()
    do while( iter%next() )
      aRadiator => this%RadiatorWareHouse_%get_radiator( iter )
      aRadiatorState = aRadiator%upDateState( GridWareHouse, ProfileWareHouse, this%radXferXsectWareHouse_ )
      RadiatorStates = [RadiatorStates,aRadiatorState]
      dscat = aRadiatorState%layer_OD_*aRadiatorState%layer_SSA_
      dscat_accum = dscat_accum + dscat
      dabs_accum  = dabs_accum + aRadiatorState%layer_OD_*(1._dk - aRadiatorState%layer_SSA_)
      select type( thesolver => this%radXferSolver_ )
        class is (delta_eddington_t)
          asym_accum  = asym_accum + aRadiatorState%layer_G_*dscat
      end select
      deallocate( aRadiator )
    enddo
    deallocate( iter )

    !> set atmosphere radiative properties
    dscat_accum = max( dscat_accum, kfloor )
    dabs_accum  = max( dabs_accum, kfloor )

    associate( OD => aRadiatorState%layer_OD_, SSA => aRadiatorState%layer_SSA_, &
               G  => aRadiatorState%layer_G_ )

    OD = dscat_accum + dabs_accum
    where( dscat_accum == kfloor )
      SSA = kfloor
    elsewhere
      SSA = dscat_accum/OD
    endwhere

    select type( thesolver => this%radXferSolver_ )
      class is (delta_eddington_t)
        G = asym_accum/dscat_accum
      class is (disord_t)
        OD  = max( OD, precis )
        SSA = max( min( SSA,1._dk - precis ),precis )
    end select

    end associate
    
    deallocate( RadiatorStates )

  end subroutine upDate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radXfer component core
  subroutine finalize( this )

    !> radXfer component
    type(radXfer_component_core_t), intent(inout) :: this

    if( associated( this%radXferXsectWareHouse_ ) ) deallocate( this%radXferXsectWareHouse_ )
    if( associated( this%RadiatorWareHouse_ ) ) deallocate( this%RadiatorWareHouse_ )
    if( associated( this%radXferSolver_ ) ) deallocate( this%radXferSolver_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module radXfer_component_core
