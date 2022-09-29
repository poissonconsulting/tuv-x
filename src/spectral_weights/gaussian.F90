! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_spectral_weight_gaussian
  ! The gaussian spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_gaussian_t

  type, extends(spectral_weight_t) :: spectral_weight_gaussian_t
    ! Calculator for Gaussian spectral weight
    real(dk) :: centroid ! The gaussian centroid
  contains
    procedure :: calculate => run
    ! Returns the number of bytes required to pack the spectral weight onto a
    ! buffer
    procedure :: pack_size
    ! Packs the spectral weight onto a character buffer
    procedure :: mpi_pack
    ! Unpacks a spectral weight from a character buffer
    procedure :: mpi_unpack
  end type spectral_weight_gaussian_t

  interface spectral_weight_gaussian_t
    module procedure constructor
  end interface spectral_weight_gaussian_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the Gaussian spectral weight

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: this   ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),            intent(inout) :: config ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'Gaussian spectral weight constructor'
    type(string_t) :: required_keys(2), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "centroid"
    optional_keys(1) = "name"
    call assert_msg( 399197574,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "gaussian spectral wght." )

    allocate( spectral_weight_gaussian_t :: this )
    select type( this )
      class is( spectral_weight_gaussian_t )
        call config%get( 'centroid', this%centroid, Iam )
    end select

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the Gaussian spectral weight

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_gaussian_t),  intent(in)     :: this ! This :f:type:`~tuvx_spectral_weight_gaussian/spectral_weight_gaussian_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    real(dk), parameter         :: rTWO = 2.0_dk
    real(dk)                    :: accum

    class(grid_t), pointer      :: lambdaGrid => null()

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    spectral_weight = exp( -( log( rTWO ) * .04_dk                            &
                      * ( lambdaGrid%mid_ - this%centroid )**2 ) )
    accum = sum( spectral_weight )
    spectral_weight = spectral_weight / accum

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the spectral weight onto a
    ! buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(spectral_weight_gaussian_t), intent(in) :: this ! spectral weight to be packed
    integer,                           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%spectral_weight_t%pack_size( comm ) +                    &
                musica_mpi_pack_size( this%centroid, comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, pos, comm )
    ! Packs the spectral weight onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(spectral_weight_gaussian_t), intent(in)    :: this      ! spectral weight to be packed
    character,                         intent(inout) :: buffer(:) ! memory buffer
    integer,                           intent(inout) :: pos       ! current buffer position
    integer,                           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = pos
    call this%spectral_weight_t%mpi_pack( buffer, pos, comm )
    call musica_mpi_pack( buffer, pos, this%centroid, comm )
    call assert( 803704767, pos - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, pos, comm )
    ! Unpacks the spectral weight onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(spectral_weight_gaussian_t), intent(out)   :: this      ! spectral weight to be unpacked
    character,                         intent(inout) :: buffer(:) ! memory buffer
    integer,                           intent(inout) :: pos       ! current buffer position
    integer,                           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = pos
    call this%spectral_weight_t%mpi_unpack( buffer, pos, comm )
    call musica_mpi_unpack( buffer, pos, this%centroid, comm )
    call assert( 284610975, pos - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_gaussian
