! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_dose_rates
  ! The dose_rates_t type and related functions

  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t
  use tuvx_spectral_weight,            only : spectral_weight_ptr
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_grid,                       only : grid_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t
  use tuvx_profile,                    only : profile_t

  implicit none

  private
  public :: dose_rates_t

  !> Photolysis rate constant calculator
  type :: dose_rates_t
    ! Spectral weights
    type(spectral_weight_ptr), allocatable :: spectral_weights_(:)
    ! Configuration label for the dose rate
    type(string_t), allocatable          :: handles_(:)
  contains
    ! Returns the dose rates for a given set of conditions
    procedure :: get
    ! Returns the names of each dose rate
    procedure :: labels
    ! Returns the number of bytes required to pack the dose rates onto a
    ! buffer
    procedure :: pack_size
    ! Packs the dose rates onto a character buffer
    procedure :: mpi_pack
    ! Unpacks dose rates from a character buffer
    procedure :: mpi_unpack
    ! Finalize the object
    final :: finalize
  end type dose_rates_t

  !> dose_rates_t constructor
  interface dose_rates_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( dose_rates )
    ! Constructor of dose_rates_t objects

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_spectral_weight_factory,  only : spectral_weight_builder

    !> Dose rate configuration
    type(config_t),            intent(inout) :: config
    !> grid warehouse
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    !> profile warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> New dose rates
    class(dose_rates_t),        pointer      :: dose_rates

    ! Local variables
    character(len=*), parameter :: Iam = "dose_rates_t constructor"

    integer        :: nRates
    type(config_t) :: wght_config, spectral_weight_config
    class(iterator_t), pointer  :: iter
    type(spectral_weight_ptr)   :: spectral_weight_ptr
    character(len=64)           :: keychar
    type(string_t)              :: netcdfFile, Object
    type(string_t)              :: wght_key
    type(string_t), allocatable :: netcdfFiles(:)

    allocate( dose_rates )

    associate( rates => dose_rates )

    allocate( string_t :: rates%handles_(0) )
    allocate( rates%spectral_weights_(0) )

    ! iterate over dose rates
    iter => config%get_iterator( )
    do while( iter%next( ) )
      keychar  = config%key( iter )
      wght_key = keychar
      rates%handles_ = [ rates%handles_, wght_key ]
      call config%get( iter, wght_config, Iam )

      ! get spectral wght
      call wght_config%get( "weights", spectral_weight_config, Iam )
      spectral_weight_ptr%val_ => &
         spectral_weight_builder( spectral_weight_config, grid_warehouse,     &
                                  profile_warehouse )
      rates%spectral_weights_ = [ rates%spectral_weights_,                    &
                                  spectral_weight_ptr ]
    end do
    deallocate( iter )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get( this, grid_warehouse, profile_warehouse, radiation_field,   &
      dose_rates, file_tag )
    ! calculates dose rate constants

    use musica_assert,                 only : die_msg
    use tuvx_constants,                only : hc
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_solver,                   only : radiation_field_t

    !> Dose rate constant calculator
    class(dose_rates_t),       intent(inout) :: this
    !> Warehouses
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> Actinic flux
    type(radiation_field_t),   intent(in)    :: radiation_field
    !> Tag used in file name of output data
    character(len=*),          intent(in)    :: file_tag
    !> Calculated dose rate constants (vertical layer, dose rate type)
    real(dk), allocatable,     intent(inout) :: dose_rates(:,:)

    ! Local variables
    character(len=*), parameter :: Iam = "dose rates calculator:"
    integer               :: wavNdx, rateNdx, nRates
    real(dk), allocatable :: spectral_weight(:)
    real(dk), allocatable :: tmp_spectral_weight(:)
    real(dk), allocatable :: sirrad(:,:)
    character(len=64), allocatable :: annotatedslabel(:)
    class(grid_t),    pointer :: zGrid => null()
    class(grid_t),    pointer :: lambdaGrid => null()
    class(profile_t), pointer :: etfl => null()

    zGrid => grid_warehouse%get_grid( "height", "km" )
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )
    etfl => profile_warehouse%get_profile( "extraterrestrial flux",           &
                                           "photon cm-2 s-1" )

    allocate( tmp_spectral_weight(0) )

    nRates = size( this%spectral_weights_ )
    if( .not. allocated( dose_rates ) ) then
      allocate( dose_rates( zGrid%ncells_ + 1, nRates ) )
    endif

    !> spectral irradiance
    sirrad = radiation_field%edr_ + radiation_field%eup_ + radiation_field%edn_
    do wavNdx = 1, lambdaGrid%ncells_
      sirrad( :, wavNdx ) = sirrad( :, wavNdx ) * etfl%mid_val_( wavNdx )
    enddo

    allocate( annotatedslabel( nRates ) )

rate_loop:                                                                    &
    do rateNdx = 1, nRates
      associate( calc_ftn => this%spectral_weights_( rateNdx )%val_ )
        spectral_weight = calc_ftn%calculate( grid_warehouse,                 &
                                              profile_warehouse )
      end associate

      tmp_spectral_weight = [ tmp_spectral_weight, spectral_weight ]
      dose_rates( :, rateNdx ) = matmul( sirrad, spectral_weight )

      if( allocated( spectral_weight ) ) deallocate( spectral_weight )
    end do rate_loop

    call diagout( 'annotatedslabels.new', this%handles_ )
    call diagout( 'sw.'//file_tag//'.new', tmp_spectral_weight )

    if( associated( zGrid ) ) deallocate( zGrid )
    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )
    if( associated( etfl ) ) deallocate( etfl )

  end subroutine get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function labels( this )
    ! Returns the names of each dose rate

    type(string_t), allocatable     :: labels(:)
    class(dose_rates_t), intent(in) :: this

    labels = this%handles_

  end function labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the dose rates onto a
    ! buffer

    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_spectral_weight_factory,  only : spectral_weight_type_name

    class(dose_rates_t), intent(in) :: this ! dose rates to be packed
    integer,             intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_elem
    type(string_t) :: type_name

    pack_size = musica_mpi_pack_size( allocated( this%spectral_weights_ ),    &
                                      comm )
    if( allocated( this%spectral_weights_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%spectral_weights_ ), comm )
      do i_elem = 1, size( this%spectral_weights_ )
      associate( weight => this%spectral_weights_( i_elem )%val_ )
        type_name = spectral_weight_type_name( weight )
        pack_size = pack_size +                                               &
                    type_name%pack_size( comm ) +                             &
                    weight%pack_size( comm )
      end associate
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%handles_  ), comm )
      do i_elem = 1, size( this%handles_ )
        pack_size = pack_size + this%handles_( i_elem )%pack_size( comm )
      end do
    end if
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the dose rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use tuvx_spectral_weight_factory,  only : spectral_weight_type_name

    class(dose_rates_t), intent(in)    :: this      ! dose rates to be packed
    character,           intent(inout) :: buffer(:) ! memory buffer
    integer,             intent(inout) :: position  ! current buffer position
    integer,             intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          allocated( this%spectral_weights_ ), comm )
    if( allocated( this%spectral_weights_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%spectral_weights_ ), &
                            comm )
      do i_elem = 1, size( this%spectral_weights_ )
      associate( weight => this%spectral_weights_( i_elem )%val_ )
        type_name = spectral_weight_type_name( weight )
        call type_name%mpi_pack( buffer, position, comm )
        call weight%mpi_pack(    buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_pack( buffer, position, allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%handles_ ), comm )
      do i_elem = 1, size( this%handles_ )
        call this%handles_( i_elem )%mpi_pack( buffer, position, comm )
      end do
    end if
    call assert( 258716172, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the dose rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_spectral_weight_factory,  only : spectral_weight_allocate

    class(dose_rates_t), intent(out)   :: this      ! dose rates to be unpacked
    character,           intent(inout) :: buffer(:) ! memory buffer
    integer,             intent(inout) :: position  ! current buffer position
    integer,             intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem, n_elems
    logical :: alloced
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%spectral_weights_( n_elems ) )
      do i_elem = 1, n_elems
      associate( weight => this%spectral_weights_( i_elem )%val_ )
        call type_name%mpi_unpack( buffer, position, comm )
        weight => spectral_weight_allocate( type_name )
        call weight%mpi_unpack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%handles_( n_elems ) )
      do i_elem = 1, n_elems
        call this%handles_( i_elem )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call assert( 143039054, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the dose rates object

    !> Dose rate object
    type(dose_rates_t), intent(inout) :: this

    integer :: ndx

    if( allocated( this%spectral_weights_ ) ) then
      do ndx = 1,size( this%spectral_weights_ )
        if( associated( this%spectral_weights_( ndx )%val_ ) ) then
          deallocate( this%spectral_weights_( ndx )%val_ )
        endif
      enddo
      deallocate( this%spectral_weights_ )
    end if

    if( allocated( this%handles_ ) ) deallocate( this%handles_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_dose_rates
