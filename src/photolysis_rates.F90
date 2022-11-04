! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_photolysis_rates
  ! The photolysis_rates_t type and related functions

  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t
  use tuvx_grid,                       only : grid_t
  use tuvx_cross_section,              only : cross_section_ptr
  use tuvx_grid_warehouse,             only : grid_warehouse_ptr
  use tuvx_profile,                    only : profile_t
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr
  use tuvx_quantum_yield,              only : quantum_yield_ptr

  implicit none

  private
  public :: photolysis_rates_t

  type :: photolysis_rates_t
    ! Photolysis rate constant calculator
    type(cross_section_ptr), allocatable :: cross_sections_(:) ! Absorption cross-sections
    type(quantum_yield_ptr), allocatable :: quantum_yields_(:) ! Quantum yields
    real(dk),                    allocatable :: scaling_factors_(:) ! Scaling factor for final rate constant
    type(string_t), allocatable              :: handles_(:) ! User-provided label for the photolysis rate constant
    logical :: enable_diagnostics_ ! Enable writing diagnostic output, defaults to false
    ! Height grid
    type(grid_warehouse_ptr) :: height_grid_
    ! Wavelength grid
    type(grid_warehouse_ptr) :: wavelength_grid_
    ! Extraterrestrial flux profile
    type(profile_warehouse_ptr) :: etfl_profile_
    ! Air density profile
    type(profile_warehouse_ptr) :: air_profile_
  contains
    ! Returns the photolysis rate constants for a given set of conditions
    procedure :: get
    ! Returns the names of each photolysis reaction
    procedure :: labels
    ! Returns the number of photolysis reactions
    procedure :: size => get_number
    ! Returns the number of bytes required to pack the rates onto a buffer
    procedure :: pack_size
    ! Packs the rates onto a character buffer
    procedure :: mpi_pack
    ! Unpacks rates from a character buffer
    procedure :: mpi_unpack
    ! Finalize the object
    final :: finalize
  end type photolysis_rates_t

  !> photolysis_rates_t constructor
  interface photolysis_rates_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of photolysis_rates_t objects
  function constructor( photolysis_config, grid_warehouse, profile_warehouse )     &
      result( photolysis_rates )

    use musica_assert,                 only : assert, assert_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_cross_section_factory,    only : cross_section_builder
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_builder
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    !> photorates rates
    !> Arguments
    type(config_t),            intent(inout) :: photolysis_config
    !> grid warehouse
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    !> profile warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> New photorates rates
    class(photolysis_rates_t),  pointer      :: photolysis_rates

    !> Local variables
    character(len=*), parameter :: Iam = "photolysis_rates_t constructor"

    real(dk)       :: rate_aliasing_factor
    type(config_t) :: reaction_set, reaction_config
    type(config_t) :: cross_section_config, quantum_yield_config
    class(iterator_t), pointer :: iter
    type(cross_section_ptr) :: a_cross_section_ptr
    type(quantum_yield_ptr) :: a_quantum_yield_ptr
    character(len=64)           :: keychar
    type(string_t)              :: reaction_key
    type(string_t)              :: required_keys(1), optional_keys(1)

    required_keys(1) = "reactions"
    optional_keys(1) = "enable diagnostics"

    call assert_msg( 425103288,                                               &
      photolysis_config%validate( required_keys, optional_keys ),             &
      "Bad configuration data format for photolysis rates." )

    allocate( photolysis_rates )

    associate( rates => photolysis_rates )

    allocate( string_t :: rates%handles_(0) )
    allocate( rates%cross_sections_(0) )
    allocate( rates%quantum_yields_(0) )
    allocate( rates%scaling_factors_(0) )

    call photolysis_config%get( "enable diagnostics",                         &
                photolysis_rates%enable_diagnostics_, Iam, default = .false. )

    rates%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    rates%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    rates%etfl_profile_ = profile_warehouse%get_ptr( "extraterrestrial flux", &
                                                     "photon cm-2 s-1" )
    rates%air_profile_ = profile_warehouse%get_ptr( "air", "molecule cm-3" )

    ! iterate over photo reactions
    call photolysis_config%get( "reactions", reaction_set, Iam )

    iter => reaction_set%get_iterator( )
    do while( iter%next( ) )
      call reaction_set%get( iter, reaction_config, Iam )

      call reaction_config%get( "name", reaction_key, Iam )
      rates%handles_ = [ rates%handles_, reaction_key ]

      ! get cross section first
      call reaction_config%get( "cross section", cross_section_config, Iam )
      a_cross_section_ptr%val_ => cross_section_builder( cross_section_config,  &
                                            grid_warehouse, profile_warehouse )
      rates%cross_sections_ = [ rates%cross_sections_, a_cross_section_ptr ]

      ! now get quantum yield
      call reaction_config%get( "quantum yield", quantum_yield_config, Iam )
      a_quantum_yield_ptr%val_ => quantum_yield_builder( quantum_yield_config,  &
                                            grid_warehouse, profile_warehouse )
      rates%quantum_yields_ = [ rates%quantum_yields_, a_quantum_yield_ptr ]

      ! finally get scaling factor factor
      call reaction_config%get( "scaling factor",                             &
                                 rate_aliasing_factor, Iam, default = 1.0_dk )
      rates%scaling_factors_ = [ rates%scaling_factors_, rate_aliasing_factor ]
    end do
    deallocate( iter )

    call assert( 613491108,                                                   &
                 size( rates%cross_sections_ )                                &
                 == size( rates%quantum_yields_ ) )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate photolysis rate constants
  subroutine get( this, la_srb, spherical_geometry, grid_warehouse,           &
      profile_warehouse, radiation_field, photolysis_rates, file_tag )

    use musica_assert,                 only : assert_msg, die_msg
    use tuvx_constants,                only : hc
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_la_sr_bands,              only : la_sr_bands_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_solver,                   only : radiation_field_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    !> Photolysis rate constant calculator
    class(photolysis_rates_t),  intent(inout) :: this
    !> Spherical geometry
    type(spherical_geometry_t), intent(inout) :: spherical_geometry
    !> Lyman Alpha, Schumann-Runge bands
    type(la_sr_bands_t),        intent(inout) :: la_srb
    !> Grid warehouse
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse
    !> Profile warehouse
    type(profile_warehouse_t),  intent(inout) :: profile_warehouse
    !> Actinic flux
    type(radiation_field_t),    intent(in)    :: radiation_field
    !> Tag used in file name of output data
    character(len=*),           intent(in)    :: file_tag
    !> Calculated photolysis rate constants
    real(dk),                   intent(inout) :: photolysis_rates(:,:)

    !> Local variables
    character(len=*), parameter :: Iam = "photolysis rates calculator"
    integer               :: vertNdx, rateNdx, nRates
    real(dk), allocatable :: airVcol(:), airScol(:)
    real(dk), allocatable :: bin_factor(:)
    real(dk), allocatable :: xsqyWrk(:)
    real(dk), allocatable :: cross_section(:,:)
    real(dk), allocatable :: quantum_yield(:,:)
    real(dk), allocatable :: xsqy(:,:)
    real(dk), allocatable :: actinicFlux(:,:)
    character(len=:),  allocatable :: annotatedRate
    character(len=64), allocatable :: annotatedjlabel(:)
    class(grid_t),    pointer :: zGrid
    class(grid_t),    pointer :: lambdaGrid
    class(profile_t), pointer :: airProfile
    class(profile_t), pointer :: etfl

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    etfl  => profile_warehouse%get_profile( this%etfl_profile_ )

    nRates = size( this%cross_sections_ )
    call assert_msg( 470014831,                                               &
                     size( photolysis_rates, 1 ) == zGrid%ncells_ + 1 .and.   &
                     size( photolysis_rates, 2 ) == nRates,                   &
                     "Bad shape for photolysis rate constant array" )

    bin_factor = ( etfl%mid_val_ * 1.e-13_dk * lambdaGrid%mid_ *              &
                   lambdaGrid%delta_ ) / hc
    actinicFlux = transpose( radiation_field%fdr_ + radiation_field%fup_ +    &
                             radiation_field%fdn_ )
    do vertNdx = 1, zGrid%ncells_ + 1
      actinicFlux( :, vertNdx ) = actinicFlux( :, vertNdx ) * bin_factor
    enddo

    if( this%enable_diagnostics_ ) then
      allocate( annotatedjlabel( nRates ) )
      allocate( xsqyWrk(0) )
    end if

rate_loop:                                                                    &
    do rateNdx = 1, nRates
      associate( calc_ftn => this%cross_sections_( rateNdx )%val_ )
        cross_section = calc_ftn%calculate( grid_warehouse, profile_warehouse )
      end associate
      associate( calc_ftn => this%quantum_yields_( rateNdx )%val_ )
        quantum_yield = calc_ftn%calculate( grid_warehouse, profile_warehouse )
      end associate

      ! O2 photolysis can have special la & srb band handling
      if( trim( this%handles_( rateNdx )%to_char( ) ) == 'O2+hv->O+O' ) then
        airProfile => profile_warehouse%get_profile( this%air_profile_ )
        allocate( airVcol( airProfile%ncells_ ),                              &
                  airScol( airProfile%ncells_ + 1 ) )
        call spherical_geometry%air_mass( airProfile%exo_layer_dens_, airVcol,&
                                          airScol )
        call la_srb%cross_section( grid_warehouse, profile_warehouse, airVcol,&
                                  airScol, cross_section )
        deallocate( airVcol, airScol )
        deallocate( airProfile )
      endif

      if( this%enable_diagnostics_ ) then
      associate( enable => this%enable_diagnostics_ )
        xsqyWrk = [ xsqyWrk, reshape( cross_section * quantum_yield,          &
                                      (/ size( cross_section ) /) ) ]
        annotatedRate = this%handles_( rateNdx )%val_//'.xsect.new'
        call diagout( trim( annotatedRate ), cross_section, enable )
        annotatedRate = this%handles_( rateNdx )%val_//'.qyld.new'
        call diagout( trim( annotatedRate ), quantum_yield, enable )
        annotatedRate = this%handles_( rateNdx )%val_//'.xsqy.new'
        call diagout( trim( annotatedRate ),                                  &
                      cross_section * quantum_yield, enable )
      end associate
      end if

      xsqy = transpose( cross_section * quantum_yield )
      do vertNdx = 1, zGrid%ncells_ + 1
        photolysis_rates( vertNdx, rateNdx ) =                                &
            dot_product( actinicFlux( :, vertNdx ), xsqy( :, vertNdx ) )
      enddo
      if( allocated( cross_section ) ) deallocate( cross_section )
      if( allocated( quantum_yield ) ) deallocate( quantum_yield )
    end do rate_loop

    call diagout( 'annotatedjlabels.new', this%handles_,                      &
      this%enable_diagnostics_ )
    call diagout( 'xsqy.'//file_tag//'.new', xsqyWrk, this%enable_diagnostics_ )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( etfl )

  end subroutine get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the names of each photolysis reaction
  function labels( this )

    !> Photolysis reaction names
    type(string_t), allocatable :: labels(:)
    !> Photolysis rate calculator
    class(photolysis_rates_t), intent(in) :: this

    labels = this%handles_

  end function labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of photolysis reactions
  integer function get_number( this )

    use musica_assert,                 only : assert_msg

    !> Photolysis rate calculator
    class(photolysis_rates_t), intent(in) :: this

    call assert_msg( 472295869, allocated( this%handles_ ),                   &
                     "Photolysis rates not initialized" )
    get_number = size( this%handles_ )

  end function get_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the rates onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_cross_section_factory,    only : cross_section_type_name
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name

    class(photolysis_rates_t), intent(in) :: this ! rates to be packed
    integer,                   intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_elem
    type(string_t) :: type_name

    pack_size = musica_mpi_pack_size( allocated( this%cross_sections_ ), comm )
    if( allocated( this%cross_sections_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%cross_sections_ ), comm )
      do i_elem = 1, size( this%cross_sections_ )
      associate( cross_section => this%cross_sections_( i_elem )%val_ )
        type_name = cross_section_type_name( cross_section )
        pack_size = pack_size +                                               &
                    type_name%pack_size( comm ) +                             &
                    cross_section%pack_size( comm )
      end associate
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( allocated( this%quantum_yields_ ), comm )
    if( allocated( this%quantum_yields_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%quantum_yields_ ), comm )
      do i_elem = 1, size( this%quantum_yields_ )
      associate( quantum_yield => this%quantum_yields_( i_elem )%val_ )
        type_name = quantum_yield_type_name( quantum_yield )
        pack_size = pack_size +                                               &
                    type_name%pack_size( comm ) +                             &
                    quantum_yield%pack_size( comm )
      end associate
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( this%scaling_factors_, comm ) +         &
                musica_mpi_pack_size( allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%handles_ ), comm )
      do i_elem = 1, size( this%handles_ )
        pack_size = pack_size + this%handles_( i_elem )%pack_size( comm )
      end do
    end if
    pack_size = pack_size +                                                   &
                musica_mpi_pack_size( this%enable_diagnostics_, comm ) +      &
                this%height_grid_%pack_size( comm ) +                         &
                this%wavelength_grid_%pack_size( comm ) +                     &
                this%etfl_profile_%pack_size( comm ) +                        &
                this%air_profile_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use tuvx_cross_section_factory,    only : cross_section_type_name
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name

    class(photolysis_rates_t), intent(in)    :: this      ! rates to be packed
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_pack( buffer, position, allocated( this%cross_sections_ ),&
                          comm )
    if( allocated( this%cross_sections_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%cross_sections_ ),   &
                            comm )
      do i_elem = 1, size( this%cross_sections_ )
      associate( cross_section => this%cross_sections_( i_elem )%val_ )
        type_name = cross_section_type_name( cross_section )
        call type_name%mpi_pack(     buffer, position, comm )
        call cross_section%mpi_pack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_pack( buffer, position, allocated( this%quantum_yields_ ),&
                          comm )
    if( allocated( this%quantum_yields_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%quantum_yields_ ),   &
                            comm )
      do i_elem = 1, size( this%quantum_yields_ )
      associate( quantum_yield => this%quantum_yields_( i_elem )%val_ )
        type_name = quantum_yield_type_name( quantum_yield )
        call type_name%mpi_pack(     buffer, position, comm )
        call quantum_yield%mpi_pack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_pack( buffer, position, this%scaling_factors_, comm )
    call musica_mpi_pack( buffer, position, allocated( this%handles_ ), comm )
    if( allocated( this%handles_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%handles_ ), comm )
      do i_elem = 1, size( this%handles_ )
        call this%handles_( i_elem )%mpi_pack( buffer, position, comm )
      end do
    end if
    call musica_mpi_pack( buffer, position, this%enable_diagnostics_, comm )
    call this%height_grid_%mpi_pack(     buffer, position, comm )
    call this%wavelength_grid_%mpi_pack( buffer, position, comm )
    call this%etfl_profile_%mpi_pack(    buffer, position, comm )
    call this%air_profile_%mpi_pack(     buffer, position, comm )
    call assert( 707537257, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Packs the rates onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_cross_section_factory,    only : cross_section_allocate
    use tuvx_quantum_yield_factory,    only : quantum_yield_allocate

    class(photolysis_rates_t), intent(out)   :: this      ! rates to be unpacked
    character,                 intent(inout) :: buffer(:) ! memory buffer
    integer,                   intent(inout) :: position  ! current buffer position
    integer,                   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem, n_elems
    logical :: alloced
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%cross_sections_( n_elems ) )
      do i_elem = 1, n_elems
      associate( cross_section => this%cross_sections_( i_elem ) )
        call type_name%mpi_unpack(     buffer, position, comm )
        cross_section%val_ => cross_section_allocate( type_name )
        call cross_section%val_%mpi_unpack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%quantum_yields_( n_elems ) )
      do i_elem = 1, n_elems
      associate( quantum_yield => this%quantum_yields_( i_elem ) )
        call type_name%mpi_unpack(     buffer, position, comm )
        quantum_yield%val_ => quantum_yield_allocate( type_name )
        call quantum_yield%val_%mpi_unpack( buffer, position, comm )
      end associate
      end do
    end if
    call musica_mpi_unpack( buffer, position, this%scaling_factors_, comm )
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elems, comm )
      allocate( this%handles_( n_elems ) )
      do i_elem = 1, n_elems
        call this%handles_( i_elem )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call musica_mpi_unpack( buffer, position, this%enable_diagnostics_, comm )
    call this%height_grid_%mpi_unpack(     buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack( buffer, position, comm )
    call this%etfl_profile_%mpi_unpack(    buffer, position, comm )
    call this%air_profile_%mpi_unpack(     buffer, position, comm )
    call assert( 534021580, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the photorates rates
  subroutine finalize( this )

    !> Photolysis rate calculator
    type(photolysis_rates_t), intent(inout) :: this

    integer :: ndx

    if( allocated( this%cross_sections_ ) ) then
      do ndx = 1,size( this%cross_sections_ )
        if( associated( this%cross_sections_( ndx )%val_ ) ) then
          deallocate( this%cross_sections_( ndx )%val_ )
        endif
      enddo
      deallocate( this%cross_sections_ )
    end if

    if( allocated( this%quantum_yields_ ) ) then
      do ndx = 1,size( this%quantum_yields_ )
        if( associated( this%quantum_yields_( ndx )%val_ ) ) then
          deallocate( this%quantum_yields_( ndx )%val_ )
        endif
      enddo
      deallocate( this%quantum_yields_ )
    end if

    if( allocated( this%scaling_factors_ ) ) then
      deallocate( this%scaling_factors_ )
    end if
    if( allocated( this%handles_ ) ) then
      deallocate( this%handles_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_photolysis_rates
