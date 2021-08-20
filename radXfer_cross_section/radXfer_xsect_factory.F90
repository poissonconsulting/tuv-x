! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_cross_section_factory module

!> Builder of cross section calculators
module micm_cross_section_factory

  use micm_abs_cross_section_type,          only : abs_cross_section_t
  use micm_base_cross_section_type,         only : base_cross_section_t
  use micm_tint_cross_section_type,         only : tint_cross_section_t
  use micm_o3_tint_cross_section_type,      only : o3_tint_cross_section_t
  use micm_no2_tint_cross_section_type,     only : no2_tint_cross_section_t
  use micm_h2o2_oh_oh_cross_section_type,   only : h2o2_oh_oh_cross_section_t
  use micm_n2o_n2_o1d_cross_section_type,   only : n2o_n2_o1d_cross_section_t
  use micm_n2o5_no2_no3_cross_section_type, only : n2o5_no2_no3_cross_section_t
  use micm_hno3_oh_no2_cross_section_type,  only : hno3_oh_no2_cross_section_t
  use micm_ch2o_cross_section_type,         only : ch2o_cross_section_t
  use micm_ch3ono2_ch3o_no2_cross_section_type, only : ch3ono2_ch3o_no2_cross_section_t
  use micm_rono2_cross_section_type,        only : rono2_cross_section_t
  use micm_nitroxy_ethanol_cross_section_type, only : nitroxy_ethanol_cross_section_t
  use micm_nitroxy_acetone_cross_section_type, only : nitroxy_acetone_cross_section_t
  use micm_t_butyl_nitrate_cross_section_type, only : t_butyl_nitrate_cross_section_t
  use micm_ch3coch3_ch3co_ch3_cross_section_type, only : ch3coch3_ch3co_ch3_cross_section_t
  use micm_cl2_cl_cl_cross_section_type,       only : cl2_cl_cl_cross_section_t
  use micm_oclo_cross_section_type,            only : oclo_cross_section_t
  use micm_clono2_cross_section_type,          only : clono2_cross_section_t
  use micm_ccl4_cross_section_type,            only : ccl4_cross_section_t
  use micm_chcl3_cross_section_type,           only : chcl3_cross_section_t
  use micm_cfc11_cross_section_type,           only : cfc11_cross_section_t
  use micm_hcfc_cross_section_type,            only : hcfc_cross_section_t
  use micm_bro_br_o_cross_section_type,        only : bro_br_o_cross_section_t
  use micm_hobr_oh_br_cross_section_type,      only : hobr_oh_br_cross_section_t
  use micm_chbr3_cross_section_type,           only : chbr3_cross_section_t

  implicit none

  private
  public :: cross_section_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cross_section_builder( config, mdlLambdaEdge ) result( new_cross_section_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use musica_constants,              only : musica_dk

    !> New rate constant calculator
    class(abs_cross_section_t), pointer :: new_cross_section_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)

    type(string_t) :: cross_section_type
    character(len=*), parameter :: Iam = 'cross section builder: '

    write(*,*) Iam,'entering'
    new_cross_section_t => null( )
    call config%get( 'cross section type', cross_section_type, Iam )

    select case( cross_section_type%to_char() )
      case( 'base cross section' )
        allocate( base_cross_section_t :: new_cross_section_t )
      case( 'tint cross section' )
        allocate( tint_cross_section_t :: new_cross_section_t )
      case( 'O3 cross section' )
        allocate( o3_tint_cross_section_t :: new_cross_section_t )
      case( 'SO2 cross section' )
        allocate( base_cross_section_t :: new_cross_section_t )
      case( 'no2 tint cross section' )
        allocate( no2_tint_cross_section_t :: new_cross_section_t )
      case( 'H2O2+hv->OH+OH cross section' )
        allocate( h2o2_oh_oh_cross_section_t :: new_cross_section_t )
      case( 'N2O+hv->N2+O(1D) cross section' )
        allocate( n2o_n2_o1d_cross_section_t :: new_cross_section_t )
      case( 'N2O5+hv->NO2+NO3 cross section' )
        allocate( n2o5_no2_no3_cross_section_t :: new_cross_section_t )
      case( 'HNO3+hv->OH+NO2 cross section' )
        allocate( hno3_oh_no2_cross_section_t :: new_cross_section_t )
      case( 'CH2O cross section' )
        allocate( ch2o_cross_section_t :: new_cross_section_t )
      case( 'CH3ONO2+hv->CH3O+NO2 cross section' )
        allocate( ch3ono2_ch3o_no2_cross_section_t :: new_cross_section_t )
      case( 'RONO2 cross section' )
        allocate( rono2_cross_section_t :: new_cross_section_t )
      case( 'nitroxy_ethanol+hv->Products cross section' )
        allocate( nitroxy_ethanol_cross_section_t :: new_cross_section_t )
      case( 'nitroxy_acetone+hv->Products cross section' )
        allocate( nitroxy_acetone_cross_section_t :: new_cross_section_t )
      case( 't_butyl_nitrate+hv->Products cross section' )
        allocate( t_butyl_nitrate_cross_section_t :: new_cross_section_t )
      case( 'CH3COCH3+hv->CH3CO+CH3 cross section' )
        allocate( ch3coch3_ch3co_ch3_cross_section_t :: new_cross_section_t )
      case( 'Cl2+hv->Cl+Cl cross section' )
        allocate( cl2_cl_cl_cross_section_t :: new_cross_section_t )
      case( 'OClO+hv->Products cross section' )
        allocate( oclo_cross_section_t :: new_cross_section_t )
      case( 'ClONO2 cross section' )
        allocate( clono2_cross_section_t :: new_cross_section_t )
      case( 'CCl4+hv->Products cross section' )
        allocate( ccl4_cross_section_t :: new_cross_section_t )
      case( 'CHCl3+hv->Products cross section' )
        allocate( chcl3_cross_section_t :: new_cross_section_t )
      case( 'CFC-11+hv->Products cross section' )
        allocate( cfc11_cross_section_t :: new_cross_section_t )
      case( 'HCFC-123+hv->Products cross section' )
        allocate( hcfc_cross_section_t :: new_cross_section_t )
      case( 'BrO+hv->Br+O cross section' )
        allocate( bro_br_o_cross_section_t :: new_cross_section_t )
      case( 'hobr+hv->oh+br cross section' )
        allocate( hobr_oh_br_cross_section_t :: new_cross_section_t )
      case( 'CHBr3+hv->Products cross section' )
        allocate( chbr3_cross_section_t :: new_cross_section_t )
      case default
        call die_msg( 450768214, "Invalid cross section type: '"//              &
                                 cross_section_type%to_char( )//"'" )
    end select
    call new_cross_section_t%initialize( config, mdlLambdaEdge )
    write(*,*) Iam,'exiting'

  end function cross_section_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_cross_section_factory
