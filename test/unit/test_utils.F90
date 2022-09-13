module tuvx_test_utils

  implicit none
  public

  interface check_values
    procedure :: check_values_1D, check_values_2D
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_1D( results, expected_results, relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    real(kind=dk), intent(in) :: results(:)
    real(kind=dk), intent(in) :: expected_results(:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_elem

    call assert( 572877200, size( results ) == size( expected_results ) )
    do i_elem = 1, size( results )
      call assert( 743034104, almost_equal(                                  &
        results( i_elem ),                                                   &
        expected_results( i_elem ),                                          &
        relative_tolerance))
    end do

  end subroutine check_values_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_2D( results, expected_results, relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    real(kind=dk), intent(in) :: results(:,:)
    real(kind=dk), intent(in) :: expected_results(:,:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_level, i_wavelength

    call assert( 577098581, &
      size( results, dim = 1 ) == size( expected_results, dim = 1) )
    call assert( 696108875, &
      size( results, dim = 2 ) == size( expected_results, dim = 2) )
    do i_wavelength = 1, size( results, dim = 2 )
      do i_level = 2, size( results, dim = 1 )
        call assert( 179372912, almost_equal(                                 &
          results( i_level, i_wavelength ),                                   &
          expected_results( i_level, i_wavelength ),                          &
          relative_tolerance))
      end do
    end do

  end subroutine check_values_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_test_utils
