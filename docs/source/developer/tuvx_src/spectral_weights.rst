.. TUV-x Spectral Weidhts

Spectral Weights
================

Spectral weights are used in the calculation of dose rates. All spectral weights are built from configuration data described here: :ref:`configuration-spectral-weights`.

To use a spectral weight, configure tuv-x for its use and then pull one
out of the :f:type:`~tuvx_spectral_weight_factory` like this: ::

  function run( config, grid_warehouse, profile_warehouse ) result( thing )
      use musica_config,                 only : config_t
      use musica_iterator,               only : iterator_t
      use musica_string,                 only : string_t
      use tuvx_grid_warehouse,           only : grid_warehouse_t
      use tuvx_profile_warehouse,        only : profile_warehouse_t
      use tuvx_spectral_weight,          only : spectral_weight_ptr
      use tuvx_spectral_weight_factory,  only : spectral_weight_builder

      character(len=*), parameter              :: Iam = "run"
      character(len=64)                        :: keychar
      class(iterator_t), pointer               :: iter
      integer                                  :: thing
      type(config_t)                           :: wght_config, spectral_weight_config
      type(config_t),            intent(inout) :: config
      type(grid_warehouse_t),    intent(inout) :: grid_warehouse
      type(profile_warehouse_t), intent(inout) :: profile_warehouse
      type(spectral_weight_ptr)                :: spectral_weight_ptr
      type(string_t)                           :: wght_key

      iter => config%get_iterator( )
      do while( iter%next( ) )
        keychar  = config%key( iter )
        wght_key = keychar
        call config%get( iter, wght_config, Iam )

      call wght_config%get( "weights", spectral_weight_config, Iam )
      spectral_weight_ptr%val_ => &
         spectral_weight_builder( spectral_weight_config, grid_warehouse,     &
                                  profile_warehouse )
      ...
    end function run

The base Spectral Weight type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight

Spectral Weight Factory
^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_factory

Eppley
^^^^^^
.. f:automodule:: tuvx_spectral_weight_eppley

Exponential Decay
^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_exp_decay

Gaussian Filter
^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_gaussian

Notch Filter
^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_notch_filter

Par
^^^
.. f:automodule:: tuvx_spectral_weight_par

Phytoplankton Boucher
^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_phytoplankton_boucher

Plant Damage Flint Caldwell Ext
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_plant_damage_flint_caldwell_ext

Plant Damage Flint Caldwell
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_plant_damage_flint_caldwell

Plant Damage
^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_plant_damage

SCUP Mice
^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_scup_mice

Standard Human Erythema
^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_standard_human_erythema

UV Index
^^^^^^^^
.. f:automodule:: tuvx_spectral_weight_uv_index
