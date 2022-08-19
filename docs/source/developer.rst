.. Instructions for TUV-x developers

Science Contributions
=====================

To submit a new feature or bug-fix, first create an issue on the
`TUV-x GitHub repository <https://github.com/NCAR/tuv-x>`_
that describes the new feature or bug. The software contribution can
then be submitted as a pull request. Please reference the issue in
the pull request description.

Software contributions must apply the
`MUSICA FORTRAN Style Guide <https://ncar.github.io/musica-core/html/coding_style.html>`_.
See also the `MUSICA Recommendations for Contributors <https://ncar.github.io/musica-core/html/contributors.html>`_.

Developers must ensure that new features have complete unit test
coverage. Pull requests that decrease the test coverage of the
TUV-x code base will not be accepted. Developers must also ensure that
new features are fully run-time configurable.

Cross Sections
--------------

The base functionality for cross section calculations is described in
``src/cross_section.F90``.
The base cross section type :f:type:`~tuvx_cross_section/cross_section_t` defined in this
module includes a constructor that reads a set of NetCDF files
specified in the configuration JSON object as an array of file
paths with the key ``netcdf files``.

Here is an exmple configuration for a base cross section:

.. code-block:: JSON

   {
     "type": "base",
     "netcdf files": [ "path/to/my/netcdf/file.nc" ]
   }


More description of the cross section configuration can
be found :ref:`here <configuration-cross-sections>`.
The reading of NetCDF file data is available to cross section
subclasses through the :f:func:`tuvx_cross_section/base_constructor` function.

Data from each NetCDF file in the array in the configuration data
will be loaded into an element of the
``cross_section_parms`` data member. If a NetCDF variable named
``cross_section_parameters`` is present, it will be used to populate
the ``array(:,:)`` data member of the ``cross_section_parms_t`` object.
This first dimension of the array is wavelength, and will be interpolated
to the native TUV-x wavelength grid if this differs from the wavelength
grid in the netcdf file (named ``wavelength``). The second dimension
is used to accomodate multiple wavelength-resolved parameters for
cross section calculations.

If a NetCDF variable named ``temperature`` is present, it will be
used to populate the ``temperature(:)`` data member of the
``cross_section_parms_t`` object.

The calculation of cross sections is done by calling the ``calculate()``
type-bound procedure on a :f:type:`~tuvx_cross_section/cross_section_t` object.

The base-class calculation of cross sections returns the
wavelength-interpolated first parameter (``array(:,1)``) from the first
NetCDF file (``cross_section_parms(1)``) specified in the configuration
data.

Cross section subclasses are located in the ``src/cross_sections/`` folder.
These each provide unique algorithms for calculating cross sections.

Before adding a new cross section class, first check to make sure there
is not an existing class that can be configured to accomodate your
needs. If you find that an existing cross section subclass could be used
by moving one or more hard-coded parameters to the configuration data, this
is preferable to adding a new subclass.

If you determine that a new cross section subclass is needed, this can be
done in four steps:

- :ref:`cs-create-subclass`
- :ref:`cs-add-to-build-scripts`
- :ref:`cs-add-to-factory`
- :ref:`cs-create-unit-test`

.. _cs-create-subclass:

Create subclass module
^^^^^^^^^^^^^^^^^^^^^^

First, choose a unique name for your cross section calculation.
Ideally, this name will describe the algorithm, rather than
the specific photolysis reaction you are applying it to.
However, many subclasses currently in TUV-x are named for
specific photolysis reactions.
For this example, we will use the name ``foo`` for our
cross section algorithm.

**Pay special attention to naming of files, modules, types, and functions
in these instructions.**

Create a file to hold your new subclass module in ``src/cross_sections/`` named
``foo.F90``. The general layout of the module will be (comments have been omitted
for this example, but should be included in an actual module):

.. code-block:: fortran

   ! Copyright (C) 2020 National Center for Atmospheric Research
   ! SPDX-License-Identifier: Apache-2.0
   !
   module tuvx_cross_section_foo

     use tuvx_cross_section,              only : cross_section_t

     implicit none

     private
     public :: cross_section_foo_t

     type, extends(cross_section_t) :: cross_section_foo_t
     contains
       procedure :: calculate
     end type cross_section_foo_t

     interface cross_section_foo_t
       module procedure constructor
     end interface cross_section_foo_t

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function constructor( config, grid_warehouse, profile_warehouse )           &
         result( this )

       use musica_assert,                 only : assert_msg
       use musica_config,                 only : config_t
       use musica_string,                 only : string_t
       use tuvx_cross_section,            only : base_constructor
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t

       class(cross_section_t),    pointer       :: this
       type(config_t),            intent(inout) :: config
       type(grid_warehouse_t),    intent(inout) :: grid_warehouse
       type(profile_warehouse_t), intent(inout) :: profile_warehouse

       type(string_t) :: required_keys(1), optional_keys(1)

       ! This block of code ensures that the configuration keys are valid for
       ! your class. These can be modified to fit your needs. The first
       ! argument to assert_msg() should be a unique integer code for this error.
       required_keys(1) = "type"
       optional_keys(1) = "name"
       call assert_msg( 465568611,                                               &
                        config%validate( required_keys, optional_keys ),         &
                        "Bad configuration data format for "//                   &
                        "foo cross section." )

       allocate( cross_section_foo_t :: this )

       ! You can call the base_constructor function to load data from NetCDF
       ! files into the `cross_section_parms(:)` data member according to the
       ! standard base class logic. Alternatively, you can perform custom
       ! initialization of the subclass object here.
       call base_constructor( this, config, grid_warehouse, profile_warehouse )

     end function constructor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
         reuslt( cross_section )

       use musica_constants,              only : dk => musica_dk
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t

       class(cross_section_foo_t), intent(in) :: this
       type(grid_warehouse_t),     intent(inout) :: grid_warehouse
       type(profile_warehouse_t),  intent(inout) :: profile_warehouse
       ! This flag indicates whether the cross-section data should be calculated
       ! at mid-points on the vertical grid. If it is false or omitted, cross-
       ! section data are calculated at interfaces on the vertical grid.
       logical, optional,          intent(in)    :: at_mid_point
       real(kind=dk), allocatable                :: cross_section(:,:)

       ! Do your calculation here

     end function calculate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end module tuvx_cross_section_foo

The constructor function is reponsible for initializing new instances of your cross
section subclass.
First, you allocate the pointer to be returned as your new type
(``cross_section_foo_t`` in this example).
Then you initialize its data members.
If you just want to use the default initialization of the base class,
you can call the ``base_constructor()`` function as shown above.
You can alternatively initialize data members of the base class
(``cross_section_parms(:)``) directly in this function or add data members to your
subclass and initialize them here (see ``src/cross_sections/o3_tint.F90`` for an example).

The ``calculate()`` function overrides the base-class ``calculate()`` function and will
be called when a user calls the ``calculate()`` type-bound procedure on an instance of
your new subclass.
You can access grid and profile data from the “warehouse” objects passed in as function
arguments, and any data in the base-class data members or in data members you’ve added
to your subclass to perform your calculations.
See the files in ``src/cross_sections/`` for examples of how to access this data in
the ``calculate()`` function.


.. _cs-add-to-build-scripts:

Add subclass module to build scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To include your new class in the build, edit the ``src/cross_sections/CMakeLists.txt`` file
and add your file name to the list saved as ``SRC``.
Files are in alphabetical order.

.. code-block:: cmake

   ################################################################################
   # Cross section source

   set(SRC acetone-ch3co_ch3.F90
           bro-br_o.F90
           ccl4.F90
           cfc-11.F90
           chbr3.F90
           chcl3.F90
           ch3ono2-ch3o_no2.F90
           ch2o.F90
           cl2-cl_cl.F90
           clono2.F90
           foo.F90
           h2o2-oh_oh.F90
           hcfc.F90
           hno3-oh_no2.F90
           hobr-oh_br.F90
           n2o-n2_o1d.F90
           n2o5-no2_no3.F90
           nitroxy_acetone.F90
           nitroxy_ethanol.F90
           no2_tint.F90
           o3_tint.F90
           oclo.F90
           rono2.F90
           t_butyl_nitrate.F90
           tint.F90
           rayliegh.F90
           )

   list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
   set(CROSS_SECTION_SRC ${SRC} PARENT_SCOPE)

   ################################################################################


.. _cs-add-to-factory:

Add subclass to factory function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to use your new subclass, you will need to add it to the
``tuvx_cross_section_factory`` module in ``src/cross_section_factory.F90``.
First, use-associate your new class at the module level:

.. code-block:: fortran

   use tuvx_cross_section_foo,            only : cross_section_foo_t

Then, inside the ``cross_section_builder()`` function, add these lines to the
``select case`` block:

.. code-block:: fortran

   case( 'foo' )
     new_cross_section => cross_section_foo_t( config, grid_warehouse,          &
                                               profile_warehouse )

Now, when you add a cross section of type ``foo`` to the configuration data,
an instance of your new subclass will be created.


.. _cs-create-unit-test:

Create unit test
^^^^^^^^^^^^^^^^

The last step to adding a cross section is to create a unit test.
This will ensure that your calculations are doing what you intended.
It will also serve as an example for how users can configure and use your
new subclass.

See :ref:`developer-add-test` for more details.

Dose Rates
----------

Dose rates apply a spectral weight to the radiation field at each
interface on the vertical grid.
The configuration for a dose rate is:


.. code-block:: JSON
   :force:

   {
     "weights": { ... }
   }

The value of ``weights`` defines the spectral weight
used to calculate the dose rate.
The standard spectral weight configuration is described
:ref:`here <configuration-spectral-weights>`.

If a new dose rate requires an algorithm for calculating the
spectral weight that TUV-x does not currently support, a new
spectral weight algorithm can be introduced in four steps:

- :ref:`dose-rate-create-subclass`
- :ref:`dose-rate-add-to-build-scripts`
- :ref:`dose-rate-add-to-factory`
- :ref:`dose-rate-create-unit-test`


.. _dose-rate-create-subclass:

Create subclass module
^^^^^^^^^^^^^^^^^^^^^^

First, choose a unique name for your spectral weight algorithm.
Ideally, this name will describe the algorithm, rather than
the specific dose rate you are applying it to.

**Pay special attention to the naming of files, modules, types, and
functions in these instructions.**

Create a file to hold your new subclass module in ``src/spectral_weights/``
named ``foo.F90``.
The general layout of the module will be (comments have been omitted
in this example, but should be included in an actual module):

.. code-block:: fortran

   ! Copyright (C) 2020 National Center for Atmospheric Research
   ! SPDX-License-Identifier: Apache-2.0
   !
   module tuvx_spectral_weight_foo

     use tuvx_spectral_weight,            only : spectral_weight_t

     implicit none

     private
     public :: spectral_weight_foo_t

     type, extends(spectral_weight_t) :: spectral_weight_foo_t
     contains
       procedure :: calculate
     end type spectral_weight_t

     interface spectral_weight_t
       module procedure :: constructor
     end interface spectral_weight_t

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function constructor( config, grid_warehouse, profile_warehouse )           &
         result ( this )

       use musica_assert,                 only : assert_msg
       use musica_config,                 only : config_t
       use musica_string,                 only : string_t
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t
       use tuvx_spectral_weight,          only : base_constructor

       class(spectral_weight_t),  pointer       :: this
       type(config_t),            intent(inout) :: config
       type(grid_warehouse_t),    intent(inout) :: grid_warehouse
       type(profile_warehouse_t), intent(inout) :: profile_warehouse

       type(string_t) :: required_keys(1), optional_keys(1)

       ! This block of code ensures that the configuration keys are valid for
       ! your class. These can be modified to fit your needs. The first
       ! argument to assert_msg() should be a unique integer code for this error.
       required_keys(1) = "type"
       optional_keys(1) = "name"
       call assert_msg( 407417332,                                               &
                        config%validate( required_keys, optional_keys ),         &
                        "Bad configuration data format for "//                   &
                        "foo spectral weight." )

       allocate( spectral_weight_foo_t :: this )

       ! You can call the base_constructor function to load data from NetCDF
       ! files into the `spectral_weight_parms(:)` data member according to the
       ! standard base class logic. Alternatively, you can perform custom
       ! initialization of the subclass object here.
       call base_constructor( this, config, grid_warehouse, profile_warehouse )

     end function constructor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine calculate( this, grid_warehouse, profile_warehouse )             &
         result( spectral_weight )

       use musica_constants,              only : dk => musica_dk
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t

       class(spectral_weight_foo_t),  intent(in)    :: this
       type(grid_warehouse_t),        intent(inout) :: grid_warehouse
       type(profile_warehouse_t),     intent(inout) :: profile_warehouse
       real(kind=dk), allocatable                   :: spectral_weight(:)

       ! do your calculations here

     end subroutine calculate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end module tuvx_spectral_weight_foo


The constructor function is reponsible for initializing new instances of your
spectral weight subclass.
First, you allocate the pointer to be returned as your new type
(``spectral_weight_foo_t`` in this example).
Then you initialize its data members.
If you just want to use the default initialization of the base class, you can
call the ``base_constructor()`` function as shown above.
You can alternatively initialize data members of the base class (``spectral_weight_parms(:)``)
directly in this function or add data members to your subclass and initialize them
here.

The ``calculate()`` function overrides the base-class ``calculate()`` function and will be
called when a user calls the ``calculate()`` type-bound procedure on an instance
of your new subclass.
You can access grid and profile data from the “warehouse” objects passed in as
function arguments, and any data in the base-class data members or in data members
you’ve added to your subclass to perform your calculations.
See the files in ``src/spectral_weights/`` for examples of how to access this data
in the ``calculate()`` function.


.. _dose-rate-add-to-build-scripts:

Add subclass module to build scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To include your new class in the build, edit the
``src/spectral_weights/CMakeLists.txt`` file and add your file name to the list
saved to ``SRC``. Files are listed in alphabetical order.

.. code-block:: cmake

   ################################################################################
   # Spectral weight source

   set(SRC notch_filter.F90
           gaussian_filter.F90
           eppley.F90
           par.F90
           exp_decay.F90
           foo.F90
           scup_mice.F90
           standard_human_erythema.F90
           UV_Index.F90
           phytoplankton_boucher.F90
           plant_damage.F90
           plant_damage_flint_caldwell.F90
           plant_damage_flint_caldwell_ext.F90
           )

   list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
   set(SPECTRAL_WGHT_SRC ${SRC} PARENT_SCOPE)

   ################################################################################


.. _dose-rate-add-to-factory:

Add subclass to factory
^^^^^^^^^^^^^^^^^^^^^^^

In order to use your new subclass, you will need to add it to the
``tuvx_spectral_weight_factory`` module in ``src/spectral_weight_factory.F90``.
First use-associate your new class at the module level:

.. code-block:: fortran

   use tuvx_spectral_weight_foo,          only : spectral_weight_foo_t


Then, inside the ``spectral_weight_builder()`` function, add these lines to the
``select case`` block:

.. code-block:: fortran

   case( 'foo' )
     new_spectral_weight => spectral_weight_foo_t( config, grid_warehouse,       &
                                                   profile_warehouse )


Now, when you add a spectral weight of type ``foo`` to the configuration data,
an instance of your new subclass will be created.



.. _dose-rate-create-unit-test:

Create unit test
^^^^^^^^^^^^^^^^

The last step to adding a spectral weight is to create a unit test.
This will ensure that your calculations are doing what you intended.
It will also serve as an example for how users can configure and use
your new subclass.

See :ref:`developer-add-test` for more details.

Quantum Yields
--------------

The base functionality for quantum yield calculations is described in
``src/quantum_yield.F90``. The base quantum yield type ``quantum_yield_t``
defined in this module includes a constructor that reads a set of
NetCDF files specified in the configuration JSON object as an
array of file paths with the key ``netcdf files`` if present, or
can set the value of the quantum yield to a constant when the
``constant value`` key is present and set to a real number.

Here is an example configuration for a quantum yield:

.. code-block:: JSON

   {
     "type": "base",
     "constant value": 1.0
   }


Data from each NetCDF file will be loaded into an element of the
``quantum_yield_parms`` data member. If a NetCDF variable named
``quantum_yield_parameters`` is present, it will be used to populate
the ``array(:,:)`` data member of the ``quantum_yield_parms_t`` object.
This first dimension of the array is wavelength, and will be interpolated
to the native TUV-x wavelength grid if this differs from the wavelength
grid in the netcdf file (named ``wavelength``). The second dimension
is used to accomodate multiple wavelength-resolved parameters for
quantum yield calculations.

If a NetCDF variable named ``temperature`` is present, it will be
used to populate the ``temperature(:)`` data member of the
``quantum_yield_parms_t`` object.

The calculation of quantum yields is done by calling the ``calculate()``
type-bound procedure on a ``quantum_yield_t`` object.

The base-class calculation of quantum yields returns the
wavelength-interpolated first parameter (``array(:,1)``) from the first
NetCDF file (``quantum_yield_parms(1)``) specified in the configuration
data.

Quantum yield subclasses are located in the ``src/quantum_yields/`` folder.
These each provide unique algorithms for calculating quantum yields.

Before adding a new quantum yield class, first check to make sure there
is not an existing class that can be configured to accomodate your
needs. If you find that an existing quantum yield subclass could be used
by moving one or more hard-coded parameters to the configuration data, this
is preferable to adding a new subclass.

If you determine that a new quantum yield subclass is needed, this can be
done in four steps:

- :ref:`qy-create-subclass`
- :ref:`qy-add-to-build-scripts`
- :ref:`qy-add-to-factory`
- :ref:`qy-create-unit-test`

.. _qy-create-subclass:

Create subclass module
^^^^^^^^^^^^^^^^^^^^^^

First, choose a unique name for your quantum yield calculation. Ideally,
this name will describe the algorithm, rather than the specific photolysis
reaction you are applying it to. However, many subclasses currently in TUV-x
are named for specific photolysis reactions. For this example, we will use
the name ``foo`` for our quantum yield algorithm.

**Pay special attention to naming of files, modules, types, and functions
in these instructions.**

Create a file to hold your new subclass module in ``src/quantum_yields/`` named
``foo.F90``. The general layout of the module will be (comments have been omitted
for this example, but should be included in an actual module):

.. code-block:: fortran

   ! Copyright (C) 2020 National Center for Atmospheric Research
   ! SPDX-License-Identifier: Apache-2.0
   !
   module tuvx_quantum_yield_foo

     use tuvx_quantum_yield,              only : quantum_yield_t

     implicit none
     private

     public :: quantum_yield_foo_t

     type, extends(quantum_yield_t) :: quantum_yield_foo_t
     contains
       procedure :: calculate
     end type quantum_yield_foo_t

     interface quantum_yield_foo_t
       module procedure constructor
     end interface

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function constructor( config, grid_warehouse, profile_warehouse )           &
         result( this )

       use musica_assert,                 only : assert_msg
       use musica_config,                 only : config_t
       use musica_string,                 only : string_t
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t
       use tuvx_quantum_yield,            only : base_constructor

       class(quantum_yield_t),    pointer       :: this
       type(config_t),            intent(inout) :: config
       type(grid_warehouse_t),    intent(inout) :: grid_warehouse
       type(profile_warehouse_t), intent(inout) :: profile_warehouse

       type(string_t) :: required_keys(1), optional_keys(1)

       ! This block of code ensures that the configuration keys are valid for
       ! your class. These can be modified to fit your needs. The first
       ! argument to assert_msg() should be a unique integer code for this error.
       required_keys(1) = "type"
       optional_keys(1) = "name"
       call assert_msg( 409635586,                                               &
                        config%validate( required_keys, optional_keys ),         &
                        "Bad configuration data format for "//                   &
                        "foo quantum yield." )

       allocate( quantum_yield_foo_t :: this )

       ! You can call the base_constructor function to load data from NetCDF
       ! files into the `quantum_yield_parms(:)` data member according to the
       ! standard base class logic. Alternatively, you can perform custom
       ! initialization of the subclass object here.
       call base_constructor( this, config, grid_warehouse, profile_warehouse )

     end function constructor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function calculate( this, grid_warehouse, profile_warehouse )               &
         result( quantum_yield )

       use musica_constants,              only : dk => musica_dk
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t

       class(quantum_yield_foo_t), intent(in)    :: this
       type(grid_warehouse_t),     intent(inout) :: grid_warehouse
       type(profile_warehouse_t),  intent(inout) :: profile_warehouse
       real(kind=dk), allocatable                :: quantum_yield(:,:)

       ! Do your calculations here

     end function calculate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end module tuvx_quantum_yield_foo


The constructor function is reponsible for initializing new instances of your
quantum yield subclass. First, you allocate the pointer to be returned as
your new type (``quantum_yield_foo_t`` in this example). Then you initialize
its data members. If you just want to use the default initialization of the
base class, you can call the ``base_constructor()`` function as shown above.
You can alternatively initialize data members of the base class
(``quantum_yield_parms(:)``) directly in this function or add data members
to your subclass and initialize them here (see
``src/quantum_yields/tint.F90`` for an example).

The ``calculate()`` function overrides the base-class ``calculate()`` function
and will be called when a user calls the ``calculate()`` type-bound procedure
on an instance of your new subclass.
You can access grid and profile data from the "warehouse" objects
passed in as function arguments, and any data in the base-class data members
or in data members you've added to your subclass to perform your calculations.
See the files in ``src/quantum_yields/`` for examples of how to access this
data in the ``calculate()`` function.

.. _qy-add-to-build-scripts:

Add subclass module to build scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To include your new class in the build, edit the ``src/quantum_yields/CMakeLists.txt``
file and add your file name to the list saved to ``SRC``. Files are listed in
alphabetical order.

.. code-block:: cmake
   :emphasize-lines: 12

   set(SRC acetone-ch3co_ch3.F90
        c2h5cho.F90
        ch2chcho.F90
        ch2o.F90
        ch3cho-ch3_hco.F90
        ch3coch2ch3-ch3co_ch2ch3.F90
        ch3cocho.F90
        clo-cl_o1d.F90
        clo-cl_o3p.F90
        clono2-clo_no2.F90
        clono2-cl_no3.F90
        foo.F90
        ho2-oh_o.F90
        mvk.F90
        no2_tint.F90
        no3_aq.F90
        o3-o2_o1d.F90
        o3-o2_o3p.F90
        tint.F90
        )

.. _qy-add-to-factory:

Add subclass to factory function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to use your new subclass, you will need to add it to the
``tuvx_quantum_yield_factory`` module in ``src/quantum_yield_factory.F90``.
First use-associate your new class at the module level:

.. code-block:: fortran

   use tuvx_quantum_yield_foo,            only : quantum_yield_foo_t

Then, inside the ``quantum_yield_builder()`` function, add these lines to the
``select case`` block:

.. code-block:: fortran

   case( 'foo' )
     quantum_yield => quantum_yield_foo_t( config, grid_warehouse,              &
                                           profile_warehouse )

Now, when you add a quantum yield of type ``foo`` to the configuration data,
an instance of your new subclass will be created.

.. _qy-create-unit-test:

Create unit test
^^^^^^^^^^^^^^^^

The last step to adding a quantum yield is to create a unit test. This will ensure
that your calculations are doing what you intended. It will also serve as an example
for how users can configure and use your new subclass.

See :ref:`developer-add-test` for more details.


Radiators
---------

Radiators are atmospheric constituents that affect the calculation of the
radiative field.
The configuration for a standard radiator is:

.. code-block:: JSON

   {
     "name": "foo",
     "type": "base",
     "cross section": "foo",
     "vertical profile": "foo",
     "vertical profile units": "molecule cm-3"
   }

A description of the components of the radiator configuration are
provided :ref:`here <configuration-radiators>`.

Most radiators can use the standard radiator configuration.
If a new algorithm for calculating the optical properties of
radiators is required, a new radiator subclass can be introduced
in four steps:

- :ref:`radiator-create-subclass`
- :ref:`radiator-add-to-build-scripts`
- :ref:`radiator-add-to-factory`
- :ref:`radiator-create-unit-test`

.. _radiator-create-subclass:

Create subclass module
^^^^^^^^^^^^^^^^^^^^^^

First, choose a unique name for your radiator algorithm.
Ideally, this name will describe the algorithm, rather than the specific
atmospheric constituent you are applying it to.
For this example, we will use the name ``foo`` for our radiator algorithm.

**Pay special attention to naming of files, modules, types, and functions
in these instructions.**

Create a file to hold your new subclass module in ``src/radiators/`` named
``foo.F90``.
The general layout of the module will be (comments have been omitted for this
example, but should be included in an actual module):

.. code-block:: fortran

   ! Copyright (C) 2020 National Center for Atmospheric Research
   ! SPDX-License-Identifier: Apache-2.0
   !
   module tuvx_radiator_foo

     use tuvx_radiator,                   only : radiator_t

     implicit none

     private
     public :: radiator_foo_t

     type, extends(radiator_t) :: radiator_foo_t
     contains
       procedure :: update_state
     end type radiator_foo_t

     interface radiator_foo_t
       module procedure :: constructor
     end interface radiator_foo_t

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     function constructor( config, grid_warehouse ) result( this )

       use musica_assert,                 only : assert_msg
       use musica_config,                 only : config_t
       use musica_string,                 only : string_t
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_radiator,                 only : base_constructor

       class(radiator_t),      pointer       :: this
       type(config_t),         intent(inout) :: config
       type(grid_warehouse_t), intent(inout) :: grid_warehouse

       type(string_t) :: required_keys(1), optional_keys(1)

       ! This block of code ensures that the configuration keys are valid for
       ! your class. These can be modified to fit your needs. The first
       ! argument to assert_msg() should be a unique integer code for this error.
       required_keys(1) = "type"
       optional_keys(1) = "name"
       call assert_msg( 302604745,                                               &
                        config%validate( required_keys, optional_keys ),         &
                        "Bad configuration data format for "//                   &
                        "foo radiator." )

       allocate( radiator_foo_t :: this )

       ! You can call the base_constructor function to load data data members
       ! with configuration data available from the standard radiator class.
       ! Alternatively, you can perform custom initialization of the subclass
       ! object here.
       call base_constructor( this, config, grid_warehouse )

     end function constructor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine update_state( this, grid_warehouse, profile_warehouse,           &
         cross_section_warehouse )

       use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
       use tuvx_grid_warehouse,           only : grid_warehouse_t
       use tuvx_profile_warehouse,        only : profile_warehouse_t

       class(radiator_foo_t),           intent(inout) :: this
       type(grid_warehouse_t),          intent(inout) :: grid_warehouse
       type(profile_warehouse_t),       intent(inout) :: profile_warehouse
       type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse

       ! Calculate optical properties (layer optical depth, layer single
       ! scattering albedo, and layer asymmetry factor) and load them into
       ! this%state_

     end subroutine update_state

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end module tuvx_radiator_foo


The ``constructor()`` function is responsible for initializing new instances of
your radiator subclass.
First, you allocate the pointer to be returned as your new type
(``radiator_foo_t`` in this example).
Then, you initialize its data members.
If you want to use the default initialization of the base class, you can
call the ``base_constructor()`` function as shown above.
You can alternatively initialize data members of the base class directly in
this function or add data members to your subclass and initialize them here.

The ``update_state()`` function overrides the base-class ``update_state()``
function and will be called when a user calls the ``update_state()`` type-bound
procedure on an instance of your new subclass.
You can access grid, profile, and cross section data from the "warehouse"
objects passed in as function arguments, and any data in the base-class data
members or in data members you've added to your subclass to perform your
calculations.
See the files in ``src/radiators/`` for examples of how to access this data
in the ``update_state()`` function.



.. _radiator-add-to-build-scripts:

Add subclass module to build scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To include your new class in the build, edit the
``src/radiators/CMakeLists.txt`` file and add your file name to the
list save to ``SRC``. Files are listed in alphabetical order.

.. code-block:: cmake

   ################################################################################
   # Radiator transfer source

   set(SRC aerosol.F90
           foo.F90
           )

   list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
   set(RADIATOR_SRC ${SRC} PARENT_SCOPE)

   ################################################################################


.. _radiator-add-to-factory:

Add subclass to factory
^^^^^^^^^^^^^^^^^^^^^^^

In order to use your new subclass, you will need to add it to the
``tuvx_radiator_factory`` module in ``src/radiator_factory.F90``.
First, use-associate your new class at the module level:

.. code-block:: fortran

   use tuvx_radiator_foo,                 only : radiator_foo_t


Then, inside the ``radiator_builder()`` function, add these lines to the
``select case`` block:

.. code-block:: fortran

   case( 'foo' )
     new_radiator => radiator_foo_t( config, grid_warehouse )


Now, when you add a radiator of type ``foo`` to the configuration data, an instance
of your new subclass will be created.



.. _radiator-create-unit-test:

Create unit test
^^^^^^^^^^^^^^^^

The last step to adding a radiator is to create a unit test.
This will ensure that your calculations are doing what you intended.
It will also serve as an example for how users can configure and use your new subclass.

See :ref:`developer-add-test` for more details.


.. _developer-add-test:

Test Creation
-------------

Unit tests are required for all new code contributions.
Source code for new unit tests should be added to the ``test/unit/`` folder
or one of its sub-folders depending on the module being tested.
Unit tests are typically Fortran programs that are linked to the ``tuv-x``
library and test the components of a single Fortran module in the ``src/``
tree.

An example of a  unit test for the fictitous ``foo`` module is shown below.

.. code-block:: fortran

   program test_foo

     implicit none

     call test_foo_t( )

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine test_foo_t( )
       ! Tests the foo_t type

       use musica_assert,              only : assert
       use tuvx_foo,                   only : foo_t

       type(foo_t) :: my_foo

       call assert( 501352581, my_foo%do_bar( ) .eq. 12.5 )
       call assert( 503258115, my_foo%do_baz( ) .eq. "qux" )

     end subroutine test_foo_t

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end program test_foo


The `musica_assert <https://ncar.github.io/musica-core/html/namespacemusica__assert.html>`_
module contains a number of functions that can be useful in
unit tests.

You will need to modify the ``CMakeLists.txt`` file in the
folder where you saved your test source code (for this example we assume the above
file is named ``test_foo.F90``) to include your new source in the build, and
your test in the test suite.
An updated ``CMakeLists.txt`` file for the ``foo`` test is shown below.


.. code-block:: cmake

   ################################################################################
   # Test utilities

   include(test_util)

   ################################################################################
   # Photo-decomp tests

   create_standard_test(NAME some_existing_test SOURCES test_bar.F90)
   create_standard_test(NAME foo SOURCES test_foo.F90)

   ################################################################################


The ``create_standard_test()`` CMake function adds your new executable to the build,
links it to the ``tuv-x`` library, and includes the test as well as a
memory check of your test to the testing suite.
The function is defined in ``cmake-modules/test_util.cmake``, but can generally used
as shown above.

If your test needs access to data files, you can place these in the ``test/data/``
folder.
By default, your test executable will be run in the build folder and can access
data files you place in this folder using a relative path: ``test/data/my_foo_data.txt``.
