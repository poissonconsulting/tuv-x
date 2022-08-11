.. TUV-x Quantum Yield

Quantum Yield
=============

A :f:type:`~tuvx_quantum_yield/tuvx_quantum_yield` object defines a calculator
used to compute the quantum yield. Various subclasses are derived from this type
for each reaction of interest to TUV-X.
The quantum yields are configured at runtime. See 
:ref:`configuration-quantum-yields` for more information.

Quantum Yield
^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield

Quantum Yield Factory
^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_factory

:math:`CH_3COCH_3+hv \rightarrow CH_3CO+CH_3`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch3coch3_ch3co_ch3
  
:math:`C_2H_5CHO+hv \rightarrow C_2H_5 + HCO`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_c2h5cho_c2h5_hco

:math:`CH_2CHCHO+hv \rightarrow Products`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch2chcho

:math:`CH_2O+hv \rightarrow H_2 + CO`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch2o_h2_co

:math:`CH_3CHO+hv \rightarrow CH_3+HCO`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch3cho_ch3_hco

:math:`CH_3COCH_2CH_3+hv \rightarrow CH_3CO+CH_2CH_3`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch3coch2ch3

:math:`CH_3COCHO+hv \rightarrow CH_3CO+HCO`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ch3cocho_ch3co_hco

:math:`ClO+hv \rightarrow Cl + O(^1D)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_clo_cl_o1d

:math:`ClO+hv \rightarrow Cl + O(^3P)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_clo_cl_o3p

:math:`ClONO_2+hv \rightarrow Cl + NO_3`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_clono2_cl_no3

:math:`ClONO_2+hv \rightarrow ClO + NO_2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_clono2_clo_no2

:math:`HO_2+hv \rightarrow OH + H`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_ho2_oh_o

:math:`MVK+hv \rightarrow Products`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_mvk

:math:`NO2` Temperature Interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_no2_tint

:math:`{NO_{3}}^-_{(aq)}+hv \rightarrow {NO_2}_{(aq)}+O^-` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_no3m_aq

:math:`O_3+hv \rightarrow O_2 + O(^1D)` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_o3_o2_o1d

:math:`O_3+hv \rightarrow O_2 + O(^3P)` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_o3_o2_o3p

Temperature Interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^
.. f:automodule:: tuvx_quantum_yield_tint