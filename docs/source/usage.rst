.. Usage information for TUV-x

Usage
=====

General
-------

Stand-alone TUV-x can be run from the command-line as:

.. code-block:: bash

   ./tuv-x configuration_file.json

The ``configuration_file.json`` contains the TUV-x configuration information described in
:doc:`configuration`.

TUV-x will generate a number of output files in a folder named ``OUTPUT/``. These fall
into several general catagories summaraized here.

+--------------------+-------------------------------------------------------+
| File name pattern  | Description                                           |
+====================+=======================================================+
| ``*->*.qyld.new``  | Quantum yield for a photolysis reaction               |
+--------------------+-------------------------------------------------------+
| ``*->*.xsect.new`` | Cross section for a photolyzing species               |
+--------------------+-------------------------------------------------------+
| ``*->*.xsqy.new``  | Product of cross-section and quantum yield for a      |
|                    | photolysis reaction                                   |
+--------------------+-------------------------------------------------------+
