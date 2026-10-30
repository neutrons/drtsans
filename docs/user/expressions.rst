.. _user.expressions:

Expression Syntax
=================

Some configuration values are mathematical expressions rather than plain numbers, for example a polarization
that varies with wavelength. The settings listed below are read with the same restricted syntax, provided by
:py:mod:`drtsans.expressions`. An invalid expression is rejected when it is read, before any data are reduced.

The ``expr`` fit constraint in the beam-center ``gaussian_centering_options`` is not covered by this page. It is
passed to the `lmfit <https://lmfit.github.io/lmfit-py/constraints.html>`_ fitting package and follows its syntax.

Where expressions are used
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 20 45

   * - Setting
     - Variable
     - Other names
   * - GPSANS attenuation coefficients file, ``formula =`` line
       (:ref:`user.corrections.attenuation`)
     - ``wavelength``, in Å
     - fitted parameters declared in the attenuator blocks
   * - GPSANS ``polarization`` settings ``polarizer.polarization``, ``polarizer.efficiency``,
       ``analyzer.polarizationZero`` and ``analyzer.polarizationPi`` (:ref:`polarization`)
     - ``x``, the wavelength in Å
     - none

The variable is optional. An expression that does not use it gives the same value at every wavelength.
Polarization settings may also be given as plain JSON numbers, for example ``0.95`` instead of ``"0.95"``.

Allowed syntax
--------------

- numbers, such as ``2``, ``0.95`` or ``1.5e-3``;
- the variable and other names listed in the table above;
- the arithmetic operators ``+``, ``-``, ``*``, ``/`` and ``**`` (power), unary ``+`` and ``-``, and parentheses;
- the functions ``abs``, ``acos``, ``asin``, ``atan``, ``cos``, ``cosh``, ``erf``, ``exp``, ``log`` (natural
  logarithm), ``log10``, ``sin``, ``sinh``, ``sqrt``, ``tan`` and ``tanh``, each called with positional
  arguments;
- the constant ``pi``.

Anything else is rejected, including undeclared names (for example ``wavelength`` in a polarization setting, which
must use ``x``), other functions, keyword arguments, attribute access such as ``np.exp``, indexing, comparisons,
the operators ``//`` and ``%``, and text in quotes. Expressions whose constant part has a complex value, such as
``sqrt(-1)``, or is infinite with undefined sign, such as ``1/0``, are also rejected, as are constant powers too
large or too small for a 64-bit floating-point number, such as ``10**400``. The error message states the reason.

When a division by zero, or a fractional power of a negative number, depends on the value of a variable or of a
fitted parameter, it cannot be detected when the expression is read. It gives an infinite or undefined value
when the expression is evaluated. The attenuation correction and the polarization correction then report an
error.

Examples
--------

Attenuation formula with fitted parameters ``A``, ``B`` and ``C``:

.. code-block:: text

    formula = A * exp(-B * wavelength) + C

Wavelength-dependent polarizer polarization in a reduction configuration:

.. code-block:: json

    "polarization": {
        "polarizer": {
            "polarization": "0.95 - 0.01*(x - 16)",
            "efficiency": 0.98
        }
    }
