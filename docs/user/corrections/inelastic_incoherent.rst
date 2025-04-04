.. _user.corrections.inelastic_incoherent:

Inelastic incoherent correction
===============================

Background
----------

In SANS, the focus is on the structural information obtained from elastic scattering, and
inelastic scattering contributions are typically assumed to be negligible, however, this assumption
does not always hold. For example, in samples containing hydrogen, which has a high
incoherent scattering cross section, there can be a significant and wavelength-dependent
contribution from inelastic scattering [Do2013]_.

Due to the wavelength dependence of inelastic scattering effects, correcting for them is of
particular concern for time-of-flight (TOF) pulsed neutron instruments like EQ-SANS, which use a
broad spectrum of neutron wavelengths. The data reduction for EQ-SANS, therefore, includes optional
wavelength dependent corrections for inelastic incoherent scattering.

Elastic reference normalization procedure
-----------------------------------------

The following steps describe the elastic reference normalization available in `drtsans`:

#. Get :math:`I(q,\lambda_i)` of the given elastic scattering reference run, making sure q-bins are
   the same for all :math:`\lambda` bins. This will be denoted :math:`I_{\text{elastic}}`.

#. Determine :math:`q_{\min}` and :math:`q_{\max}` that exist in all :math:`I(q,\lambda_i)` for the
   fitting (minimization) process.

#. Find scale factor, :math:`K(\lambda_i)`, that minimizes following

   .. math::
      \sum_{q_k=q_{\min}}^{q_{\max}} \lvert I_{\text{elastic}}(q_k,\lambda_{ref})-K(\lambda_i)
      I_{\text{elastic}}(q_k,\lambda_i) \rvert^2

   Here, :math:`\lambda_{ref}` is the shortest wavelength bin of :math:`\lambda_{i}`.
   The resulting K, with the subscript "elastic" is omitted for convenience, is

   .. math::
      K(\lambda_i) = \frac{\sum_{q_k} I(q_k, \lambda_{ref})
      I(q_k, \lambda_i)}{\sum_{q_k} \left( I(q_k, \lambda_i) \right)^2}

#. Update the data (sample data) by

   .. math::
      I'(q,\lambda_i) = K(\lambda_i) I(q,\lambda_i)

#. Update the error.

   .. math::
      \delta I'(q, \lambda_i)^2 = \delta K(\lambda_i)^2 I(q, \lambda_i)^2 +
       K(\lambda_i)^2 \delta I(q, \lambda_i)^2

   The above scaling is applied for :math:`I(q,\lambda_i)` except for the reference wavelength,
   :math:`I(q, \lambda_{ref})`.


Inelastic incoherent correction procedure
-----------------------------------------

The following steps describe the inelastic incoherent correction available in `drtsans`:

#. Get :math:`I(q,\lambda_i)` of the sample run, making sure q-bins are same for all :math:`\lambda`
   bins.

#. Determine :math:`q_{\min}` and :math:`q_{\max}` that exist in all :math:`I(q,\lambda_i)` for the
   fitting (minimization) process.

#. Calculate the inelastic incoherent correction factor, :math:`b(\lambda)`. Here,
   :math:`\lambda_{ref}` is the shortest wavelength bin and :math:`N` is the number of :math:`q`
   points between :math:`q_{\min}` and :math:`q_{\max}` inclusive.

   #. If ``"incohfit_intensityweighted"`` is ``True``:

      .. math::
         b(\lambda_i) = -\frac{1}{N \sum_{q_k=q_{\min}}^{q_{\max}} \frac{1}{I(q_k,\lambda_{ref})}}
         \sum_{q_k=q_{\min}}^{q_{\max}} \frac{I(q_k,\lambda_{ref})-I(q_k,\lambda_i)}{I(q_k,\lambda_{ref})}

   #. If ``"incohfit_intensityweighted"`` is ``False``:

      .. math::
         b(\lambda_i)=-{\frac{1}{N}}{\sum_{q_k=q_{\min}}^{q_{\max}} (I(q_k,\lambda_{ref})-( I(q_k,\lambda_i) )}

#. If JSON parameter ``"selectMinIncoh"`` is true, find :math:`\lambda_i` that has smallest :math:`b`,
   and choose that :math:`\lambda_i` as the new :math:`\lambda_{ref}`. :math:`b(\lambda)` is then
   recalculated using the new :math:`\lambda_{ref}`, after which all :math:`b(\lambda)` should be
   greater than zero.

#. Update data for :math:`\lambda_i` except :math:`\lambda_{ref}`:

   .. math::
      I'(q,\lambda_i) = I(q,\lambda_i) - b(\lambda_i)

#. Update errors for :math:`\lambda_i` except :math:`\lambda_{ref}`. In the following, we use
   :math:`I_k^{\lambda_i} = I(q_k,\lambda_i)`.

   If :math:`q_k` is within :math:`q_1...q_N` defined above use:

   .. math::
      \left( {\delta I'}_k^{\lambda_i} \right)^2 = \left( \delta I_k^{\lambda_i} \right)^2
      \left(1-\frac{2}{N}\right) + \frac{1}{N^2} \sum_{k=1}^{k=N} \left[\left(\delta I_k^{\lambda_{ref}}\right)^2 +
      \left(\delta I_k^{\lambda_i}\right)^2
      \right]

   If :math:`q_k` is outside :math:`q_1...q_N`, then use:

   .. math::
      \left( {\delta I'}_k^{\lambda_i} \right)^2 = \left( \delta I_k^{\lambda_i} \right)^2
      + \frac{1}{N^2} \sum_{k=1}^{k=N} \left[\left(\delta I_k^{\lambda_{ref}}\right)^2 +
      \left(\delta I_k^{\lambda_i}\right)^2
      \right]




Parameters
----------

.. note::
   The following parameters are only available for EQ-SANS data reduction.


.. list-table::
   :widths: 25 65 10
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``"fitInelasticIncoh"``
     - If ``"true"``, inelastic incoherent correction will be applied.
     - ``false``
   * - ``"incohfit_intensityweighted"``
     - If ``"true"``, the intensity weighted method is used in the inelastic incoherent correction.
     - ``false``
   * - ``"incohfit_qmin"``
     - :math:`q_{\min}` for the inelastic incoherent correction. If ``null``, the minimum valid
       :math:`q_{\min}` will be used.
     - ``null``
   * - ``"incohfit_qmax"``
     - :math:`q_{\max}` for the inelastic incoherent correction. If ``null``, the maximum valid
       :math:`q_{\max}` will be used.
     - ``null``
   * - ``"incohfit_factor"``
     -
     - ``null``
   * - ``"outputWavelengthDependentProfile"``
     - If ``"true"``, output intensity profiles for each wavelength before and after elastic
       reference normalization and inelastic incoherent correction.
     - ``false``
   * - ``"elasticReference"``
     - Elastic reference run. If empty, the elastic reference normalization will be skipped.
     -
   * - ``"elasticReferenceBkgd"``
     - Background run for the elastic reference run.
     -
   * - ``"selectMinIncoh"``
     - If ``"true"``, use the smallest wavelength as reference wavelength.
     - ``false``


Example
-------

These are the relevant parameters in the JSON schema with their default values (the corrections are
disabled by default).

.. code-block:: json

    {
        "fitInelasticIncoh": false,
        "incohfit_intensityweighted": false,
        "incohfit_qmin": null,
        "incohfit_qmax": null,
        "incohfit_factor": null,
        "outputWavelengthDependentProfile": false,
        "elasticReference": {
          "runNumber": null,
          "thickness": 1.0,
          "transmission": {
            "runNumber": null,
            "value": null
          }
        },
        "elasticReferenceBkgd": {
          "runNumber": null,
          "transmission": {
            "runNumber": null,
            "value": null
          }
        },
        "selectMinIncoh": false
    }

References
----------

.. [Do2013]	C. Do, W. T. Heller, C. Stanley, F. X. Gallmeier, M. Doucet, and G. S. Smith,
   “Understanding inelastically scattered neutrons from water on a time-of-flight small-angle neutron
   scattering (SANS) instrument,” Nucl. Instruments Methods Phys. Res. Sect. A Accel. Spectrometers,
   Detect. Assoc. Equip. 737 42–46 (2013),
   `doi:10.1016/j.nima.2013.11.030 <https://doi.org/10.1016/j.nima.2013.11.030>`_.
