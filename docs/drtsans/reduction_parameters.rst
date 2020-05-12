====================
Reduction Parameters
====================

.. contents::

.. exec::
    from drtsans.redparms import _instrument_json_generator
    # render the default values of the reduction parameter for each instrument suitable to restructured text
    docs = [default_json.to_rest() for _, default_json in _instrument_json_generator()]
    print(r'{}'.format(''.join(docs)))  # we require a raw string
