====================
Reduction Parameters
====================

.. contents::

Default Values
--------------

.. exec::
    from drtsans.redparams import _instrument_json_generator
    # render the default values of the reduction parameter for each instrument suitable to restructured text
    docs = [default_json.to_rest() for _, default_json in _instrument_json_generator()]
    print(r'{}'.format(''.join(docs)))  # we require a raw string

Descriptions
------------

.. exec::
    from drtsans.redparams import _instrument_json_generator
    # render the descriptions of the reduction parameters for each instrument as a table
    docs = ''
    for name, default_json in _instrument_json_generator(field='description'):
        docs += default_json.to_rest_table(name)
    print(r'{}'.format(docs))  # we require a raw string
