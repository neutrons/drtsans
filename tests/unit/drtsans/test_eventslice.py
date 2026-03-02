import pytest

from drtsans.eventslice import resolve_slicing


def test_resolve_slicing():
    options = {"configuration": {"useTimeSlice": True, "useLogSlice": False}, "sample": {"runNumber": "12345"}}
    assert resolve_slicing(options) == (True, False)

    with pytest.raises(ValueError) as except_info:
        options["configuration"]["useLogSlice"] = True
        resolve_slicing(options)
    assert "Can't do both time and log slicing" in str(except_info.value)

    options = {"configuration": {"useTimeSlice": True, "useLogSlice": False}, "sample": {"runNumber": "1,2,3"}}
    with pytest.raises(ValueError) as except_info:
        resolve_slicing(options)
    assert "Can't do slicing on summed data sets" in str(except_info.value)
