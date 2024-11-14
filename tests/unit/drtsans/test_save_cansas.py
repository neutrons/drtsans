import os
import pytest
import tempfile

from drtsans.save_cansas import save_cansas_nx

from mantid.simpleapi import CreateSampleWorkspace, ConvertUnits, CreateWorkspace, LoadInstrument


def test_save_cansas_nx():
    ws = CreateSampleWorkspace("Histogram", NumBanks=1, BankPixelWidth=1)
    ws = ConvertUnits(ws, Target="MomentumTransfer")
    LoadInstrument(ws, False, InstrumentName="SANS2D")

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
        save_cansas_nx(ws, tf.name)

        assert os.path.exists(tf.name)
        assert os.path.getsize(tf.name) > 0

        os.remove(tf.name)


def test_save_cansas_nx_ValueError_uncommonBinBoundaries():
    ws = CreateWorkspace(
        DataX=[1, 2, 3],
        UnitX="MomentumTransfer",
        DataY=[4, 5, 6],
        VerticalAxisValues=[7, 8, 9],
        VerticalAxisUnit="MomentumTransfer",
        NSpec=3,
    )

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
        save_cansas_nx(ws, tf.name)

        assert os.path.exists(tf.name)
        assert os.path.getsize(tf.name) == 0

        os.remove(tf.name)


if __name__ == "__main__":
    pytest.main([__file__])
