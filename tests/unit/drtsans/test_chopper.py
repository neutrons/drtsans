import pytest
from numpy.testing import assert_almost_equal

from drtsans.chopper import DiskChopper, DiskChopperConfiguration, DiskChopperConfigurationParsingError
from drtsans.frame_mode import FrameMode


class TestDiskChopper:
    ch = DiskChopper(1.0, 45, 60, 2000, 850)

    def test_pulse_width(self):
        assert self.ch.pulse_width == DiskChopper._pulse_width
        self.ch.pulse_width = 0
        assert self.ch.pulse_width == 0
        assert self.ch.pulse_width != DiskChopper._pulse_width
        self.ch.pulse_width = DiskChopper._pulse_width  # restore state

    def test_cutoff_wl(self):
        assert self.ch.cutoff_wl == DiskChopper._cutoff_wl
        self.ch.cutoff_wl = 0
        assert self.ch.cutoff_wl == 0
        assert self.ch.cutoff_wl != DiskChopper._cutoff_wl
        self.ch.cutoff_wl = DiskChopper._cutoff_wl  # restore state

    def test_phase(self):
        assert self.ch.phase == 1150

    def test_period(self):
        assert_almost_equal(self.ch.period, 16666, decimal=0)

    def test_transmission_duration(self):
        assert_almost_equal(self.ch.transmission_duration, 2083, decimal=0)

    def test_opening_phase(self):
        assert_almost_equal(self.ch.opening_phase, 108, decimal=0)

    def test_closing_phase(self):
        assert_almost_equal(self.ch.closing_phase, 2191, decimal=0)

    def test_rewind(self):
        assert_almost_equal(self.ch.rewind, 108, decimal=0)
        self.ch.offset += 109
        assert_almost_equal(self.ch.rewind, -1, decimal=0)
        self.ch.offset -= 109  # restore state

    def test_wavelength(self):
        assert_almost_equal(self.ch.wavelength(1200), 4.7, decimal=1)
        assert_almost_equal(self.ch.wavelength(1200, pulsed=True), 4.3, decimal=1)

    def test_tof(self):
        assert_almost_equal(self.ch.tof(4.747), 1200, decimal=0)
        assert_almost_equal(self.ch.tof(4.399, pulsed=True), 1200, decimal=0)

    def test_transmission_bands(self):
        wb = self.ch.transmission_bands()
        assert len(wb) == 1
        assert_almost_equal((wb[0].min, wb[0].max), (0.42, 8.67), decimal=2)
        ch = DiskChopper(1.0, 30, 240, 0, 0)
        wb = ch.transmission_bands()
        assert len(wb) == 3
        assert_almost_equal((wb[0].min, wb[0].max), (0, 0.687), decimal=2)


class TestDiskChopperConfiguration:
    """Unit tests for DiskChopperConfiguration class."""

    def _json_data(self):
        """Helper method to provide sample JSON data for testing."""
        return [
            {
                "daystamp": 20240101,
                "n_choppers": 4,
                "aperture": [45, 45, 45, 45],
                "to_source": [1.0, 2.0, 3.0, 4.0],
                "offsets": {"not_skip": [0, 0, 0, 0], "skip": [10, 10, 10, 10]},
            }
        ]

    def test_initialization(self):
        """Test basic initialization of DiskChopperConfiguration."""
        config = DiskChopperConfiguration(
            n_choppers=4,
            aperture=[45, 45, 45, 45],
            to_source=[1.0, 2.0, 3.0, 4.0],
            offsets={FrameMode.skip: [0, 0, 0, 0], FrameMode.not_skip: [10, 10, 10, 10]},
        )
        assert config.n_choppers == 4
        assert config.aperture == [45, 45, 45, 45]
        assert config.to_source == [1.0, 2.0, 3.0, 4.0]

    def test_from_json(self):
        """Test loading configuration from JSON string."""
        json_data = [
            {
                "daystamp": 20250101,
                "n_choppers": 4,
                "aperture": [45, 45, 45, 45],
                "to_source": [1.0, 2.0, 3.0, 4.0],
                "offsets": {"not_skip": [0, 0, 0, 0], "skip": [10, 10, 10, 10]},
            },
            {
                "daystamp": 20260101,
                "n_choppers": 6,
                "aperture": [45, 30, 45, 30, 45, 30],
                "to_source": [1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
                "offsets": {"not_skip": [0, 0, 0, 0, 0, 0], "skip": [15, 15, 15, 15, 15, 15]},
            },
        ]

        # Test loading pre-2026 configuration
        config = DiskChopperConfiguration.from_json(json_data, 20250601)
        assert config.n_choppers == 4
        assert config.aperture == [45, 45, 45, 45]
        assert config.to_source == [1.0, 2.0, 3.0, 4.0]

        # Test loading post-2026 configuration
        config = DiskChopperConfiguration.from_json(json_data, 20260201)
        assert config.n_choppers == 6
        assert config.aperture == [45, 30, 45, 30, 45, 30]
        assert config.to_source == [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

    def test_from_json_invalid_date(self):
        """Test from_json raises error for date before any configuration."""
        json_data = self._json_data()
        with pytest.raises(DiskChopperConfigurationParsingError):
            DiskChopperConfiguration.from_json(json_data, 20230101)

    @pytest.mark.parametrize("key_str", ["aperture", "to_source"])
    def test_from_json_wrong_length(self, key_str):
        """Test from_json raises error for mismatched list lengths."""
        json_data = self._json_data()
        json_data[0][key_str] = [45, 45, 45]  # Incorrect length
        with pytest.raises(DiskChopperConfigurationParsingError, match="list length"):
            DiskChopperConfiguration.from_json(json_data, 20250601)

    @pytest.mark.parametrize("key_str", ["aperture", "to_source"])
    def test_from_json_non_float(self, key_str):
        """Test from_json raises error for non-float values."""
        json_data = self._json_data()
        json_data[0][key_str] = [1.0, 2.0, 3.0, "bad data"]
        with pytest.raises(DiskChopperConfigurationParsingError, match="Invalid"):
            DiskChopperConfiguration.from_json(json_data, 20250601)

    @pytest.mark.parametrize(
        "key_str,value,expected_error,expected_msg",
        [
            ("skip", [1.0, 2.0, 3.0], DiskChopperConfigurationParsingError, "length"),
            ("bogus_mode", [1.0, 2.0, 3.0, 4.0], DiskChopperConfigurationParsingError, "Invalid frame mode"),
            ("not_skip", [1.0, 2.0, 3.0, "bad_data"], DiskChopperConfigurationParsingError, "Invalid offsets value"),
        ],
    )
    def test_from_json_offsets(self, key_str, value, expected_error, expected_msg):
        json_data = self._json_data()
        json_data[0]["offsets"][key_str] = value
        with pytest.raises(expected_error, match=expected_msg):
            DiskChopperConfiguration.from_json(json_data, 20250601)


if __name__ == "__main__":
    pytest.main([__file__])
