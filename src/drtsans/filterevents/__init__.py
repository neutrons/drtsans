"""
Event filtering strategies for time, log, and spin-based slicing.

This package provides a strategy pattern implementation for filtering event data
in various ways, including time-based intervals, sample log values, and polarization
(spin) states.

The main entry points are:
- :class:`FilterStrategy`: Abstract base class for all filtering strategies
- :func:`create_filter_strategy`: Factory function to create appropriate filter instances
- :func:`resolve_slicing`: Determine which filtering strategy to use from reduction config
"""

from drtsans.filterevents.basefilter import FilterStrategy, create_filter_strategy, resolve_slicing

__all__ = ["FilterStrategy", "create_filter_strategy", "resolve_slicing"]
