.. release_notes

=============
Release Notes
=============

<Next Release>
--------------
(date of release)

**Of interest to the User**:

- MR #XYZ: one-liner description
- MR 1165: update Mantid version to 6.11
- MR 1162: When reducing `gpsans` data with `direct_beam` scaling, the `direct_beam_scaling` parameter is now logged during
  the reduction process and stored in the output Nexus file at `reduction_information/special_parameters/direct_beam_scaling/value`.
- MR #7718: Add a parameters  removeAlgorithmHistory to write less data and speed up I/O during reduction

**Of interest to the Developer:**

- MR #XYZ: one-liner description

?.??.?
------
????-??-??

**Of interest to the User**:

- MR #1159:  Separate configuration for elastic normalization and inelastic incoherence correction

1.11.0
------
2024-04-05

**Of interest to the User**:

- Pixel calibration and reduction of data collected at the BIOSANS midrange detector.
- Python notebooks for pixel calibration of the BIOSANS midrange detector.
- Overlap wedge stitching parameters for the BIOSANS panel detectors.
- The script to calculate pixel sensitivities for the monochromatic instruments has been splited into one script
  for BIOSANS and one script for GPSANS.
- Ability to split and reduce data according to a periodic process-variable. For instance, if the sample
  was subjected to a sinusoidal external magnetic field with a period `T` of 10 seconds,
  User can split the run into 0.1 seconds segments and add segments separated by `T` seconds,
  yielding 100 different `I(Q)` profiles.
  Each profile is associated to a particular value of the external magnetic field.
- Fix a defect preventing plotting after importing `drtsans` in Mantid's workbench.
- Option to export the wavelength-dependent `I(Q)` profiles before and after applying the incoherent-inelastic and
  the coherent elastic corrections.
  Files are exported to different subdirectories when time or log slicing is requested, and also when in
  skipped-frame mode (EQSANS only).
- Standard reduction scripts for each instrument available and described in User's online documentation.
- For EQSANS, default files have been removed for dark run, mask, and sensitivity corrections.
  Associated parameters in the JSON configuration must provide paths to these files or be left emtpy if
  the corrections are not desired.
- It is highly likely that in the future, new parameters will be introduced into the instrument JSON schema.
  Still, User wants to run the current `drtsans` with a JSON file compliant with these future schema.
  User now has the option to set parameter `permissible=True` in API reduction functions
  `validate_reduction_parameters` and `update_reduction_parameters`.
  This setting permits using current `drtsans` to reduce "futuristic" JSON files.

**Of interest to the Developer:**

- Move from `versioneer` to `versioningit`.
- Functionality to generate a fake set of TOF events to generate typical intensity patterns
  (e.g. concentric rings, a flood pattern) in the detector panels.
  Useful for testing without using real event Nexus files.
- Store large data files into dedicated Git LFS repository `drtsans-data` for purposes of integration testing.
  Developer's docs explain how to use this repository and how to run tests that make use of this dataset.
- Exclude from GitLab CI testing those tests requiring files stored in the `/SNS` or `/HFIR` file systems.
  These tests may be run manually from a machine where these file systems are reachable.
  These tests can be selected with pytest marker `eqsans_mount`. Developer's documentation explains how-to.
- Results from the elastic-correction of `I(Q)` are now reused for the elastic-correction of `I(Qx, Qy)`.
- Python wheel/Conda package process modernization and documentation of the process

1.10.2
------
2023-10-19

**Of interest to the User**:

- MR 1060: New option in configuration file to split/sum with periodic log profiles in EQSANS
- MR 1059: Make stitch_profiles backwards compatible by allowing two formats for parameter overlaps


**Of interest to the Developer:**

- Mr 1085: use mamba instead of conda to build the docs
- MR 1074: EQSANS integration test with simulated TOF scattering
- MR 1058: Add the test repository as a git submodule
