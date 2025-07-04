.. release_notes

=============
Release Notes
=============
..
  Use the following template to add a new release note.

  <Next Release>
  --------------
  (date of release, format YYYY-MM-DD)

  **Of interest to the User**:
  - PR #XYZ: one-liner description

  **Of interest to the Developer:**
  - PR #XYZ: one-liner description
..

..
  1.17.0
  --------------
  XXXX-XX-XX

  **Of interest to the User**:
  - PR #XYZ: one-liner description

  **Of interest to the Developer:**
  - PR #XYZ: one-liner description
..

1.16.0
--------------
2025-07-08

**Of interest to the User**:
- PR #1034: Now able to select reference detector for stitching

**Of interest to the Developer:**
- PR #1032: Add splitter of polarized runs

1.15.0
-------
2025-06-10

**Of interest to the User**:
- PR #1030: Introduces command grasp_cg2 to export gpsans reduction runs to GRASP format
- PR #1029: Calculation of the I(Q2D) error is simplified for the incoherent inelastic correction
- PR #1025: Includes v1.14.0 release notes
- PR #1016: User documentation for inelastic incoherent correction and elastic reference normalization

**Of interest to the Developer:**
- PR #1031: Pins muparser to get around an issue in mantid

1.14.0
-------
2025-04-29

**Of interest to the User**:
- PR #1018: Add initial user documentation on polarization
- PR #1015: Upgrades to Mantid 6.12 and numpy 2

**Of interest to the Developer:**
- PR #1023: Github actions - conda build versiongit information from tags
- PR #1022: Generator of simulated runs for testing
- PR #1021: stub for the polarization module
- PR #1019: Fixes an edge case bug in creating the sample logs periodic index log where the number of entries is less than the number of times

1.13.0
------
2025-04-01

**Of interest to the User**:

- PR #1005: Annular binning output is now annotated with phi rather than Q and plotted with a linear x axis

**Of interest to the Developer:**

1.12.0
------
2025-03-17

Notable changes in this release include migrating the repository from GitLab to GitHub,
fixing a bug in the auto wedge finding function, and saving I(Qx, Qy) in NXCANSAS format.
Additionally, several improvements were made to the JSON schema, documentation,
and configuration options.

**Of interest to the User**:

- PR #1010: Fixes rendering of README.rst and adds some badges
- PR #1008: Add support for annular binning with no wavelength binning.
- PR #999: Adds support for flexible time slicing.
- PR #997: Fix bug causing symmetric auto wedge finding to fail. Add mirrored wedge to auto wedge fit function plot.
- PR #998: remove the TOF offset that is done by the data acquisition system
- PR #996: Iq.dat files are now written even if they fail the assumption check
- PR #994: Remove unused module `drtsans/tof/eqsans/reduce.py`
- PR #993: Skip slices with too high transmission error when using time sliced sample transmission run
- PR #325: Migrates repository from GitLab to GitHub
- MR #1185: Rename output directory for coherent and incoherent-inelastic corrections as "info"
- MR #1184: Document wedge binning
- MR #1183: Added a script generate_report to generate a summary report from an hdf5 log file
- MR #1180: User documentation on time and log slicing
- MR #1179: Add validation for timeslice parameters, period must be integer multiple of interval
- MR #1177: Correction output now goes to `outputDir/debug/elastic_norm` and `outputDir/debug/incoherent_inelastic`
- MR #1176: The directory /SNS/EQSANS/IPTS-XXXX/nexus/ has a priority in file search algorithm
- MR #1175: Input parameters JSON field `sample::loadOptions::additionalProperties` now accepts a boolean.
- MR #1169: I(Qx, Qy) is saved in NXCANSAS format for biosans, gpsans, and eqsans workflows
- MR #1168: Document scaling the detector panels in the user's guide
- MR #1167: Allow separate configurations for inelastic incoherence correction per frame in frame skipping mode
- MR #1166: Added option to _fitSpectrum to auto-find one wedge and mirror it
- MR #1162: When reducing `gpsans` data with `direct_beam` scaling, the `direct_beam_scaling` parameter is now logged during
  the reduction process and stored in the output Nexus file at `reduction_information/special_parameters/direct_beam_scaling/value`.
- MR #1161: Add a parameters removeAlgorithmHistory to write less data and speed up I/O during reduction
- MR #1160: Expose pixel detector rescalings to the instrument API's
- MR #1159: Separate configuration for elastic normalization and inelastic incoherence correction

**Of interest to the Developer:**

- MR #1171: update environment.yml to replace build in favor of python-build
- MR #1169: creates src/drtsans/save_cansas.py to define CANSAS format file handling methods
- MR #1165: update Mantid version to 6.11
- MR #1158: added options in the JSON schemae to rescale whole detector panels

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
