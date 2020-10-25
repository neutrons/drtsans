import os


def clear_buffer():
    # Body of Reduction - DO NOT CHANGE....
    from mantid.simpleapi import mtd  # noqa: E401
    mtd.clear()


def reduce_biosans_nexus(IPTS_Number, EXPERIMENT_NUMBER, sample_names,
                         samples, backgrounds, samples_trans, backgrounds_trans, empty_trans,
                         beam_center,
                         sample_thick, sample_identifier, overWrite,
                         dark_mfname, dark_wfname,
                         base_output_directory, sens_mfname, sens_wfname,
                         scalefac, scaling_beam_radius, Lin1DQbins_Main, Lin1DQbins_Wing,
                         Lin2DQxy_Main, Lin2DQxy_Wing, Plot_type, Plot_binning,
                         LogQbinsPerDecade_Main, LogQbinsPerDecade_Wing, q_range_main, q_range_wing,
                         OL_range,    flexible_pixelsizes, wedge_min_angles, wedge_max_angles, Qmin_TDW, Qmax_TDW,
                         Qdelta_TDW, PeakWidth_TDW, AziDelta_TDW, BkgWidth_TDW, MinSigtoNoise_TDW,
                         q_range_main_wedge0, q_range_wing_wedge0, OL_range_wedge0,
                         q_range_main_wedge1, q_range_wing_wedge1, OL_range_wedge1, refreshCycle
                         ):
    import numpy as np  # noqa: E401
    import warnings  # noqa: E401
    warnings.filterwarnings('ignore')

    # get_ipython().run_line_magic('matplotlib', 'inline')
    import time  # noqa: E401
    from drtsans.mono.spice_data import map_to_nexus  # noqa: E401

    # Convert SPICE scan-pt tuple to NeXus files
    CG3 = 'CG3'
    samples = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, samples, nexus_dir=None)
    samples_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, samples_trans, nexus_dir=None)
    backgrounds = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, backgrounds, nexus_dir=None)
    backgrounds_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, backgrounds_trans, nexus_dir=None)
    beam_center = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [beam_center], nexus_dir=None)[0]
    empty_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [empty_trans], nexus_dir=None)[0]
    dark_mfname = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [dark_mfname], nexus_dir=None)[0]
    dark_wfname = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [dark_wfname], nexus_dir=None)[0]

    from drtsans.mono.biosans import (load_all_files, reduce_single_configuration, plot_reduction_output,
                                      reduction_parameters, update_reduction_parameters)  # noqa: E401

    # reduction parameters common to all the reduction runs to be carried out in this notebook
    common_configuration = {
        "iptsNumber": IPTS_Number,
        "beamCenter": {"runNumber": beam_center},
        "emptyTransmission": {"runNumber": empty_trans},
        "configuration": {
            "outputDir": base_output_directory,
            "darkMainFileName": dark_mfname,
            "darkWingFileName": dark_wfname,
            "sensitivityMainFileName": sens_mfname,
            "sensitivityWingFileName": sens_wfname,
            "defaultMask": [{'Pixel': '1-18,239-256'}, {'Bank': '18-24,42-48'}, {'Bank': '49', 'Tube': '1'}],
            'StandardAbsoluteScale': scalefac,
            "DBScalingBeamRadius": scaling_beam_radius,
            "mmRadiusForTransmission": "",
            "absoluteScaleMethod": "standard",
            "numMainQBins": Lin1DQbins_Main,
            "numWingQBins": Lin1DQbins_Wing,
            "numMainQxQyBins": Lin2DQxy_Main,
            "numWingQxQyBins": Lin2DQxy_Wing,
            "1DQbinType": Plot_type,
            "QbinType": Plot_binning,
            "LogQBinsPerDecadeMain": LogQbinsPerDecade_Main,
            "LogQBinsPerDecadeWing": LogQbinsPerDecade_Wing,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "sampleApertureSize": 14,
            "QminMain": q_range_main[0],
            "QmaxMain": q_range_main[1],
            "QminWing": q_range_wing[0],
            "QmaxWing": q_range_wing[1],
            "overlapStitchQmin": OL_range[0],
            "overlapStitchQmax": OL_range[1],
            "usePixelCalibration": flexible_pixelsizes,
            "useTimeSlice": False,
            "timeSliceInterval": 60,
            "WedgeMinAngles": wedge_min_angles,
            "WedgeMaxAngles": wedge_max_angles,
            "autoWedgeQmin": Qmin_TDW,
            "autoWedgeQmax": Qmax_TDW,
            "autoWedgeQdelta": Qdelta_TDW,
            "autoWedgePeakWidth": PeakWidth_TDW,
            "autoWedgeAzimuthalDelta": AziDelta_TDW,
            "autoWedgeBackgroundWidth": BkgWidth_TDW,
            "autoWedgeSignalToNoiseMin": MinSigtoNoise_TDW,
            "wedge1QminMain": q_range_main_wedge0[0],
            "wedge1QmaxMain": q_range_main_wedge0[1],
            "wedge1QminWing": q_range_wing_wedge0[0],
            "wedge1QmaxWing": q_range_wing_wedge0[1],
            "wedge1overlapStitchQmin": OL_range_wedge0[0],
            "wedge1overlapStitchQmax": OL_range_wedge0[1],
            "wedge2QminMain": q_range_main_wedge1[0],
            "wedge2QmaxMain": q_range_main_wedge1[1],
            "wedge2QminWing": q_range_wing_wedge1[0],
            "wedge2QmaxWing": q_range_wing_wedge1[1],
            "wedge2overlapStitchQmin": OL_range_wedge1[0],
            "wedge2overlapStitchQmax": OL_range_wedge1[1],
        }
    }

    common_configuration_full = reduction_parameters(common_configuration, 'BIOSANS', validate=False)
    # pretty_print(common_configuration_full)

    if len(backgrounds) == 1 and len(samples) > len(backgrounds):
        backgrounds = backgrounds * len(samples)
    if len(backgrounds_trans) == 1 and len(samples_trans) > len(backgrounds_trans):
        backgrounds_trans = backgrounds_trans * len(samples_trans)
    if len(sample_thick) == 1 and len(samples) > len(sample_thick):
        sample_thick = sample_thick * len(samples)

    # Checking if output directory exists, if it doesn't, creates the folder
    # Also, if do not overwrite, then makes sure the directory does not exists.
    output_dir = base_output_directory
    if not overWrite:
        suffix = 0
        while os.path.exists(output_dir):
            output_dir = base_output_directory[0, len(base_output_directory) - 2] + "_" + str(suffix) + "/"
            suffix += 1

    if sample_identifier is not '':
        if sample_identifier is not '':
            output_dir = base_output_directory + str(sample_identifier) + "/"
            change_outputdir = {
                'configuration': {
                    'outputDir': output_dir,
                },
            }
            common_configuration_full = update_reduction_parameters(common_configuration_full,
                                                                    change_outputdir,
                                                                    validate=False)
        for subfolder in ['1D', '2D']:
            output_folder = os.path.join(output_dir, subfolder)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

    start_time = time.time()
    # Loop for samples
    for i, sample_name in enumerate(sample_names):
        start_time_loop = time.time()

        # form base output file name
        if isinstance(samples[i], str) and os.path.exists(samples[i]):
            # a full path to NeXus is given
            part1 = os.path.basename(samples[i]).split('.')[0]
        else:
            part1 = samples[i]
        output_file_name = f'r{part1}_{sample_names[i]}'

        run_data = {
            'sample': {
                'runNumber': samples[i],
                'thickness': sample_thick[i],
                'transmission': {'runNumber': samples_trans[i]}
            },
            'background': {
                'runNumber': backgrounds[i],
                'transmission': {'runNumber': backgrounds_trans[i]}
            },
            'outputFileName': output_file_name,
        }

        # Update our common settings with the particulars of the current reduction
        reduction_input = update_reduction_parameters(common_configuration_full, run_data, validate=True)
        # pretty_print(reduction_input)
        reduction_input['configuration']['WedgeMinAngles'] = wedge_min_angles
        reduction_input['configuration']['WedgeMaxAngles'] = wedge_max_angles

        # Load all files
        loaded = load_all_files(reduction_input, use_nexus_idf=True)

        # Reduced from workspaces loaded from NeXus files
        out = reduce_single_configuration(loaded, reduction_input)

        plot_reduction_output(out, reduction_input)

        print('\nloop_' + str(i + 1) + ": ", time.time() - start_time_loop)

        if np.remainder(i, refreshCycle) == 0 and i > 0:
            # mtd.clear()
            clear_buffer()

    print('Total Time : ', time.time() - start_time)

    def generate_output_log_file(output_directory, file_base_name, file_suffix):
        filename = os.path.join(output_directory,  f'{file_base_name}_reduction_log{file_suffix}.hdf')
        return filename

    # samples shall be list of file names
    test_log_files = list()
    for i_s, sample in enumerate(samples):

        if isinstance(samples[i_s], str) and os.path.exists(samples[i_s]):
            # a full path to NeXus is given
            part1 = os.path.basename(samples[i_s]).split('.')[0]
        else:
            part1 = samples[i_s]
        output_file_name = f'r{part1}_{sample_names[i_s]}'

        test_log_file = generate_output_log_file(base_output_directory, output_file_name, '')
        print(test_log_file)
        assert os.path.exists(test_log_file), f'Output log file {test_log_file} cannot be found.'
        test_log_files.append(test_log_file)

    return test_log_files
