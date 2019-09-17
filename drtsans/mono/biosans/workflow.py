from drtsans.mono.biosans.beam_finder import center_detector


def prepare_data(input_workspace,
                 x_center=None, y_center=None, y_center_gravity=None,
                 dark_current=None,
                 flux_method=None, flux=None,
                 mask=None, panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None,
                 output_workspace=None):
    ##########
    # TODO

    center_detector(input_workspace, x_center, y_center, y_center_gravity)
