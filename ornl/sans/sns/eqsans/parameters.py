import os
from configparser import RawConfigParser
from itertools import chain
from glob import glob

from ornl.settings import MultiOrderedDict

from dotenv import load_dotenv
load_dotenv()

CONFIG_DIRECTORY = os.getenv('EQSANS_CONFIG_DIRECTORY')
CONFIG_FILE_PREFIX = os.getenv('EQSANS_CONFIG_FILE_PREFIX')


def _get_config_file(run_number):
    '''
    Given a run number get the respective configuration file
    The numbers are for this run or when the run starts
    '''
    files = glob(os.path.join(CONFIG_DIRECTORY,
                              CONFIG_FILE_PREFIX+"[0-9]*[0-9]"))
    extensions = [int(os.path.splitext(f)[-1][1:]) for f in files]
    extensions = sorted(extensions)
    run_to_use = None
    if run_number in extensions:
        run_to_use = run_number
    else:
        for ext in extensions:
            if ext < run_number:
                run_to_use = ext
            else:
                break
    return os.path.join(CONFIG_DIRECTORY,
                        CONFIG_FILE_PREFIX+str(run_to_use))


def get_parameters(run_number):
    '''
    Get the parameters from the configuration file
    If the same key exist, the value is appended with '\n'
    Returns a dictionary
    '''
    conf_file = _get_config_file(run_number)
    parser = RawConfigParser(
        dict_type=MultiOrderedDict,
        strict=False,
        inline_comment_prefixes="#")
    with open(conf_file, 'r') as f:
        f = chain(("[DEFAULT]",), f)  # This line does the trick.
        parser.read_file(f)
    return {i[0]: i[1] for i in parser['DEFAULT'].items()}
