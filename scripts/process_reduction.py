#!/usr/bin/env python3
import json
import os
import subprocess
import sys

# call as:
# python3 process_reduction.py filename.json
# the json content is {"run_number": 88980, "instrument": "EQSANS",
#                      "output_dir": "./SANS_output", "fail": false}

if __name__ == '__main__':
    if len(sys.argv) == 2:
        json_file = sys.argv[1]
        with open(json_file, 'r') as jsf:
            json_parameters = json.load(jsf)
        instrument = json_parameters['instrumentName']
        reduction_script = 'sans_reduction_test.py'
        if instrument == 'EQSANS':
            reduction_script = 'eqsans_reduction.py'
        elif instrument == 'BIOSANS':
            reduction_script = 'biosans_reduction.py'
        filename_string = json_parameters['outputFilename']
        output_folder = json_parameters['configuration']['outputDir']
        out_log = os.path.join(output_folder,
                               filename_string+'.out')
        out_err = os.path.join(output_folder,
                               filename_string+'.err')
        logFile = open(out_log, "w")
        errFile = open(out_err, "w")
        json_string = json.dumps(json_parameters)
        script_folder = os.path.abspath(os.path.dirname(__file__))
        cmd = "python3 {}".format(os.path.join(script_folder, reduction_script))
        cmd += " '{}'".format(json_string)
        proc = subprocess.Popen(cmd,
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=logFile,
                                stderr=errFile,
                                universal_newlines=True)
        proc.communicate()
        logFile.close()
        errFile.close()
        rc = proc.returncode
        if rc:
            exit(rc)
        if os.path.isfile(out_err) and os.stat(out_err).st_size == 0:
            os.remove(out_err)
        exit(os.path.isfile(out_err))
    else:
        raise RuntimeError("A parameter json string is required as input")
