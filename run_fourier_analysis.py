import src.analysis as analysis
from shutil import copyfile
import src.util as util
import ROOT as r
import json
import sys
import re

# do not show the ROOT INFO messages (e.g. Canvas being saved to EPS)
r.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")


def main():

    # read in configuration file
    with open(sys.argv[1]) as config_file:
        # Remove '#'-denoted comments from the configuration file.
        comment_pattern = r"(?m)^\s*#.*\n?"
        config_text = re.sub(comment_pattern, "", config_file.read())
        # Parse the configuration string.
        config = json.loads(config_text)

    # manage results directory
    util.manage_results_directory('results/' + config['tag'])

    # copy configuration file to results directoy
    copyfile(sys.argv[1], 'results/' + config['tag'] + '/config.json')

    # run fourier analysis
    return analysis.run_fourier(config)


if __name__ == '__main__':

    # print welcome message
    util.print_welcome_message()

    # check that the configuration file was provided
    # exit if not provided
    if len(sys.argv) < 2:
        print('must provide config!')
        sys.exit(0)

    # run the analysis for associated configuration file
    sys.exit(main())
