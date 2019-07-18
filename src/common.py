import src.configparser as configparser
import matplotlib.pyplot as plt
import src.style as style
import ROOT as r

class Init(configparser.ParseConfig):

    def __init__(self, config):
        super(Init, self).__init__(config)

        # matplotlib configuration: set params for all 
        # the subsequent calls to matplotlib
        plt.rcParams['figure.figsize'] = [9, 6]
        plt.rcParams.update({'font.size': 17})
