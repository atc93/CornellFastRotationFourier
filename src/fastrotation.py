import src.statfluctuation as statfluctuation
import src.configparser as configparser
import src.constants as constants
import src.plotting as plotting
import matplotlib.pyplot as plt
import src.style as style
import src.util as util
import numpy as np
import ROOT as r
import array
import json


class FastRotation(configparser.ParseConfig):

    def __init__(self, config):
        super(FastRotation, self).__init__(config)
        self.times_to_plot = [1, 10, 50, 100, 200, self.tM]
        self.canvas = r.TCanvas('c', 'c', 900, 600)
        self.out_file = r.TFile('results/' + self.tag + '/fastrotation.root','RECREATE')

    def optimize_tS_tM(self):

        # reset range user
        self.histogram.GetXaxis().SetRangeUser(1, self.histogram.GetNbinsX())

        # refine ts: we want ts corresponding to a signal intensity of 1
        # to match the intensity of 1 at late time (modulo the noise)
        min_bin = self.histogram.FindBin(self.tS-0.030)
        max_bin = self.histogram.FindBin(self.tS+0.030)
        delta = 999
        opt_tS = -1
        for bin_idx in range(min_bin, max_bin+1):
            if (abs(self.histogram.GetBinContent(bin_idx)-1) < delta):
                delta = abs(self.histogram.GetBinContent(bin_idx)-1)
                opt_tS = self.histogram.GetBinCenter(bin_idx)

        min_bin = self.histogram.FindBin(self.tM-0.100)
        max_bin = self.histogram.FindBin(self.tM+0.100)
        delta = 999
        opt_tM = -1
        for bin_idx in range(min_bin, max_bin+1):
            if (abs(self.histogram.GetBinContent(bin_idx)-1) < delta):
                delta = abs(self.histogram.GetBinContent(bin_idx)-1)
                opt_tM = self.histogram.GetBinCenter(bin_idx)

        return opt_tS, opt_tM

    def apply_stat_fluctutation(self):

        self.histogram = statfluctuation.poisson(self.histogram)

        return self.histogram

    def return_frs_np_array(self, frs_histogram):

        bin_center, bin_content = util.rootHistToNumpArray(
            frs_histogram, self.tS, self.tM)

        return (bin_center, bin_content)

    def fit_wiggle(self, wiggle_histogram):

        if (self.n_fit_param == 2):
            return self.fit_wiggle_2p(wiggle_histogram)
        elif (self.n_fit_param == 5):
            return self.fit_wiggle_5p(wiggle_histogram)
        elif (self.n_fit_param == 9):
            return self.fit_wiggle_9p(wiggle_histogram)

    def fit_wiggle_2p(self, wiggle_histogram):

        two_param_fit = r.TF1(
            "two_param_fit", "[0]*exp(-x/[1])", self.start_fit_time, self.tM)
        two_param_fit.SetParameters(wiggle_histogram.GetBinContent(
            wiggle_histogram.FindBin(self.start_fit_time))*1.3, 64.4)
        wiggle_fit = two_param_fit
        wiggle_histogram.Fit('two_param_fit', 'REMQ')

        if (self.verbose > 1):
            print('\n=== 2-parameter fit ===\n')
            print('N0 = ', two_param_fit.GetParameter(0))
            print('Lifetime = ', two_param_fit.GetParameter(1))

        return wiggle_fit

    def fit_wiggle_5p(self, wiggle_histogram):

        five_param_fit = r.TF1(
            "five_param_fit", "[0]*exp(-x/[1])*(1+[2]*cos(2*TMath::Pi()*[3]*x)+[4]*sin(2*TMath::Pi()*[3]*x))", self.start_fit_time, self.tM)
        five_param_fit.SetParameters(wiggle_histogram.GetBinContent(
            wiggle_histogram.FindBin(self.start_fit_time))*1.3, 64.4, -0.1, 0.229, -0.1)
        wiggle_fit = five_param_fit
        wiggle_histogram.Fit('five_param_fit', 'REMQ')

        if (self.verbose > 1):
            print('\n=== 5-parameter fit ===\n')
            print('N0 = ', five_param_fit.GetParameter(0))
            print('Lifetime = ', five_param_fit.GetParameter(1))
            print('Acos = ', five_param_fit.GetParameter(2))
            print('Asin = ', five_param_fit.GetParameter(4))
            print('Omega_a = ', five_param_fit.GetParameter(3))

        return wiggle_fit

    def fit_wiggle_9p(self, wiggle_histogram):

        five_param_fit = self.fit_wiggle_5p(wiggle_histogram)

        cbo_fit = r.TF1(
            "cbo", "1+exp(-x/[0])*([1]*cos(2*TMath::Pi()*[2]*x)+[3]*sin(2*TMath::Pi()*[2]*x))", self.start_fit_time, self.tM)
        nine_param_fit = r.TF1(
            "nine_param_fit", "five_param_fit*cbo", self.start_fit_time, self.tM)
        nine_param_fit.SetParameter(0, five_param_fit.GetParameter(0))
        nine_param_fit.SetParameter(1, five_param_fit.GetParameter(1))
        nine_param_fit.SetParameter(2, five_param_fit.GetParameter(2))
        nine_param_fit.SetParameter(3, five_param_fit.GetParameter(3))
        nine_param_fit.SetParameter(4, five_param_fit.GetParameter(4))
        nine_param_fit.SetParameter(5, 150)
        nine_param_fit.SetParameter(6, 0.004)
        nine_param_fit.SetParameter(7, 0.370)
        nine_param_fit.SetParameter(8, 0.004)
        nine_param_fit.SetNpx(10000)

        wiggle_fit = nine_param_fit
        wiggle_histogram.Fit('nine_param_fit', 'REMQ')

        if (self.verbose > 1):
            print('\n=== 9-parameter fit ===\n')
            print('N0 = ', nine_param_fit.GetParameter(0))
            print('Lifetime muon = ', nine_param_fit.GetParameter(1))
            print('Acos = ', nine_param_fit.GetParameter(2))
            print('Asin = ', nine_param_fit.GetParameter(4))
            print('Omega_a = ', nine_param_fit.GetParameter(3))
            print('Lifetime cbo = ', nine_param_fit.GetParameter(5))
            print('Acos cbo = ', nine_param_fit.GetParameter(6))
            print('Asin cbo = ', nine_param_fit.GetParameter(8))
            print('Omega cbo = ', nine_param_fit.GetParameter(7))

        return wiggle_fit

    def produce(self):

        print(' ### Step 2/4: produce fast rotation signal\n')

        # open input ROOT file
        in_file = r.TFile(self.root_file)

        # retrive and style input histogram
        self.histogram = in_file.Get(self.histo_name)
        style.setTH1Style(self.histogram, '', 'Time [#mus]', 'Intensity')

        # define and style canvas
        style.setTCanvasStyle(self.canvas)

        # vary statistics in each bin of the input histogram if option specified
        if (self.stat_fluctuation):
            self.histogram = self.apply_stat_fluctutation()

        # rebin fast rotation histogram
        if (self.n_fit_param == 0):
            self.histogram.Rebin(self.rebin_frs_factor)

        # plot input histogram if option specified
        if (self.print_plot):
            for time in self.times_to_plot:
                plotting.plot(self.canvas, self.histogram, self.tag +
                              '/Intensity', self.tS, self.tS + time, self.tM)

        # return histogram if no need for fit
        if (self.n_fit_param == 0):

            opt_tS, opt_tM = self.optimize_tS_tM()
            bin_center, bin_content = self.return_frs_np_array(self.histogram)

            return opt_tS, opt_tM, bin_center, bin_content

        # clone input histogram to produce wiggle plot
        wiggle_histogram = self.histogram.Clone()

        # rebin cloned histograms
        wiggle_histogram.Rebin(self.rebin_wiggle_factor)

        # plot fitted wiggle plot if option specified
        if (self.print_plot):
            for time in self.times_to_plot:
                plotting.plot(self.canvas, wiggle_histogram, self.tag +
                              '/Wiggle', self.tS, self.start_fit_time + time, self.tM)

        # fit wiggle plot
        wiggle_fit = self.fit_wiggle(wiggle_histogram)

        # create histogram of the fit residuals
        residuals = wiggle_histogram
        residuals.GetYaxis().SetTitle('Residual [%]')
        for i in range(wiggle_histogram.GetNbinsX()):
            if (wiggle_histogram.GetBinContent(i+1) != 0):
                residuals.SetBinContent(i+1, (wiggle_histogram.GetBinContent(i+1)-wiggle_fit.Eval(
                    wiggle_histogram.GetBinCenter(i+1)))/wiggle_histogram.GetBinContent(i+1)*100)

        # plot histogram of residuals
        if (self.print_plot):
            for time in self.times_to_plot:
                plotting.plot(self.canvas, residuals, self.tag + '/Residuals',
                              self.tS, self.start_fit_time + time, self.tM)

        # create fast rotation histogram
        for i in range(self.histogram.GetNbinsX()):
            self.histogram.SetBinContent(i+1, self.histogram.GetBinContent(i+1)/(
                wiggle_fit.Eval(self.histogram.GetBinCenter(i+1))/self.rebin_wiggle_factor))

        # rebin fast rotation histogram
        self.histogram.Rebin(self.rebin_frs_factor)

        opt_tS, opt_tM = self.optimize_tS_tM()

        # plot frs plot if option specified
        if (self.print_plot):
            plotting.plot(self.canvas, self.histogram, self.tag + '/FRS', 0, 10, self.tM)
            for time in self.times_to_plot:
                plotting.plot(self.canvas, self.histogram, self.tag +
                              '/FRS', round(opt_tS, 6), round(opt_tS + time, 6), round(opt_tM, 6))

        # save signal to ROOT file
        self.out_file.cd()
        self.histogram.Write('fr')

        bin_center, bin_content = self.return_frs_np_array(self.histogram)

        return opt_tS, opt_tM, bin_center, bin_content
