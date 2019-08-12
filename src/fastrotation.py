import src.statfluctuation as statfluctuation
import src.constants as constants
import src.plotting as plotting
import matplotlib.pyplot as plt
import src.common as common
import src.style as style
import src.util as util
import numpy as np
import ROOT as r
import array
import json
import sys

# class that produces the fast rotation signal
# it inherits from the Common class that itself
# inherits from the ParseConfig class


class FastRotation(common.Init):

    # class initialization
    def __init__(self, config):
        super(FastRotation, self).__init__(config)

        # define the time ranges to plot the various signal from tS t0 tS+self.times_to_plot in micro-seconds
        self.times_to_plot = [1, 10, 50, 100, 200, self.tM]

        # create the TCanvas for all the plots
        canvas_fastrotation = r.TCanvas(
            'c_fastrotation', 'c_fastrotation', 900, 600)
        style.setTCanvasStyle(canvas_fastrotation)

        # create output ROOT file to store the fast rotation signal
        self.out_file = r.TFile(
            'results/' + self.tag + '/fastrotation.root', 'RECREATE')

    # method that optimizes tS and tM : we want tS corresponding to the fast rotation signal intensity
    # of 1 to match the intensity of 1 at tM (not exactly 1 because of statistical noise)
    # see docdb...
    def optimize_tS_tM(self):

        # reset range user in case it was changed
        self.histogram.GetXaxis().SetRangeUser(1, self.histogram.GetNbinsX())

        # use a 60 micro-seconds time window (less than half a cyclotron period)
        # to find the signal intensity of 1
        min_bin = self.histogram.FindBin(self.tS-0.030)
        max_bin = self.histogram.FindBin(self.tS+0.030)

        # initialize required variables
        delta = 999
        opt_tS = -1

        for bin_idx in range(min_bin, max_bin+1):
            if (abs(self.histogram.GetBinContent(bin_idx)-1) < delta):
                delta = abs(self.histogram.GetBinContent(bin_idx)-1)
                opt_tS = self.histogram.GetBinCenter(bin_idx)

        # use a 200 micro-secs time window to find tM such as
        # the signal intensity is as close to 1 as possible
        min_bin = self.histogram.FindBin(self.tM-0.100)
        max_bin = self.histogram.FindBin(self.tM+0.100)

        # initialize required variables
        delta = 999
        opt_tM = -1

        for bin_idx in range(min_bin, max_bin+1):
            if (abs(self.histogram.GetBinContent(bin_idx)-1) < delta):
                delta = abs(self.histogram.GetBinContent(bin_idx)-1)
                opt_tM = self.histogram.GetBinCenter(bin_idx)

        # in case the optimization did not complete succesfully
        if (opt_tM == -1):
            opt_tM = self.tM
            print(' Warning -- could not optimize tM')

        # returns results to the main method 'produce'
        return opt_tS, opt_tM

    # method that applies statistical fluctuation to the input
    # positron counts spectrum
    def apply_stat_fluctutation(self):

        self.histogram = statfluctuation.poisson(self.histogram)

        return self.histogram

    # method that converts a ROOT histogram into a numpy array
    def return_frs_np_array(self, frs_histogram):

        bin_center, bin_content = util.rootHistToNumpArray(
            frs_histogram, self.tS, self.tM)

        return (bin_center, bin_content)

    # method that perform the wiggle fit (with different parameter numbers)
    def fit_wiggle(self, wiggle_histogram):

        # for 2-parameter fit (exponential decay)
        if (self.n_fit_param == 2):
            return self.fit_wiggle_2p(wiggle_histogram)
        # for 5-parameter fit (exponential decay+omega_a)
        elif (self.n_fit_param == 5):
            return self.fit_wiggle_5p(wiggle_histogram)
        # for 9-parameter fit (exponential decay+omega_a+omega_cbo)
        elif (self.n_fit_param == 9):
            return self.fit_wiggle_9p(wiggle_histogram)

    # method that fits the wiggle plot with a 2-parameter function (exponential decay)
    def fit_wiggle_2p(self, wiggle_histogram):

        # define ROOT TF1 fit function
        two_param_fit = r.TF1(
            "two_param_fit", "[0]*exp(-x/[1])", self.start_fit_time, self.tM)
        # provide initial guesses for the 2 fit parameters
        two_param_fit.SetParameters(wiggle_histogram.GetBinContent(
            wiggle_histogram.FindBin(self.start_fit_time))*1.3, 64.4)
        # perform the fit
        wiggle_histogram.Fit('two_param_fit', 'REMQ')

        if (self.verbose > 1):
            print('\n=== 2-parameter fit ===\n')
            print('N0 = ', two_param_fit.GetParameter(0))
            print('Lifetime = ', two_param_fit.GetParameter(1))

        return two_param_fit

    # method that fits the wiggle plot with a 5-parameter function (exponential decay+omega_a)
    def fit_wiggle_5p(self, wiggle_histogram):

        # define ROOT TF1 fit function
        five_param_fit = r.TF1(
            "five_param_fit", "[0]*exp(-x/[1])*(1+[2]*cos(2*TMath::Pi()*[3]*x)+[4]*sin(2*TMath::Pi()*[3]*x))", self.start_fit_time, self.tM)
        # provide initial guesses for the 5 fit parameters
        five_param_fit.SetParameters(wiggle_histogram.GetBinContent(
            wiggle_histogram.FindBin(self.start_fit_time))*1.3, 64.4, -0.1, 0.229, -0.1)
        # perform the fit
        wiggle_histogram.Fit('five_param_fit', 'REMQ')

        if (self.verbose > 1):
            print('\n=== 5-parameter fit ===\n')
            print('N0 = ', five_param_fit.GetParameter(0))
            print('Lifetime = ', five_param_fit.GetParameter(1))
            print('Acos = ', five_param_fit.GetParameter(2))
            print('Asin = ', five_param_fit.GetParameter(4))
            print('Omega_a = ', five_param_fit.GetParameter(3))

        return five_param_fit

    # method that fits the wiggle plot with a 9-parameter function (exponential decay+omega_a+omega_cbo)
    def fit_wiggle_9p(self, wiggle_histogram):

        # perform first the 5-parameter fit
        five_param_fit = self.fit_wiggle_5p(wiggle_histogram)

        # define the CB0 function
        cbo_fit = r.TF1(
            "cbo", "1+exp(-x/[0])*([1]*cos(2*TMath::Pi()*[2]*x)+[3]*sin(2*TMath::Pi()*[2]*x))", self.start_fit_time, self.tM)
        # define the 9-parameter fit as the 5-parameter fit + the CB0 function
        nine_param_fit = r.TF1(
            "nine_param_fit", "five_param_fit*cbo", self.start_fit_time, self.tM)
        # provide initial guesses for the 9 fit parameters
        nine_param_fit.SetParameter(0, five_param_fit.GetParameter(0))
        nine_param_fit.SetParameter(1, five_param_fit.GetParameter(1))
        nine_param_fit.SetParameter(2, five_param_fit.GetParameter(2))
        nine_param_fit.SetParameter(3, five_param_fit.GetParameter(3))
        nine_param_fit.SetParameter(4, five_param_fit.GetParameter(4))
        nine_param_fit.SetParameter(5, 150)
        nine_param_fit.SetParameter(6, 0.003)
        nine_param_fit.SetParameter(7, self.cbo_freq/1000)
        nine_param_fit.SetParameter(8, 0.004)
        nine_param_fit.SetNpx(10000)

        # perform the fit
        wiggle_histogram.Fit('nine_param_fit', 'REMQ', 'same')

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

        return nine_param_fit

    # main method that produces the fast rotation signal
    def produce(self):

        print(' ### Step 2/4: produce fast rotation signal\n')

        # open input ROOT file containing the input positron counts histogram
        in_file = r.TFile(self.root_file)

        # retrive and style input histogram
        self.histogram = in_file.Get(self.histo_name)
        style.setTH1Style(self.histogram, '', 'Time [#mus]', 'Intensity')

        # save input histogram to ROOT file
        self.out_file.cd()
        self.histogram.Write('positron_spectrum')

        # define and style canvas (canvas local to 'produce' method because not needed elsewhere)
        canvas_fastrotation = r.TCanvas(
            'c_fastrotation', 'c_fastrotation', 900, 600)
        style.setTCanvasStyle(canvas_fastrotation)

        # vary statistics in each bin of the input histogram if option specified
        if (self.stat_fluctuation):
            self.histogram = self.apply_stat_fluctutation()

        # rebin fast rotation histogram
        if (self.n_fit_param == 0):
            self.histogram.Rebin(self.rebin_frs_factor)

        # plot input histogram if option specified
        if (self.print_plot >= 2):
            for time in self.times_to_plot:
                plotting.plot(canvas_fastrotation, self.histogram, self.tag +
                              '/Intensity', self.tS, self.tS + time, self.tM)

        # check number of hits in positron spectrum if option specified
        if (self.check_positron_hits):
            n_positron_hits = self.histogram.Integral(
                self.histogram.FindBin(self.tS), self.histogram.FindBin(self.tM))
            if (self.verbose > 1):
                print(' Number of positron analyzed: ', n_positron_hits)
            if (n_positron_hits < self.positron_hits_threshold):
                print(' Warning -- low positron statistics: program exiting')
                sys.exit(1)
        else:
            n_positron_hits = -1

        # return histogram if no need for fit (typically for simulated fast rotation signals)
        if (self.n_fit_param == 0):
            # optimize tS and tM
            opt_tS, opt_tM = self.optimize_tS_tM()
            bin_center, bin_content = self.return_frs_np_array(self.histogram)

            # save signal to ROOT file
            self.out_file.cd()
            self.histogram.Write('fast_rotation')

            return opt_tS, opt_tM, n_positron_hits, bin_center, bin_content

        # clone input histogram to produce wiggle plot (histogram local to 'produce' method because not needed elsewhere)
        wiggle_histogram = self.histogram.Clone()

        # rebin cloned histograms
        wiggle_histogram.Rebin(self.rebin_wiggle_factor)

        # plot wiggle plot if option specified
        if (self.print_plot >= 2):
            for time in self.times_to_plot:
                plotting.plot(canvas_fastrotation, wiggle_histogram, self.tag +
                              '/Wiggle', self.tS, self.start_fit_time + time, self.tM)

        # fit wiggle plot (fit local to 'produce' method because not needed elsewhere)
        wiggle_fit = self.fit_wiggle(wiggle_histogram)

        # plot fitted wiggle plot if option specified
        if (self.print_plot >= 2):
            for time in self.times_to_plot:
                plotting.plot(canvas_fastrotation, wiggle_histogram, self.tag +
                              '/FittedWiggle', self.tS, self.start_fit_time + time, self.tM, wiggle_fit)

        # create histogram of the fit residuals
        residuals = wiggle_histogram
        residuals.GetYaxis().SetTitle('Residual [%]')
        for i in range(wiggle_histogram.GetNbinsX()):
            if (wiggle_histogram.GetBinContent(i+1) != 0):
                residuals.SetBinContent(i+1, (wiggle_histogram.GetBinContent(i+1)-wiggle_fit.Eval(
                    wiggle_histogram.GetBinCenter(i+1)))/wiggle_histogram.GetBinContent(i+1)*100)

        # plot histogram of residuals if option specified
        if (self.print_plot >= 2):
            for time in self.times_to_plot:
                plotting.plot(canvas_fastrotation, residuals, self.tag + '/Residuals',
                              self.tS, self.start_fit_time + time, self.tM)
                plotting.plot(canvas_fastrotation, residuals, self.tag + '/Residuals',
                              self.start_fit_time, self.start_fit_time + time, self.tM)

        # create fast rotation histogram
        for i in range(self.histogram.GetNbinsX()):
            self.histogram.SetBinContent(i+1, self.histogram.GetBinContent(i+1)/(
                wiggle_fit.Eval(self.histogram.GetBinCenter(i+1))/self.rebin_wiggle_factor))

        # rebin fast rotation histogram if option specified
        self.histogram.Rebin(self.rebin_frs_factor)

        # optimize tS and tM
        opt_tS, opt_tM = self.optimize_tS_tM()

        # plot frs plot if option specified
        if (self.print_plot >= 1):
            plotting.plot(canvas_fastrotation, self.histogram,
                          self.tag + '/FRS', 0, 10, self.tM)
            for time in self.times_to_plot:
                plotting.plot(canvas_fastrotation, self.histogram, self.tag +
                              '/FRS', round(opt_tS, 6), round(opt_tS + time, 6), round(opt_tM, 6))

        # save signal to ROOT file
        self.out_file.cd()
        self.histogram.Write('fast_rotation')

        # retrieve arrays of bin contents and bin centers
        bin_center, bin_content = self.return_frs_np_array(self.histogram)

        # return results
        return opt_tS, opt_tM, n_positron_hits, bin_center, bin_content
