from scipy.optimize import curve_fit
import src.constants as constants
import src.plotting as plotting
import matplotlib.pyplot as plt
import src.common as common
import src.style as style
import src.util as util
import numpy as np
import ROOT as r
import array

# class that runs the Fourier analysis given
# the optimized t0 provided by the previous
# step


class Fourier(common.Init):

    # class initialization
    def __init__(self, config, bin_center, bin_content, opt_t0, bkg_stat_noise, fit_boundary1, fit_boundary2, n_positron_hits):
        super(Fourier, self).__init__(config)
        self.bin_center = bin_center
        self.bin_content = bin_content
        self.opt_t0 = opt_t0
        self.fit_boundary1 = fit_boundary1
        self.fit_boundary2 = fit_boundary2
        self.bkg_stat_noise = bkg_stat_noise
        self.n_positron_hits = n_positron_hits

        # create the TCanvas for all the plots
        self.canvas_fourier = r.TCanvas('c_fourier', 'c_fourier', 900, 600)
        style.setTCanvasStyle(self.canvas_fourier)

        # create the TH1 to store the cosine and sine Fourier transforms
        self.cosine_histogram = r.TH1D(
            'cosine', 'cosine', self.n_freq_step, self.lower_freq, self.upper_freq)
        self.sine_histogram = r.TH1D(
            'sine', 'sine', self.n_freq_step, self.lower_freq, self.upper_freq)

        # create the TFile that will store the various results
        self.out_file2 = r.TFile(
            'results/' + self.tag + '/results.root', 'UPDATE')

        # when performing scan and wanting to keep in separate ROOT files
        # the results from each scan point
        if (self.print_plot >= 2):
            self.out_file = r.TFile('results/' + self.tag + '/results_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.root'.format(
                self.opt_t0, self.tS, self.tM, self.freq_step_size), 'RECREATE')

    # method that does a little bit of printing to the terminal (input info)
    def print_info(self):

        print('\n ### Step 4/4: produce results\n')

        print('    -- input parameters --')
        print('    t0 = ', round(self.opt_t0, 6), ' us')
        print('    tS = ', self.tS, 'us')
        print('    tM = ', self.tM, ' us')
        print('    df = ', self.freq_step_size, ' kHz')
        print('    n  = ', self.field_index)

    # method that does a little bit of printing to the terminal (results)
    def print_results(self, xe, width, ce, mean_freq, std_freq):

        print('\n    -- results --')
        print('    mean freq   =', round(mean_freq, 2), ' kHz')
        print('    width freq  =', round(std_freq, 2), ' kHz')
        print('    xe          =', round(xe, 2), ' mm')
        print('    width       =', round(width, 2), ' mm')
        print('    ce          =', round(ce,), ' ppb')

    # method that computes the cosine and sine Fourier transform
    def produce_cosine_transform(self):

        # compute cosine Fourier transform
        freq, intensity = util.calc_cosine_transform(
            self.opt_t0, self.bin_content, self.bin_center, self.freq_step_size, self.n_freq_step, self.lower_freq)

        # fill ROOT TH1 histogram
        for x, y in zip(freq, intensity):
            iBin = self.cosine_histogram.FindBin(x)
            self.cosine_histogram.SetBinContent(iBin, y)

        # compute sine Fourier transform if specified in config file
        if (self.calc_sine):
            freq, intensity = util.calc_sine_transform(
                self.opt_t0, self.bin_content, self.bin_center, self.freq_step_size, self.n_freq_step, self.lower_freq)

            # fill ROOT TH1 histogram
            for x, y in zip(freq, intensity):
                iBin = self.sine_histogram.FindBin(x)
                self.sine_histogram.SetBinContent(iBin, y)

        # define histogram(s) name
        histoName = ['Cosine transform: t_{0} = ' + '{0:.2f} ns'.format(
            self.opt_t0*1000), 'Sine transform: t_{0} = ' + '{0:.2f} ns'.format(self.opt_t0)]

        # clone cosine histogram for plotting/styling puporses
        cosine_clone = self.cosine_histogram.Clone()

        # clone sine histogram for plotting/styling puporses
        sine_clone = self.sine_histogram.Clone()

        # create a list containing both cosine and sine histograms
        clone_hist_list = [cosine_clone, sine_clone]

        # styling and plotting the cosine and sine Fourier transform
        for idx in range(0, 2):

            # style the Fourier transform histogram ==#
            style.setTH1Style(
                clone_hist_list[idx], histoName[idx], 'Frequency [kHz]', 'Arbitrary units', 1.2, 1.1)

            # define lines to be drawn at collimator apertures (frequency space) ==#
            inner_line, outer_line = style.set_collimator_aperture_tline(
                clone_hist_list[idx].GetMinimum(), clone_hist_list[idx].GetMaximum(), 'frequency')

            # define pave text to go along the collimator apertures lines ==#
            pt, pt2 = style.setCollimatorAperturePaveText(clone_hist_list[idx].GetMaximum(
            )*0.38, clone_hist_list[idx].GetMaximum()*0.52, 'frequency')

            # draw it all
            list_to_draw = [clone_hist_list[idx],
                            inner_line, outer_line, pt, pt2]
            plotting.plotMultipleObjects('', list_to_draw)
            self.canvas_fourier.Draw()

            # save TH1 histogram to output ROOT file
            self.out_file2.cd()
            clone_hist_list[idx].Write()

            # save TH1 histogram to output ROOT file if option specified
            # for scan purposes
            if (self.print_plot >= 2):
                self.out_file.cd()
                clone_hist_list[idx].Write()

            # print plot if option provided
            if (idx == 0 and self.print_plot >= 2):
                self.canvas_fourier.Print('results/' + self.tag + '/Cosine_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                    self.opt_t0, self.tS, self.tM, self.freq_step_size))

            # limit the frequency (x-axis) to the collimator aperture
            clone_hist_list[idx].GetXaxis().SetRangeUser(
                constants.lowerCollimatorFreq, constants.upperCollimatorFreq)
            clone_hist_list[idx].SetTitle('')
            clone_hist_list[idx].Draw()
            self.canvas_fourier.Draw()

            # print plot if option provided
            if (idx == 0 and self.print_plot >= 1):
                self.canvas_fourier.Print('results/' + self.tag + '/CosineCollAperture_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                    self.opt_t0, self.tS, self.tM, self.freq_step_size))

            # print plot if option provided and if sine transform was performed
            if (idx == 1 and self.print_plot >= 2 and self.calc_sine == 1):
                self.canvas_fourier.Print('results/' + self.tag + '/Sine_t0_{0:.3f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                    self.opt_t0, self.tS, self.tM, self.freq_step_size))

    # method that performs the background fit
    def fit_background(self):

        # create list to store the background of the cosine Fourier transform
        a = []
        b = []

        # append the background data point according to the optimized fit
        # foundary in the previous step as well as the option regarding
        # physical and/or non-physical frequencies allowed in the background
        for bin_idx in range(1, self.cosine_histogram.GetNbinsX()+1):
            if (self.background_frequencies == 'physical'):
                if ((self.cosine_histogram.GetBinCenter(bin_idx) <= self.fit_boundary1 and self.cosine_histogram.GetBinCenter(bin_idx) >= constants.lowerCollimatorFreq) or
                        (self.cosine_histogram.GetBinCenter(bin_idx) >= self.fit_boundary2 and self.cosine_histogram.GetBinCenter(bin_idx) <= constants.upperCollimatorFreq)):
                    a.append(self.cosine_histogram.GetBinCenter(bin_idx))
                    b.append(self.cosine_histogram.GetBinContent(bin_idx))
            elif (self.background_frequencies == 'all'):
                if (self.cosine_histogram.GetBinCenter(bin_idx) <= self.fit_boundary1 or self.cosine_histogram.GetBinCenter(bin_idx) >= self.fit_boundary2):
                    a.append(self.cosine_histogram.GetBinCenter(bin_idx))
                    b.append(self.cosine_histogram.GetBinContent(bin_idx))

        # create list with background statistical noise
        err = [self.bkg_stat_noise]*len(a)

        # create list to store fit parameters
        popt = []

        # fit the background according to the fit function
        # specified in the config file
        if (self.background_fit == 'pol'):
            fit_fn, chi2, sigma = util.fit_pol(
                a, b, err, self.poly_order, self.bkg_stat_noise)

        if (self.background_fit == 'sinc'):
            fit_fn, chi2, sigma, popt, pcov = util.fit_sinc(
                a, b, err, self.rebin_frs_factor, self.verbose)

        if (self.background_fit == 'erfi'):
            fit_fn, chi2, sigma, popt, pcov = util.fit_erfi(
                a, b, err, self.tS, self.verbose)

        if (self.background_fit == 'triangle'):
            fit_fn, chi2, sigma, popt, pcov = util.fit_triangle(
                a, b, err, self.tS, self.verbose)

        # plot the background
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        plt.subplots_adjust(left=0.125, bottom=0.125,
                            right=0.95, top=0.95, wspace=0, hspace=0)
        plt.errorbar(a, b, yerr=err, marker='o', ms=5, markerfacecolor='black',
                     ecolor='black', markeredgecolor='black', linestyle='', label='background', zorder=4)

        a = []
        b = []

        # append the entire cosine Fourier transform
        for bin_idx in range(1, self.cosine_histogram.GetNbinsX()+1):
            a.append(self.cosine_histogram.GetBinCenter(bin_idx))
            b.append(self.cosine_histogram.GetBinContent(bin_idx))

        # plot frequency distribution alongside the background and its fit (lines above)
        plt.plot(a, b, marker='o', ms=5, zorder=1)
        plt.xlabel('Frequency [kHz]')
        plt.ylabel('Arbitrary units')
        plt.plot(a, fit_fn(a, *popt), label='background fit',
                 linewidth=3, zorder=3, color='green')
        plt.legend(loc="upper right", frameon=False)

        # print plot if enabled by config file
        if (self.print_plot >= 1):
            plt.savefig('results/' + self.tag + '/Background_fit_t0_{0:.6f}_tS_{1}_tM_{2}_df_{3}_{4}_{5}.eps'.format(
                        self.opt_t0, self.tS, self.tM, self.freq_step_size, self.fit_boundary1, self.fit_boundary2), format='eps')

        # returns results
        return fit_fn(a, *popt), chi2

    # method that corrects the cosine Fourier transform
    # for its background using the background fit results
    def correct_transform(self, fit):

        max_amplitude_bin_idx = self.cosine_histogram.GetMaximumBin()

        # loop over the bin entries and correct the cosine Fourier transform by the fit
        # function and/or set the entry to zero if the background removal is allowed
        # the fit function is currently an array corresponding to the frequency window
        for bin_idx in range(max_amplitude_bin_idx, 0, -1):
            if (self.remove_background and abs(self.cosine_histogram.GetBinContent(bin_idx) - fit[int(bin_idx)-1]) < self.background_removal_threshold*self.bkg_stat_noise):
                self.cosine_histogram.SetBinContent(bin_idx, 0)
            else:
                self.cosine_histogram.SetBinContent(bin_idx, self.cosine_histogram.GetBinContent(
                    bin_idx)-fit[int(bin_idx)-1])

        for bin_idx in range(max_amplitude_bin_idx+1, self.cosine_histogram.GetNbinsX()+1, +1):
            if (self.remove_background and abs(self.cosine_histogram.GetBinContent(bin_idx) - fit[int(bin_idx)-1]) < self.background_removal_threshold*self.bkg_stat_noise):
                self.cosine_histogram.SetBinContent(bin_idx, 0)
            else:
                self.cosine_histogram.SetBinContent(bin_idx, self.cosine_histogram.GetBinContent(
                    bin_idx)-fit[int(bin_idx)-1])

        # style the corrected frequency distribution histogram
        style.setTH1Style(self.cosine_histogram, '',
                          'Frequency [kHz]', 'Arbitrary units', 1, 1.1)
        self.cosine_histogram.SetMinimum(-0.5)

        # define lines to be drawn at collimator apertures (frequency space)
        inner_line, outer_line = style.set_collimator_aperture_tline(
            self.cosine_histogram.GetMinimum(), self.cosine_histogram.GetMaximum(), 'frequency')

        # define pave text to go along the collimator apertures lines
        pt, pt2 = style.setCollimatorAperturePaveText(self.cosine_histogram.GetMaximum(
        )*0.38, self.cosine_histogram.GetMaximum()*0.52, 'frequency')
        pt3 = style.setRecFrequenciesPaveText(self.cosine_histogram.GetMean(), self.cosine_histogram.GetRMS(
        ), self.cosine_histogram.GetMaximum()*0.38, self.cosine_histogram.GetMaximum()*0.52)

        # draw it all
        list_to_draw = [self.cosine_histogram, inner_line, outer_line, pt, pt2]
        plotting.plotMultipleObjects('', list_to_draw)
        self.canvas_fourier.Draw()

        # print plot if enabled by config file
        if (self.print_plot >= 1):
            self.canvas_fourier.Print('results/' + self.tag + '/CorrectedCosine_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                self.opt_t0, self.tS, self.tM, self.freq_step_size))

        # save TH1 histogram to output ROOT file
        self.out_file2.cd()
        self.cosine_histogram.Write('corrected_cosine')

        # save TH1 histogram to output ROOT file if option specified
        # for scan purposes
        if (self.print_plot >= 2):
            self.out_file.cd()
            self.cosine_histogram.Write('corrected_cosine')

        # plot the corrected cosine Fourier transform limited to the collimator aperture
        self.cosine_histogram.GetXaxis().SetRangeUser(
            constants.lowerCollimatorFreq, constants.upperCollimatorFreq)
        self.cosine_histogram.SetTitle('')
        self.cosine_histogram.Draw('hist')
        self.canvas_fourier.Draw()

        # print plot if enabled by config file
        if (self.print_plot >= 1):
            self.canvas_fourier.Print('results/' + self.tag + '/CorrectedCosineCollAperture_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                self.opt_t0, self.tS, self.tM, self.freq_step_size))

        # compare frequency distributions at truth/reconstructed levels if optiion specified for simulated data
        if (self.compare_with_truth):

            # normalize the integral of the reco frequency distribution to 1
            self.cosine_histogram.Scale(1/self.cosine_histogram.Integral())

            # retrieve the truth level frequency distribution
            truth_file = r.TFile(self.truth_root_file)
            truth_histogram = truth_file.Get(self.truth_freq_histo_name)

            # rebin the truth level distribution if needed to adjust to the reconstructed one
            if (self.freq_step_size > truth_histogram.GetBinWidth(1)):
                truth_histogram.Rebin(
                    int(self.freq_step_size / truth_histogram.GetBinWidth(1)))

            # normalize the integral of the truth frequency distribution to 1
            truth_histogram.Scale(1/truth_histogram.Integral())

            # draw reconstructed level distribution
            self.cosine_histogram.Draw("hist")

            # style the truth level distribution
            truth_histogram.SetMarkerColor(2)
            truth_histogram.SetLineColor(2)
            truth_histogram.SetLineStyle(1)
            truth_histogram.SetLineWidth(0)
            truth_histogram.SetMarkerStyle(20)
            truth_histogram.SetMarkerSize(0.8+0.4*(self.freq_step_size-1))

            # draw the markers of the truth level distribution
            truth_histogram.Draw("samehistP0")

            # draw the line of the truth level distribution
            truth_histogram.Draw("samehist")

            # define pave text to display truth/reco level mean/width frequencies
            pt = style.setRecTruthFrequenciesPaveText(self.cosine_histogram.GetMean(), truth_histogram.GetMean(),
                                                      self.cosine_histogram.GetRMS(), truth_histogram.GetRMS(),
                                                      self.cosine_histogram.GetMaximum()*1., self.cosine_histogram.GetMaximum()*0.8)

            # draw the TPaveText
            pt.Draw("same")

            # draw TCanvas
            self.canvas_fourier.Draw()

            # print plot if enabled by config file
            if (self.print_plot >= 1):
                self.canvas_fourier.Print('results/' + self.tag + '/CorrectedCosineCollApertureWithTruth_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(
                    self.opt_t0, self.tS, self.tM, self.freq_step_size))

        # return results: x_e and width
        return (self.cosine_histogram.GetMean(), self.cosine_histogram.GetRMS())

    # method that converts the frequency distribution into its radial counter part
    def produce_radial(self, chi2, mean_freq, std_freq):

        # define arrays containing radius and intensity information
        # the array will be used to produce a TGraph (required since non equidistance radial points)
        intensity, radius = array.array('d'), array.array('d')

        # fill radius and intensity arrays from the frequency distribution
        util.convert_freq_to_radius(
            self.cosine_histogram, radius, intensity, self.n_freq_step)

        # extract equilibirum radius (average radius)
        eq_radius = util.compute_radial_mean(radius, intensity)

        # extract maximum intensity for normalization purpose
        max_intensity = np.amax(intensity)

        # normalize intensity to 1
        intensity /= max_intensity

        # create TGraph to store radial distribution
        graph = r.TGraph(self.n_freq_step, radius, intensity)
        graph_min = -0.05
        graph_max = 1.05
        style.setTGraphStyle(
            graph, '', 'Radius [mm]', 'Arbitrary units', graph_min, graph_max)

        # define lines to be drawn at collimator apertures (radial space)
        inner_line, outer_line = style.set_collimator_aperture_tline(
            graph_min, graph_max, 'radial')

        # define pave text to go along the radial collimator apertures lines
        pt, pt2 = style.setCollimatorAperturePaveText(
            graph_max*0.45, graph_max*0.55, 'radial')

        # define line to be draw at the magic radius
        magic_line = style.setMagicRadiusTLine(graph_min, graph_max)

        # define pave text to go along the magic radius line
        magic_radius_pave_text = style.setMagicRadiusPaveText(
            graph_max*0.04, graph_max*0.11)

        # compute standard deviation of the radial distribution (within the collimator aperture)
        try:
            std = util.compute_radial_std(radius, intensity, eq_radius, 'ring')
        except:
            print('\n Error -- calculation of the width of the radial distribution failed due to negative number(s) under the square root --> exiting program before completion\n')
            return(0)

        # compute E-field correction
        c_e = util.compute_efield_correction(
            self.field_index, eq_radius-constants.magic_r, std)

        # define pave text with x_e, width and c_e information
        results_pave_text = style.set_radial_results_pave_text(
            eq_radius, std, c_e, graph_max*0.7, graph_max*0.96, 7070, 7095, 'ring')

        # draw it all
        list_to_draw = [graph, magic_line, inner_line, outer_line,
                        pt, pt2, magic_radius_pave_text, results_pave_text]
        plotting.plotMultipleObjects('APL', list_to_draw)
        self.canvas_fourier.Draw()

        # print plot if option specified
        if (self.print_plot >= 2):
            self.canvas_fourier.Print(
                'results/' + self.tag + '/Radial_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(self.opt_t0, self.tS, self.tM, self.freq_step_size))

        # limit radial distribution to collimator aperture
        graph.GetXaxis().SetRangeUser(
            constants.lowerCollimatorRad, constants.upperCollimatorRad)

        # define pave text with x_e, width and CE information
        results_pave_text = style.set_radial_results_pave_text(
            eq_radius, std, c_e, graph_max*0.667, graph_max*0.952, 7070, 7090, 'ring')

        # redraw radial distribution within collimator aperture
        list_to_draw = [graph, magic_line,
                        magic_radius_pave_text, results_pave_text]
        plotting.plotMultipleObjects('APL', list_to_draw)

        # print plot if option specified
        if (self.print_plot >= 2):
            self.canvas_fourier.Print(
                'results/' + self.tag + '/RadialCollAperture_t0_{0:.5f}_tS_{1}_tM_{2}.eps'.format(self.opt_t0, self.tS, self.tM))
            self.canvas_fourier.Print(
                'results/' + self.tag + '/RadialCollAperture_t0_{0:.5f}_tS_{1}_tM_{2}.C'.format(self.opt_t0, self.tS, self.tM))

        # convert from ring global coordinate to beam local coordinate
        util.global_to_local_radial_coordinate(graph)
        graph.GetXaxis().SetRangeUser(-45, +45)
        results_pave_text = style.set_radial_results_pave_text(
            eq_radius-constants.magic_r, std, c_e, graph_max*0.7, graph_max*0.96, -40, -20, 'beam')
        graph.Draw('APL')
        results_pave_text.Draw("same")

        # compare reco/truth level distribution if option specified for simulated data
        if (self.compare_with_truth):

            # retrieve truth radial TGraph
            truth_file = r.TFile(self.truth_root_file)
            truth_graph = truth_file.Get(self.truth_rad_histo_name)

            # convert from ring global coordinate to local beam coordinate
            util.global_to_local_radial_coordinate(truth_graph)
            truth_graph.GetXaxis().SetRangeUser(-45, +45)

            # normalize truth TGraph intensity to 1
            n_point = truth_graph.GetN()
            max_amp = max(truth_graph.GetY())
            truth_intensity, truth_radius = array.array('d'), array.array('d')
            for i in range(1, n_point):
                x, y = r.Double(), r.Double()
                truth_graph.GetPoint(i, x, y)
                truth_graph.SetPoint(i, x, y/max_amp)
                truth_radius.append(x)
                truth_intensity.append(y)

            truth_graph.SetMarkerColor(2)
            truth_graph.SetLineColor(2)
            truth_graph.SetMarkerStyle(20)
            truth_graph.Draw("sameP")

            # compute x_e, width, C_E
            truth_eq_radius = np.average(
                truth_radius, axis=0, weights=truth_intensity)
            truth_std = util.compute_radial_std(
                truth_radius, truth_intensity, truth_eq_radius, 'beam')
            truth_c_e = util.compute_efield_correction(
                self.field_index, truth_eq_radius, truth_std)

            # define and plot the TPaveText that hols truth-level results
            pt5 = style.set_radial_results_pave_text(
                truth_eq_radius, truth_std, truth_c_e, graph_max*0.7, graph_max*0.96, 20, 40, 'truth')
            pt5.Draw("same")

            self.out_file2.cd()
            truth_graph.Write('truth_rad')

        if (self.print_plot >= 1):
            self.canvas_fourier.Print(
                'results/' + self.tag + '/RadialBeamCoordinate_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(self.opt_t0, self.tS, self.tM, self.freq_step_size))

        if (self.print_plot >= 2):
            self.out_file.cd()
            graph.Write('rad')

        # save TH1 histogram to output ROOT file
        self.out_file2.cd()
        graph.Write('rad')

        # convert equilibrium radius from ring to beam coordinate
        eq_radius -= constants.magic_r

        # print results to terminal
        self.print_results(eq_radius, std, c_e, mean_freq, std_freq)

        # save results to a text file
        if (self.append_results):
            results_text_file = open(
                str('results/' + self.tag + '/results.txt'), "a")
        else:
            results_text_file = open(
                str('results/' + self.tag + '/results.txt'), "w")

        results_text_file.write('t0 %f chi2 %f noise %f t0_bkg_th %f bkg_rev_t %f tS %f tM %f df %f fieldIndex %f eq_radius %f std %f c_e %f n_hits %f \n' %
                                (self.opt_t0, chi2, self.bkg_stat_noise, self.t0_background_threshold, self.background_removal_threshold,
                                 self.tS, self.tM, self.freq_step_size, self.field_index, eq_radius, std, c_e, self.n_positron_hits))

    def run(self):

        # print to terminal relevant information
        self.print_info()

        # produce the cosine (sine) transform
        self.produce_cosine_transform()

        # perform the background fit in order to correct for it
        if (self.background_correction == 'fit'):
            # fit he cosine transform background
            try:
                bkg_fit, chi2 = self.fit_background()
            except:
                print(
                    '\n !!ERROR!! Failure to fit the background --> program exiting before completion\n')
                return(0)

            # correct the cosine transform using the fit to the background
            mean_freq, std_freq = self.correct_transform(bkg_fit)

        ''' INTEGRAL BACKGROUND CORRECTION: NOT RELEASED
        if (self.background_correction == 'integral'):
            truth_file = r.TFile(self.truth_root_file)
            approx = truth_file.Get(self.truth_freq_histo_name)
            parabola = self.cosine_histogram.Clone()
            util.calc_parabola( self.opt_t0, self.tS, approx, parabola, self.n_freq_step, self.freq_step_size, self.lower_freq)
            parabola.Draw()
            self.canvas_fourier.Draw()
            self.canvas_fourier.Print('results/' + self.tag + '/parabola.eps')
            a, b = util.minimization(parabola, self.cosine_histogram, self.n_freq_step, self.fit_boundary1, self.fit_boundary2)

            for bin_idx in range(1, self.cosine_histogram.GetNbinsX()+1):
                parabola.SetBinContent( bin_idx, -1*( a*parabola.GetBinContent(bin_idx)+b) )
            parabola.Draw()
            parabola.SetTitle('Integral correction')
            self.canvas_fourier.Draw()
            self.canvas_fourier.Print('results/' + self.tag + '/parabola2.eps')

            for bin_idx in range(1, self.cosine_histogram.GetNbinsX()+1):
                self.cosine_histogram.SetBinContent( bin_idx, self.cosine_histogram.GetBinContent(bin_idx)  + parabola.GetBinContent(bin_idx) )

            chi2 = -1'''

        # convert from frequency to radial distribution
        self.produce_radial(chi2, mean_freq, std_freq)
