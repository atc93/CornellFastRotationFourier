import src.configparser as configparser
import src.constants as constants
import src.plotting as plotting
import src.style as style
import src.util as util
import numpy as np
import ROOT as r
import array

class Fourier(configparser.ParseConfig):

    def __init__(self, config, bin_center, bin_content, opt_t0, noise_sigma, fit_boundary1, fit_boundary2):
        super(Fourier, self).__init__(config)
        self.bin_center = bin_center
        self.bin_content = bin_content
        self.opt_t0 = opt_t0
        self.fit_boundary1 = fit_boundary1
        self.fit_boundary2 = fit_boundary2
        self.noise_sigma = noise_sigma

        self.canvas = r.TCanvas( 'c_results', 'c_results', 900, 600 )
        self.cosine_histogram = r.TH1D( 'cosine', 'cosine', self.n_freq_step, self.lower_freq, self.upper_freq )
        self.sine_histogram = r.TH1D( 'sine', 'sine', self.n_freq_step, self.lower_freq, self.upper_freq )

        style.setTCanvasStyle(self.canvas)

    def print_info(self):

        print('\n ### Step 4/4: produce results\n')
        print ( '    t0 = ', self.opt_t0 )
        print ( '    tS = ', self.tS )
        print ( '    tM = ', self.tM )
        print ( '    df = ', self.freq_step_size )
        print ( '    n  = ', self.field_index )

    def produce_cosine_transform(self):

        # compute cosine Fourier transform
        freq, intensity = util.calc_cosine_transform(self.opt_t0, self.bin_content, self.bin_center, self.freq_step_size, self.n_freq_step, self.lower_freq)

        # fill ROOT histogram
        for x, y in zip(freq, intensity):
            iBin = self.cosine_histogram.FindBin(x)
            self.cosine_histogram.SetBinContent(iBin, y)

        # compute sine Fourier transform
        if (self.calc_sine):
            freq, intensity = util.calc_sine_transform(self.opt_t0, self.bin_content, self.bin_center)

            for x, y in zip(freq, intensity):
                iBin = self.sine_histogram.FindBin(x)
                self.sine_histogram.SetBinContent(iBin, y)

        # define histogram(s) name
        histoName = [ 'Cosine transform: t_{0} = ' + '{0:.2f} ns'.format(self.opt_t0*1000), 'Sine transform: t_{0} = ' + '{0:.2f} ns'.format(self.opt_t0) ]

        # clone cosine histogram for plotting/styling puporses
        cosine_clone = self.cosine_histogram.Clone()

        # clone sine histogram for plotting/styling puporses
        sine_clone = self.sine_histogram.Clone()

        # create a list containing both cosine and sine histograms
        clone_hist_list = [ cosine_clone, sine_clone ]

        # styling and plotting for cosine and sine Fourier transform
        for idx in range ( 0, 2 ):

            # style the Fourier transform histogram ==#
            style.setTH1Style( clone_hist_list[ idx ], histoName[ idx ], 'Frequency [kHz]', 'Arbitrary', 1.2, 1.3 )

            # define lines to be drawn at collimator apertures (frequency space) ==#
            inner_line, outer_line = style.setCollimatorApertureTLine ( clone_hist_list[ idx ].GetMinimum(), clone_hist_list[ idx ].GetMaximum(), 'frequency' )

            # define pave text to go along the collimator apertures lines ==#
            pt, pt2 = style.setCollimatorAperturePaveText( clone_hist_list[ idx ].GetMaximum()*0.38, clone_hist_list[ idx ].GetMaximum()*0.52, 'frequency' )

            # draw it all
            list_to_draw = [ clone_hist_list[ idx ], inner_line, outer_line, pt, pt2 ]
            plotting.plotMultipleObjects( '', list_to_draw )
            self.canvas.Draw()

            # save plot if option provided
            if ( idx == 0 and self.print_plot == 1 ):
                self.canvas.Print('results/'+ self.tag + '/Cosine_t0_{0:.5f}_tS_{1}_tM_{2}_df_{3}.eps'.format(self.opt_t0, self.tS, self.tM, self.freq_step_size))

            # save plot if option provided and if sine transform was performed
            if ( idx == 1 and self.print_plot == 1 and self.calc_sine == 1 ):
                self.canvas.Print('results/' + self.tag + '/Sine_t0_{0:.3f}_tS_{1}_tM_{2}_df_{3}.eps'.format(self.opt_t0, self.tS, self.tM, self.freq_step_size))


    def fit_background(self):

        a = []
        b = []

        for bin_idx in range (1, self.cosine_histogram.GetNbinsX()+1):
            if ( self.cosine_histogram.GetBinCenter(bin_idx) < self.fit_boundary1 or self.cosine_histogram.GetBinCenter(bin_idx) > self.fit_boundary2 ):
                a.append( self.cosine_histogram.GetBinCenter(bin_idx))
                b.append( self.cosine_histogram.GetBinContent(bin_idx))

        # create list with error values for the intensity
        err = [self.noise_sigma]*len(a)

        fit = np.polyfit(a, b, self.poly_order, w=err)
        
        # compute chi2
        chi2 = np.sum((np.polyval(fit, a) - b) ** 2 / self.noise_sigma ** 2)/(len(a)-self.poly_order)

        return fit, chi2

    def correct_transform(self, fit):

        max_amplitude_bin_idx = self.cosine_histogram.GetMaximumBin()

        zero_out = False
        for bin_idx in range(max_amplitude_bin_idx, 0, -1):
            if ( ( self.cosine_histogram.GetBinContent(bin_idx) - np.polyval(fit, self.cosine_histogram.GetBinCenter(bin_idx)) ) < constants.noiseThreshold*self.noise_sigma  ):
                zero_out = True
            if (zero_out):
                self.cosine_histogram.SetBinContent(bin_idx, 0)
            else:
                self.cosine_histogram.SetBinContent(bin_idx, self.cosine_histogram.GetBinContent(bin_idx)-np.polyval(fit, self.cosine_histogram.GetBinCenter(bin_idx)))

        zero_out = False
        for bin_idx in range(max_amplitude_bin_idx+1, self.cosine_histogram.GetNbinsX()+1, +1):
            if ( ( self.cosine_histogram.GetBinContent(bin_idx) - np.polyval(fit, self.cosine_histogram.GetBinCenter(bin_idx)) ) < constants.noiseThreshold*self.noise_sigma  ):
                zero_out = True
            if (zero_out):
                self.cosine_histogram.SetBinContent(bin_idx, 0)
            else:
                self.cosine_histogram.SetBinContent(bin_idx, self.cosine_histogram.GetBinContent(bin_idx)-np.polyval(fit, self.cosine_histogram.GetBinCenter(bin_idx)))

        #== Style the approximated frequency distribution histogram ==#
        style.setTH1Style( self.cosine_histogram, '', 'Frequency [kHz]', 'Arbitrary units', 1, 1.3 )
        self.cosine_histogram.SetMinimum( -0.5 )

        #== Define lines to be drawn at collimator apertures (frequency space) ==#
        inner_line, outer_line = style.setCollimatorApertureTLine ( self.cosine_histogram.GetMinimum(), self.cosine_histogram.GetMaximum(), 'frequency' )
        #== Define pave text to go along the collimator apertures lines ==#
        pt, pt2 = style.setCollimatorAperturePaveText( self.cosine_histogram.GetMaximum()*0.38, self.cosine_histogram.GetMaximum()*0.52, 'frequency' )

        #== Draw it all ==#
        list_to_draw = [ self.cosine_histogram, inner_line, outer_line, pt, pt2 ]
        plotting.plotMultipleObjects( '', list_to_draw )
        self.canvas.Draw()

        if ( self.print_plot == 1 ):
            self.canvas.Print('results/' + self.tag + '/CorrectedCosine_t0_{0:.5f}_tS_{1}_tM_{2}df_{3}.eps'.format(self.opt_t0, self.tS, self.tM, self.freq_step_size))

    def produce_radial(self, chi2):

        #== Define arrays containing radius and intensity information ==#
        #== The array will be used to produce a TGraph (required since non equidistance radial points) ==#
        intensity, radius = array.array( 'd' ), array.array( 'd' )

        #== Fill radius and intensity arrays from the frequency distribution ==#
        util.convert_freq_to_radius( self.cosine_histogram, radius, intensity, self.n_freq_step)

        #== Extract equilibirum radius (average radius) ==#
        eq_radius = util.computeRadialMean( radius, intensity )

        #== Extract maximum intensity for normalization purpose ==#
        max_intensity = np.amax(intensity)

        #== Normalize intensity to 1 ==#
        intensity /= max_intensity

        #== Create TGraph to store radial distribution ==#
        graph    = r.TGraph( self.n_freq_step, radius, intensity )
        graph_min = -0.05
        graph_max = 1.05
        style.setTGraphStyle( graph, '', 'Radius [mm]', 'Arbitrary units', graph_min, graph_max )

        #== Define lines to be drawn at collimator apertures (frequency space) ==#
        inner_line, outer_line = style.setCollimatorApertureTLine( graph_min, graph_max, 'radial' )

        #== Define pave text to go along the radial collimator apertures lines ==#
        pt, pt2 = style.setCollimatorAperturePaveText( graph_max*0.45, graph_max*0.55, 'radial' )

        #== Define line to be draw at the magic radius ==#
        magic_line = style.setMagicRadiusTLine( graph_min, graph_max)

        #== Define pave text to go along the magic radius line ==#
        magic_radius_pave_text = style.setMagicRadiusPaveText( graph_max*0.04, graph_max*0.11 );

        #== Compute Standard Deviation of the radial distribution (within the collimator aperture) ==#
        std = util.computeRadialSTD( radius, intensity, eq_radius )

        #== Compute E-field correction
        c_e = util.computeEfieldCorrection( self.field_index, eq_radius-constants.magicR, std )

        #== Define pave text with x_e, width and c_e information ==#
        results_pave_text = style.setRadialResultsPaveText( eq_radius, std, c_e, graph_max*0.7, graph_max*0.96, 7070, 7095, 'ring' )

        #== Draw it all ==#
        list_to_draw = [ graph, magic_line, inner_line, outer_line, pt, pt2, magic_radius_pave_text, results_pave_text ]
        plotting.plotMultipleObjects( 'APL', list_to_draw )

        self.canvas.Draw()

        if ( self.print_plot == 1 ):
            self.canvas.Print('results/' + self.tag + '/Radial_t0_{0:.5f}_tS_{1}_tM_{2}.eps'.format(self.opt_t0, self.tS, self.tM))
            self.canvas.Print('results/' + self.tag + '/Radial_t0_{0:.5f}_tS_{1}_tM_{2}.C'.format(self.opt_t0, self.tS, self.tM))

        #== Limit radial distribution to collimator aperture ==#
        graph.GetXaxis().SetRangeUser( constants.lowerCollimatorRad, constants.upperCollimatorRad )

        #== Define pave text with x_e, width and CE information ==#
        results_pave_text = style.setRadialResultsPaveText( eq_radius, std, c_e, graph_max*0.667, graph_max*0.952, 7070, 7090, 'ring' )

        #== Redraw radial distribution within collimator aperture ==#
        list_to_draw = [ graph, magic_line, magic_radius_pave_text, results_pave_text ]
        plotting.plotMultipleObjects( 'APL', list_to_draw )

        if ( self.print_plot == 1 ):
            self.canvas.Print('results/' + self.tag + '/RadialInAperture_t0_{0:.5f}_tS_{1}_tM_{2}.eps'.format(self.opt_t0, self.tS, self.tM))
            self.canvas.Print('results/' + self.tag + '/RadialInAperture_t0_{0:.5f}_tS_{1}_tM_{2}.C'.format(self.opt_t0, self.tS, self.tM))

        #== Convert from ring global coordinate to beam local coordinate ==#
        util.globalToLocalRadialCoordinate( graph )
        graph.GetXaxis().SetRangeUser( -45, +45 )
        results_pave_text = style.setRadialResultsPaveText( eq_radius-constants.magicR, std, c_e, graph_max*0.7, graph_max*0.96, -40, -20, 'beam' )
        graph.Draw('APL')
        results_pave_text.Draw("same")
        if ( self.print_plot == 1 ):
            self.canvas.Print('results/' + self.tag + '/RadialBeamCoordinate_t0_{0:.5f}_tS_{1}_tM_{2}.eps'.format(self.opt_t0, self.tS, self.tM))

        # convert equilibrium radius from ring to beam coordinate
        eq_radius -= constants.magicR

        # save results to a text file
        if (self.append_results):
            results_text_file = open(str('results/' + self.tag + '/results.txt'), "a")
        else:
            results_text_file = open(str('results/' + self.tag + '/results.txt'), "w")

        results_text_file.write('t0 %f chi2 %f noise %f tS %f tM %f df %f fieldIndex %f eq_radius %f std %f c_e %f \n' %
            (self.opt_t0, chi2, self.noise_sigma, self.tS, self.tM, self.freq_step_size, self.field_index, eq_radius, std, c_e) )


    def run(self):

        # print to terminal relevant information
        self.print_info()

        # produce the cosine (sine) transform
        self.produce_cosine_transform()

        # fit he cosine transform background
        bkg_fit, chi2 = self.fit_background()

        # correct the cosine transform using the fit to the background
        self.correct_transform(bkg_fit)

        # convert from frequency to radial distribution
        self.produce_radial(chi2)

