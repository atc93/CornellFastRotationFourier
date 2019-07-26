from scipy.optimize import curve_fit
import src.configparser as configparser
import src.constants as constants
import matplotlib.pyplot as plt
import src.util as util
import numpy as np

# class that optimized the t0 parameter
# it inherits from the ParseConfig class


class Optimize_t0(configparser.ParseConfig):

    # class initialization
    def __init__(self, config, bin_center, bin_content):
        super(Optimize_t0, self).__init__(config)
        self.bin_center = bin_center
        self.bin_content = bin_content

    # method that extracts for a given iteration the optimized t0 value
    # from fitting the chi2 distribution as a function of t0
    def extract_optimum(self, x, y, label):

        # scale from micro-s to ns because ns is the
        # natural unit for t0, will return the t0
        # value in micro-s which is the default unit
        # of the analysis code
        x = [i * 1000 for i in x]

        # fit chi2 distribution with a 4th order polynomial
        # using order higher than 2 is better when scanning
        # over large t0 window that leads to non-quadratic
        # behaviors

        fit = np.polyfit(x, y, 4)

        # create function from fit result
        fit_fn = np.poly1d(fit)

        # extract optimum chi2 by finding the derivative's
        # root constrained to be lies within the t0 window
        # being scanned
        #
        derivative = fit_fn.deriv().r
        t0_roots = derivative[derivative.imag == 0].real
        # factor 1000 for micro-s/ns conversion
        t0_roots = np.ma.masked_where(t0_roots > self.upper_t0*1000, t0_roots)
        # factor 1000 for micro-s/ns conversion
        t0_roots = np.ma.masked_where(t0_roots < self.lower_t0*1000, t0_roots)
        t0_roots = np.ma.compressed(t0_roots)

        # if the 4th order minimzation of the chi2 failed,
        # meaning more than one root was found, go back to
        # a classic 2nd order minization (less accurate for
        # a large t0 window)
        if (len(t0_roots) != 1):
            print('\n    Warning -- problem extracting the optimized t0 from fitting the chi2 vs t0 distribution --> switching to back-up method')

            # fit chi2 distribution with second order polynomial
            fit = np.polyfit(x, y, 2)

            # create function from fit result
            fit_fn = np.poly1d(fit)

            # extract the optimized t0 corresponding to the minimum chi2
            opt_t0 = np.real(fit_fn.r[0])

            # extract corresponding chi2
            opt_chi2 = fit_fn(opt_t0)

        # if the 4th order polynomial fit worked properly
        elif (len(t0_roots) == 1):
            # extract the optimized t0 corresponding to the minimum chi2
            opt_t0 = t0_roots[0]
            # extract corresponding chi2
            opt_chi2 = fit_fn(opt_t0)

        # plot the results
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        plt.subplots_adjust(left=0.125, bottom=0.125,
                            right=0.95, top=0.95, wspace=0, hspace=0)
        plt.plot(x, y, 'ro', ms=8, zorder=2)
        smoothed_x = np.linspace(min(x), max(x), 100)
        plt.plot(smoothed_x, fit_fn(smoothed_x),
                 color='black', linewidth=2.5, zorder=1)
        plt.xlabel('$\mathregular{t_{0}}$ [ns]')
        plt.ylabel('$\mathregular{{\chi}^2}}$/d.o.f.')

        # save plot if specified to do so in config file
        if (self.print_plot >= 1):
            plt.savefig('results/' + self.tag + '/t0_optimization/t0Opt_' + label +
                        '_fit_tS_{0}_tM_{1}.eps'.format(self.tS, self.tM))

        # close the plot
        plt.close()

        # return results (convert t0 back to mciro-s)
        return opt_t0/1000, opt_chi2

    # method that performs the t0 optimization by
    # producing cosine Fourier transforms for many
    # t0 values and fitting their background in
    # order to extract the chi2 and residuals
    def optimization_loop(self, fit_bound1, fit_bound2, bkg_stat_noise):

        # compute how many loops required to scan t0
        n_t0_step = round((self.upper_t0-self.lower_t0)/self.t0_step_size)+1

        # create list to store the results
        t0_list = []
        chi2_list = []
        bkg_stat_noise_list = []
        bound1_list = []
        bound2_list = []

        # print outs
        if (self.verbose > 0):
            print('    scanning over t0 values:\n')

        # loop over the t0 parameters
        for idx_t0 in range(0, n_t0_step):

            # compute the t0 value to use for this loop iteration
            t0 = self.lower_t0+idx_t0*self.t0_step_size

            # compute cosine Fourier transfom and return the results to a list
            freq, intensity = util.calc_cosine_transform(
                t0, self.bin_content, self.bin_center, self.freq_step_size, self.n_freq_step, self.lower_freq)

            # create empty list to store cosine Fourier transform values
            a = []
            b = []

            # set the boundaries of the background window depending
            # if physical and/or non-physical frequencies are allowed
            # via the config file ('physical' or 'all)
            for x, y in zip(freq, intensity):

                # set the background boundaries in case of only
                # physical frequencies allowed
                if (self.background_frequencies == 'physical'):
                    if (fit_bound1 > constants.lowerCollimatorFreq):
                        lower_fit_bound1 = constants.lowerCollimatorFreq
                    else:
                        lower_fit_bound1 = self.lower_freq

                    if (fit_bound2 < constants.upperCollimatorFreq):
                        upper_fit_bound2 = constants.upperCollimatorFreq
                    else:
                        upper_fit_bound2 = self.upper_freq

                # set the background boundaries in case of allowing
                # physical and non-physical frequencies
                elif (self.background_frequencies == 'all'):
                    lower_fit_bound1 = self.lower_freq
                    upper_fit_bound2 = self.upper_freq

                # append background data point to the list according
                # to the aboved defined background boundaries
                if ((x > lower_fit_bound1 and x < fit_bound1) or (x > fit_bound2 and x < upper_fit_bound2)):
                    a.append(x)
                    b.append(y)

            # create list containing the background statistical noise
            err = [bkg_stat_noise]*len(a)

            # create list to store fit parameters
            popt = []

            # fit the background according to the fit function
            # specified in the config file
            if (self.background_fit == 'pol'):
                fit_fn, chi2, sigma = util.fit_pol(
                    a, b, err, self.poly_order, bkg_stat_noise)

            if (self.background_fit == 'sinc'):
                fit_fn, chi2, sigma, popt, pcov = util.fit_sinc(
                    a, b, err, self.rebin_frs_factor, self.verbose)

            if (self.background_fit == 'erfi'):
                fit_fn, chi2, sigma, popt, pcov = util.fit_erfi(
                    a, b, err, self.tS, self.verbose)

            if (self.background_fit == 'triangle'):
                fit_fn, chi2, sigma, popt, pcov = util.fit_triangle(
                    a, b, err, self.tS, self.verbose)

            # skip this lopp iteration because of failed fit
            if (chi2 == -1):
                continue

            # plot frequency distribution alongside the background and its fit
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            plt.subplots_adjust(left=0.125, bottom=0.125,
                                right=0.95, top=0.95, wspace=0, hspace=0)
            plt.plot(freq, intensity, marker='o', ms=5, zorder=1)
            plt.xlabel('Frequency [kHz]')
            plt.ylabel('Arbitrary units')
            plt.errorbar(a, b, yerr=bkg_stat_noise, marker='o', ms=5, markerfacecolor='black',
                         ecolor='black', markeredgecolor='black', linestyle='', label='background', zorder=4)
            fit_label = 'bkg ' + self.background_fit + ' fit'
            plt.plot(freq, fit_fn(freq, *popt),
                     label=fit_label, linewidth=3, zorder=2)
            plt.legend(loc="upper right", frameon=False)

            # show plot if enabled by config file
            if (self.print_plot >= 2):
                plt.savefig('results/' + self.tag + '/t0_optimization/Cosine_t0_{0:.6f}_tS_{1}_tM_{2}_{3}_{4}.eps'.format(
                    t0, self.tS, self.tM, fit_bound1, fit_bound2), format='eps')

            # close plot
            plt.close()

            # optimize the lower fit boundary using background statistical noise and the noise threshold
            for i in range(int(len(freq)/2-1), 0, -1):
                if (abs(intensity[i]-fit_fn(freq[i], *popt)) < self.t0_background_threshold*bkg_stat_noise):
                    opt_bound1 = freq[i]
                    break
                else:
                    opt_bound1 = constants.lowerCollimatorFreq

            # optimize the upper fit boundary using background statistical noise and the noise threshold
            for i in range(int(len(freq)/2), int(len(freq)), +1):
                if (abs(intensity[i]-fit_fn(freq[i], *popt)) < self.t0_background_threshold*bkg_stat_noise):
                    opt_bound2 = freq[i]
                    break
                else:
                    opt_bound2 = constants.upperCollimatorFreq

            # set the boundaries if specified in the config file
            if (self.fix_fit_bound):
                opt_bound1 = self.fit_lower_bound
                opt_bound2 = self.fit_upper_bound

            # print outs
            if (self.verbose > 0):
                print('    t0: {0:.3f}'.format(t0*1000), 'ns\tchi2/dof: {0:.3f}'.format(chi2), '\tnoise: {0:.5f}'.format(sigma),
                      '\topt bound1: {0:.0f}'.format(opt_bound1), 'kHz \topt bound2: {0:.0f}'.format(opt_bound2), 'kHz')

            # append results to lists
            t0_list.append(t0)
            chi2_list.append(chi2)
            bkg_stat_noise_list.append(sigma)
            # rounding important to avoid random last digit like 6733.000000000001
            bound1_list.append(round(opt_bound1, 8))
            # rounding important to avoid random last digit like 6733.000000000001
            bound2_list.append(round(opt_bound2, 8))

        # return results
        return (t0_list, chi2_list, bkg_stat_noise_list, bound1_list, bound2_list)

    # main method that runs the t0 optimization
    # with as many iterations as specified in the
    # config file
    def run_t0_optimization(self):

        print(' ### Step 3/4: optimize t0 (' +
              str(self.n_t0_opt) + ' iterations)\n')

        # iteration the t0 optimization
        for iOpt in range(self.n_t0_opt):

            # first iteration with an inital guess of the background stat noise
            if (iOpt == 0):

                t0_list, chi2_list, bkg_stat_noise_list, bound1_list, bound2_list = self.optimization_loop(
                    constants.lowerCollimatorFreq, constants.upperCollimatorFreq, 0.02)
                opt_t0, opt_chi2 = self.extract_optimum(
                    t0_list, chi2_list, 'Opt_1')
                opt_idx = min(range(len(t0_list)),
                              key=lambda i: abs(t0_list[i]-opt_t0))

                if (self.verbose >= 1):
                    print('\n    --> t0 optimization, iteration #' + str(iOpt + 1) + ' done, t0 = ' +
                          '{0:.3f}'.format(opt_t0*1000) + ' ns\n')

            # consecutive iterations using the measured background stat noise
            else:

                t0_list, chi2_list, bkg_stat_noise_list, bound1_list, bound2_list = self.optimization_loop(
                    bound1_list[opt_idx], bound2_list[opt_idx], bkg_stat_noise_list[opt_idx])
                opt_t0, opt_chi2 = self.extract_optimum(
                    t0_list, chi2_list, 'Opt_' + str(iOpt+1))
                opt_idx = min(range(len(t0_list)),
                              key=lambda i: abs(t0_list[i]-opt_t0))
                if (self.verbose >= 1):
                    print('\n    --> t0 optimization, iteration #' + str(iOpt+1) + ' done, t0 = ' +
                          '{0:.3f}'.format(opt_t0*1000) + ' ns\n')

        # return results
        return opt_t0, opt_chi2, bkg_stat_noise_list[opt_idx], bound1_list[opt_idx], bound2_list[opt_idx]
