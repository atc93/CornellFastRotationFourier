import src.configparser as configparser
import scipy.optimize
import src.constants as constants
import matplotlib.pyplot as plt
import src.util as util
import numpy as np


class Optimize_t0(configparser.ParseConfig):

    def __init__(self, config, bin_center, bin_content):
        super(Optimize_t0, self).__init__(config)
        self.bin_center = bin_center
        self.bin_content = bin_content

    def extract_optimum(self, x, y, label):

        # scale from mu-s to ns
        x = [i * 1000 for i in x]

        # fit chi2 distribution with 5 order polynomial
        # using order higher than 2 in case of scanning
        # over large t0 window leading to non-quadratic
        # behaviors
        fit = np.polyfit(x, y, 5)

        # create function from fit result
        fit_fn = np.poly1d(fit)

        # extract optimum chi2 by finding the derivative
        # root that lies within the t0 window
        derivative = fit_fn.deriv().r
        t0_roots = derivative[derivative.imag==0].real
        t0_roots = np.ma.masked_where(t0_roots > self.upper_t0*1000, t0_roots)
        t0_roots = np.ma.masked_where(t0_roots < self.lower_t0*1000, t0_roots)
        t0_roots = np.ma.compressed(t0_roots)

        # if the 5th order minimzation of the chi2 failed
        # go back to a classic 2nd order minization (less
        # accurate for a large t0 window)
        if (len(t0_roots) != 1):
            print('   Error in extracting optimized t0 from fitting chi2 vs t0 distribution --> switching to back-up method')

            # fit chi2 distribution with second order polynomial
            fit = np.polyfit(x, y, 2)

            # create function from fit result
            fit_fn = np.poly1d(fit)

            # extract the optimized t0 from the chi2 minimization
            opt_t0 = np.real(fit_fn.r[0])

            # extract corresponding optimum chi2
            opt_chi2 = fit_fn(opt_t0)

        elif (len(t0_roots) == 1):
            opt_t0 = t0_roots[0]
            opt_chi2 = fit_fn(opt_t0)

        # plot the results
        plt.plot(x, y, 'ro')
        smoothed_x = np.linspace(min(x), max(x), 100)
        plt.plot(smoothed_x, fit_fn(smoothed_x), color='black')
        plt.xlabel('$\mathregular{t_{0}}$ [ns]')
        plt.ylabel('$\mathregular{{\chi}^2}}$/d.o.f.')

        # save plot if specified to do so in config file
        if (self.print_plot):
            plt.savefig('results/' + self.tag + '/t0_optimization/t0Opt_' + label +
                        '_fit_tS_{0}_tM_{1}.eps'.format(self.tS, self.tM))

        # close the plot
        plt.close()

        return opt_t0/1000, opt_chi2 # convert back to mu-s

    def optimization_loop(self, fit_bound1, fit_bound2, noise_sigma):

        # compute how many loops required to scan t0
        n_t0_step = round((self.upper_t0-self.lower_t0)/self.t0_step_size)+1

        # create list to store the results
        t0_list = []
        chi2_list = []
        noise_list = []
        bound1_list = []
        bound2_list = []

        # print outs
        if (self.verbose > 0):
            print('    scanning over t0 values:\n')

        # loop over the t0 parameters
        for idx_t0 in range(0, n_t0_step):

            # compute the t0 value to use
            t0 = self.lower_t0+idx_t0*self.t0_step_size

            # compute cosine Fourier transfom and return the results to a list
            freq, intensity = util.calc_cosine_transform(
                t0, self.bin_content, self.bin_center, self.freq_step_size, self.n_freq_step, self.lower_freq)

            a = []
            b = []

            for x, y in zip(freq, intensity):
                if (x < fit_bound1 or x > fit_bound2):
                    a.append(x)
                    b.append(y)

            # create list with error values for the intensity
            err = [noise_sigma]*len(a)

            # fit the background
            fit = np.polyfit(a, b, self.poly_order, w=err)

            # compute chi2
            chi2 = np.sum((np.polyval(fit, a) - b) ** 2 /
                          noise_sigma ** 2)/(len(a)-self.poly_order)

            # create function from fit results
            fit_fn = np.poly1d(fit)

            # plot frequency distribution alongside the background and its fit
            plt.plot(freq, intensity, marker='o', ms=2)
            plt.plot(freq, fit_fn(freq))
            plt.errorbar(a, b, yerr=noise_sigma, marker='o', ms=2, markerfacecolor='black',
                         ecolor='black', markeredgecolor='black', linestyle='')

            # show plot if enabled by config file
            if (self.print_plot):
                plt.savefig('results/' + self.tag + '/t0_optimization/Cosine_t0_{0:.4f}_tS_{1}_tM_{2}_{3}_{4}.eps'.format(
                    t0, self.tS, self.tM, fit_bound1, fit_bound2), format='eps')

            # close plot
            plt.close()

            # estimate the noise from fit residuals
            residuals = []
            for i in range(len(a)):
                residuals.append(np.polyval(fit, a[i])-b[i])

            # optimize the fit boundaries using noise_sigma and the noise threshold
            for i in range(int(len(freq)/2-1), 0, -1):
                if (intensity[i]-np.polyval(fit, freq[i]) < constants.noiseThreshold*noise_sigma):
                    opt_bound1 = freq[i]
                    break
            for i in range(int(len(freq)/2), int(len(freq)), +1):
                if (intensity[i]-np.polyval(fit, freq[i]) < constants.noiseThreshold*noise_sigma):
                    opt_bound2 = freq[i]
                    break

            # print outs
            if (self.verbose > 0):
                print('    t0: {0:.3f}'.format(t0*1000), 'ns\tchi2/dof: {0:.3f}'.format(chi2), '\tnoise: {0:.5f}'.format(np.std(residuals)),
                      '\topt bound1: {0:.0f}'.format(opt_bound1), 'kHz \topt bound2: {0:.0f}'.format(opt_bound2), 'kHz')

            # append results to lists
            t0_list.append(t0)
            chi2_list.append(chi2)
            noise_list.append(np.std(residuals))
            bound1_list.append(opt_bound1)
            bound2_list.append(opt_bound2)

        return (t0_list, chi2_list, noise_list, bound1_list, bound2_list)

    def run_t0_optimization(self):

        print(' ### Step 3/4: optimize t0 (' + str(self.n_t0_opt) + ' iterations)\n')

        for iOpt in range(self.n_t0_opt):

            if (iOpt == 0):

                t0_list, chi2_list, noise_list, bound1_list, bound2_list = self.optimization_loop(
                        constants.lowerCollimatorFreq, constants.upperCollimatorFreq, 0.02)
                opt_t0, opt_chi2 = self.extract_optimum(t0_list, chi2_list, 'Opt_1')
                opt_idx = min(range(len(t0_list)),
                key=lambda i: abs(t0_list[i]-opt_t0))
                if ( self.verbose > 0 ):
                    print('')
                print('    t0 optimization, iteration #' + str(iOpt + 1) + ' done, t0 = ' +
            '{0:.3f}'.format(opt_t0*1000) + ' ns')
                if ( self.verbose > 0 ):
                    print('')

            else:

                t0_list, chi2_list, noise_list, bound1_list, bound2_list = self.optimization_loop(
                bound1_list[opt_idx], bound2_list[opt_idx], noise_list[opt_idx])
                opt_t0, opt_chi2 = self.extract_optimum(t0_list, chi2_list, 'Opt_' + str(iOpt+1))
                opt_idx = min(range(len(t0_list)),
                key=lambda i: abs(t0_list[i]-opt_t0))
                if ( self.verbose > 0 ):
                    print('')
                print('    t0 optimization, iteration #' + str(iOpt+1) + ' done, t0 = ' +
            '{0:.3f}'.format(opt_t0*1000) + ' ns')
                if ( self.verbose > 0 ):
                    print('')

        return opt_t0, opt_chi2, noise_list[opt_idx], bound1_list[opt_idx], bound2_list[opt_idx]