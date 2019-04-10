import src.t0optimization as t0optimization
import src.fastrotation as fastrotation
import src.constants as constants
import src.fourier as fourier
import src.util as util
import ROOT as r


def run_tS_scan(config):

    # compute number of iterations for tS scan
    n_tS_step = round(
        (config['upper_tS']-config['lower_tS'])/config['tS_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_tS in range(0, n_tS_step):

        # set tS
        config['tS'] = round(config['lower_tS']+idx_tS *
                             config['tS_step_size'], 3)

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        bin_center, bin_content = fr.produce()

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2)
        results.run()

        del fr, t0, results


def run_tM_scan(config):

    # compute number of iterations for tM scan
    n_tM_step = round(
        (config['upper_tM']-config['lower_tM'])/config['tM_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_tM in range(0, n_tM_step):

        # set tM
        config['tM'] = round(config['lower_tM']+idx_tM *
                             config['tM_step_size'], 3)

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        bin_center, bin_content = fr.produce()

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2)
        results.run()

        del fr, t0, results


def run_t0_scan(config):

    # compute number of iterations for t0 scan
    n_t0_step = round(
        (config['upper_t0']-config['lower_t0'])/config['t0_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_t0 in range(0, n_t0_step):

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        bin_center, bin_content = fr.produce()

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(config, bin_center, bin_content, round(
            config['lower_t0']+idx_t0 * config['t0_step_size'], 6), noise, fit_boundary1, fit_boundary2)
        results.run()

        del fr, t0, results


def run_freq_step_scan(config):

    # compute number of iterations for frequency step size scan
    n_freq_step = round(
        (config['upper_freq_step_size']-config['lower_freq_step_size'])/config['freq_step_size_increment'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_freq in range(0, n_freq_step):

        # set tM
        config['freq_step_size'] = round(config['lower_freq_step_size']+idx_freq *
                             config['freq_step_size_increment'], 3)

        print(config['freq_step_size'])

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        bin_center, bin_content = fr.produce()

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2)
        results.run()

        del fr, t0, results


def run_default(config):

    # instantiate fast rotation class
    fr = fastrotation.FastRotation(config)
    # class method returns numpy arrays of the fast rotation signal
    bin_center, bin_content = fr.produce()

    # instantiate t0 optimization class
    t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
    # iterative optimization of t0 (2 iterations are usually enough)
    opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()
    # print(opt_t0, ' ', chi2, ' ', noise, ' ', fit_boundary1, ' ', fit_boundary2)

    results = fourier.Fourier(
        config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2)
    results.run()


def run_fourier(config):

    if (config['run_tS_scan']):
        run_tS_scan(config)
    elif (config['run_tM_scan']):
        run_tM_scan(config)
    elif (config['run_t0_scan']):
        run_t0_scan(config)
    elif (config['run_freq_step_scan']):
        run_freq_step_scan(config)
    else:
        run_default(config)
