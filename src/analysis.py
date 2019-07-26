import src.t0optimization as t0optimization
import src.fastrotation as fastrotation
import src.constants as constants
import src.fourier as fourier
import src.util as util
import ROOT as r
import sys


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
                             config['tS_step_size'], 6)

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

        save_default_tM = config['tM']

        config['tS'] = round(opt_tS, 6)
        config['tM'] = round(opt_tM, 6)

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        if (config['fix_t0']):
            opt_t0 = config['fixed_t0_value']

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del fr, t0, results

        config['tM'] = save_default_tM


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
        opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

        save_default_tM = config['tM']
        save_default_tS = config['tS']

        config['tS'] = round(opt_tS, 6)
        config['tM'] = round(opt_tM, 6)

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        if (config['fix_t0']):
            opt_t0 = config['fixed_t0_value']

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del fr, t0, results

        config['tM'] = save_default_tM
        config['tS'] = save_default_tS


def run_t0_scan(config):

    # compute number of iterations for t0 scan
    n_t0_step = round(
        (config['upper_t0']-config['lower_t0'])/config['t0_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    # disable saving plots
    config['print_plot'] = False

    # instantiate fast rotation class
    fr = fastrotation.FastRotation(config)
    # class method returns numpy arrays of the fast rotation signal
    opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

    config['tS'] = round(opt_tS, 6)
    config['tM'] = round(opt_tM, 6)

    # instantiate t0 optimization class
    t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
    # iterative optimization of t0 (2 iterations are usually enough)
    opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

    if (config['fix_t0']):
        opt_t0 = config['fixed_t0_value']

    # re-enable saving plots
    config['print_plot'] = True

    del fr, t0

    # loop over the t0 parameters
    for idx_t0 in range(0, n_t0_step):

        # produce results
        results = fourier.Fourier(config, bin_center, bin_content, round(
            config['lower_t0']+idx_t0 * config['t0_step_size'], 6), noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del results


def run_freq_step_scan(config):

    # compute number of iterations for frequency step size scan
    n_freq_step = round(
        (config['upper_freq_step_size']-config['lower_freq_step_size'])/config['freq_step_size_increment'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_freq in range(0, n_freq_step):

        # set freq step size
        config['freq_step_size'] = round(config['lower_freq_step_size']+idx_freq *
                                         config['freq_step_size_increment'], 3)

        # check that the number of frequency bins is an integer factor of the freq step size
        freq_window = int(config['upper_freq']-config['lower_freq'])
        if not ((freq_window/config['freq_step_size']).is_integer()):
            continue

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

        save_default_tM = config['tM']
        save_default_tS = config['tS']

        config['tS'] = round(opt_tS, 6)
        config['tM'] = round(opt_tM, 6)

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        if (config['fix_t0']):
            opt_t0 = config['fixed_t0_value']

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del fr, t0, results

        config['tM'] = save_default_tM
        config['tS'] = save_default_tS


def run_background_threshold_scan(config):

    # compute number of iterations for frequency step size scan
    n_noise_step = round(
        (config['upper_background_threshold']-config['lower_background_threshold'])/config['background_threshold_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    # loop over the t0 parameters
    for idx_noise in range(0, n_noise_step):

        # set noise threshold
        config['t0_background_threshold'] = round(config['lower_background_threshold']+idx_noise *
                                                  config['background_threshold_step_size'], 3)

        # disable saving plots
        config['print_plot'] = False

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

        save_default_tM = config['tM']
        save_default_tS = config['tS']

        config['tS'] = round(opt_tS, 6)
        config['tM'] = round(opt_tM, 6)

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        if (config['fix_t0']):
            opt_t0 = config['fixed_t0_value']

        # re-enable saving plots
        config['print_plot'] = True

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del fr, t0, results

        config['tM'] = save_default_tM
        config['tS'] = save_default_tS


def run_background_removal_threshold_scan(config):

    # compute number of iterations for frequency step size scan
    n_step = round(
        (config['upper_background_removal_threshold']-config['lower_background_removal_threshold'])/config['background_removal_threshold_step_size'])+1

    # append results in output text file
    config['append_results'] = True

    config['remove_background'] = True

    # disable saving plots
    config['print_plot'] = False

    # instantiate fast rotation class
    fr = fastrotation.FastRotation(config)
    # class method returns numpy arrays of the fast rotation signal
    opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

    save_default_tM = config['tM']
    save_default_tS = config['tS']

    config['tS'] = round(opt_tS, 6)
    config['tM'] = round(opt_tM, 6)

    # instantiate t0 optimization class
    t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
    # iterative optimization of t0 (2 iterations are usually enough)
    opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

    if (config['fix_t0']):
        opt_t0 = config['fixed_t0_value']

    # re-enable saving plots
    config['print_plot'] = True

    # loop over the t0 parameters
    for idx_noise in range(0, n_step):

        # set noise threshold
        config['background_removal_threshold'] = round(
            config['lower_background_removal_threshold']+idx_noise * config['background_removal_threshold_step_size'], 3)

        # produce results
        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        del results

    del fr, t0


def run_stat_fluc(config):

    # append results in output text file
    config['append_results'] = True

    # disable saving plots
    config['print_plot'] = False

    for i in range(config['n_stat_fluctuation']):

        # instantiate fast rotation class
        fr = fastrotation.FastRotation(config)
        # class method returns numpy arrays of the fast rotation signal
        opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

        save_default_tM = config['tM']
        save_default_tS = config['tS']

        config['tS'] = round(opt_tS, 6)
        config['tM'] = round(opt_tM, 6)

        # instantiate t0 optimization class
        t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
        # iterative optimization of t0 (2 iterations are usually enough)
        opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

        if (config['fix_t0']):
            opt_t0 = config['fixed_t0_value']

        config['print_plot'] = True

        results = fourier.Fourier(
            config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
        results.run()

        config['print_plot'] = False

        del fr, t0, results

        config['tM'] = save_default_tM
        config['tS'] = save_default_tS


def run_default(config):

    # instantiate fast rotation class
    fr = fastrotation.FastRotation(config)
    # class method returns numpy arrays of the fast rotation signal
    opt_tS, opt_tM, n_positron_hits, bin_center, bin_content = fr.produce()

    config['tS'] = round(opt_tS, 6)
    config['tM'] = round(opt_tM, 6)

    # instantiate t0 optimization class
    t0 = t0optimization.Optimize_t0(config, bin_center, bin_content)
    # iterative optimization of t0 (2 iterations are usually enough)
    opt_t0, chi2, noise, fit_boundary1, fit_boundary2 = t0.run_t0_optimization()

    if (config['fix_t0']):
        opt_t0 = config['fixed_t0_value']

    results = fourier.Fourier(
        config, bin_center, bin_content, opt_t0, noise, fit_boundary1, fit_boundary2, n_positron_hits)
    results.run()


def run_fourier(config):

    # return error if more than one scan is set to true
    scans = (config['run_tS_scan'], config['run_tM_scan'], config['run_t0_scan'], config['run_freq_step_scan'],
             config['run_background_threshold_scan'], config['run_background_removal_threshold_scan'], config['stat_fluctuation'])
    count_scan = sum(1 for count in scans if count)

    if (count_scan > 1):
        print(' ERROR: more than one scan is configured to run. Only one scan at a time is allowed --> exiting')
        sys.exit(0)

    if (config['run_tS_scan']):
        run_tS_scan(config)
    elif (config['run_tM_scan']):
        run_tM_scan(config)
    elif (config['run_t0_scan']):
        run_t0_scan(config)
    elif (config['run_freq_step_scan']):
        run_freq_step_scan(config)
    elif (config['run_background_threshold_scan']):
        run_background_threshold_scan(config)
    elif (config['run_background_removal_threshold_scan']):
        run_background_removal_threshold_scan(config)
    elif (config['stat_fluctuation']):
        run_stat_fluc(config)
    else:
        run_default(config)
