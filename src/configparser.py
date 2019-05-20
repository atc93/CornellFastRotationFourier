class ParseConfig:

    def __init__(self, config):
        self.lower_t0 = config["lower_t0"]
        self.upper_t0 = config["upper_t0"]
        self.t0_step_size = config["t0_step_size"]
        self.tS = config["tS"]
        self.tM = config["tM"]
        self.poly_order = config["poly_order"]
        self.print_plot = config["print_plot"]
        self.tag = config["tag"]
        self.verbose = config["verbose"]
        self.start_fit_time = config["start_fit_time"]
        self.n_fit_param = config["n_fit_param"]
        self.rebin_wiggle_factor = config["rebin_wiggle_factor"]
        self.rebin_frs_factor = config["rebin_frs_factor"]
        self.histo_name = config['histo_name']
        self.truth_histo_name = config['truth_histo_name']
        self.stat_fluctuation = config['stat_fluctuation']
        self.n_t0_opt = config['n_t0_opt']
        self.field_index = config['field_index']
        self.calc_sine = config['calc_sine']
        self.append_results = config['append_results']
        self.root_file = config['root_file']
        self.truth_root_file = config['truth_root_file']
        self.freq_step_size = config['freq_step_size']
        self.lower_freq = config['lower_freq']
        self.upper_freq = config['upper_freq']
        self.compare_with_truth = config['compare_with_truth']
        self.n_freq_step = int(
            (self.upper_freq - self.lower_freq) / self.freq_step_size)
        self.t0_background_threshold = config['t0_background_threshold']
        self.background_removal_threshold = config['background_removal_threshold']
        self.run_background_threshold_scan = config['run_background_threshold_scan']
        self.lower_background_threshold = config['lower_background_threshold']
        self.upper_background_threshold = config['upper_background_threshold']
        self.background_threshold_step_size = config['background_threshold_step_size']
        self.n_stat_fluc = config['n_stat_fluctuation']
        self.fix_t0 = config['fix_t0']
        self.fixed_t0_value = config['fixed_t0_value']
        self.background_correction = config['background_correction']
        self.remove_background = config['remove_background']
        self.background_frequencies = config['background_frequencies']
        self.background_fit = config['background_fit']
        self.fit_lower_bound = config['fit_lower_bound']
        self.fit_upper_bound = config['fit_upper_bound']
        self.fix_fit_bound = config['fix_fit_bound']
