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
        self.rebin_factor = config["rebin_factor"]
        self.histo_name = config['histo_name']
        self.stat_fluctuation = config['stat_fluctuation']
        self.n_t0_opt = config['n_t0_opt']
        self.field_index = config['field_index']
        self.calc_sine = config['calc_sine']
        self.append_results = config['append_results']
        self.root_file = config['root_file']
        self.freq_step_size = config['freq_step_size']
        self.lower_freq = config['lower_freq']
        self.upper_freq = config['upper_freq']
        self.n_freq_step = int(
            (self.upper_freq - self.lower_freq) / self.freq_step_size)
