
cd('/Users/Alex/Dropbox/KET-PMP-GABZOL/OLD_FreqVer_NotRVE/example');

% Set up config
%--------------------------------------------------------
cfg.Dname = '110511-51_Standard_20120620_VisualCut.ds';
cfg.NC    = [];
cfg.UL    = [];
cfg.fname = 'ica';
cfg.bonf  = 0;
cfg.writelog = 0;
cfg.window   = 20;
cfg.bpf_data = [1 100];
cfg.bpf_eogs = [1 100];

% functions to perform spatial correlations on
%f1 = @(x) max(x')';
%f2 = @(x) min(x')';
%f3 = @(x) PEig(x);
%f4 = @(x) mean(abs(hilbert(x))')';
%funcs = {f1 f2 f3 f4};

funcs = {@(x) PEig(x) @(x) mean(abs(hilbert(x))')'};

cfg.funcs  = funcs;
cfg.dfuncs = {@(x) abs(hilbert(x))};

% Set up thresholds
%--------------------------------------------------------
thr.threshold_raw   = 0;
thr.threshold_comps = 1;
thr.samp_std_over_trials     = 0;
thr.chan_std_over_samples    = 1;
thr.std_per_chan_over_trials = 1;

% Set up plots
%--------------------------------------------------------
plots.review_topo = 0; 
plots.plot_change = 1;
plots.review_spatial = 1;


figure('position',[1000          72         788         906])
FastICA_CTF_SpatioTemp(cfg,thr,plots);