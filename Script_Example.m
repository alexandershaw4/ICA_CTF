
cd('/cubric/scratch/sapas10/WAND/Visual/314_09540/');

% Set up config
%--------------------------------------------------------
cfg.Dname = '314_09540_Visual.ds'; % name of the raw (unepoiched) CTF ds
cfg.NC    = 128;             % number of components to look for
cfg.UL    = [];             % upper limit on num to reject     
cfg.fname = '_icanew';      % output filename appendix
cfg.bonf  = 0;              % Bonferonni correct correlations
cfg.writelog = 0;           % flag: write a logfile or to command window
cfg.window   = 4;          % length of windows to compute ICA
cfg.overlap  = cfg.window;  % window overlap (set to len window == no overlap)
cfg.bpf_data = [1 250];     % bandpass data
cfg.bpf_eogs = [1 20];     % bandpass EOG / peripheral

% functions to perform spatial correlations on
%f1 = @(x) max(x')';
%f2 = @(x) min(x')';
%f3 = @(x) PEig(x);
%f4 = @(x) mean(abs(hilbert(x))')';
%funcs = {f1 f2 f3 f4};

funcs = {@(x) sum(PEig(x)')' @(x) mean(abs(hilbert(x))')'};
cfg.funcs  = funcs;

cfg.dfuncs = {@(x) real(x)  @(x) abs(hilbert(x))};

% Set up plots
%--------------------------------------------------------
plots.review_topo = 0; 
plots.plot_change = 1;
plots.review_spatial = 0;


figure('position',[1000          72         788         906])
New_FICA_CTF_ST(cfg,plots);