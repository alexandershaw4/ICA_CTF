
cd('/cubric/scratch/sapas10/WAND/Visual/314_09540/');

% Set up config
%--------------------------------------------------------------------------
cfg.Dname = '314_09540_Visual.ds'; % name of the raw (unepoiched) CTF ds
cfg.NC    = 274;             % number of components to look for
cfg.UL    = [];             % upper limit on num to reject     
cfg.fname = '_icanew';      % output filename appendix
cfg.bonf  = 1;              % Bonferonni correct correlations
cfg.writelog = 0;           % flag: write a logfile or to command window
cfg.window   = 4;          % length of windows to compute ICA
cfg.overlap  = cfg.window;  % window overlap (set to len window == no overlap)
cfg.bpf_data = [1 300];     % bandpass data
cfg.bpf_eogs = [1 20];     % bandpass EOG / peripheral

% feature selection functions to perform spatial correlations on
%--------------------------------------------------------------------------
%f1 = @(x) max(x')';
%f2 = @(x) min(x')';
%f3 = @(x) PEig(x);
%f4 = @(x) mean(abs(hilbert(x))')';
%funcs = {f1 f2 f3 f4};

cfg.funcs = {@(x) sum(PEig(x)')' @(x) mean(abs(hilbert(x))')'};

% feature selection functions to perform temporal correlations on
%--------------------------------------------------------------------------
%cfg.dfuncs = {@(x) real(x)  @(x) abs(hilbert(x))};
cfg.dfuncs = {@(x) real(x)  @(x) abs(hilbert(x)) @(x) unwrap(angle(x)) };

% Set up plots
%--------------------------------------------------------------------------
plots.review_topo = 0; 
plots.plot_change = 0;
plots.review_spatial = 0;


%figure('position',[1000          72         788         906])
New_FICA_CTF_ST(cfg,plots);