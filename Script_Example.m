
cd('/cubric/scratch/sapas10/WAND/Visual/314_09540/');

% Set up config
%--------------------------------------------------------
cfg.Dname = '314_09540_Visual.ds';
cfg.NC    = [];
cfg.UL    = [];
cfg.fname = 'ica';
cfg.bonf  = 0;
cfg.writelog = 0;
cfg.window   = 4;
cfg.overlap  = cfg.window;
cfg.bpf_data = [1 100];
cfg.bpf_eogs = [1 100];

% functions to perform spatial correlations on
%f1 = @(x) max(x')';
%f2 = @(x) min(x')';
%f3 = @(x) PEig(x);
%f4 = @(x) mean(abs(hilbert(x))')';
%funcs = {f1 f2 f3 f4};

funcs = {@(x) sum(PEig(x)')' @(x) mean(abs(hilbert(x))')'};

cfg.funcs  = funcs;
cfg.dfuncs = {@(x) abs(hilbert(x))};

% Set up plots
%--------------------------------------------------------
plots.review_topo = 0; 
plots.plot_change = 1;
plots.review_spatial = 0;


figure('position',[1000          72         788         906])
New_FICA_CTF_ST(cfg,plots);