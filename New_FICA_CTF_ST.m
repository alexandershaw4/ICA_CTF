function New_FICA_CTF_ST(cfg,plots)
% FastICA for *RAW* CTF MEG dataset - removal verison.
%
% Generates 'windows' of n-seconds
%
% Inputs
%
% cfg.Dname    = 'MEG.ds'; % the CTF MEG dataset
% cfg.NC       = [];       % number of components to find [def num chans]
% cfg.UL       = [];       % max number of components to reject [def none]
% cfg.fname    = 'ica';    % appendix for output file [def 'ica']
% cfg.bonf     = 0;        % whether to Bonferonni correct correls [def 0]
% cfg.writelog = 0;        % whether to write logfile (1) or to console (0)
% cfg.window   = 20;       % window length for ICA computation (def 20 seconds)
% cfg.overlap  = cfg.window; % set to same to use non-overlapping windows
% cfg.funcs    = [];       % *functions to apply to spatial data for correlation
% cfg.dfuncs   = [];       % *functions to apply to component (time) data (e.g. hilbert)
% cfg.bpf_data = [1 100];  % bandpass edges for the MEG channel data
% cfg.bpf_eogs = [1 100];  % bandpass edges for EOG + peripheral channels
% 
% plots.review_topo    = 0; % flag to plot topography for each window
% plots.plot_change    = 1;
% plots.review_spatial = 1;

% Parse input structures: config
if ~isfield(cfg,'Dname');   return;           else; Dname    = cfg.Dname;end
if ~isfield(cfg,'NC');      NC       = [];    else; NC       = cfg.NC;   end
if ~isfield(cfg,'UL');      UL       = [];    else; UL       = cfg.UL;   end
if ~isfield(cfg,'fname');   fname    = 'ica'; else; fname    = cfg.fname;end
if ~isfield(cfg,'bonf');    bonf     = 0;     else; bonf     = cfg.bonf; end
if ~isfield(cfg,'writelog');writelog = 0;     else; writelog = 1;        end
if ~isfield(cfg,'window');  window   = 20;    else; window   = cfg.window; end
if ~isfield(cfg,'overlap'); overlap  = window;else; overlap = cfg.overlap; end
if ~isfield(cfg,'funcs');   funcs    = defs;  else; funcs    = cfg.funcs;  end
if ~isfield(cfg,'dfuncs');  dfuncs   = defs_d;  else; dfuncs = cfg.dfuncs; end
if ~isfield(cfg,'bpf_data');data_bpf = [1 100]; else; data_bpf = cfg.bpf_data;end
if ~isfield(cfg,'bpf_eogs');peri_bpf = [1 100]; else; peri_bpf = cfg.bpf_eogs;end

% Plots
if ~isfield(plots,'review_topo');   review_topo   = 0; else; review_topo    = plots.review_topo;   end
if ~isfield(plots,'plot_change');   plot_change   = 1; else; plot_change    = plots.plot_change;   end
if ~isfield(plots,'review_spatial');review_spatial= 1; else; review_spatial = plots.review_spatial;end

% Paths
addpath('/home/sapas10/spm12/'); 
addpath('/home/sapas10/FastICA_25/');
addpath(genpath('/home/sapas10/code/ICA_CTF'));
if ~exist('fasticag'); getfastica;    end

% write log, or to screen
persistent loc;
if nargin < 6; writelog = 0; end
if writelog;
    [fp,fn,ne] = fileparts(Dname);
    loc        = fopen([fn '.txt'],'w');
else
    loc = 1;
end

% read ctf - using fieldtrip
fprintf(loc,'Processing dataset: %s\n',Dname);
ds = ft_read_data(Dname);
D  =  ft_read_header(Dname);

Periph = strmatch('eeg',D.chantype);
MEGid  = strmatch('meggrad',D.chantype); 

% Remove periph from data matrix
ii   = find(~ismember(1:size(ds,2),Periph));
Orig = ds;
Data = ds(MEGid,:);

% for newer cardiff datasets
if isempty(Periph)
    Periph = [Periph strmatch('hEOG', D.res4.chanNames)];
    Periph = [Periph strmatch('vEOG', D.res4.chanNames)];
    Periph = [Periph strmatch('ECG', D.res4.chanNames)];
    Periph = [Periph strmatch('EMG', D.res4.chanNames)];
end
E    = Orig(Periph,:);

% How many components
if isempty(NC); NC = size(Data,1); end

% Bonferroni or not
if ~bonf
     thrp = .05;
     fprintf(loc,'Not Bonferroni correcting...\n');
else thrp = .05 / NC;
     fprintf(loc,'Bonferroni correcting...\n');
end

% compute (1) channel covariance and (2) distances between grads
grad_i    = strncmp('M',D.grad.label,1) ;
chan_pos  = D.grad.chanpos(grad_i,:);
chan_dist = cdist(chan_pos,chan_pos);
chan_cov  = cov(Data');


% time
t = ((0:D.nSamples)/D.Fs ) - (D.nSamplesPre/D.Fs);
time = t;
onset = t(1);
offset = t(end);

NS = D.nSamples;
SR = D.Fs;
NT = D.nTrials;

%fprintf(loc,'Using trial times %d to %d (@ sr %d Hz)\n',onset,offset,D.res4.sample_rate);
fprintf(loc,'Using trial times %d to %d (s) (@ sr %d Hz)\n',onset,offset,D.Fs);

% Generate a sliding window of length 'window' that will glide through the
prog   = overlap * SR;
wlen   = window  * SR;

if prog == wlen
    % if aiming for non-overlapping windows
    wlen = wlen - 1;
end

clear block
for i  = 1:inf
    i0 = i*prog;
    i1 = i0 + wlen;
    
    if i1 < (NS+wlen+1)
        block(i,:) = [i0 i1];
    else
        break
    end
end
block = block - wlen   

Nwind = size(block,1);


Data  = Data * 10^15;                     % convert T to fT
scale = norm((Data*Data')./size(Data,2)); % scale the data
scale = sqrt(scale);

% Start loop over trials
%-------------------------
for t = 1:Nwind
    clear p_temp p Q
    
    fprintf(loc,'\n\n+Analysing window %d of %d\n',t,Nwind);
    Win  = block(t,1):block(t,2);
    
    cD  = squeeze(Data(:,Win)) / scale;  
    e   = E(:,Win);
        
    % some descriptives to print for comparison after ica
    vcD  = var(cD');
    ivcD = find(vcD~=0);
       
    % do the ica 
    %---------------------------------------------------------------
    
    % bandpass filter
    fprintf(loc,' -Bandpass filtering\n');
    cD = bandpassfilter(cD,SR,data_bpf);
    e  = bandpassfilter(e ,SR,peri_bpf);
    
    % reduce dims of peripheral channels
    fprintf(loc,' -Taking prinipcal component from peripheral (EOG etc) channel matrix\n');
    try
        e = sum( PEig(e') ,2);
    catch
        if t == 1
            fprintf(loc,[' -EOG / peripheral channel is too big to do an ' ...
                'eigen decomposition\nso using mean channel data instead\n']);
        end
        e = mean(e,1);
        opt.eig = 0;
    end
        
    try
        fprintf(loc,' -Trying ICA\n');
        if nargin < 2 || isempty(NC)
            [C , A , W]   = fastica(cD,'verbose','off');
        else [C , A , W]  = fastica(cD,'numOfIC',NC,'verbose','off');
        % 'whiteSig',ws,'whiteMat',wm,...
           %     'dewhiteMat',dwm,
        end
        fprintf(' -Components in ica: %d\n',size(C,1));
    catch
        fprintf(loc,' -Could not compute ICA: skipping this trial...\n\n');
        continue;
    end
    
    % temporal correlates of the coomponents
    p_temp = zeros(size(C,1),1);
    for i  = 1:size(C,1)
        %c = C  (i,:);
        
        % functions of C
        for i0 = 1:length(dfuncs)
            f0 = dfuncs{i0};
            c  = C (i,:);
            c  = f0(c);
            
            try
                [Q,p] = corr(c(:), e(:));
                if any(p < thrp)
                    p_temp(i) = (p);
                end
            catch
                fprintf(loc,' -Could not run correlations for component %d\n',i);
                p = [];
            end
        end
    end
    
    try  p_temp; catch p_temp = []; end
    
    p_temp(p_temp==0)=1;
    p_temp = (p_temp < .05);
    
    if   isempty(p_temp); continue;
    else fprintf(loc,' -Found %d temporally correlated components...\n',length(find(p_temp)));
         fprintf(loc,' -Constructing topographies from these components...\n');
    end
    
    % only include as bad if all periph ncomp correlate
    tmp = find(p_temp);
    iW  = pinv(W);
    
    % make topographies for these
    TOPOi        = iW*0;
    TOPOi(:,tmp) = iW(:,tmp);
    tC           = C*0;
    tC(tmp,:)    = C(tmp,:);
    TOPO         = (tC'*TOPOi')';
    
    % use max, min and peig
    for i0 = 1:length(funcs)
        f  = funcs{i0};
        
        Y{i0}    = f(TOPO);
        [Q,ps]   = corr(iW,Y{i0});
        p_sp{i0} = find(ps<thrp);
        %p_all(i0,:) = ps;
    end
    
    %p_all(p_all==0)=1;
    %test=p_all*inv(chan_dist.*chan_cov);
    %p_spat = unique([ find(test(1,:)<.05) find(test(2,:)<.05) ])
    
    p_spat  = unique(cat(1,p_sp{:}));
    
    
    if review_spatial
        try thelay; catch; load('ctflay'); end
        x   = thelay.pos(1:length(MEGid),1);
        y   = thelay.pos(1:length(MEGid),2);
        tri = delaunay(x,y);
        npl = length(funcs);
        
        for i0 = 1:npl
            subplot(1,npl,i0);
            trisurf(tri,x,y,Y{i0}); title(char(funcs{i0})); view([0 90]);
            hold on;shading interp;set(gca,'visible','off');
        end
        drawnow;
    end
    
    
    % temporally and sptially bad
    common = tmp(ismember(tmp,p_spat));
    
    % log num comps rejected per trial
    if   isempty(common); continue;
    else fprintf(loc,' -A total of %d components correlate both temporally & spatially\n',length(common));
    end
    try  AllGone(t,:) = length(common); end
    
    
    if review_topo
        try thelay; catch; load('ctflay'); end
        
        x   = thelay.pos(1:length(MEGid),1);
        y   = thelay.pos(1:length(MEGid),2);
        tri = delaunay(x,y);
        
        bC   = C(common,:);
        biW  = iW(:,common);
        BADS = (bC(:,:)'*biW(:,:)')' ;
        
        SCL =       spm_vec( mean( (C'*iW')' ,2));
        SCL = [SCL; spm_vec( mean(BADS,2))];
        ylm(1) = min(SCL(:))*1.1;
        ylm(2) = max(SCL(:))*1.1;
        
        % MEG Orig
        ORIG = ((C'*iW')');
        subplot(131),...
            trisurf(tri,x,y,mean(ORIG,2));
        shading interp; view(0,90);caxis([ylm(1) ylm(2)]);
        set(gca,'visible','off');
        title('ORIG','fontsize',20);
        set(findall(gca, 'type', 'text'), 'visible', 'on');
        
        % MEG BAD Comps
        subplot(132),...
            trisurf(tri,x,y,mean(BADS,2));
        shading interp; view(0,90);caxis([ylm(1) ylm(2)]);
        set(gca,'visible','off');
        title('BAD','fontsize',20);
        set(findall(gca, 'type', 'text'), 'visible', 'on');
        
        % remove bad comps from MEG
        C(common,:)   = 0;
        iW(:,common)  = 0;
        
        % MEG Cleaned Comps
        subplot(133),...
            trisurf(tri,x,y,mean( (C'*iW')' ,2));
        shading interp; view(0,90);caxis([ylm(1) ylm(2)]);
        set(gca,'visible','off');
        title('CLEANED','fontsize',20);
        set(findall(gca, 'type', 'text'), 'visible', 'on');
        drawnow;
    else
        % grab a copy of the BAD components, too
        nc       = size(C,1);
        goodi    = find(~ismember(1:nc,common));
        
        badcomp          = C;
        badcomp(goodi,:) = 0;
        badW             = iW;
        badW(:,goodi)    = 0;
        
        C(common,:)  = 0;
        iW(:,common) = 0;
    end
    
    % store
    if ~isempty(common)
        this = ( C'*iW')';
       % this = reshape(this,size(cD));
        Data(:,Win) = this * scale;
        
        % compute post V and nV==0
        vThis  = var(this');
        ivThis = find(vThis~=0);
        
        fprintf(' +Mean prior variance: %d (%d non-zero)\n',mean(vcD),sum(ivcD));
        fprintf(' +Mean post variance: %d (%d non-zero)\n',mean(vThis),sum(ivThis));
        
    end
    
    if plot_change && ~review_topo
        thiscat = (this);
        %timecat = linspace(time(1),time(end)*LWind,size(thiscat,2));
        subplot(421),imagesc(cD); title('Orig Chan x Time Data');
        subplot(422),imagesc(thiscat);title('Cleaned Chan x Time Data');
        
        subplot(4,2,[3 4]), plot(cD'     ,'color',[.4 .4 .4]);hold on;
        eplot = TSNorm(e,5)*max(this(:));
        plot(eplot,'r','LineWidth',3);hold off
        title('Original Time x Channel Data');
        
        subplot(4,2,[5 6]), plot(thiscat','color',[.4 .4 .4]);hold on;
        plot(eplot,'r','LineWidth',3);hold off
        title('Cleaned time x Channel Data');
        
        thebad = (badcomp'*badW')';
        subplot(4,2,[7 8]), plot(thebad','color',[.4 .4 .4]);hold on;
        plot(thebad','b','LineWidth',2);
        plot(eplot','r','LineWidth',0.5);alpha .3; hold off;
        title('Bad components');
        drawnow;
    end
    
end

fprintf(loc,'Removed an average of %i components per trial\n',round(mean(AllGone)));

% return
Orig(MEGid,:)   = Data / 10^15; % fT -> T

% writeCTFds
[fp,fn,fe] = fileparts(Dname);
if isempty(fp); fp = evalinContext('pwd'); end

% write the dataset including bad markers
cd(fp);
fname = [fn fname];
%D.TrialClass.trial = allbad;

DatHead = readCTFds(Dname);
%DataTest   = getCTFdata(DatHead);
%Orig = Orig';
writeCTFds([fp '/' fname fe],DatHead,Orig');

% close log
if writelog;
    fclose(loc);
end

end

function getfastica()

url= 'https://github.com/davidkun/FastICA/archive/master.zip';

fprintf('Installing fastica algorithm from github\n');
str = (url);
urlwrite(str,'fastica.zip');
unzip('fastica.zip');

addpath('FastICA-master');
!rm fastica.zip

end


function funcs = defs
fprintf('Using default feature functions for spatial correlations\n');

% functions to perform spatial correlations on
f1 = @(x) max(x')';
f2 = @(x) min(x')';
f3 = @(x) PEig(x);
f4 = @(x) mean(abs(hilbert(x))')';

funcs = {f1 f2 f3 f4};
end

function dfuncs = defs_d
fprintf('Using default feature functions for temporal correlations\n');

dfuncs = {@(x) real(x)};

end