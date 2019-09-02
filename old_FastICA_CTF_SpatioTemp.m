function FastICA_CTF_SpatioTemp(cfg,thr,plots)
% FastICA for CTF MEG dataset - removal verison.
%
% Generates 'windows' of n-seconds by concatenating m-trials.
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
% cfg.funcs    = [];       % *functions to apply to spatial data for correlation
% cfg.dfuncs   = [];       % *functions to apply to component (time) data (e.g. hilbert)
% cfg.bpf_data = [1 100];  % bandpass edges for the MEG channel data
% cfg.bpf_eogs = [1 100];  % bandpass edges for EOG + peripheral channels
% 
% 
% thr.threshold_raw            = 0; % flag to threshold the raw data [def 0]
% thr.threshold_comps          = 1; % flag to threshold components [def 1]
% thr.samp_std_over_trials     = 0; % what to threshold: samp_std_over_trials
% thr.chan_std_over_samples    = 1; %                    chan_std_over_samples
% thr.std_per_chan_over_trials = 1; %                    std_per_chan_over_trials
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
if ~isfield(cfg,'funcs');   funcs    = defs;  else; funcs    = cfg.funcs;  end
if ~isfield(cfg,'dfuncs');  dfuncs   = defs_d;  else; dfuncs = cfg.dfuncs; end
if ~isfield(cfg,'bpf_data');data_bpf = [1 100]; else; data_bpf = cfg.bpf_data;end
if ~isfield(cfg,'bpf_eogs');peri_bpf = [1 100]; else; peri_bpf = cfg.bpf_eogs;end

% Thresholds
if ~isfield(thr,'threshold_raw');           thr.threshold_raw            = 0;end
if ~isfield(thr,'threshold_comps');         thr.threshold_comps          = 1;end
if ~isfield(thr,'samp_std_over_trials');    thr.samp_std_over_trials     = 0;end
if ~isfield(thr,'chan_std_over_samples');   thr.chan_std_over_samples    = 1;end
if ~isfield(thr,'std_per_chan_over_trials');thr.std_per_chan_over_trials = 1;end

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

% read ctf
fprintf(loc,'Processing dataset: %s\n',Dname);
D      = readCTFds(Dname);
Data   = getCTFdata(D);
Periph = strmatch('EEG', D.res4.chanNames);
chans  = cellstr(D.res4.chanNames);
MEGid  = cat(1,D.res4.senres.sensorTypeIndex);
MEGid  = find(MEGid==5);

% Remove periph from data matrix
ii   = find(~ismember(1:size(Data,2),Periph));
Orig = Data;
Data = Data(:,MEGid,:);
E    = Orig(:,Periph,:);


% How many components
if isempty(NC); NC = size(Data,2); end

% Bonferroni or not
if ~bonf
     thrp = .05;
     fprintf(loc,'Not Bonferroni correcting...\n');
else thrp = .05 / NC;
     fprintf(loc,'Bonferroni correcting...\n');
end


% Read epoch info from hist
T{1} = [' regexp(D.hist, ''StartTime: *'', ''match'') '];
T{2} = [' regexp(D.hist, ''EndTime: *'', ''match'') '];

[o1,o2] = eval(T{1});
onset   = str2num(D.hist(o2+11:o2+16));
[o1,o2] = eval(T{2});
offset  = str2num(D.hist(o2+9:o2+13));
t       = onset:(1/D.res4.sample_rate):offset;
time    = t(2:end);

if isempty(t) % resting
    NS = D.res4.no_samples;
    SR = D.res4.sample_rate;
    NT = D.res4.no_trials;
    onset  = 0;
    offset = NS/SR;
    t = linspace(onset,offset,NS);
    time = t;
end
   
fprintf(loc,'Using trial times %d to %d (@ sr %d Hz)\n',onset,offset,D.res4.sample_rate);


% Permute so dim1 = channels
Data = permute(Data,[2 1 3]);
E    = permute(E,[2 1 3]);
nt   = size(Data,3);
BAD  = zeros(1,nt);

% make some thresholds
thr_chan = 3*squeeze(std(Data,[],3)); % std over trials (chans * samps) [3]
thr_samp = 3*squeeze(std(Data,[],2)); % std over samples (chans * trials)[3]

samp_thr = 1/6; % proportion of noisy* samples for trial rejection (1/6)
chan_thr = 1/7; % proportion of channels to simultaneously have noise* (1/7)
                % *Noise as defined by thr_chan & thr_samp

% Generate 'block' window of n-trials
window;
LWind = window / (NS/SR);
LWind = round(LWind);
W1    = 1:LWind:(NT-LWind);
W2    = LWind:LWind:(NT); 
Nwind = length(W1);
fprintf('Reconstructing %d second windows, concatenating %d trials\n',window,LWind);

% Reshape thresholds for block, not trials
thr_chan = repmat(thr_chan,[1,LWind]);


% Start loop over trials
%-------------------------
for t = 1:Nwind
    clear p_temp p Q
    
    fprintf(loc,'\n\nAnalysing window %d of %d\n',t,Nwind);
    Win  = W1(t):W2(t);
    
    cD1  = squeeze(Data(:,:,Win)); 
    cD   = reshape(cD1,[size(cD1,1),size(cD1,2)*size(cD1,3)]);
    e1   = E(:,:,Win);
    e    = reshape(e1,[size(e1,1),size(e1,2)*size(e1,3)]);
    
    %thresh check trial first
    %--------------------------
    if thr.threshold_raw
        [BAD(Win),STR] = threshold(cD,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,Win,1,thr);
    
        if any(BAD(Win));
            fprintf(loc,STR);
        end
    end
    
   
    % do the ica if trial not bad
    %-----------------------------
    if any(~BAD(Win))
    
        % bandpass filter
        fprintf(loc,'Bandpass filtering\n');
        cD = bandpassfilter(cD,SR,data_bpf);
        e  = bandpassfilter(e ,SR,peri_bpf);

        % reduce dims of peripheral channels
        fprintf(loc,'Taking prinipcal component from peripheral (EOG etc) channel matrix\n');
        try
            e = PEig(e');
        catch
            fprintf(loc,['EOG / peripheral channel is too big to do an ' ...
                         'eigen decomposition\nso using mean channel data instead\n']);
            e = mean(e,1);
            opt.eig = 0;
        end

        % the ica
        try
            fprintf(loc,'Trying ICA..........\n');
            if nargin < 2 || isempty(NC)
                 [C , A , W]  = fastica(cD,'verbose','off');
            else [C , A , W]  = fastica(cD,'numOfIC',NC,'verbose','off');
            end  
        catch
            fprintf(loc,'Could not compute ICA: skipping this trial...\n\n');
            BAD(t) = 1;
            continue;
        end
        
        % do same thresholding on each component
        if thr.threshold_comps
            for ch = 1:size(C,1)
                iw = pinv(W);
                xx = C(ch,:);
                yy = iw(:,ch);
                thiscomp = ( xx'* yy' )';            
                bad_comp = threshold(thiscomp,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,ch,0,thr);

                if bad_comp
                    fprintf(loc,'Rejecting component %d based on noise\n',ch);
                    C(ch,:) = 0;
                end
            end
        end
        

        % temporal correlates of the coomponents
        for i = 1:size(C,1)        
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
                    fprintf(loc,'Could not run correlations for component %d\n',i);
                    p = [];
                end
            end
        end  

        try  p_temp; catch p_temp = []; end
        if   isempty(p_temp); continue;
        else fprintf(loc,'Found %d temporally correlated components...\n',length(find(p_temp)));
             fprintf(loc,'Constructing topographies from these components...\n');  
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
        end
        
        p_spat  = unique(cat(1,p_sp{:}));
        
        if review_spatial
            try thelay; catch; load('ctflay'); end
            x   = thelay.pos(1:length(MEGid),1);
            y   = thelay.pos(1:length(MEGid),2);
            tri = delaunay(x,y);  
            npl = length(funcs);
            
            for i0 = 1:npl
                subplot(1,npl,i0);
                trisurf(tri,x,y,Y{i0}); title('max'); view([0 90]);
                hold on;shading interp;
            end
        end

        
        % temporally and sptially bad
        common = tmp(ismember(tmp,p_spat));

        % log num comps rejected per trial
        if   isempty(common); continue; 
        else fprintf(loc,'A total of %d components correlate both temporally & spatially\n',length(common));
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
            this = reshape(this,size(cD1));
            Data(:,:,Win) = this;
        end

        if plot_change && ~review_topo
            thiscat = reshape(this,[size(this,1),size(this,2)*size(this,3)]);
            timecat = linspace(time(1),time(end)*LWind,size(thiscat,2));
            subplot(421),imagesc(cD); title('Orig Chan x Time Data');
            subplot(422),imagesc(thiscat);title('Cleaned Chan x Time Data');
            
            subplot(4,2,[3 4]), plot(timecat,cD     ,'color',[.4 .4 .4]);hold on;
            eplot = TSNorm(e,5)*max(this(:));
            plot(timecat,eplot,'r','LineWidth',3);hold off
            title('Original Time x Channel Data');
            
            subplot(4,2,[5 6]), plot(timecat,thiscat,'color',[.4 .4 .4]);hold on;
            plot(timecat,eplot,'r','LineWidth',3);hold off
            title('Cleaned time x Channel Data');
            
            
            thebad = (badcomp'*badW')';
            subplot(4,2,[7 8]), plot(timecat,thebad,'color',[.4 .4 .4]);hold on;
            plot(timecat,thebad,'b','LineWidth',2);
            plot(timecat,eplot,'r','LineWidth',0.5);alpha .3; hold off;
            title('Bad components');
            drawnow;
        end
    
    end
    
end

% Put the dimensions back: ( samp * chan * trial)
Data = permute(Data,[2 1 3]);

fprintf(loc,'Removed an average of %i components per trial\n',round(mean(AllGone)));


% return
Orig(:,MEGid,:) = Data(:,:,1:size(Orig,3));

% writeCTFds
[fp,fn,fe] = fileparts(Dname);
if isempty(fp); fp = evalinContext('pwd'); end

% Bad trials
if isempty(UL)
    UL = NT*(8*(t/NT));
end
allbad = find(BAD);
try
    allbad = allbad(1:UL);
end

% write the dataset including bad markers
cd(fp);
fname = [fn fname];
D.TrialClass.trial = allbad;

writeCTFds([fp '/' fname fe],D,Orig);

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


function [BAD,STR] = threshold(cD,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,t,v,thr)

    % threshold check first, just reject trial
    %-----------------------------------------
    BAD = 0;
    STR = [];
    noise  = cD > thr_chan;
    if thr.samp_std_over_trials
        for ch = 1:size(cD,1) 
            if sum(noise(ch,:)) > round(NS*samp_thr); 
                if v; STR=sprintf('Rejecting trial(s) %d on channel noise (samp std over trials)\n',t);end
                BAD = 1; 
                return;%break;continue
            end
        end
    end
    
    if thr.chan_std_over_samples
        for ch = 1:size(cD,2)
            if sum(noise(:,ch)) > round(size(MEGid,1)*chan_thr);
                if v; STR=sprintf('Rejecting trial(s) %d on channel noise (chan std over samples)\n',t);end
                BAD = 1;
                return;%break;continue
            end
        end
    end
    
    if thr.std_per_chan_over_trials
        for ch = 1:NS
            if noise(:,ch) > (mean(thr_samp,2)/NS);
                if v; STR=sprintf('Rejecting trial(s) %d on sample noise (std per chan over trials)\n',t);end
                BAD = 1;
                return;%break;continue;
            end
        end
    end
end

function funcs = defs

% functions to perform spatial correlations on
f1 = @(x) max(x')';
f2 = @(x) min(x')';
f3 = @(x) PEig(x);
f4 = @(x) mean(abs(hilbert(x))')';

funcs = {f1 f2 f3 f4};
end

function dfuncs = defs_d

dfuncs = @(x) real(x);

end