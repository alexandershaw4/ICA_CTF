function FastICA_CTF_SpatioTemp(Dname,NC,UL,fname,bonf,writelog)
% FastICA for CTF MEG dataset
%
% Thresholds each trial of epoched data, checking that no value exceeds
% more than 3 * standard deviation for:
%
% - sample standard deviation of over trials
% - channel standard deviation over samples
% - standard deviation per channel over trials
%
% otherwise trial is marked as bad (written to ClassFile.cls of output .ds)
%
% For 'good' trials it implements:
%
% - ICA of bandpass filtered data, per trial.
% - Find correlations between component time-series and principal component 
%   of peripheral (EOG) channels. 
% - Project spatial components of these correlated components and robust
%   (weighted) average, forming a subject specific topography. 
% - Find correlations between this topo and spatial ICA components
% - Remove ICA components that are both temporally & spatially correlated.
% - Also check the same thresholds per-component as done for real trial
%   data. If surpases 3SDs of mean, switch off component. 
% - After all trials, write new dataset (.ds) and write BAD to ClassFile.
%
% Example usage:
%               FastICA_CTF_SpatioTemp(MEG_Cut.ds,[],[],'_ica',0)
% Inputs:
%         Dname = (epoched) CTF .ds dataset
%         NC    = number of components in data (optional, def = num chans)
%         UL    = upper limit on number of comps to reject (optional)
%         fname = output dataset appendix, e.g. 'ICA' for MEG_Cut_ICA.ds
%         bonf  = flag for bonferonni correction (def 0)
%         
%         writelog = flag (0/1) to write to console (default) or a logfile
%
% Dependencies: (matlab) 
%   - fieldtrip (for ft_topoplotter)
%   - spm (for spm_vec & spm_robust_average) [included]
%   - fastica algorithm
%   - the other tools in this github package
%
% AS2017


review_topo = 0; 
plot_change = 0;

% Paths
addpath('/home/sapas10/spm12/'); 
addpath('/home/sapas10/FastICA_25/');
%addpath('/home/sapas10/fieldtrip-20170509/');
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

% Max number rejected components
if nargin < 3
    UL = [];
end

% New file name [clone]
if nargin < 4 || isempty(fname)
    fname = 'ica';
end


if nargin == 5 && isempty(NC)
    NC = size(Data,2);
end

% Bonferroni or not
if nargin < 5; bonf = 0; end
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
time = t(2:end);

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
thr_chan = 3*squeeze(std(Data,[],3)); % std over trials (chans * samps)
thr_samp = 3*squeeze(std(Data,[],2)); % std over samples (chans * trials)

samp_thr = 1/6; % proportion of noisy* samples for trial rejection
chan_thr = 1/7; % proportion of channels to simultaneously have noise*
                % *Noise as defined by thr_chan & thr_samp


% Start loop over trials
%-------------------------
for t = 1:nt
    clear p_temp p Q
    
    fprintf(loc,'\n\nAnalysing trial %d of %d\n',t,nt);
    
    cD  = squeeze(Data(:,:,t)); 
    e   = E(:,:,t);
    
    %thresh check trial first
    %--------------------------
    BAD(t) = threshold(cD,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,t,1);
    
    
%     % threshold check first, just reject trial
%     %-----------------------------------------
%     noise  = cD > thr_chan;
%     for ch = 1:size(cD,1) 
%         if sum(noise(ch,:)) > round(NS*samp_thr); 
%             fprintf('Rejecting trial %d on channel noise (samp std over trials)\n',t);
%             BAD(t) = 1; 
%             break;continue
%         end
%     end
%     for ch = 1:size(cD,2)
%         if sum(noise(:,ch)) > round(size(MEGid,1)*chan_thr);
%             fprintf('Rejecting trial %d on channel noise (chan std over samples)\n',t);
%             BAD(t) = 1;
%             break;continue
%         end
%     end
%     for ch = 1:NS
%         if noise(:,ch) > (mean(thr_samp,2)/NS);
%             fprintf('Rejecting trial %d on sample noise (std per chan over trials)\n',t);
%             BAD(t) = 1;
%             break;continue;
%         end
%     end
        
    
    % do the ica if trial not bad
    %-----------------------------
    if ~BAD(t)
    
        % bandpass filter
        fprintf(loc,'Bandpass filtering\n');
        cD = bandpassfilter(cD,SR,[1 100]);
        e  = bandpassfilter(e ,SR,[1 20]);

        % Sort regressors
        for ej = 1:size(e,1)
            e(ej,:) = HighResMeanFilt(e(ej,:),1,16);
        end

        % reduce dims of peripheral channels
        fprintf(loc,'Taking prinipcal component from peripheral channel matrix\n');
        e = PEig(e');

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
        for ch = 1:size(C,1)
            iw = pinv(W);
            xx = C(ch,:);
            yy = iw(:,ch);
            thiscomp = ( xx'* yy' )';            
            bad_comp = threshold(thiscomp,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,ch,0);
            
            if bad_comp
                fprintf(loc,'Rejecting component %d based on noise\n',ch);
                C(ch,:) = 0;
            end
        end
        

        % temporal correlates
        for i = 1:size(C,1)        
            c = C  (i,:);
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
        % project spatial portion of these comps
        TOPO         = (C'*TOPOi')'; 
        topo         = spm_robust_average(TOPO,2);

        % correlate this trial-specific topo with mixing matrix
        [Q,ps] = corr(iW,topo);
        p_spat = find(ps<thrp);

        % temporally and sptially bad
        common = tmp(ismember(tmp,p_spat));

        % log num comps rejected per trial
        if   isempty(common); continue; 
        else fprintf(loc,'A total of %d components correlate both temporally & spatially\n',length(common));
        end
        try  AllGone(t,:) = length(common); end


        if review_topo
            cfg.layout = 'CTF275.lay';
            thelay     = ft_prepare_layout(cfg);

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
            C(common,:)  = 0;
            iW(:,common) = 0;
        end

        % store
        if ~isempty(common)
            this = ( C'*iW')';
            Data(:,:,t) = this;
        end

        if t > 50
            % debugging breakpoint
            %disp('test');
        end


        if plot_change
            subplot(221),imagesc(cD);
            subplot(222),imagesc(this);
            subplot(2,2,[3 4]), plot(time,this,'color',[.4 .4 .4]);hold on;
            eplot = TSNorm(e,5)*max(this(:));
            plot(time,eplot,'r','LineWidth',3);hold off
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

%newname = [Dname(1:end-3) fname];
cd(fp);
%fname = [newname];
fname = [fn fname];
writeCTFds([fp '/' fname fe],D,Orig);

% close log
if writelog;
    fclose(loc);
end

% write bad trials
ctf_write_BadTrials((find(BAD)),[fp '/' fname fe],t)

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


function [BAD] = threshold(cD,thr_chan,thr_samp,samp_thr,chan_thr,NS,MEGid,t,v)
persistent loc;

    % threshold check first, just reject trial
    %-----------------------------------------
    BAD = 0;
    noise  = cD > thr_chan;
    for ch = 1:size(cD,1) 
        if sum(noise(ch,:)) > round(NS*samp_thr); 
            if v; fprintf(loc,'Rejecting trial %d on channel noise (samp std over trials)\n',t);end
            BAD = 1; 
            return;%break;continue
        end
    end
    for ch = 1:size(cD,2)
        if sum(noise(:,ch)) > round(size(MEGid,1)*chan_thr);
            if v; fprintf(loc,'Rejecting trial %d on channel noise (chan std over samples)\n',t);end
            BAD = 1;
            return;%break;continue
        end
    end
    for ch = 1:NS
        if noise(:,ch) > (mean(thr_samp,2)/NS);
            if v; fprintf(loc,'Rejecting trial %d on sample noise (std per chan over trials)\n',t);end
            BAD = 1;
            return;%break;continue;
        end
    end



end