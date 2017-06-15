function FastICA_CTF_SpatioTemp(Dname,NC,UL,fname,bonf)
% FastICA for CTF MEG dataset
%
%
%
% AS


review_topo = 0; 
plot_change = 0;

% Paths
addpath('/home/sapas10/spm12/'); 
addpath('/home/sapas10/FastICA_25/');
addpath('/home/sapas10/fieldtrip-20170509/');
if ~exist('fasticag'); getfastica;    end

% read ctf
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
if nargin < 5
     thrp = .05;
     fprintf('Not Bonferroni correcting...\n');
else thrp = .05 / NC;
     fprintf('Bonferroni correcting...\n');
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
   
fprintf('Using trial times %d to %d (@ sr %d Hz)\n',onset,offset,D.res4.sample_rate);


% Permute so dim1 = channels
Data = permute(Data,[2 1 3]);
E    = permute(E,[2 1 3]);
nt   = size(Data,3);
BAD  = zeros(1,nt);

% make some thresholds
thr_chan = 3*squeeze(std(Data,[],3)); % std over trials (chans * samps)
thr_samp = 3*squeeze(std(Data,[],2)); % std over samples (chans * trials)

samp_thr = 1/8; % proportion of noisy* samples for trial rejection
chan_thr = 1/8; % proportion of channels to simultaneously have noise*
                % *Noise as defined by thr_chan & thr_samp


% Start loop over trials
%-------------------------
for t = 1:nt
    clear p_temp p Q
    
    fprintf('Finding components in trial %d of %d\n',t,nt);
    
    cD  = squeeze(Data(:,:,t)); 
    e   = E(:,:,t);
    
    
    % threshold check first, just reject trial
    %-----------------------------------------
    noise  = cD > thr_chan;
    for ch = 1:size(cD,1) 
        if sum(noise(ch,:)) > round(NS*samp_thr); 
            fprintf('Rejecting trial %d\n',t);pause(1);
            BAD(t) = 1; 
            break;continue
        end
    end
    for ch = 1:size(cD,2)
        if sum(noise(:,ch)) > round(size(MEGid,1)*chan_thr);
            fprintf('Rejecting trial %d\n',t);pause(1);
            BAD(t) = 1;
            break;continue
        end
    end
    for ch = 1:NS
        if noise(:,ch) > (PEig(thr_samp)/NS);
            fprintf('Rejecting trial %d\n',t);pause(1);
            BAD(t) = 1;
            break;continue;
        end
    end
        
    
    % do the ica if trial not bad
    %-----------------------------
    if ~BAD(t)
    
        % bandpass filter
        fprintf('Bandpass filtering\');
        cD = bandpassfilter(cD,SR,[1 100]);

        % Sort regressors
        for ej = 1:size(e,1)
            e(ej,:) = HighResMeanFilt(e(ej,:),1,16);
        end

        % reduce dims of peripheral channels
        [junk,ncomp] = PEig90(e);
        e            = PEig(e',1:ncomp);
        fprintf('Reduced dimensionality of regressor matrix to %d\n',ncomp);

        % the ica
        if nargin < 2 || isempty(NC)
             [C , A , W]  = fastica(cD,'verbose','off');
        else [C , A , W]  = fastica(cD,'numOfIC',NC,'verbose','off');
        end    

        % temporal correlates
        for i = 1:size(C,1)        
            c = C  (i,:);

            try
                [Q,p] = corr(c(:), e(:));
                if any(p < thrp)
                    p_temp(i,:) = p;
                end
            catch
                fprintf('Could not run correlations for component %d\n',i);
                p = [];
            end
        end  

        try  p_temp; catch p_temp = []; end
        if   isempty(p_temp); continue;
        else fprintf('Found %d temporally correlated components...\n',length(find(p_temp)));
             fprintf('Constructing topographies from these components...\n');  
        end

        % only include as bad if all periph ncomp correlate 
        tmp = sum(logical(p_temp),2);
        tmp = find(tmp==ncomp);
        iW  = pinv(W);

        % make topographies for these
        TOPOi = iW*0;
        TOPOi(tmp,:) = iW(tmp,:); 
        TOPO = (C'*TOPOi')';      % project spatial portion of these comps
        topo = spm_robust_average(TOPO,2);

        % correlate this trial-specific topo with all mixing matrix
        [Q,ps] = corr(iW,topo);
        p_spat = find(ps<thrp);

        % temporally and sptially bad
        common = tmp(ismember(tmp,p_spat));

        % log num comps rejected per trial
        if   isempty(common); continue; 
        else fprintf('A total of %d components correlate both temporally & spatially\n',length(common));
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




        if plot_change
            subplot(221),imagesc(cD);
            subplot(222),imagesc(this);
            subplot(2,2,[3 4]), plot(time,this);hold on;
            plot(time,e*100,'r','LineWidth',3);hold off
            drawnow;
        end
        if t > 50;
            disp('rest');
        end
    
    end
    
end

% Put the dimensions back: ( samp * chan * trial)
Data = permute(Data,[2 1 3]);

fprintf('Removed an average of %i components per trial\n',round(mean(AllGone)));


% return
Orig(:,MEGid,:) = Data(:,:,1:size(Orig,3));

% writeCTFds
[fp,fn,fe] = fileparts(Dname);
if isempty(fp); fp = evalinContext('pwd'); end
newname = [Dname(1:end-3) fname];
cd(fp);
fname = [newname];
writeCTFds([fp '/' fname fe],D,Orig);

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
