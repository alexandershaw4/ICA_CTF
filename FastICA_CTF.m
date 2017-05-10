function FastICA_CTF(Dname,NC,UL,fname,time,bonf)
% FastICA for CTF MEG dataset
% - Temporal only ica [at the moment]
% - Dname is a CTF MEG dataset [.ds]
% - NC is maximum number of components [optional]
% - UL is upper limit on number to remove [optional]
% - fname is the prefix to the filename, eg. 'ica_'
% - time is widow size in seconds [optional, def 10]
% - bonf is a flag for bonferroni correction [0/1] [optional]
%
% - This version works on a specified time window [time input, in s]
%   Concatenates n trials / epochs needed to make time window, runs ica & correls,
%   adjusts data, then places back into epoch / trial positions in data
%
% AS17 basod on my SPM {FastICA_MEEG_AS_5}


review_topo = 1; 
plotcor     = 0;

% Paths
addpath('/home/sapas10/spm12/'); 
addpath('/home/sapas10/FastICA_25/');
addpath('/home/sapas10/fieldtrip-20170509/');
if ~exist('fasticag'); getfastica;    end

D      = readCTFds(Dname);
Data   = getCTFdata(D);
Periph = strmatch('EEG', D.res4.chanNames);
chans  = cellstr(D.res4.chanNames);
MEGid  = cat(1,D.res4.senres.sensorTypeIndex);
MEGid  = find(MEGid==5);

for i = 1:length(Periph)
    if ~ any(any(Data(:,Periph(i),:))); Periph(i) = []; end
end

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
    fname = 'nica_';
end

% Time bin length
if nargin < 5 || isempty(time)
    time = 10 ;
end

if nargin == 6 && isempty(NC)
    NC = size(Data,2);
end

% Bonferroni or not
if nargin < 6
     thrp = .05;
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
fprintf('Using trial times %d to %d (@ sr %d Hz)\n',onset,offset,D.res4.sample_rate);

% sort num trials to concat
tocat = round(time/( t(end)-t(1) ));
nc    = ceil(size(Data,3)/tocat);
win   = 1;

% Permute so dim1 = channels
Data = permute(Data,[2 1 3]);
E    = permute(E,[2 1 3]);

for t = 1:nc
    thewin{t} = [win:win+tocat-1];
    win       = win + tocat;
end

% Start loop over windows
for t = 1:nc
    clear p_temp p Q
    
    fprintf('Finding components in time window %d\n',t);
    
    if t == length(thewin); 
        thewin{t} = thewin{t}(1):size(Data,3);
        tocat     = length(thewin{t});
    end
    
    cD  = squeeze(Data(:,:,thewin{t})); 
    cD  = reshape(cD,[size(cD,1) size(cD,2)*size(cD,3)]);
    e   = E(:,:,thewin{t});
    e   = reshape(e,[size(e,1) size(e,2)*size(e,3)]);
    
    % Sort regressors
    for ej = 1:size(e,1)
        e(ej,:) = HighResMeanFilt(e(ej,:),1,16);
    end
    
    % reduce dim
    [junk,ncomp] = PEig90(e);
    e            = PEig(e',1:ncomp);
    fprintf('Reduced dimensionality of regressor matrix to %d\n',ncomp);
    
    
    fprintf('Including trials %d to %d in this window\n',thewin{t}(1),thewin{t}(end));

    if nargin < 2 || isempty(NC)
         [C , A , W]  = fastica(cD);
    else [C , A , W]  = fastica(cD,'numOfIC',NC);
    end    
    
    for i = 1:size(C,1)        
        c = C  (i,:);
        
        try
            [Q,p] = corr(c', e');
            
            if any(p < thrp)
                fprintf('found peripherally correlated component: %d ... \n',i);
                p_temp(i,:) = p;
            end
        catch
            fprintf('Could not run correlations for component %d\n',i);
            p = [];
        end
    end  
    
    try p_temp; catch p_temp = []; end
    if isempty(p_temp); continue; end
    
    % only include as bad if all ncomp correlate 
    tmp = sum(logical(p_temp),2);
    tmp = find(tmp==ncomp);
    
    forkill = (tmp);
    iW      = pinv(W);
    
    if isempty(forkill); continue; end
    try AllGone(t,:) = length(forkill); end
    
    if plotcor
        
        bC  = C(forkill,:);
        biW = iW(:,forkill);
        BADS = (bC(:,:)'*biW(:,:)')' ;
        
        subplot(411), plot(e'); title('Peripherals [EOG, EMG] components');
        subplot(412), plot(mean(cD,1)),title('Original MEG [Mean]');
        subplot(413), plot(mean(BADS,1)); title('Correlated components [Mean]');
        
        C(forkill,:)   = 0;
        iW(:,forkill)  = 0;  
        
        subplot(414), plot(mean( (C'*iW')' ,1)); title('Corrected MEG');
        
        
        drawnow 
        
        
    elseif review_topo
        cfg.layout = 'CTF275.lay';
        thelay     = ft_prepare_layout(cfg);
        
        x   = thelay.pos(1:length(MEGid),1);
        y   = thelay.pos(1:length(MEGid),2);
        tri = delaunay(x,y);  
        
        bC   = C(forkill,:);
        biW  = iW(:,forkill);
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
        C(forkill,:)   = 0;
        iW(:,forkill)  = 0;  
        
        % MEG Cleaned Comps
        subplot(133),...
        trisurf(tri,x,y,mean( (C'*iW')' ,2));    
        shading interp; view(0,90);caxis([ylm(1) ylm(2)]);
        set(gca,'visible','off');        
        title('CLEANED','fontsize',20);
        set(findall(gca, 'type', 'text'), 'visible', 'on');
        drawnow;
    else
        C(forkill,:)  = 0;
        iW(:,forkill) = 0;
    end
    
    % store
    if ~isempty(forkill)
        mnew = reshape( ( C'*iW')'   ,[size(Data,1) size(Data,2) (tocat)]);
        Data(:,:,thewin{t}) = mnew;
    end
    
end

% Put the dimensions back: ( samp * chan * trial)
Data = permute(Data,[2 1 3]);

fprintf('Removed an average of %d components per trial\n',mean(AllGone));

% return
Orig(:,MEGid,:) = Data(:,:,1:size(Orig,3));

% writeCTFds
[fp,fn,fe] = fileparts(Dname);
if isempty(fp); fp = evalinContext('pwd'); end
fname = [fp '/' fname];
writeCTFds(fname,D,Orig);


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








