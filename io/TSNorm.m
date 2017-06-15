function y = TSNorm(x,varargin)
% amplitude scaling function for timeseries, x.
% second argument is method [optional]: 
% 1 = scale between 0-1. [default]
% 2 = divide by maximum.
% 3 = mean correct and divide by std. 
% 4 = mean correct.
% 5 = scale approx [-1 1]
% 6 = scale between [-n and n]*
% 7 = scale by unit length (./norm(x))
%
% works on most datatypes (vect,matrix,cell etc) by calling SPMs
% [un]vectorising functions.
%
% if x is a matrix, function will either:
% - vectorise, normalise along whole vector and reshape [default], or
% - normalise each row separately: input forth argument = 1 for this. 
%
% *update: 3rd input n = scale between [-n n]
% *update: now appropriately handles sparse inputs 
%
% AS2016 [util]

% methods:
%--------------------------
T = {'Def','Max','Std','Mean','minusone','between','norm'};

if nargin > 1; t = T{varargin{1}};
else           t = T{1}; 
end

if nargin > 2
    n = varargin{2};
end

% matrix treatment:
%--------------------------
try varargin{3}; catch varargin{3} = 0; end

if  varargin{3} ~= 1   
    % vectorise and normalise against one constant
    s = x;
    x = spm_vec(x);
else
    % normalise each row of matrix 
    for i = 1:size(x,1)
        y(i,:) = TSNorm(squeeze(x(i,:)),varargin{1},varargin{2},0);
    end
    return
end


% handle sparse inputs
if issparse(x)
     dx = full(x(find(x)));
     xi = find(x);
     x  = dx;
     Sp = 1;
else Sp = 0;
end


switch t
    case 'Def';
        
    % default: 0-1 scaling norm
    y = (x - min(x)) / ( max(x) - min(x) );

    case 'Max';
        
    % divide by the maximum value
    y = x / max(x);
     
    case 'Std';
        
    % minus mean, divided by std
    y = (x - mean(x) ) / std(x);
    
    case 'Mean';
    
    % mean correct
    y = x - mean(x);
    
    case 'minusone';
        
    % scale [-1 1]
    y = (x-mean(x))/max(x);
    
    case 'between';
    
    % scale between -n and n
    y = -n + (abs(n)*2).*(x - min(x))./(max(x) - min(x));  
    
    case 'norm';
        
    % scale by unit length
    y = x./norm(x);
    
end


% re-instate sparse & index
if Sp;
   s(xi) = y;
   y     = s;
end
    

y = spm_unvec(y,s);