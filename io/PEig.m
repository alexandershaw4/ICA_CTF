function y = PEig(x,varargin)
% Retain principal eigenmode from a vector or matrix.
% [or specify which mode in second argument]
%
% AS2016 [util]


if nargin > 1; 
      n = varargin{1};
else  n = 1;
end

if ndims(x) > 2
    x = VecRetainDim(x);
end

if length(n) > 1
    for i = 1:length(n)
        y(i,:) = PEig(x,n(i));
    end
    return
end



[u,s,v] = svd(x);
y       = u(:,n)*s(n,n)*mean(v(:,n));