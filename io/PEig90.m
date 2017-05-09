function [y,n] = PEig90(x,varargin)
% SVD based PCA, retaining n-modes needed to explain 90% variance
%
% AS2016 [util]


[u,s,v] = svd(x);
eigVals = diag(s);

for i = 1:length(eigVals)
    energy(i) = sum(eigVals(1:i));
end

propEnergy  = energy./energy(end);
n           = min(find(propEnergy > 0.9));

fprintf('%d components explain 90 percnt variance\n',n);

%for i = 1:n
%    y(i,:)       = u(:,(i))*s((i),(i))*mean(v(:,(i)));
%end

y       = u(:,1:n)*s(1:n,1:n)*mean(v(:,1:n))';