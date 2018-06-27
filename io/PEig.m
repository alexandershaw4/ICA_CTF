function [y,coeff,score,latent,tsquared,explained,mu] = PEig(x,varargin);
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


% the pca function
[coeff, score, latent, tsquared, explained, mu] = pca(x);

% % find num to explain 90% variance
% tot_exp = 0;
% n       = 0;
% 
% while tot_exp < 90;
%     n       = n + 1;
%     tot_exp = tot_exp + explained(n);
% end
% 
% fprintf('%d components explain 90 percnt variance\n',n);

y = score(:,1:n) * coeff(:,1:n)' + repmat(mu, size(x,1), 1);



% [u,s,v] = svd(x);
% y       = u(:,n)*s(n,n)*mean(v(:,n));