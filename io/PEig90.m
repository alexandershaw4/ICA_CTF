function [y,n] = PEig90(x,varargin)
% SVD based PCA, retaining n-modes needed to explain 90% variance
%
% AS2016 [util]


% if exist('pca')
%     fprintf('Using installed pca function\n');
%     
    % use the pca function from stats toolbox/built in
    [coeff, score, latent, tsquared, explained, mu] = pca(x);

    % find num to explain 90% variance
    tot_exp = 0; 
    n       = 0;

    while tot_exp < 90;
          n       = n + 1;
          tot_exp = tot_exp + explained(n);
    end
    
    fprintf('%d components explain 90 percnt variance\n',n);

    y = score(:,1:n) * coeff(:,1:n)' + repmat(mu, size(x,1), 1);
%     
% 
% else
%     fprintf('PCA function not installed: trying SVD routine\n');
%     % otherwise try this crude svd routine
%     
%     % de mean
%     mu = mean(x);
%     x  = x - repmat(mu,[size(x,1),1]);
%     
%     % svd
%     [u,s,v] = svd(x);
%     eigVals = diag(s);
% 
%     % loop to get eigenval prop
%     for i = 1:length(eigVals)
%         energy(i) = sum(eigVals(1:i));
%     end
% 
%     propEnergy  = energy./energy(end);
%     n           = min(find(propEnergy > 0.9));
% 
%     fprintf('%d components explain 90 percnt variance\n',n);
% 
%     %for i = 1:n
%     %    y(i,:)       = u(:,(i))*s((i),(i))*mean(v(:,(i)));
%     %end
% 
%     y       = u(:,1:n)*s(1:n,1:n)*mean(v(:,1:n))';
%     y       = y + repmat(mu,[size(x,1),1]);
% end