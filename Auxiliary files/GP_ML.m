function NLL = GP_ML(x,X,y,cov)
h1 = size(X,2)+1 ; % number of hyperparameters from covariance 
h2 = 1;            % number of hyperparameters from likelihood

hyp.cov = x(1:h1);
hyp.lik = x(end);
hyp.mean = [];

% calc K
n = size(y,1);
K       = cov(hyp.cov,X)+exp(hyp.lik*2)*eye(n) + 1e-4*eye(n);
K       = (K+transpose(K))/2;
CH      = chol(K);
invK    = solve_chol(CH,eye(n));

logDetK = 2*sum(log(abs(diag(CH))));
NLL = logDetK + y'*invK*y ;
end