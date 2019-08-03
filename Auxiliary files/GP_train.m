function [hyp,meany,invK] = GP_train(X,y)

%% Calculate hyperparameters
cov = @covSEard;
h1  = size(X,2)+1; % number of hyperparameters from covariance 
h2  = 1;           % number of hyperparameters from likelihood 

% Normalize y to zero-mean
meany = mean(y); y = y - meany;

% Definition of lower and upper bound of optimization algorithm
stdX = std(X); stdy = std(y);
lb = zeros(h1+h2,1);
ub = zeros(h1+h2,1);
lb(1:h1-1) =  log(stdX'/4);
ub(1:h1-1) =  log(stdX'*4);
lb(h1) = log(stdy'/4);
ub(h1) = log(stdy'*4);
lb(h1+h2) = log(10^(-6));
ub(h1+h2) = log(10^(2));

% x0 = (lb+ub)/2
x0 = [log(stdX');log(stdy');0];

% Defintion of options for fmincon solver
%options.Algorithm = 'sqp';
options.DerivativeCheck = 'off';
options.TolCon = 10^-12;
options.Display = 'off';
options.Hessian = 'lbfgs';
options.TolFun = 10^-12;
options.PlotFcns = [];
options.GradConstr = 'off';
options.GradObj = 'off';
options.TolX = 10^-10;
options.UseParallel = 0;

%keyboard
hyp_opt = fmincon(@(x) GP_ML(x,X,y,cov),...
    x0,[],[],[],[], lb, ub,[],options);
hyp.cov = hyp_opt(1:h1);
hyp.lik = hyp_opt(h1+1:h2+h1);

%% Calculate invK
n = length(y);
K       = cov(hyp.cov,X) + exp(hyp.lik*2)*eye(n) + 1e-6*eye(n)*(stdy^2);
K       = (K+transpose(K))/2;
try
    CH      = chol(K);
    invK    = solve_chol(CH,eye(n));
catch 
    warning(' used pseudo-inverse since Erics code failed ')
    invK = pinv(K,1e-12);
end