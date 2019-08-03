function u_opt = ModifierAdapt_GP(lb,ub,GP)
options = optimset('disp','iter');
u_opt = fmincon(@(u)objective(u),[mean(GP.U_GP(:,:),2);0],...
[],[],[],[],[lb;-inf],[ub;inf],@(u)constraints(u,GP),options);
u_opt = u_opt(1:end-1);
end

function obj =  objective(u)
    obj = u(end);
end

function [c,ceq] =  constraints(u,GP)

for i = 1:length(GP.Y_GP_norm)
    [f_pred(i),f_var(i)] = GP_predict(GP.invK{i},GP.hyp{i},GP.U_GP,GP.Y_GP_norm{i},u(1:end-1),GP.meany{i});
end
F = Approx_data(u(1:end-1));

% objective
ceq = u(end) - (f_pred(1) + F(1));

% constraints g
for i = 2:length(GP.Y_GP_norm)
c(i) = f_pred(i) + F(i);
end
function [mean,var] = GP_predict(invK,hyp,X,y,x,meany)
    k    = covSEard(hyp.cov,X',x'); 
    mean = meany + k'*invK*y;
    var  = exp(2*hyp.cov(end)) - k'*invK*k;
end
end