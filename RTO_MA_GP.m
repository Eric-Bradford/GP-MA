% GP-MA

%% Specifications
Noise_mat     = diag([5e-1,5e-4,5e-4].^2); % Noise covariance
Tot_iter      = 15;                        % Number of iterations
lb            = [4,70]';                   % Lower bounds of inputs
ub            = [7,100]';                  % Upper bounds of inputs
u0            = [6.9,83]';                 % Initial inputs
filter_par    = 0.6;                       % Filter parameter
data_set_size = 3;                         % Initial data set size

%% Initial dataset  
U_GP0 = lhsdesign(data_set_size,size(u0,1));
U_GP0 = U_GP0'.*repmat(ub-lb,1,data_set_size) + repmat(lb,1,data_set_size); 
U_GP  = U_GP0;
GP.U_GP = U_GP;
for i = 1:size(U_GP,2)
    F_GP(i,:) = Plant_data(U_GP(:,i));
    Y_GP(i,:) = Plant_data(U_GP(:,i)) - Approx_data(U_GP(:,i))...
    + mvnrnd(zeros(1,length(Noise_mat)),Noise_mat);  
end

%% building GP models
for i = 1:size(Y_GP,2)
[GP.hyp{i},GP.meany{i},GP.invK{i}] = GP_train(U_GP',Y_GP(:,i));
GP.Y_GP_norm{i}                    = Y_GP(:,i) - GP.meany{i};  
end

%% GP real-time optimization
for rto_loop = 1:Tot_iter
        
        % Update input dataset
        U_GP    = [U_GP,u0]; GP.U_GP = U_GP;
        
        % Measurement
        f_new   = Plant_data(u0);
        F_GP    = [F_GP;f_new]; 
        y_new   = f_new - Approx_data(u0) + mvnrnd(zeros(1,length(Noise_mat)),Noise_mat);
        Y_GP    = [Y_GP;y_new];  
        for i = 1:size(Y_GP,2)
            [GP.hyp{i},GP.meany{i},GP.invK{i}] = GP_train(U_GP',Y_GP(:,i));
            GP.Y_GP_norm{i}                    = Y_GP(:,i) - GP.meany{i};  
        end
        
        % Determine optimal input
        u_opt = ModifierAdapt_GP(lb,ub,GP); 
        u0    = u0 + filter_par*(u_opt-u0);    
end

%% Plot graphs
u_opt_true      = fmincon(@(u) Plant_obj(u) ,[lb(1);lb(2)],[],[],[],[],lb,ub,[]);
u_opt_approx    = fmincon(@(u) Approx_obj(u),[lb(1);lb(2)],[],[],[],[],lb,ub,[]);
u1cont          = linspace(lb(1),ub(1),50);
u2cont          = linspace(lb(2),ub(2),50);
[U1cont,U2cont] = meshgrid(u1cont,u2cont); PROFIT = zeros(50,50);
i               = 0;
for i = 1:50
    for j = 1:50
        F           = Plant_data([u1cont(i),u2cont(j)]);
        PROFIT(i,j) = F(1);
    end
end

u1cont = linspace(lb(1),ub(1),100);
u2cont = linspace(lb(2),ub(2),100);
G      = [u1cont;u2cont]';
GC     = []; GU = [];
for i = 1:100
    for j = 1:100
        F  = Plant_data([G(i,1),G(j,2)]);     
        GC = [GC,F(2:end)'];
        GU = [GU,[G(i,1),G(j,2)]'];
    end
end
GC = GC';

Gplot = {};
for n_c = 1:(length(F)-1)
    [B,index]  = sortrows(abs(GC),n_c);
    GUsorted   = GU(:,index);        
    Gplot{n_c} = GUsorted(:,1:50);
    Gplot{n_c} = sortrows(Gplot{n_c}',1)';
end
        
figure
xlabel('Fb')
ylabel('Temperature (°C)')
zlabel('Profit')
contour(U1cont,U2cont,PROFIT',40)
hold on
plot(u_opt_true(1),u_opt_true(2),'rx','LineWidth',2)
plot(u_opt_approx(1),u_opt_approx(2),'rx','LineWidth',2)
for i = 1:n_c
    Gplota = Gplot{i};
    plot(Gplota(1,:),Gplota(2,:),'k--','LineWidth',2)
end
plot(U_GP0(1,:),U_GP0(2,:),'ro','LineWidth',2)
plot(U_GP(1,1+size(U_GP0,2):end),U_GP(2,1+size(U_GP0,2):end),'b-o','LineWidth',2)
xlim([lb(1),ub(1)]); ylim([lb(2),ub(2)]);  