function F = Plant_data(u)

Xa0  = 0.2;   % Initial condition Mass fraction
Xb0  = 0.2;   % Initial condition Mass fraction
Xe0  = 0.2;   % Initial condition Mass fraction
Xp0  = 0.2;   % Initial condition Mass fraction
Xg0  = 0.2;   % Initial condition Mass fraction
Xs0  = 0.2;   % Initial condition Mass fraction
x0 =[Xa0;Xb0;Xe0;Xp0;Xg0;Xs0];

options = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
x       = fsolve(@(x) Plant_model(x,u),x0,options);                                                            
Fa      = 1.8275; Fb = u(1); Fr = Fa + Fb;

Obj  = -(1043.38*x(5,1)*Fr + 20.92*x(4,1)*Fr - 79.23*Fa - 118.34*u(1));
G(1) = x(1,1) - 0.12;
G(2) = x(6,1) - 0.08;
F    = [Obj,G];
end