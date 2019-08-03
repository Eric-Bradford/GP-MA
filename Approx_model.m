function sys = Approx_model(x,u)

    % Control inputs
    Fb = u(1);
    Tr = u(2);
    
    % Parameters
    phi1 = -3;
    psi1 = -17;
    phi2 = -4;
    psi2 = -29;
    Fa = 1.8275;
    Mt = 2105.2;
    
    % Referency temperature 
     Tref = 110 + 273.15; %[=] K.
     
    % Estados do sistema
     Xa= x(1); % Mass fraction
     Xb= x(2); % Mass fraction
     Xe= x(3); % Mass fraction
     Xp= x(4); % Mass fraction
     Xg= x(5); % Mass fraction
    
    %Reparametrization
     k1 = exp(phi1) * exp( (Tref/(Tr+273.15)-1) * psi1);
     k2 = exp(phi2) * exp( (Tref/(Tr+273.15)-1) * psi2);

    %taxas de reação
     r1 = k1*Xa*Xb*Xb*Mt;
     r2 = k2*Xa*Xb*Xp*Mt;
    
    %Mass balance
     Fr = Fa + Fb;
    
     sys(1,1) = Fa -   r1 -   r2 - Fr*Xa;
     sys(2,1) = Fb - 2*r1 -   r2 - Fr*Xb;
     sys(3,1) =    + 2*r1        - Fr*Xe;
     sys(4,1) =    +   r1 -   r2 - Fr*Xp;
     sys(5,1) =           + 3*r2 - Fr*Xg;
end