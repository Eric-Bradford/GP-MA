function sys = Plant_model(x,u)
    
    % Control inputs
    Fb = u(1);
    Tr = u(2);

    % Parameters
    Fa = 1.8275;
    Mt = 2105.2;
    
    % Estados do sistema
    Xa= x(1); % Mass fraction
    Xb= x(2); % Mass fraction
    Xc= x(3); % Mass fraction
    Xe= x(4); % Mass fraction
    Xp= x(5); % Mass fraction
    Xg= x(6); % Mass fraction

    % Constantes cinéticas
    k1 = 1.6599e6 * exp(-6666.7/(Tr+273.15));
    k2 = 7.2117e8 * exp(-8333.3/(Tr+273.15));
    k3 = 2.6745e12* exp(-11111 /(Tr+273.15));
    
    % taxas de reação
    r1 = k1*Xa*Xb*Mt;
    r2 = k2*Xb*Xc*Mt;
    r3 = k3*Xc*Xp*Mt;
    
    % Banço mássico
    Fr = Fa + Fb;
    
    sys(1,1) = (Fa -   r1                 - Fr*Xa )/Mt;
    sys(2,1) = (Fb -   r1 -   r2          - Fr*Xb )/Mt;
    sys(3,1) = (   + 2*r1 - 2*r2 -     r3 - Fr*Xc )/Mt;
    sys(4,1) = (          + 2*r2          - Fr*Xe )/Mt;
    sys(5,1) = (          +   r2 - 0.5*r3 - Fr*Xp )/Mt;
    sys(6,1) = (                 + 1.5*r3 - Fr*Xg )/Mt;
end