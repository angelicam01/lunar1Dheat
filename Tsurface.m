function T0 = Tsurface(T,rho_s,Qs,K,epsxsigm,four_epsxsigm,twodz0)
% Newton's root finding method for calculating surface temperature. %
% See Apendix A of Hayne et al. (2017) for a complete description %

T0          = T(1);
[dKdT]      = derivk(T0,rho_s);
f           = epsxsigm*(T0^4) - Qs - K(1)*(-3*T0 + 4*T(2) - T(3))/(twodz0);
fprime      = four_epsxsigm*(T0^4) - (dKdT*((4*T(2) - 3*T0 - T(3))/(twodz0))) + ((3/(twodz0))*K(1));
dt          = -f/fprime;
T0          = T0 + dt; 
counter     =0;

while abs(f/fprime) >= 0.001    
    counter = counter + 1;
    K           = updateRK(T0,rho_s); 
    [dKdT]      = derivk(T0,rho_s);
    f           = epsxsigm*(T0^4) - Qs - K(1)*(-3*T0 + 4*T(2) - T(3))/(twodz0);
    fprime      = four_epsxsigm*(T0^4) - (dKdT*((4*T(2) - 3*T0 - T(3))/(twodz0))) + ((3/(twodz0))*K(1));
    dT          = -f/fprime; 
    T0          = T0 + dT;
    if counter > 50
        break
                
    end
end
