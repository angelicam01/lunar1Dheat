function T0 = TsurfaceHayne(T,kc,Qs,epsxsigm,four_epsxsigm,twodz0)
% Newton's root finding method for calculating surface temperature %
% See Apendix A of Hayne et al. (2017) for a complete description %

chi         = 2.7;  
B1          = (kc(1)*chi)/(350^3); 
T0          = T(1);
f           = epsxsigm*(T0^4) - Qs - (kc(1)+B1*(T0^3))*(-3*T0 + 4*T(2) - T(3))/(twodz0);
fprime      = four_epsxsigm*(T0^3) - (3*B1*(T0^2)*((4*T(2) - 3*T0 - T(3))/(twodz0))) + ((3/(twodz0))*(kc(1)+B1*(T0^3)));
dt          = -f/fprime;
T0          = T0 + dt;
counter     = 0; 

while abs(f/fprime) >= 0.001    
    counter = counter + 1;
    f       = epsxsigm*(T0^4) - Qs - (kc(1)+B1*(T0^3))*(-3*T0 + 4*T(2) - T(3))/(twodz0);
    fprime  = four_epsxsigm*(T0^3) - (3*B1*(T0^2)*((4*T(2) - 3*T0 - T(3))/(twodz0))) + ((3/(twodz0))*(kc(1)+B1*(T0^3)));
    dT      = -f/fprime;
    T0      = T0 + dT;
    if counter > 50
        break
    end
end

end
