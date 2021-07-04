function [Qs,solF,albedo]  = insolationFeng2(dt,delta,phi,P,S,surface)


% Calculate surface illumination conditions
%%% Calculates albedo using the derived constants from Feng
%%% et al. (2020). A0 = 0.12 for highland and 0.07 for Mare 

%%%%%%% surface = 1 for highlands; 
%%%%%%% surface = 0 for mare; 

if surface == 1
    A0 = 0.12; 
else
    A0 = 0.07; 
end

t           = 0:dt:P;
h           = 2*pi*t/P;  % Hour angle. Has size of total amounts of steps.
costh       = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h); 
th          = acos(costh); % Solar incidence angle 
clipfunc    = 0.5*(costh + abs(costh));
albedo      = A0 + (1 - (costh).^(0.2752)); 
solF        = S*clipfunc; 
Qs          = (1-albedo).*solF;

end
