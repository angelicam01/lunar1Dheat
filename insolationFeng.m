function [Qs,solF,albedo]  = insolationFeng(dt,delta,phi,P,S,avgLola)
 
% Calculate surface illumination conditions 
% See Apendix A of Hayne et al. (2017) and Feng et al. (2020) for bolometric Bond
% albedo model

t           = 0:dt:P;
h           = 2*pi*t/P;  % Hour angle. Has size of total amounts of steps.
costh       = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h); 
th          = acos(costh); % Solar incidence angle 
clipfunc    = 0.5*(cos(th) + abs(costh));
albedo      = 0.39*avgLola+ (1 - (costh).^(0.2752)); 
solF        = S*clipfunc; 
Qs          = (1-albedo).*solF;

end
