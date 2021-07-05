function [Qs]  = insolationcrater(dt,delta,P,S,latitude,D)

phi         = (pi/180)*latitude; % Latitude [radians]
t           = 0:dt:P;
A           = 0.12; % Bond albedo at normal solar incidence
h          = 2*pi*t/P; 
emissivity = 0.95;
costh      = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h); 
Qs         = (S*(4*emissivity*(1-A))/D^2)*(1+(A/emissivity))*(1/2*(costh + abs(costh)));

end
