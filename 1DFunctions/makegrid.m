function [z,dz,d3z,g1,g2,rho,kc] = grid(zmax,m,n,H,rho_s,rho_d,ks,kd,skinDepth)

z(1)  = 0; 
dz(1) = skinDepth / m; 
%zmax  = skinDepth * bn; 
%zmax  = 1.4;
% Or use 

i = 1;
while z(i) < zmax 
    i           = i+1;
    h           = dz(i-1)*(1+1/n);
    dz(i)       = h; 
    z(i)        = z(i-1) + dz(i);
end

dz = diff(z); 
d3z = dz(2:end).*dz(1:end-1).*(dz(2:end)+ dz(1:end-1));
g1 = 2*dz(2:end)./d3z; 
g2 = 2*dz(1:end-1)./d3z; 

% Turn row into column vector
dz  = dz';
d3z = d3z';
g1  = g1'; 
g2  = g2';
z   = z';

rho = zeros(length(z),1);
kc = zeros(length(z),1);
rho(1) = rho_s;
kc(1) = ks; 

for j=2:length(z)
    rho(j) = rho_d - (rho_d - rho_s)*exp(-z(j)/H);
    kc(j) = kd - (kd-ks)*(((rho_d - rho(j))/(rho_d-rho_s)));
end

disp("Layer variables:")
disp("     Model depth:  " + max(z) + " m")
disp("     Model layers: " + (length(z)-1))
end
