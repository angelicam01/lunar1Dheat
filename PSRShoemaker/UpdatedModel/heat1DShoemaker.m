
load('shoemakerIllumination.mat','IRillumination','visibleillumination','juliandate'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Model Parameters %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latitude    = 0; % Latitude [degrees]
phi         = (pi/180)*latitude; % Latitude [radians]
delta       = (pi/180)*0.0; % Solar declination angle 
dt          = 50; % Time step [s]
P           = 2.55024e6;   % Diurnal period [s]
S           = 1361.0; % Solar constant [W/m^2]
Q           = 0.018; % Interior heat flow [W/m^2]
zmax        = 2.5; % Maximum depth of grid
m           = 20; % Grid parameter
n           = 30; % Grid parameter
jd          = juliandate(1):dt/86400:juliandate(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Thermophysical parameters %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_s      = 1100.0; % Surface layer density [kg/m^3]
rho_d      = 1800.0; % Deep layer density [kg/m^3]
kd         = 3.4e-3; % Deep layer conductivity [W/(m*K)] 
ks         = 7.4e-4; % Surface layer conductivity [W/(m*K)]
H          = 0.056; % H-parameter, density scale factor [m]
epsilon    = 0.95; %Infrared emissivity 
sigma      = 5.67051196e-8; % Stefan-Boltzman constant 



%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1D Heat Flow %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
skinDepth               = skDepth(delta,S,rho_s,ks,P,latitude);

% Initialize grid
[z,dz,d3z,g1,g2,rho,kc] = makegrid(zmax,m,n,H,rho_s,rho_d,ks,kd,skinDepth);

% Define surface illumination conditions calculated through ray tracing
% algorithm
Qs                      = IRillumination+visibleillumination;
Qsnew                   = spline(juliandate,Qs,jd); % for interpolation


%%%%%% Temperature initialization %%%%%% 
load('1DShoemakerUpdated_30DayProfile.mat','temperature');
T                = temperature(:,1); 
Tnew             = T; 
temperature      = zeros(length(z),length(Qsnew));  
temperature(:,1) = T; 


% Variables used for efficeincy
epsxsigm        = epsilon*sigma;     
four_epsxsigm   = 4*epsilon*sigma;   
three_over2dz0  = 3/(2*dz(1));
twodz0          = 2*dz(1);


% Variables for equilibrium condition
equiltemp1       = temperature(:,1);
equiltemp2       = equiltemp1*1000;
count            = 0;
equilcondition   = 0.0001; % Equilibrium condition to be used in while loop

% Perform one-dimensional heat transfer calculation until model results between to consecutive runs varies by less than the number defined in equilcondition.


while any(abs(equiltemp1-equiltemp2) >= equilcondition)
equiltemp1       = temperature(:,1); 
  for i=1:length(Qsnew) 
       K                = updateRK(T,rho);
       cp               = updateC(T);  
       Tnew(2:end-1)    = Tlayers(T,dt,rho,cp,K,g1,g2);
       Tnew(1)          = Tsurface(T(1:3),rho_s,Qsnew(i),K,epsxsigm,four_epsxsigm,twodz0);
       Tnew(end)        = T(end-1) + Q/K(end-1)*dz(end);
       T                = Tnew; 
       temperature(:,i) = Tnew;  
  end 
  count            = count + 1;  
  equiltemp2       = temperature(:,1);
  disp("Day:" + count)     
end


filename = sprintf("1DShoemakerUpdated_Fullrun.mat");
save(filename)







