% One-dimensional heat transfer model for shadowed near-polar bowl-shaped craters
% with diameter-to-depth ratio, D.
% This script uses the standard Hayne et al. (2017) thermal conductivity model


function [temperature,P,totalsteps,z,D] = heat1DcraterStandard(latitude,D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi         = (pi/180)*latitude; % Latitude [radians]
delta       = (pi/180)*0.0; % Solar declination angle 
dt          = 150; % Time step [s]
P           = 2.55024e6;   % Diurnal period [s]
S           = 1361.0; % Solar constant [W/m^2]
Q           = 0.018; % Interior heat flow [W/m^2]
totalsteps  = P/dt; % Amount of steps to take
zmax        = 2.0; % Maximum depth of grid 
m           = 20; % Grid parameter
n           = 30; % Grid parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Thermophysical parameters %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_s      = 1100.0; % Surface layer density [kg/m^3]
rho_d      = 1800.0; % Deep layer density [kg/m^3]
kd         = 3.8e-3; % Deep layer conductivity [W/(m*K)] 
ks         = 8.0e-4;  % Surface layer conductivity [W/(m*K)]
H          = 0.06; % H-parameter, density scale factor [m]
epsilon    = 0.95; %Infrared emissivity 
sigma      = 5.67051196e-8; % Stefan-Boltzman constant 



%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1D Heat Flow %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
skinDepth               = skDepth(delta,S,rho_s,ks,P,latitude);

% Initialize grid
[z,dz,d3z,g1,g2,rho,kc] = makegrid(zmax,m,n,H,rho_s,rho_d,ks,kd,skinDepth);

% Calculate surface illumination conditions
[Qs]                    = insolationcrater(dt,delta,P,S,latitude,D);

%%%%%% Temperature initialization %%%%%% 
T          = z*0; 
Tnew       = T; 
T(1:end)   = 250;

% Variables used for efficeincy
epsxsigm        = epsilon*sigma;     
four_epsxsigm   = 4*epsilon*sigma;   
twodz0          = 2*dz(1);

% Intialize temperature matrix
temperature = zeros(length(z),round(totalsteps)-1); 
temperature(:,1) = T;


% Variables for equilibrium condition
equiltemp1 = temperature(:,1);
equiltemp2 = equiltemp1*1000;
count = 0;
equilcondition = 0.000001; % Equilibrium condition to be used in while loop

% Perform one-dimensional heat transfer calculation until model results between to consecutive runs varies by less than the number defined in equilcondition (currently 1e-6 [K]).
% This value may be increased, as this current number is chosen conservatively.

while any(abs(equiltemp1-equiltemp2) >= equilcondition)
equiltemp1       = temperature(:,1); 
  for i=1:totalsteps 
       K                = updateK(T,kc);
       cp               = updateC(T); 
       Tnew(2:end-1)    = Tlayers(T,dt,rho,cp,K,g1,g2);
       Tnew(1)          = TsurfaceHayne(T(1:3),kc,Qs(i),epsxsigm,four_epsxsigm,twodz0);
       Tnew(end)        = T(end-1) + Q/K(end-1)*dz(end);
       T                = Tnew; 
       temperature(:,i) = Tnew; 
      
   end
    count            = count + 1; 
    equiltemp2       = temperature(:,1);
    disp("Day:" + count)
       
end
end
