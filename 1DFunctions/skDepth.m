function skinDepth  = skDepth(delta,S,rho_s,ks,P,latitude)
phi         = (pi/180)*latitude; % Conversion to radians 
A           = 0.12; % Bond albedo at normal solar incidence 
a           = 0.06; % Constant for Bond albedo 
b           = 0.25; % Contstant for Bond albedo 
epsilon     = 0.95; %Infrared emissivity 
sigma       = 5.67051196e-8; % Stefan-Boltzman constant 

%%%% Initial temperature based on solar constant %%%%
h1          = 0;
sola        = acos((sin(phi)*sin(delta))+(cos(phi)*cos(delta)*cos(h1)));
cossa1      = (sin(phi)*sin(delta)) + (cos(phi)*cos(delta)*cos(h1));
solFinit    = S*cos(phi);
albedoinit  = A + a*(sola/(pi/4))^3 + b*(sola/(pi/2))^8;
initTemp    = (((1-albedoinit)/(epsilon*sigma))*solFinit)^(1.0/4.0);

%%%% Initial specific heat capacity %%%%
c0      = -3.6125; % [J/(kg*K)]
c1      = 2.7431; % [J/(kg*K)]
c2      = 2.3616e-3; % [J/(kg*K)]
c3      = -1.2340e-5; % [J/(kg*K)]
c4      = 8.9093e-9; % [J/(kg*K)]
initCap = c0 + c1*initTemp + c2*(initTemp^2) + c3*(initTemp^3) + c4*(initTemp^4);


%%%%%% Surface phonon conductivity needs work %%%%%%%
%%%% Calculate surface phonon conductivity %%%%
%rhobulk1    = rho_s;
%rhosolid1   = 2800.0;
%p1           = 1 - (rhobulk1/rhosolid1);
%B1a         = 9.9e-4;
%C1a         = 9.2e-10;
%Aam1        = -2.03297e-1;
%Bam1        = -11.472; 
%Cam1        = 22.5793;
%Dam1        = -14.3084;
%Eam1        = 3.41742;
%Fam1        = 0.01101;
%Gam1        = -2.80491e-5;
%Ham1        = 3.35837e-8;
%Iam1        = -1.40021e-11;
%kam1        = Aam1 + Bam1*(initTemp^(-4.0)) + Cam1 *(initTemp^(-3.0)) + Dam1*(initTemp^(-2.0)) + Eam1*(initTemp^(-1.0)) + Fam1*initTemp + Gam1*(initTemp^(2.0)) + Ham1*(initTemp^(3.0)) + Iam1*(initTemp^(4.0));
%kc_surf     = (B1a + C1a*(rhobulk1^2))*(1-p1)*kam1;

%%%% Final calculation 
kappa = ks/(1100*initCap); % Thermal diffusivity 
skinDepth = sqrt((kappa*P)/pi); % Skin depth. Before each model run, skin depth is calculated based on the expected temperature range



end

