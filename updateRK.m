function [rkout] = updateRK(T,rho)

% Constants for temperature-dependent amorphous solid conductivity matrix found in Woods-Robinson et el. (2019)
Aam = -2.03297e-1;
Bam = -11.472; 
Cam = 22.5793;
Dam = -14.3084;
Eam = 3.41742;
Fam = 0.01101;
Gam = -2.80491e-5;
Ham = 3.35837e-8;
Iam = -1.40021e-11;
kam = Aam + Bam*(T.^-4.0) + Cam*(T.^-3.0) + Dam*(T.^-2.0) + Eam*(T.^-1.0) + Fam*(T) + Gam*(T.^2.0) + Ham*(T.^3.0) + Iam*(T.^4.0);

% Derived constants for density-dependent functions
A1 = (5.0821e-6);
A2 = (-0.0051);
B1 = 2.022e-13;
B2 = -1.953e-10;

rkout = ((A1*rho + A2)).*kam + (B1*rho + B2).*T.^3;

end
