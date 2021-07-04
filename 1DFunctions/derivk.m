function [dKdT] = derivk(T,rho)
A1 = (5.0821e-6);
A2 = (-0.0051);
B1 = 0.944*(2.121e-13); 
B2 = -1.953e-10;

Aam = -2.03297e-1;
Bam = -11.472; 
Cam = 22.5793;
Dam = -14.3084;
Eam = 3.41742;
Fam = 0.01101;
Gam = -2.80491e-5;
Ham = 3.35837e-8;
Iam = -1.40021e-11;
dkam =  - 4*Bam*(T.^-5) - 3*Cam*(T.^-4) - 2*Dam*(T.^-3) - Eam*(T.^-2) + Fam + 2*Gam*T + 3*Ham*(T.^2) + 4*Iam*(T.^3); 
dKdT = ((A1*rho + A2)).*dkam + 3*(B1*rho + B2).*T.^2; 


end
