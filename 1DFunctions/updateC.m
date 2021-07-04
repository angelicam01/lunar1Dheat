function [cOut]= updateC(temp)
% Heat capacity model found in Appendix A of Hayne et al. (2017)
c0 = -3.6125; % [J/(kg*K)]
c1 = 2.7431; % [J/(kg*K)]
c2 = 2.3616e-3; % [J/(kg*K)]
c3 = -1.2340e-5; % [J/(kg*K)]
c4 = 8.9093e-9; % [J/(kg*K)]
cOut=c0+c1*temp+c2*temp.^2+c3*temp.^3+ c4*temp.^4;

end
