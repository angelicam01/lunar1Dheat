function [kout] = updateK(temp,kc)
% Thermal conductivity model found in Appendix A of Hayne et sl. (2017)
chi     = 2.7;
kout    = kc.*(1+chi.*(temp/350).^3);
end
