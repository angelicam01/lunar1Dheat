function Tnew = Tlayers(T,dt,rho,cp,K,g1,g2)

alpha            = g1.*K(1:end-2);
beta             = g2.*K(2:end-1);
Tnew             = T(2:end-1) + (dt ./ (rho(2:end-1) .*cp(2:end-1))).* (alpha.*T(1:end-2) - ((alpha + beta).* T(2:end-1)) + (beta.* T(3:end)));    

end
