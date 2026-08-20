function F = odefun_CV(Y,t,gamma,zeta,Cs,ss,Cm,sm,parameters) %,kappa

Nmodes = parameters.NmodesSax + parameters.NmodesCV ;
pks = Y(1:parameters.NmodesSax) ;
pkm = Y(parameters.NmodesSax+1:Nmodes) ;
ps = 2*sum(real(pks)) ; pm = 2*sum(real(pkm)) ;
Dp = ps-pm ;
F = zeros(length(Y),1) ;

if parameters.dyn == 0 % sans dynamique :
    x = Dp-gamma ;  
    u = zeta/2*(x+1+sqrt((x+1)^2+parameters.epsi))*(gamma-Dp)/(((gamma-Dp)^2+parameters.epsi)^(1/4)); % - kappa*xdot/om_r +
else 
    x = Y(Nmodes + 1) ;
    xdot = Y(Nmodes + 2) ; 
    u = zeta/2*(x+1+sqrt((x+1)^2+parameters.epsi))*(gamma-Dp)/(((gamma-Dp)^2+parameters.epsi)^(1/4)); % - kappa*xdot/om_r +
    F(Nmodes + 1) = xdot ;
    F(Nmodes + 2) = - parameters.q_r*parameters.om_r*xdot+parameters.om_r^2.*(Dp-gamma-x) ; % Fc = 0, ghost reed
end

F(1:parameters.NmodesSax) = Cs.*u + ss.*pks ;
F(parameters.NmodesSax+1:Nmodes)= -parameters.Zms*Cm.*u + sm.*pkm ;

end
