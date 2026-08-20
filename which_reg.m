function reg = which_reg(pk)

% à ajouter : calcul du ZCR pour régime quasi-périodique et périodicité

epsi1 = 1e-1 ;                             % seuil equilibre/oscillation
[~,kmax] = max(rms(pk)) ;
pmax = pk(:,kmax) ; 
pRMS = rms(pmax) ;

% NTs = diff(find(diff(sign(sum(real(pk),2)))>0));
% disp(var(NTs)) ; 

if pRMS < epsi1
   reg = 0 ;
else 
   reg = kmax ;
end

end