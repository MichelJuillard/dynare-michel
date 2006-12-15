function forcs = mcforecast3(cL,H,mcValue,shocks,forcs,T,R,mv,mu)
% stephane.adjemian@ens.fr [06-11-2006]

if cL
    e = zeros(size(mcValue,1),cL);
    for t = 1:cL
        e(:,t) = inv(mv*R*mu)*(mcValue(:,t)-mv*T*forcs(:,t)-mv*R*shocks(:,t));
        forcs(:,t+1) = T*forcs(:,t)+R*(mu*e(:,t)+shocks(:,t));
    end
end
for t = cL+1:H
    forcs(:,t+1) = T*forcs(:,t)+R*shocks(:,t);
end