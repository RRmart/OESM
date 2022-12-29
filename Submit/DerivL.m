function [df] = DerivL(u,dens,p,nelx,nely,connect,KE,alpha)

global l1 l2 l3 c1 c2 c3

for iel = 1:nelx*nely

        dc(iel,1) = u(connect(iel,:))'*((-p*dens(iel)^(p-1).*KE))*u(connect(iel,:));
end

df = dc + max(0, l1+c1*sum(dens-alpha));