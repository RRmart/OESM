function [df] = Deriv(u,dens,p,nelx,nely,connect,KE)
E = 200*10^3; nu = 0.3;

% Transform to plane stress properties
E1 = E/(1-nu^2);
E2=nu*E1;
G=E/(2*(1+nu));

%KE=K_TP4(E,nu,L1,L2);
for iel = 1:nelx*nely
%         nodesx = connect(iel,1:2:end);    % for element iel keep vertices nodes
%         nodesy = connect(iel,2:2:end);
%         idx = floor((nodesx+1)/2);
%         ex     = coord(idx,1);       % get the x-coord from the above vertices
%         ey     = coord(idx,2);       % get the y-coord from the above vertices
%     
%         C=[E1 E2 0; E2 E1 0; 0 0 G];
        %Ke=KE.iel;
        dc(iel,1) = u(connect(iel,:))'*((-p*dens(iel)^(p-1).*KE))*u(connect(iel,:));
end
df = dc;