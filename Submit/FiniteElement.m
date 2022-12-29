function [U,K] = FiniteElement(nelx,nely,rho,connect,penalty,KE)

K = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
F = sparse(2*(nely+1)*(nelx+1),1);

for iel=1:nely*nelx       

    K(connect(iel,:),connect(iel,:)) = K(connect(iel,:),connect(iel,:)) + rho(iel)^penalty*KE;

end

%vector of applied concentrated forces - [node number, force magnitude]
AppF=[2*(nely+1)*(nelx)+2, -1];
%add concentrated forces
for i=1:size(AppF,1)
    F(AppF(i,1))=F(AppF(i,1))+AppF(i,2);
end

%Displacement boundary conditions
U = zeros(2*(nely+1)*(nelx+1),1);

%Displacement boundary conditions
DispNodes=[1 0; 2 0; 2*(nely+1) 0; 2*(nely+1)-1 0];
DispValues=0;

Fixeddofs=DispNodes(:,1);
alldofs     = (1:2*(nely+1)*(nelx+1));
freedofs    = setdiff(alldofs,Fixeddofs);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(Fixeddofs,:)= DispValues;

end