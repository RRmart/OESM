function KE=KeQ9(ex,ey,C)
% element stiffness matrix - calculated using Gauss Quadrature rule with 3*3 points

KE = zeros(18,18);

% Gauss quadrature points: 
csi = [-sqrt(3/5) 0 sqrt(3/5) -sqrt(3/5) 0 sqrt(3/5) -sqrt(3/5) 0 sqrt(3/5)];
eta = [-sqrt(3/5) -sqrt(3/5) -sqrt(3/5) 0 0 0 sqrt(3/5) sqrt(3/5) sqrt(3/5)];

pg = [25/81 40/81 25/81 40/81 64/81 40/81 25/81 40/81 25/81]; %quadrature weights

for i=(1:9)

%element shape functions and derivatives calculated at Quadrature points
[N,dNdxi]=shape_funcQ9(csi(i),eta(i));

J0=[ex ey]'*dNdxi; %Jacobian Matrix
dNdx  = J0\dNdxi';

% N1=[N(1) 0 N(2) 0 N(3) 0 N(4) 0 ;...
%     0 N(1) 0 N(2) 0 N(3) 0 N(4)];

Ne=zeros(2,18);
Ne(1,1:2:end)=N(:);
Ne(2,2:2:end)=N(:);

% Be=[dNdx(1,1) 0 dNdx(1,2) 0 dNdx(1,3) 0 dNdx(1,4) 0; ...
%     0 dNdx(2,1) 0 dNdx(2,2) 0 dNdx(2,3) 0 dNdx(2,4); ...
%     dNdx(2,1) dNdx(1,1) dNdx(2,2) dNdx(1,2) dNdx(2,3) dNdx(1,3) dNdx(2,4) dNdx(1,4)];

Be = zeros(3,2*9); 
Be(1,1:2:2*9) = dNdx(1,:) ; 
Be(2,2:2:2*9) = dNdx(2,:) ; 
Be(3,1:2:2*9) = dNdx(2,:) ; 
Be(3,2:2:2*9) = dNdx(1,:) ;

Jacdet=det(J0); 

% Approximate the integrals that are non zero at triangle K:

KE=KE+pg(i)*(Be'*C*Be)*(Jacdet);
end

end