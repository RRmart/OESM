%%% OESM PROJECT - P2
%%Rita Martins nº86804
% -----------------------
% OC METHOD WITH MODEL 1
% -----------------------

clear all
%% Input Parameters

%design space dimensions
L1=30;
L2=30;

%number of elements
nelx=50;
nely=50;
nel=nelx*nely;

%Area upper bound constraint value
alpha=0.5;

%% Defining the mesh
[coord,connect]=elementsQ4(nelx,nely,L1,L2);

% plot the mesh (makes the code slower)
% figure()
% for i=(1:size(connect,1))
%     nodesx=connect(i,1:2:end);
%     nodesy=connect(i,2:2:end);
%     ord=[1,2,3,4,1];
%     x=zeros(size(ord,2),1);
%     y=zeros(size(ord,2),1);
%     for j=1:size(ord,2)
%         n=ord(j);
%         x(j)=coord(floor((nodesx(n)+1)/2),1);
%         y(j)=coord(floor((nodesy(n)+1)/2),2);
%     end
%     inputnodes=[];
%     plot(x,y)
%     hold on
%     text((x(1)+((x(4)-x(1))/2)), (y(1)+((y(3)-y(1))/2)), num2str(i),'FontSize',6);
%     for i=1:length(x)-1
%         if ~ismember(nodesx(i),inputnodes)
%         text(x(i)-1,y(i),num2str(nodesx(i)),'FontSize',10,'Color','b');
%         inputnodes=[nodesx(i),inputnodes];
%         end
%         if ~ismember(nodesy(i),inputnodes)
%         text(x(i)+1,y(i),num2str(nodesy(i)),'FontSize',10,'Color','r');
%         inputnodes=[nodesy(i),inputnodes];
%         end
%     end
% end
% xlabel('x')
% ylabel('y')
% title('mesh')
% hold off

%% Element Material property
%start with the same density in all elements
rho(1:nelx*nely,1)=alpha;

iter=1;
penalty=3;
rmin=3;

E = 1; nu = 0.3;

% Transform to plane stress properties
E1 = E/(1-nu^2);
E2=nu*E1;
G=E/(2*(1+nu));

while iter<100
%% Finite Element Analysis

K = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
F = sparse(2*(nely+1)*(nelx+1),1);

%Assemble Stiffness Matrix
KE=K_Q4(E,nu,L1,L2);
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

DispNodes=[1 0; 2 0; 2*(nely+1) 0; 2*(nely+1)-1 0];
DispValues=0;

Fixeddofs=DispNodes(:,1);
alldofs     = (1:2*(nely+1)*(nelx+1));
freedofs    = setdiff(alldofs,Fixeddofs);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(Fixeddofs,:)= DispValues;

%% Sensitivity Analysis
[df{1}] = Deriv(U,rho,penalty,nelx,nely,connect,KE);

%% Filtering
[df{1}] = filtering(nelx,nely,rmin,rho,df{1});

%% Optimatlity Criteria
rho_2 = OC(nelx,nely,rho,alpha,df{1});

c=U'*K*U;
dx = L1/nelx;
dy = L2/nely;
vol=sum(((dx*dy)/(L1*L2)).*rho);
vol1=sum(rho)/(nelx*nely);

epsilon = max(max(abs(rho_2-rho)));
if epsilon < 0.01
    break
else
    rho = rho_2;
end

%%% outputs
disp([' N.iter.: ' sprintf('%3i',iter)...
' Func.obj.: ' sprintf('%10.4f',c)...
' Vol. ' sprintf('%6.3f', vol)]);
iter=iter+1;


dens=zeros(nely,nelx);
for i=1:nelx
    dens(:,i)=rho_2((i-1)*nely+1:i*nely);
end
    %dens=flip(dens,1);
    dx = L1/nelx;
    dy = L2/nely;
    center1=[dx/2 L1-dx/2];
    center2=[L2-dy/2 dy/2];
    colormap(gray);
    imagesc(center1,center2,-dens);
    axis equal;
    axis tight;
    set(gca,'YDir','normal')
    pause(1e-6);
end