%%% OESM PROJECT - P2
%%Rita Martins nº86804

% -------------------
% Augmented Lagrangean model 1
% -------------------

clear all
%% Input Parameters

%design space dimensions
L1=30;
L2=30;

%number of elements
nelx=30;
nely=30;
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
%         text(x(i)-0.02,y(i),num2str(nodesx(i)),'FontSize',10,'Color','b');
%         inputnodes=[nodesx(i),inputnodes];
%         end
%         if ~ismember(nodesy(i),inputnodes)
%         text(x(i)+0.02,y(i),num2str(nodesy(i)),'FontSize',10,'Color','r');
%         inputnodes=[nodesy(i),inputnodes];
%         end
%     end
% end
% xlabel('x')
% ylabel('y')
% title('mesh')
% hold on

%% Element Material property
%start with same density in all elements
rho(1:nelx*nely,1)=alpha;

iter=1;
penalty=3;
rmin=1.5;

E =1; nu = 0.3;

% Transform to plane stress properties
E1 = E/(1-nu^2);
E2=nu*E1;
G=E/(2*(1+nu));

%% Augmented Lagrangean parameters
global l1 l2 l3 c1 c2 c3
l1=5;
l2=5;
l3=5;
c1=1;
c2=100;
c3=100;

delta = 0.5*10^-1;
Kval=10^12;

while iter<100
%% Finite Element Analysis

%Assemble Stiffness Matrix
KE=K_Q4(E,nu,L1,L2);
[U,K]=FiniteElement(nelx,nely,rho,connect,penalty,KE);

%% Sensitivity Analysis
[df{1}] = DerivL(U,rho,penalty,nelx,nely,connect,KE,alpha);

%% Filtering
[df{1}] = filtering(nelx,nely,1.2,rho,df{1});

%% Augmented Lagrangean

c=U'*K*U;
rho_CG=rho;

% Conjugate Gradient
counter2=1;
j=counter2;
n=size(rho,1); %number of design variables
while counter2<100
    d=-df{1};
    d=d/norm(d);
    
    %line search
    [rho_2,count]=Bisection((rho_CG-d),(rho_CG+d),delta,penalty,nelx,nely,connect,KE,alpha);
    

    [U,K]=FiniteElement(nelx,nely,rho_2,connect,penalty,KE);
    df{2}=DerivL(U,rho_2,penalty,nelx,nely,connect,KE,alpha);
    [df{2}] = filtering(nelx,nely,1.2,rho,df{2});
    
    if norm(df{2}) < 0.01 %Critério de paragem
        break
    end
    
    if j<n+1
        beta=dot(df{1},df{1})/dot(df{2},df{2});
        df{1}=(df{1}+beta.*d);
    
    rho_CG = rho_2;
    clear rho_2;
    else
        j=1;
    end
    
    counter2=counter2+1;
    j=j+1;
end

    Kbar=max([abs(max(sum(rho-alpha),-l1/c1)),abs(max(sum(rho-1),-l2/c2)),abs(max(sum(-rho+0.001),-l3/c3))]);

c_new=U'*K*U;
dx = L1/nelx;
dy = L2/nely;
vol=sum(((dx*dy)/(L1*L2)).*rho_CG);
vol1=sum(rho_CG)/(nelx*nely); 

% Lagrange multipliers update
        l1=max(l1+c1*sum(rho-alpha),0);
        l2=max(l2+c2*sum(rho-1),0);
        l3=max(l3+c3*sum(-rho+0.001),-l3/c3);
% penalty update
    if min(sum(rho-alpha),0)>0.25*min(sum(rho_CG-alpha),0)
        c1=max(5*c1,iter^2);
    end

epsilon = abs(c_new-c)/c;
if epsilon < 0.01 && Kbar<0.01
    break
else
    rho = rho_CG;
end
    
    %%% outputs
disp([' N.iter.: ' sprintf('%3i',iter)...
' Compliance.: ' sprintf('%10.4f',c)...
' Vol. ' sprintf('%6.3f', vol)]);

iter=iter+1;

%plot density as image intensity
dens=zeros(nely,nelx);
for i=1:nelx
    dens(:,i)=rho((i-1)*nely+1:i*nely);
end
    %dens=flip(dens);
    dx = L1/nelx;
    dy = L2/nely;
    center1=[dx L2-dy];
    center2=[L1-dx dy];
    colormap(gray);
    imagesc(center1,center2,-dens);
    axis equal;
    axis tight;
    set(gca,'YDir','normal')
    pause(1e-6);
end