function [coord,connect]=elementsQ4(nelx,nely,L1,L2)
%design the mesh according to inputs

%number of nodes on each boundary
nx = nelx+1; 
ny = nely+1;
dx = L1/nelx;
dy = L2/nely;

coord=zeros(nx*ny,2);
connect=zeros(nelx*nely,8);

% coordinates of the vertices
for j=1:nx
    i=1:ny;
    k = (j-1)*ny+i;
   coord(k,1) = (j-1)*dx;
   coord(k,2) = L2-(i-1)*dy;
end
% for i=nx*ny
%     
% end
% connectivity matrix (defining the elements)
% for j=1:nely
%     for i=1:nelx
%         k=(j-1)*nelx+i;
%         connect(k,:) = (j-1)*2*nx + [2*i-1 2*i 2*i+1 2*i+2 2*(i+nx)+1 2*(i+nx)+2 2*(i+nx)-1 2*(i+nx)];
%     end
% end
j=1;
for i=1:nely*nelx
    connect(i,:)=2*i+[2*j-1 2*j (2*j-3)+2*(nely+1)+2 2*j-3+2*(nely+1)+3 ...
         2*j-3+2*(nely+1) 2*j-3+2*(nely+1)+1 2*j-3 2*j-2];
    j=1+floor(i/nely);
end