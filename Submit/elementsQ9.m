function [coord,connect]=elementsQ9(nelx,nely,L1,L2)
%design the mesh according to inputs

%number of nodes on each boundary
nx = 2*nelx+1; 
ny = 2*nely+1;
dx = L1/(2*nelx);
dy = L2/(2*nely);

coord=zeros(nx*ny,2);
connect=zeros(nelx*nely,18);

% coordinates of the vertices
for j=1:ny
    i=1:nx;
    k = (i-1)*ny+j;
   coord(k,1) = (i-1)*dx;
   coord(k,2) = L2-(j-1)*dy;
end

% connectivity matrix (defining the elements)
for j=1:nelx
    for i=1:nely
        k=(j-1)*nely+i;
        connect(k,:) = (j-1)*4*ny + [4*i+1 4*i+2 4*(i+ny)+1 4*(i+ny)+2 4*(i+ny)-3 4*(i+ny)-2 4*i-3 4*i-2 ...
                                     2*(2*i+ny)+1 2*(2*i+ny)+2 4*(i+ny)-1 4*(i+ny) 2*(2*i+ny)-3 2*(2*i+ny)-2 4*i-1 4*i ...
                                     2*(2*i+ny)-1 2*(2*i+ny)];
    end
end
% j=1;
% for i=1:nely*nelx
%     connect(i,:)= 4*i+[4*j-3 4*j-2 (4*j-3)+6*(nely+1)+5 4*j-2+6*(nely+1)+5 ...
%         4*j-7+6*(nely+1)+6 4*j-6+6*(nely+1)+6 4*j-7 4*j-6 (4*j-5)+3*(nely+1)+5 (4*j-4)+3*(nely+1)+5 ...
%         (4*j-5)+6*(nely+1)+6 (4*j-4)+6*(nely+1)+6 4*j-7+3*(nely+1)+3 4*j-6+3*(nely+1)+3 4*j-5 4*j-4 ...
%           (4*j-5)+3*(nely+1)+3 (4*j-4)+3*(nely+1)+3];
%     j=2+4*floor(i/nely);
% end