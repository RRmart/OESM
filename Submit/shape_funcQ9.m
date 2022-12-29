function [N,dNdxi]=shape_funcQ9(xi,eta)

N=[(1/4)*(xi^2-xi)*(eta^2-eta);
    (1/4)*(xi^2+xi)*(eta^2-eta);
    (1/4)*(xi^2+xi)*(eta^2+eta);
    (1/4)*(xi^2-xi)*(eta^2+eta);
    (1/2)*(1-xi^2)*(eta^2-eta);
    (1/2)*(xi^2+xi)*(1-eta^2);
    (1/2)*(1-xi^2)*(eta^2-eta);
    (1/2)*(xi^2-xi)*(1-eta^2);
    (1-xi^2)*(1-eta^2)];
dNdxi=(1/4)*[eta*(2*xi-1)*(eta-1), xi*(xi-1)*(2*eta-1);
            eta*(2*xi+1)*(eta-1), xi*(xi+1)*(2*eta-1);
            eta*(2*xi+1)*(eta+1), xi*(xi+1)*(2*eta+1);
            eta*(2*xi-1)*(eta+1), xi*(xi-1)*(2*eta+1);
            -4*xi*eta*(eta-1), -2*(xi+1)*(xi-1)*(2*eta-1);
            -2*(2*xi+1)*(eta+1)*(eta-1), -4*xi*eta*(xi+1); 
            -4*xi*eta*(eta+1), -2*(xi+1)*(xi-1)*(2*eta+1);
            -2*(2*xi-1)*(eta+1)*(eta-1), -4*xi*eta*(xi-1); 
            8*xi*(eta^2-1),8*eta*(xi^2-1)];
end