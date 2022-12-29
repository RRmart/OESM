function [N,dNdxi]=shape_funcQ4(xi,eta)

N=[(1/4)*(1-xi)*(1-eta);
    (1/4)*(1+xi)*(1-eta);
    (1/4)*(1+xi)*(1+eta);
    (1/4)*(1-xi)*(1+eta);];
dNdxi=[(1/4)*(-1+eta), (1/4)*(-1+xi);
    (1/4)*(1-eta), (1/4)*(-1-xi); 
    (1/4)*(1+eta), (1/4)*(1+xi);
    (1/4)*(-1-eta), (1/4)*(1-xi);];
end