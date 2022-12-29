function [Ke] = K_Q4(E,v,L1,L2)
gamma=L1/L2;
K1=(1+v)*gamma;
K2=(1-3*v)*gamma;
K3=2+(1-v)*gamma^2;
K4=2*gamma^2+(1-v);
K5=(1-v)*gamma^2-4;
K6=(1-v)*gamma^2-1;
K7=4*gamma^2-(1-v);
K8=gamma^2-(1-v);
Ke= E/(24*gamma*(1-v^2)).* [4*K3 3*K1 2*K5 -3*K2 -2*K3 -3*K1 -4*K6 3*K2;...
    3*K1 4*K4 3*K2 4*K8 -3*K1 -2*K4 -3*K2 -2*K7; ...
    2*K5 3*K2 4*K3 -3*K1 -4*K6 -3*K2 -2*K3 3*K1; ...
    -3*K2 4*K8 -3*K1 4*K4 3*K2 -2*K7 3*K1 -2*K4; ...
    -2*K3 -3*K1 -4*K6 3*K2 4*K3 3*K1 2*K5 -3*K2; ...
    -3*K1 -2*K4 -3*K2 -2*K7 3*K1 4*K4 3*K2 4*K8; ...
    -4*K6 -3*K2 -2*K3 3*K1 2*K5 3*K2 4*K3 -3*K1; ...
    3*K2 -2*K7 3*K1 -2*K4 -3*K2 4*K8 -3*K1 4*K4];
end