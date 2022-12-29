function [x, i] = Bisection(a, b, delta, penalty,nelx,nely,connect,KE,alpha)
global c1 c2 c3
%% Inicialização
i = 1; 
d = b - a; %Search direction
xold = (a+b)/2;

move = 0.2;
while norm(b-a) > delta && i < 100

    x = (a+b)/2; %Ponto médio
    dens = max(0.001,max(xold-move,min(1,min(xold+move,x))));
    
    u=FiniteElement(nelx,nely,dens,connect,penalty,KE);
    [df] = DerivL(u,dens,penalty,nelx,nely,connect,KE,alpha);
    deriv = dot(df,d);

    if deriv > 0 %Se a derivada é positiva restringe-se o intervalo à esquerda
        b = dens;
    else % Se é negativa restringe-se à direita
        a = dens;
    end
    i = i+1;

    xold = dens;
    
end

    x = xold;
end