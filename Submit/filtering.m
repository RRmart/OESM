function [dcn]=filtering(nelx,nely,rmin,x,dc)
dcn=zeros(nely*nelx,1);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j+(i-1)*(nely),1) = dcn(j+(i-1)*(nely),1) + max(0,fac)*x(l+(k-1)*(nely),1)*dc(l+(k-1)*(nely),1);
      end
    end      
    dcn(j+(i-1)*(nely)) = dcn(j+(i-1)*(nely))/(x(j+(i-1)*(nely))*sum);
  end
end