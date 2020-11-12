function [gss] = Gauss_IC(tt,dx,dy,x0,y0,Nx,Ny,U_MEAN,Dd)

gss=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        gss(i,j)=(1/sqrt(4*pi*Dd*tt))*exp(-(((i*dx-x0-U_MEAN*tt)^2)+(j*dy-y0)^2)/(4*Dd*tt));
    end
end