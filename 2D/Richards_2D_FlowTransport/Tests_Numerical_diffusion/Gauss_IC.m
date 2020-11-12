function [gss] = Gauss_IC(t,dx,dy,x0,y0,Nx,Ny,U_MEAN,D1,D2)

gss=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        gss(i,j)=(1/(4*pi*t*sqrt(D1*D2)))*exp(-(((i*dx-x0)^2)/(4*D1*t)+((j*dy-y0-U_MEAN*t)^2)/(4*D2*t)));
    end
end