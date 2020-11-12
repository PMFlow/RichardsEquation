%% Contours of solutions to the coupled flow-transport in the (x,z) plane
% (random Ksat case)

close all; clear all;
itest=0; % 0=loam; 1=clay; 
levels=12;
if itest == 1
    disp('Beit Netofa clay')
    load solution_clay
else
    disp('slit loam')
    load solution_loam
end

figure;
contourf(x,y,p,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\psi(x,z,t)$','Interpreter','latex'); 

figure;
contourf(x,y,c,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$c(x,z,t)$','Interpreter','latex'); 

figure;
contourf(x,y,tht,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\theta(x,z,t)$','Interpreter','latex'); 

figure
contourf(x,y,Vx,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_x(x,z)$','Interpreter','latex'); 

figure
contourf(x,y,Vy,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_z(x,z)$','Interpreter','latex'); 
