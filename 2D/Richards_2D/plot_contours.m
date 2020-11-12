function plot_contours(I,J,x,y,p,tht,Vx,Vy,levels)

figure;
contourf(x,y,p,levels); colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\psi(x,z,t)$','Interpreter','latex'); 

figure;
contourf(x,y,tht,levels); colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\theta(x,z,t)$','Interpreter','latex'); 

xx=x(2:I-1);
yy=y(2:J-1);
figure
contourf(xx,yy,Vx,levels);  colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_x(x,z)$','Interpreter','latex'); 

figure
contourf(xx,yy,Vy,levels);  colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_z(x,z)$','Interpreter','latex'); 
