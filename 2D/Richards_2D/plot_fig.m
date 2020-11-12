function plot_fig(I,J,x,y,p,tht,Vx,Vy)

figure;
mesh(x,y,p);
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;

figure;
mesh(x,y,tht);
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\theta(x,z,t)$','Interpreter','latex'); view(115,15); %view(-45,25); %(20,50);
grid on ;

figure
mesh(x(2:I-1),y(2:J-1),Vx)
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
zlabel('$q_x(x,z)$','Interpreter','latex'); view(115,15);

figure
mesh(x(2:I-1),y(2:J-1),Vy)
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
zlabel('$q_z(x,z)$','Interpreter','latex'); view(115,15);
