function plot_fig(x,p,c,tht,q)

figure;
plot(p,x); 
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); 

figure;
plot(c,x); 
ylabel('$z$','Interpreter','latex');
xlabel('$c(z,t)$','Interpreter','latex'); 

figure;
plot(tht,x); 
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 

figure
plot(q,x)
ylabel('$z$','Interpreter','latex'); 
xlabel('$q(z,t)$','Interpreter','latex'); 
