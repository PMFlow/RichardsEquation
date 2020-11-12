function plot_fig(x,p,c,tht,q)

figure;
plot(x,p); 
xlabel('$z$','Interpreter','latex');
ylabel('$\psi(z,t)$','Interpreter','latex'); 

figure;
plot(x,c); 
xlabel('$z$','Interpreter','latex');
ylabel('$c(z,t)$','Interpreter','latex'); 

figure;
plot(x,tht); 
xlabel('$z$','Interpreter','latex');
ylabel('$\theta(z,t)$','Interpreter','latex'); 

figure
plot(x,q)
xlabel('$z$','Interpreter','latex'); 
ylabel('$q(z,t)$','Interpreter','latex'); 

load test_dt
figure
plot(test_dt(:,1),test_dt(:,2));
xlabel('$t$','Interpreter','latex'); 
ylabel('$\Delta t(t)$','Interpreter','latex'); 
