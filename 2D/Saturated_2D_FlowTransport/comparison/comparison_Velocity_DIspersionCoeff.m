%% GRW "ensemble" velocities and dispersion coefficients 
%% compared to 1-st order (linear) approximation results
%  close all; clear all;
 
D=0.01; 
%% mean velocity of center of mass 
figure;
load ensemble_coefficients
plot(t*dt,Mx,'r',t*dt,My,':b','LineWidth',1.5); hold all; ylim([-0.4 1.2]);
load linear_approximation
plot(t*dt,Mx,'.k',t*dt,My,'.k','MarkerSize',12);  
legend('$V_x(t)$','$V_y(t)$','Interpreter','latex','Location','best'); legend('boxoff');
xlabel('$t$','Interpreter','latex'); ylabel('mean velocity components'); set(gca, 'XTick', [0:2:10])
%% effective dispersion coefficients
figure
load ensemble_coefficients
plot(t*dt,Dx/D,'r',t*dt,Dy/D,':b','LineWidth',1.5); hold all;
load linear_approximation
plot(t*dt,Dx/D,'.k',t*dt,Dy/D,'.k','MarkerSize',12);
xlabel('$t$','Interpreter','latex'); ylabel('Dispersion coefficients'); set(gca, 'XTick', [0:2:10])
legend('$D_x(t) / D$','$D_y(t)/ D$','Interpreter','latex','Location','best'); legend('boxoff');



