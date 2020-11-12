close all;

load Richy_Pressure_scenario1
Rp1=Richy_Pressure_scenario1;
load Richy_Pressure_scenario2;
Rp2=Richy_Pressure_scenario2;

load Richy_WaterContent_scenario1;
Rt1=Richy_WaterContent_scenario1;
load Richy_WaterContent_scenario2;
Rt2=Richy_WaterContent_scenario2;

load Richy_WaterFlux_scenario1
Rf1=Richy_WaterFlux_scenario1;
load Richy_WaterFlux_scenario2
Rf2=Richy_WaterFlux_scenario2;

load solution_scenario1.mat
p1=p; tht1=tht; q1=q;
load solution_scenario2.mat
p2=p; tht2=tht; q2=q;

figure
plot(p1,x,'+c'); hold all; plot(p2,x,'.m');
plot(Rp1(:,2),(Rp1(:,1)+1),'--k'); plot(Rp2(:,2),(Rp2(:,1)+1),'.-k'); 
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); set(gca, 'XTick', [-0.8:0.2:0.6]); 
legend('GRW (1)','GRW (2)','MFEM (1)','MFEM (2)','Location','best'); legend('boxoff');

figure
plot(tht1,x,'+c'); hold all; plot(tht2,x,'.m');
plot(Rt1(:,2),(Rt1(:,1)+1),'--k'); plot(Rt2(:,2),(Rt2(:,1)+1),'.-k'); 
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex');  set(gca, 'XTick', [0.05:0.05:0.4]); 
legend('GRW (1)','GRW (2)','MFEM (1)','MFEM (2)','Location','best'); legend('boxoff');

figure
plot(q1,x2,'+c'); hold all; plot(q2,x2,'.m');
plot(Rf1(:,2),(Rf1(:,1)+1),'--k'); plot(Rf2(:,2),(Rf2(:,1)+1),'.-k'); 
ylabel('$z$','Interpreter','latex'); 
xlabel('$q(z,t)$','Interpreter','latex')
legend('GRW (1)','GRW (2)','MFEM (1)','MFEM (2)','Location','best'); legend('boxoff');
 
figure 
load dt_scenario1
tdt=dt_scenario1;
plot(tdt(:,1),tdt(:,2),'+c'); hold;
load dt_scenario2
tdt=dt_scenario2;
plot(tdt(:,1),tdt(:,2),'.m');
set(gca,'yscale','log'); set(gca,'xscale','log'); ylim([0.1 1e2]); xlim([0 1e4])
xlabel('$t$','Interpreter','latex'); ylabel('$\Delta t$','Interpreter','latex')
legend('Scenario (1)','Scenario (2)','Location','best'); legend('boxoff');

format shortE
norm_c1=norm(p1'-Rp1(:,2))/norm(Rp1(:,2))
norm_theta1=norm(tht1(1:end-1)'-Rt1(1:2:end,2))/norm(Rt1(1:2:end,2))
norm_q1=norm(q1'-Rf1(2:end,2))/norm(Rf1(2:end,2))

norm_c2=norm(p2'-Rp2(:,2))/norm(Rp2(:,2))
norm_theta2=norm(tht2(1:end-1)'-Rt2(1:2:end,2))/norm(Rt2(1:2:end,2))
norm_q2=norm(q2'-Rf2(2:end,2))/norm(Rf2(2:end,2))
