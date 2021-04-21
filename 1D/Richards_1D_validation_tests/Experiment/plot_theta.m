%% Comparison Experiment-Hydrus 1D-GRW solution
%% GRW solution
I= 61; 
a=0; b=6; 
dx=(b-a)/(I-1);
x=a:dx:b;
xb=(x-b)*100;
load ('thT','thT')
load ('pT','pT')
%% HYDRUS 1D solution for  drainage in a large caisson
xH=0:5:600; 
xH=-xH;
load ('thH','thH')
load ('pH','pH')
%% theta - experimantal data from [Zambra et al., 2012, Fig. 2]
xE=[-345 -270 -190 -115 -45];
xE2=[-270 -190 -115 -45];
thE1=[0.32 0.32 0.312 0.3 0.282];
thE2=[0.28 0.262 0.257 0.24];
thE3=[0.24 0.23 0.215 0.21 0.198];
thE4=[0.202 0.190 0.178 0.177 0.165];
%% p - experimantal data from [Zadeh, 2011, Fig. 4]
pE1=[0 0 -20 -25 -26];
pE2=[-45 -60 -80 -90];
pE3=[-105 -125 -130 -160 -190];
pE4=[-195 -230 -228 -265 -320];

%% theta:

figure; hold all;
for k=1:jT_plot
P(k)=plot(thT(k,:),xb,'MarkerSize',4); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); xlim([0.1 0.35]);
legend(strvecT,'Location','best','AutoUpdate','off'); legend('boxoff');
grid on
for k=1:jT_plot
plot(thH(k,:),xH,'-k');
end

plot(thE1,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(thE2,xE2,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(thE3,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(thE4,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)

%% presure:

figure; hold all;
for k=1:jT_plot
P(k)=plot(pT(k,:)*100,xb,'MarkerSize',4); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best','AutoUpdate','off'); legend('boxoff'); 
grid on
for k=1:jT_plot
plot(pH(k,:),xH,'-k');
end

plot(pE1,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(pE2,xE2,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(pE3,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(pE4,xE,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
