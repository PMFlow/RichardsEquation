function plot_T(jT_plot,xb,strvecT)

load ('pT','pT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(pT(k,:)*100,xb); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best'); legend('boxoff');
grid on

load ('thT','thT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(thT(k,:),xb); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best'); legend('boxoff');
grid on

load ('qT','qT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(qT(k,:),xb); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$q(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best'); legend('boxoff');
grid on
