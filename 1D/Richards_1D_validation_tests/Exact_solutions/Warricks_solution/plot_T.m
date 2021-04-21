function plot_T(jT_plot,xb,strvecT)

load ('pT','pT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(pT(k,:),xb); 
end
NameArray = {'Marker'}; ValueArray = {'o'};
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best'); legend('boxoff');
grid on

load ('WX_n1_5','thtW','xW')
load ('thT','thT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(thT(k,:),xb); 
end
NameArray = {'Marker'}; ValueArray = {'o'};
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best','AutoUpdate','off'); legend('boxoff');
grid on
for k=1:jT_plot
plot(thtW,-xW(k,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
end

load ('qT','qT')
figure; hold all;
for k=1:jT_plot
P(k)=plot(qT(k,:),xb); 
end
NameArray = {'Marker'}; ValueArray = {'o'};
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$q(z,t)$','Interpreter','latex'); 
legend(strvecT,'Location','best'); legend('boxoff');
grid on

