function plot_conv_Alt(kt_plot,S,strvect)

load convf
iS=1:S;
figure; hold all;
for k=1:kt_plot
P(k)=plot(iS,convf(k,:)); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');

load convt
figure; hold all;
for k=1:kt_plot
P(k)=plot(iS,convt(k,:));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|c^s - c^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');

