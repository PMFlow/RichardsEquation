%% plots of L-scheme convergence for (p,c) - coupled flow transport 2D
% (random Ksat case)

close all; clear all;
itest=0; % 0=loam; 1=clay; 
S= 50000; iS=1:S;
kt_plot=3;
strvect=['t=1'; 't=2'; 't=3'];
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
if itest == 1
    disp('Beit Netofa clay')
    load convf_clay
    load convt_clay
    dk=150;
else
    disp('slit loam')
    load convf_loam
    load convt_loam
    dk=20;
end
figure; hold all;
for k=1:kt_plot
P(k)=plot(iS(1:dk:end),convf(k,1:dk:end));
end
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff'); 

figure; hold all;
for k=1:kt_plot
P(k)=plot(iS(1:dk:end),convt(k,1:dk:end));
end
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|c^s - c^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');
