% plot comparison GRW-FD
close all

load pGRWe3.mat
pGRWe3=p;
load pGRWe6.mat
pGRWe6=p;
load pGRWe10.mat
pGRWe10=p;
load pGRWe18.mat
pGRWe18=p;
load pGRWe24.mat
pGRWe24=p;

load pFD.mat
pFD=p;

figure; hold all;
plot(pGRWe3,x,'rd','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(pGRWe6,x,'b+','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(pGRWe10,x,'m*','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(pGRWe18,x,'gs','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(pGRWe24,x,'yo','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(pFD,x,'k.','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
ylabel('$z$','Interpreter','latex'); xlabel('$\psi(z)$','Interpreter','latex');
legend('N=1e3','N=1e6','N=1e10','N=1e18','N=1e24','FD','Interpreter','latex'); box on;

figure; hold all;
plot(abs(pFD-pGRWe3)./abs(pFD),x,'rd','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(abs(pFD-pGRWe6)./abs(pFD),x,'b+','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(abs(pFD-pGRWe10)./abs(pFD),x,'m*','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(abs(pFD-pGRWe18)./abs(pFD),x,'gs','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
plot(abs(pFD-pGRWe24)./abs(pFD),x,'k.','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
ylabel('$z$','Interpreter','latex'); xlabel('$(|\psi^{GRW}(z)-\psi^{FD}(z)|/|\psi^{FD}(z)|$','Interpreter','latex')
legend('N=1e3','N=1e6','N=1e10','N=1e18','N=1e24','Interpreter','latex','location','best'); xlim([1e-20 1e5])
set(gca,'xscale','log'); box on;
