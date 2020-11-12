close all; 
irand=1; % 0=determinist; 1=random;
dx=0.05;
D0 = 1e-03;
if irand==0
    load loam_det
    xx=x; p1=p; c1=c; tht1=tht; q1=q; D1=D;
    load clay_det
    p2=p; c2=c; tht2=tht; q2=q; D2=D;
else
    load loam
    xx=x; p1=p; c1=c; tht1=tht; q1=q; D1=D;
    load clay
    p2=p; c2=c; tht2=tht; q2=q; D2=D;
end
figure;
plot(p1,xx,p2,xx,'LineWidth',1.5); 
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t)$','Interpreter','latex'); 
legend('loam','clay','Location','best'); legend('boxoff');
figure;
plot(c1,xx,c2,xx,'LineWidth',1.5); 
ylabel('$z$','Interpreter','latex');
xlabel('$c(z,t)$','Interpreter','latex'); 
legend('loam','clay','Location','best'); legend('boxoff');
figure;
plot(tht1,xx,tht2,xx,'LineWidth',1.5); 
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 
legend('loam','clay','Location','south'); legend('boxoff');
figure
plot(q1,xx,q2,xx,'LineWidth',1.5)
ylabel('$z$','Interpreter','latex'); 
xlabel('$q(z,t)$','Interpreter','latex'); 
legend('loam','clay','Location','southeast'); legend('boxoff');
%% Peclet numbers:
for i=1:2
    if i==1
        disp('Pe_loam');
        D=D1; q=q1;
        disp('p_Peclet:')
        Dk=mean(mean(D))
        pV=mean(mean(abs(diff(D))))
        pPe=pV/Dk % Pe=(diff(K)/dx)*dx/Dk=diff(K)/Dk
        maxpV=max(max(abs(diff(D))))
        maxpPe=maxpV/Dk
        disp('c_Peclet:')
        V=mean(mean(abs(q)))
        Pe=V*dx/D0
        maxV=max(max(abs(q)))
        maxPe=maxV*dx/D0        
    else
        disp('Pe_clay');
        D=D2; q=q2;
        disp('p_Peclet:')
        Dk=mean(mean(D))
        pV=mean(mean(abs(diff(D))))
        pPe=pV/Dk % Pe=(diff(K)/dx)*dx/Dk=diff(K)/Dk        
        maxpV=max(max(abs(diff(D))))
        maxpPe=maxpV/Dk        
        disp('c_Peclet:')
        V=mean(mean(abs(q)))
        Pe=V*dx/D0
        maxV=max(max(abs(q)))
        maxPe=maxV*dx/D0
    end    
end