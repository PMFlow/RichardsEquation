%% Warrick's solution for n=1.5 according to [Warrick et al., 1985, Table 3]

clear all;
close all;
Ksat = 0.5184; 
theta_res=0.1;
theta_sat=0.45;
alpha_W=1; 

W=[0.4 0.6 0.8]; % = (theta-theta_res)/(theta_sat-theta_res) 
%                  corresponding to W* in [Warrick et al., 1985, Table 3]

nt=4;
X=zeros(nt,3); it=zeros(1,nt); IT=zeros(1,nt);
DT=0.5/24; % =1/2 hours
T=0;
for k=1:nt
    Tk=k*DT;
    strTW=['T=',num2str(Tk*24)];
    strvecTW(k,1:length(strTW))=strTW;  
    it(k)=Tk;
    T=Tk*alpha_W*Ksat/(theta_sat-theta_res);
    IT(k)=T;
    X(k,1)=0.828*T^(0.5)+0.278*T+0.159*T^(1.5);
    X(k,2)=0.791*T^(0.5)+0.302*T+0.165*T^(1.5);
    X(k,3)=0.691*T^(0.5)+0.364*T+0.186*T^(1.5);
end

figure; hold all;
for k=1:nt
P(k)=plot(W,-X(k,:));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$Z$','Interpreter','latex');
xlabel('$\Theta(z,t)$','Interpreter','latex'); 
legend(strvecTW,'Location','best','AutoUpdate','off'); legend('boxoff');
title('Dimensionless variables [Warrick et al., 1985]');
grid on

thtW=theta_res+(theta_sat-theta_res)*W;
xW=X/alpha_W; 
figure; hold all;
for k=1:4
P(k)=plot(thtW,-xW(k,:));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 
legend(strvecTW,'Location','best','AutoUpdate','off'); legend('boxoff');
title('Actual (rescaled) variables');
grid on

% save('WX_n1_5','thtW','xW','strvecTW')

% it = [0.0208    0.0417    0.0625    0.0833] % days
% IT = [0.0309    0.0617    0.0926    0.1234] % dimensionless time
% xW =
% 
%     0.1549    0.1492    0.1336
%     0.2253    0.2177    0.1970
%     0.2821    0.2733    0.2492
%     0.3321    0.3223    0.2958
% W=[0.4 0.6 0.8];
% thtW =
% 
%     0.2568    0.2922    0.3276
