%% errors of z(\theta) relative to Warrick's solution [Warrick et al., 1985]
clear all
close all
I= 101;
a=0; b=1;
dx=(b-a)/(I-1);
xthT=a:dx:b;

%
load thT;
load WX_n1_5
% load thT2;
% load WX_n2
xq=thtW;

nt=4; dx=1e-2;
err=zeros(nt,3);
figure; hold all;
for k=1:nt
    y=xthT(find(thT(k,:)>0.175));
    x=thT(k,find(thT(k,:)>0.175)); 
    P(k)=plot(x,y);
end
NameArray = {'Marker'}; ValueArray = {'o','+','x','*'}';
set(P,NameArray,ValueArray);
legend(strvecTW,'Location','best','AutoUpdate','off'); legend('boxoff');

for k=1:nt
    y=xthT(find(thT(k,:)>0.175));
    x=thT(k,find(thT(k,:)>0.175)); 
    p = pchip(x,y);
    pp = ppval(p,xq);
    xthtW=1-xW(k,:);
    err(k,:)=(pp-xthtW)./thtW;
    plot(xq,pp,'-s')   
    plot(thtW,1-xW(k,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4); % '-*')
    xlabel('$\theta(z,t)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); %legend('off')
    fprintf('k = %2d err = %0.2e  %0.2e %0.2e\n',k,err(k,:)) ;
end

% k =  1 err = 5.31e-02  5.31e-02 5.69e-02
% k =  2 err = -2.46e-03  2.175e-02 5.10e-02
% k =  3 err = -5.70e-02  -1.41e-02 4.55e-02
% k =  4 err = -9.69e-02  -3.88e-02 4.91e-02

