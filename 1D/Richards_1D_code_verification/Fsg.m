function [fs,fg]=Fsg(t,x,p)

    % % with Matlab:
%     syms p(t,x)
%     p = -t*x*(1-x)+x/4;
%     syms c(t,x)
%     c = t*x*(1-x)+1;
%     theta = 1./(3.4333-p-1/10.*c); 
%     % K= p/.p==1;
%     q = -(diff(p+x))
%     sol_q = @(t,x) - t*x - t*(x - 1) - 5/4;
%     F1 = diff(theta,t)-diff(diff(p+x,x),x); % q=-K*diff(p+x,x)=-diff(p+x,x)
    F1 = @(t,x) (9*x.*(x - 1))./(10*(x/4 + (9*t*x.*(x - 1))/10 - 33333/10000).^2) - 2*t;
%     F2 = -diff(diff(p+x,x),x); 
    F2 = @(t,x) - 2*t;
%     G1 = diff(theta*c,t)-diff(diff(p+x,x)*c,x)-D1*diff(diff(c,x),x);
    G1 = @(t,x) 2*t + (t*x + t*(x - 1)).*(t*x + t*(x - 1) + 5/4) + 2*t*(t*x.*(x - 1) - 1) ...
        + (x.*(x - 1))./(x/4 + (9*t*x.*(x - 1))/10 - 33333/10000) ...
        - (9*x.*(t*x.*(x - 1) - 1).*(x - 1))./(10*(x/4 + (9*t*x.*(x - 1))/10 - 33333/10000).^2);
%     G2 = diff(0.3*c,t)-diff(diff(p+x,x)*c,x)-D1*diff(diff(c,x),x);   %  theta=0.3, saturated
    G2 = @(t,x) 2*t + (t*x + t*(x - 1)).*(t*x + t*(x - 1) + 5/4) + 2*t*(t*x.*(x - 1) - 1) - (3*x.*(x - 1))/10;

    fs_neg=F1(t,x); fs_pos=F2(t,x);
    fg_neg=G1(t,x); fg_pos=G2(t,x);
    neg=p<0;
    fs=fs_neg.*neg+fs_pos.*~neg;
    fg=fg_neg.*neg+fg_pos.*~neg;
