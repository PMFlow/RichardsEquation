function [tht]=theta_GM(theta_res,theta_sat,p,c,alpha,nGM,ng,aGM,bGM)

theta_neg = theta_res+(theta_sat-theta_res)*(1./(1+real((-alpha*(1./(1-bGM*log(c/aGM+1)).*p)).^nGM))).^ng;
theta_pos = theta_sat;
neg = p < 0; 
tht = theta_neg.*neg + theta_pos.*~neg;

