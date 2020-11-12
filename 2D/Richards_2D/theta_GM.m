function [tht]=theta_GM(theta_res,theta_sat,c,alpha,nGM,ng)
           
theta_neg = theta_res+(theta_sat-theta_res)*(1./(1+(-alpha*c).^nGM)).^ng;
theta_pos = theta_sat;
neg = c < 0; 
tht = theta_neg.*neg + theta_pos.*~neg;

