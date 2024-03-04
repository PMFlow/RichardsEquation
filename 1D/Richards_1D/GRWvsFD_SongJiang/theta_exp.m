function [tht]=theta_exp(theta_res,theta_sat,alpha,p)

theta_neg = theta_res+(theta_sat-theta_res)*exp(alpha*p);
theta_pos = theta_sat;
neg = p < 0; 
tht = theta_neg.*neg + theta_pos.*~neg;

