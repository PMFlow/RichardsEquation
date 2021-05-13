function [tht]=theta(p,c)

theta_neg = 1./(3.4333-p-1/10.*c);
theta_pos = 0.3;
neg = p < 0; 
tht = theta_neg.*neg + theta_pos.*~neg;

