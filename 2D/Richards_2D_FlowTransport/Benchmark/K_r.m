function sol = K_r(x,y,Nmod,KMean,varK,C1,C2,phi)    

x = x(:)'; y=y(:)';
coeff = sqrt(varK*2/Nmod) ;       % amplitude


s = sum(cos( (C1.*repmat(x,size(C1,1),1) + C2.*repmat(y,size(C2,1),1))*(2*pi) + phi)) ;
ak = coeff*s ;

sol = KMean * exp(-varK/2) * exp(ak) ;
sol = sol(:);

end

