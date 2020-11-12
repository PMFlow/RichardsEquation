function sol = K_r(x,Nmod,KMean,varK,C1,C2,phi)    

x = x(:)';
coeff = sqrt(varK*2/Nmod) ;       % amplitude

s = sum(cos( (C1.*repmat(x,size(C1,1),1) + C2)*(2*pi) + phi)) ;
ak = coeff*s ;

sol = KMean * exp(-varK/2) * exp(ak) ;
%sol = sol(:);

end

