function [p_t] = bvp(x,t,v,s,b)
var = s.^2;
vt = v.*t;
denom = 2.*var.*t;
e1 = exp(-((x-vt).^2)./denom);
e2 = exp(1).^((2.*v.*b)./var);
e3 = exp(-((x-2.*b-vt).^2)./denom);
c = 1./sqrt(pi.*denom);

p_t = c.*(e1-e2.*e3);
end