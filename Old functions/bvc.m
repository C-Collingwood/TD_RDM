function [p_t] = bvc(t,v,s,b)
var = s.^2;
denom = s.*sqrt(t);
e1 = (b-v.*t)./denom;
e2 = exp((2.*v.*b)./var);
e3 = (-b-v.*t)./denom;

p_t = normcdf(e1)-e2.*normcdf(e3);
end