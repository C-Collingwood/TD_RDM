function [CDF_out] = twodriftc(rt,t2,v1,v2,s,b)
%%% Calculates CDF of second accumulation period

iv1 = b./v1; 
is = (b./s).^2;
ftq = cdf_ig(t2,iv1,is);
t_all = rt-t2;
fun_cdf = @(x)(cdf_ig(t_all,(b-x)./v2,((b-x)./s).^2).*bvp(x,t2,v1,s,b));
x = -b:0.01:b;
min_lim = find(~isnan(fun_cdf(x)),1,'first');
if ~isempty(min_lim)
    CDF_out = integral(fun_cdf,x(min_lim),b)+ftq;
else
    CDF_out = ftq;
end
if any(isnan(CDF_out))
    CDF_out(isnan(CDF_out))=1;
end
end