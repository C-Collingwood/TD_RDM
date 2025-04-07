function [CDF_out] = twodriftc2(t,tq,v1,v2,s,b,a)

fun1 = @(b) cdf_ig2(tq,b,a,v1,s);
ftq = integral(fun1,b-a,b);
t2 = t-tq;
fun_cdf = @(b,x)(cdf_ig2(t2,b,a,v2,s).*bvp(x,tq,v1,s,b));
x = 0:0.01:b;
min_lim = find(~isnan(fun_cdf(x)),1,'first');
if ~isempty(min_lim)
    CDF_out = integral2(fun_cdf,b-a,b,x(min_lim),b)+ftq;
else
    CDF_out = ftq;
end
if any(isnan(CDF_out))
    CDF_out(isnan(CDF_out))=1;
end
%CDF_out = conv(ftq,ptq);
end