%%% calculate cost for free trials
function [LL] = free_cost(rt,t1,mu1,mu2,chose,non_c,theta,s)
mu1(mu1<=0)=0.00001;


iv = theta./mu1;
is = (theta./s).^2;
nstim = length(non_c)+1;
pdf_bin = 0.01;

if rt <t1
    PDF_out = pdf_ig(rt,iv(chose),is);
    CDF_out = cdf_ig(rt,iv,is);
    if any(isnan(CDF_out))
        CDF_out(isnan(CDF_out))=0;
    end
else
    for c = 1:nstim
        CDF_out(c)=twodriftc(rt,t1,mu1(c),mu2(c),s,theta);
    end
    cdf_up = twodriftc(rt+pdf_bin,t1,mu1(chose),mu2(chose),s,theta);
    if rt - pdf_bin <=0 %%% t_start and t_switch very close
        cdf_down = nan; 
    elseif rt-pdf_bin <=t1 %%% rt on boundary between rates
        cdf_down = cdf_ig(rt-pdf_bin,iv(chose),is);
    else
        cdf_down = twodriftc(rt-pdf_bin,t1,mu1(chose),mu2(chose),s,theta);
    end
    PDF_out =  (cdf_up-cdf_down)./(2*pdf_bin);
end

if PDF_out <=0
    PDF_out = 0.00000000000001;
end
CDF_out(CDF_out>0.9999)=0.9999;
LL = log(PDF_out*prod(1-CDF_out(non_c)));

end