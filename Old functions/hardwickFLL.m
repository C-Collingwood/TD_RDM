function [LL] = hardwickFLL(Q,H,par,rt,chose,nlr)
%%% Calculates the likelihood of a given RT and choice for an interogation
%%% process
%%% Flag determines whether a constant or changing drift rate is used.
%%% Flag = 1 : Constant drift of Q+H
%%% Flag = 2: Drift change from H to Q+H at t_q.


%%%%%%%%%%%
bq= par(2);
bh = par(4);
s = par(7);

t_h = par(5); %minimum RT minus the non-decision limit
t_q = par(6); % Fixed parameter
t_2 = t_q-t_h;
t_1 = rt-t_h;
dt = 0.001;

nc = 1:length(Q); %idx of non chosen options
if ~isnan(chose) & chose~=0
    nc(chose)=[];
end
if rt<=t_h+dt
    mu = zeros(1,length(Q));
    P=nan;
else
    %v1=repmat(bh*H,1,max_rt);
    %v2=repmat(bh*H+bq*Q,1,max_rt);

    if nlr ==1
        mu = bq.*Q.*(t_1);
    elseif rt<t_q
        mu = bh.*H.*(t_1);
    else%if rt>=t_q
        mu = (bh.*H).*t_1 + (bq.*Q).*(t_1-t_2);% (bh.*H).*t_1 + (bq.*Q).*(t_1-t_2);%bh.*H.*t_2 + (bq.*Q + bh.*H).*(rt-t_q);
        
    end
    if isnan(chose)
        mu(end+1)=0;
        chose = 5;
    end
    %%
    s_t = s.*sqrt(t_1);
    mu_nc = mu(nc);
    if chose~=5
        fun = @(x)(normpdf(x,mu(chose),s_t).*(normcdf(x,mu(nc(1)),s_t)).*(normcdf(x,mu(nc(2)),s_t)).*(normcdf(x,mu(nc(3)),s_t)));
    else
        fun = @(x)(normpdf(x,mu(chose),s_t).*(normcdf(x,mu_nc(1),s_t)).*(normcdf(x,mu_nc(2),s_t)).*(normcdf(x,mu_nc(3),s_t)).*(normcdf(x,mu_nc(4),s_t)));
    end

    P = integral(fun,-Inf,Inf);
end
if P<=0
    P=0.0000000001;
end

if length(par) == 11 %%% if we are including contaminants
    p0=par(9);
    uf = (par(11)+par(10))./2;
    LL = log((1-p0).*P + p0.*uf);
else 
    LL = log(P);
end






end