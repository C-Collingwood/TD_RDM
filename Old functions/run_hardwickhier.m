function [FR_NLL,fr_par,FD_NLL,fd_par,trial_num] = run_hardwickhier(model,fr_free,starting_par,nlr,mini,maxi,data)
%HIERARCHICAL Running hardwick trials with hierarchical independent fitting
%   Runs fitting separately on free response versus fixed duration trials,
%   then uses hierarchical fitting to find best fitting parameters
fd_free = fr_free(fr_free~=8);
fr_par = starting_par;
fd_par = starting_par;
options = optimoptions('fmincon', 'Display','notify');

mfunFR = @(x)run_hardwicktrials(model,x,fr_par,data,nlr,fr_free,1); %% free only
mfunFD = @(x)run_hardwicktrials(model,x,fd_par,data,nlr,fd_free,0); %% forced only
init_fr = fr_par(fr_free);
[~,LL_1] = mfunFR(init_fr);
trial_num(1) = sum(~isnan(LL_1));

init_fd = fd_par(fd_free);
[~,LL_2] = mfunFD(init_fd);
%trial_num(2) = sum(~isnan(LL_2));


%%% Run initial fitting
[p1,FR_NLL(1),~,~,~,~,~]=fmincon(mfunFR,init_fr,[],[],[],[],mini(fr_free),maxi(fr_free),[],options);
fr_par(1,:)=fr_par(1,:);
fr_par(1,fr_free)=p1;

[p2,FD_NLL(1),~,~,~,~,~]=fmincon(mfunFD,init_fd,[],[],[],[],mini(fd_free),maxi(fd_free),[],options);
fd_par(1,:)=fd_par(1,:);
fd_par(1,fd_free)=p2;


%%% Run 10 rounds:

mfunFR = @(x,fr_free,av_par,s_par)run_hardwicktrialsH(model,x,fr_par,data,nlr,fr_free,1,av_par,s_par); %% free only
mfunFD = @(x,fd_free,av_par,s_par)run_hardwicktrialsH(model,x,fd_par,data,nlr,fd_free,0,av_par,s_par); %% forced only


%%% Run initial fitting
for i = 1:4
av_par = mean([p1(1:end-1);p2]);
s_par = std([p1(1:end-1);p2]);
idx = find(s_par<10^-4);
av_par(idx)=[];
s_par(idx)=[];
free = fr_free;
free(idx)=[];


if isempty(s_par)
    fr_par(i+1,:)=fr_par(i,:);
    FR_NLL(i)=nan;
    fd_par(i+1,:)=fd_par(i,:);
    FD_NLL(i)=nan;
else
%%% Run free
    funfr = @(x)mfunFR(x,free,av_par,s_par);
    funfr(fr_par(i,free));
[p1,FR_NLL(i),~,~,~,~,~]=fmincon(funfr,fr_par(i,free),[],[],[],[],mini(free),maxi(free),[],options);
fr_par(i+1,:)=fr_par(i,:);
fr_par(i+1,free)=p1;

%%% Run fixed
free = fd_free;
free(idx)=[];
funfd = @(x)mfunFD(x,free,av_par,s_par);

[p2,FD_NLL(i),~,~,~,~,~]=fmincon(funfd,fd_par(i,free),[],[],[],[],mini(free),maxi(free),[],options);
fd_par(i+1,:)=fd_par(i,:);
fd_par(i+1,free)=p2;

end
end



end