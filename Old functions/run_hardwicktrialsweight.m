function [NLL,LL,Q,H] = run_hardwicktrialsweight(model,init,par_all,data,nlr,free,w)
%% This function runs a simulation of the Hardwick experiment based on data provided and returns NLL.
%%% Input:
%%% model: string name for the model used to determine choices
%%% par: the parameters to be fit for the model (aq,bq,ah,bh)
%%% t_par: theta, t_r and t_h (Not fitted)
%%% data: the subject data

%data=data(data.Exp==2,:); %% delete minimal training data
%%% UF: the uniform distribution f(t) (1/(maxRT-minRT)) par(9)=p0. par(10)=uf
%% set up parameters:

fun = str2func(model);
% Set parameters to correct places [aq,bq,ah,bh,theta]
par=par_all;
par(free)=init;
par(5) = round(par(5),4);
par(6) = round(par(6),4);

% par_ft = par; %% Make forced choice tq a different parameter
% par_ft(6) = par(9);
% par_ft(9)=[]; par(9)=[];

%%%Transform alpha
% % if alpha ==1
% % par(1) = 1./(1+exp(-par(1)));
% % par(3) = 1./(1+exp(-par(3)));
% % end

%%%%%par = par.*[0.1,10,0.01,1,0.1,0.1,1,1];
% Set internal parameters
chose = data.Response;
RT = data.RT;
Q=0.5*ones(4,4,height(data)); %Action, Stimulus, Time
H = zeros(4,4,height(data));
r = data.Correct;
LL=nan(1,length(r));

f = data.First;
bl = data.Block;
stim = data.Stimulus;
if ismember('Exp', fieldnames(data))
    e2_start = find(data.Exp==2,1,'first');
else
    e2_start=1;
end

trial_num = 0;
if sum(bl==3 & f == 1)>200
    NLL = 1*10^9;
  %  NLL = nan;
else
    %%% Learning and Observation model
    for t = 1:length(bl)
        if t == e2_start
            %%% This trial is new experiment, so reset Q and H
            Q(:,:,t)=0.5;
            H(:,:,t)=0;
        end
        Q(:,:,t+1)=Q(:,:,t);
        H(:,:,t+1)=H(:,:,t);
        variable{1} = Q(:,stim(t),t);
        variable{2} = H(:,stim(t),t);
        if ~isnan(chose(t)) % & f(t)==1 %%% We learn from all trials
            [Q(chose(t),stim(t),t+1),H(:,stim(t),t+1)] = fun(par,variable,"learn",[r(t),chose(t)]);
        end

        if f(t)~=1 | isnan(chose(t)) | isnan(RT(t)) %%If it isn't first response made, or a response wasnt made
            LL(t)=nan; % We don't calculate likelihood
        elseif bl(t)==4 %%%%& exist('fr_only','var')==0
            if RT(t)>=par(5)
                par_ft = par;
                LL(t) = hardwickFLL(Q(:,stim(t),t),H(:,stim(t),t),par_ft,(RT(t)),chose(t),nlr);
            else
                LL(t) = nan;
            end
        elseif  t<=50 | (bl(t)==3) | (t>=e2_start & t<=e2_start+49)%%%bl(t)~=4  %%
             LL(t) = (hardwickLL(Q(:,stim(t),t),H(:,stim(t),t),par,(RT(t)),chose(t),nlr));
        else
            LL(t)=nan;
        end

    end

    NLL = -(1-w)*sum(LL(bl~=4),'omitnan')-w*sum(LL(bl==4),'omitnan');%-sum(LL,'omitnan');%%
end
end

