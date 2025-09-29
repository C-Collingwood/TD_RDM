function [NLL,LL,MODEL] = runfitting(init,init_names,fit,MODEL,new_cond)
%% Runs fitting for free and forced trials using two-drift race model.
%%% Input
% init: Initial values of free parameters
% init_names: the parameter names ordering of the init vector
% fit: Vector flagging trials for inclusion in fitting (default = all)
% MODEL: model object containing parameter and data information.
%%% Output
% NLL = Negative loglikelihood as calculated by two-drift race model cost
%       function.
% LL = Loglikelihood for each trial, as calculated by two-drift race model
%       cost function.


%%% Extract relevant data
fun = str2func(MODEL.type);
data = MODEL.data;
LL=nan(1,height(data));
forced_flag = MODEL.data.Time_Controlled;
w = MODEL.weight;

%%% Assign new parameters

for i = 1:length(init_names)
    MODEL.par.(init_names{i})=init(i);
end

%%% Extract time parameters

switch MODEL.type
    case {"Habit1_Race", "RL2_Race", "Habit2_Race"}
        t_start = MODEL.par.t1;
        t_switch = MODEL.par.t2;
    case "RL_Race"
        t_start = MODEL.par.t1;
        t_switch = t_start+0.15; % Irrelevant as before and after are the same.
end

CHOICE_OPT = unique(MODEL.map(:,2:3));
s = MODEL.par.s;

for t = 1:height(data)
    %%% Updating values
    if any(new_cond==t)
        f = fieldnames(MODEL.values);
        for i = 1:length(f)
            if contains(f{i},"Q")
                MODEL.values.(f{i})(:,:,t:end)=0.5;
            elseif contains(f{i},"d") % don't need to reset delta
            else
                MODEL.values.(f{i})(:,:,t:end)=0;
            end
        end

    end

    [MODEL] = fun(MODEL,"learn",t,CHOICE_OPT);

    %%% calculating likelihood
    if fit(t)~=0 
        %%% Calculate drift rates
        [mu1,mu2] = fun(MODEL,"mu",t);
        
        %%% Extract RT
        RT = MODEL.data.Reaction_Time(t);
        CHOSE = MODEL.data.Choice(t);
        NON_C = CHOICE_OPT;
        NON_C(NON_C==CHOSE)=[];

        rt = RT-t_start;
        t12 = t_switch-t_start;
        t2e = RT-t_switch;

        if forced_flag(t)==1
                if rt<t12
                    mu = mu1.*rt;
                else
                    mu = mu1.*t12+mu2.*(t2e);
                end
                LL(t)=forced_cost(mu(CHOSE),mu(NON_C),s.*sqrt(rt));
        elseif forced_flag(t)==0
                theta = MODEL.par.theta;
                LL(t)= free_cost(rt, t12,mu1,mu2,CHOSE,NON_C,theta,s);
                
        end
    end
end

NLL = -(1-w).*sum(LL(forced_flag==0),'omitnan') - w*sum(LL(forced_flag==1),'omitnan');

end