function [output1,output2] = RL_Race(obj,step,T,additional)
%% This function holds all the observation model and learning rules for the RL RACE model
%%% Input:
% obj: Model Object containing data and variables
% step: which process is occurring ("learn","observe","mu")
% T: trial number
% additional: If "learn", contains available choices, 
%             If "observe" and time-controlled trial, contains RT. 

%%% Output:
% Depends on step.
%   "learn"; 
%   "obs";
%   "mu";

%%% Extract Parameters and Variables
par = obj.par;
aq = par.aq;
bq = par.bq;
t_1 = par.t1;
theta = par.theta;
s=par.s;


Q = obj.values.Q;
num_stim = size(Q,2);
num_choice = size(Q,1);


%%% Extract trial data
STIM = obj.data.Stimulus(T);
CHOSE = obj.data.Choice(T);
OUTCOME = obj.data.Outcome(T);

if step == "learn"
    %%% Calculate prediction error and update variables
    Q(:,:,T+1)=Q(:,:,T);
    dq =OUTCOME-Q(CHOSE,STIM,T);
    Q(CHOSE,STIM,T+1)= Q(CHOSE,STIM,T)+aq.*dq;

    obj.values.Q = Q;
    obj.values.dq(CHOSE,STIM,T)=dq;
    output1 = obj;

    output2 = [];

end

if step == "mu"
     %%% Extract drift-rate
    output1 = bq.*Q(:,STIM,T);
    output2 = output1;
end




if step == "obs" % observation model
    dt = 0.001;
    
     if ~isempty(additional{1})
        RT = additional{1};
        force = 1;
    else
        force = 0;
    end
    
    
    
    % %         Adapted to free v forced choice

    t0 = round(t_1/dt);
    noise = randn(num_stim,2/dt);
    noise(:,1:t0)=0;
    
    [mu,~] = RL_Race(obj,"mu",T);

    t1 = round((t_1-dt)/dt:2/dt);
    V = zeros(num_stim,2/dt);
    if t1(1)==0
        t1(1)=[];
    end
    V(:,t1)=repmat(mu,1,length(t1));
    X=cumsum(dt.*V+sqrt(dt).*noise,2);

     %%% Free or Time-Controlled, calculate choice and RT
    if force == 0
        for n=1:num_stim
            bound = find(X(n,:)>theta,1,'first');
            if ~isempty(bound)
                RT(n)=bound;
            else
                RT(n)=nan;
            end
        end
        RT=RT*dt;
        if all(isnan(RT))
            obj.data.Reaction_Time(T) = nan;
            obj.data.Choice(T)=nan;
        else
            [obj.data.Reaction_Time(T),obj.data.Choice(T)] = min(RT);
        end

    elseif force ==1 %% t is predetermined
        if RT<=t_1
            obj.data.Choice(T) = randi(num_choice);
        else
            X_force = X(:,round(RT/dt));
            [~,obj.data.Choice(T)] = max((X_force == max(X_force)));
        end
        obj.data.Reaction_Time(T)=RT;
    end
    output1=obj; output2 = [];
end





end