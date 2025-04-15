function [output1,output2] = RL_race(obj,step,t,additional)
%% This function holds all the observation model and learning rules for the RL RACE model
%%% Input:
% obj: Model Object containing data and variables
% step: which process is occurring ("learn","observe","mu")
% t: trial number
%%% Output:
% Depends on step.
%   "learn"; 
%   "observe";
%   "mu";

par = obj.par;

aq = par.aq;
bq = par.bq;
t_q = par.tq;
theta = par.theta;
s=par.s;


Q = obj.values.Q;
num_stim = size(Q,2);
num_choice = size(Q,1);

STIM = obj.data.Stimulus(t);
CHOSE = obj.data.Choice(t);
OUTCOME = obj.data.Outcome(t);

if step == "learn"
    Q(:,:,t+1)=Q(:,:,t);
    dq =OUTCOME-Q(CHOSE,STIM,t);
    Q(CHOSE,STIM,t+1)= Q(CHOSE,STIM,t)+aq.*dq;

    obj.values.Q = Q;
    obj.values.dq(CHOSE,STIM,t)=dq;
    output1 = obj;

    output2 = [];

end

if step == "mu"
    output1 = bq.*Q(:,STIM,t);
    output2 = output1;
end




if step == "obs" % observation model
    dt = 0.001;
    
     if ~isempty(additional)
        RT = additional{1};
        force = 1;
    else
        force = 0;
    end
    
    
    
    
    
    
    
    % %         Adapted to free v forced choice

    t0 = round(t_q/dt);
    noise = randn(num_stim,2/dt);
    noise(:,1:t0)=0;
    t = additional{1};
    mu = bq.*Q(:,STIM,t);
    t1 = round((t_q-dt)/dt:2/dt);
    V = zeros(num_stim,2/dt);
    if t1(1)==0
        t1(1)=[];
    end
    V(:,t1)=repmat(mu,1,length(t1));
    X=cumsum(dt.*V+sqrt(dt).*noise,2);

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
            obj.data.Reaction_Time = nan;
            obj.data.chose(T)=nan;
        else
            [obj.data.Reaction_Time(T),obj.data.chose(T)] = min(RT);
        end

    else %% t is predetermined
	    if t<=t_q
	        obj.data.chose(T) = randi(num_choice);
        else
            X_force = X(:,round(t/dt));
            [~,obj.data.chose(T)] = max((X_force == max(X_force)));
        end
        obj.data.Reaction_Time(T)=RT;
    end
end





end