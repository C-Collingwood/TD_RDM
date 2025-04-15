function [output1,output2] = RL2_race(obj,step,t,additional)
%% This function holds all the observation model and learning rules for the RL2 RACE model
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
aq2 = par.aq2;
aq1 = par.aq1;
bq2 = par.bq2;
bq1 = par.bq1;
t_q1 = par.tq1;
t_q2= par.tq2;
theta = par.theta;
s=par.s;


Q1 = obj.values.Q1;
Q2 = obj.values.Q2;
num_stim = size(Q1,2);
num_choice = size(Q1,1);


STIM = obj.data.Stimulus(t);
CHOSE = obj.data.Choice(t);
OUTCOME = obj.data.Outcome(t);


if step == "learn"
    Q1(:,:,t+1)=Q1(:,:,t);
    Q2(:,:,t+1)=Q2(:,:,t);
    
    dq1 =OUTCOME-Q1(CHOSE,STIM,t);
    Q1(CHOSE,STIM,t+1)=Q1(CHOSE,STIM,t)+aq1.*dq1;

    dq2 =OUTCOME-Q2(CHOSE,STIM,t);
    Q2(CHOSE,STIM,t+1)=Q2(CHOSE,STIM,t)+aq2.*dq2;
    
    obj.values.Q1 = Q1;
    obj.values.Q2 = Q2;
    obj.values.dq1(CHOSE,STIM,t)=dq1;
    obj.values.dq2(CHOSE,STIM,t)=dq2;
    output1 = obj;
    output2 = [];
end


if step == "mu"
    output1 = bq1.*Q1(:,STIM,t);
    output2 = bq1.*Q1(:,STIM,t)+bq2.*Q2(:,STIM,t);
end


if step == "obs" % observation model
    dt = 0.001;
    if ~isempty(additional)
        RT = additional{1};
        force = 1;
    else
        force = 0;
    end
    
    %%% Accumulation 
     t0 = round(t_q1/dt);
     noise =  randn(num_stim,2/dt);
     noise(:,1:t0)=0;
     mu1 = bq1.*Q1(:,STIM,T);
     mu2 =  bq1.*Q1(:,STIM,T)+bq2.*Q2(:,STIM,T);
    t1 = round((t_q1-dt)/dt:t_q2/dt);
    t2 = round((t_q2+dt)/dt:2/dt);
     V = zeros(num_stim,2/dt);
    if t1(1)==0
        t1(1)=[];
    end
    V(:,t1)=repmat(mu1,1,length(t1));
    V(:,t2)=repmat(mu2,1,length(t2));
    X=cumsum(dt.*V+sqrt(dt).*noise,2);


    if force ==0 % free choice
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
    elseif force == 1
        if RT<=t_q1
    	    obj.data.chose(T) = randi(num_choice);
        else
            X_force = X(:,round((T)/dt));
            [~,obj.data.chose(T)] = max((X_force == max(X_force)));
        end
        obj.data.Reaction_Time(T)=RT;
    end

    output1=obj; output2 = [];
end

end