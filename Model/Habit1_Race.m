function [output1,output2] = Habit1_Race(obj,step,T,additional)
%% This function holds all the observation model and learning rules for the HABIT RACE_1b model
%%% Input:
% obj: Model Object containing data and variables
% step: which process is occurring ("learn","observe","mu")
% T: trial number
% additional: If "learn", contains available choices, 
%             If "obs" and time-controlled trial, contains RT. 

%%% Output:
% Depends on step.
%   "learn"; 
%   "obs";
%   "mu";

%%% Extract Parameters and Variables
par = obj.par;
aq = par.aq;
ah = par.ah;
bq = par.bq;
bh = par.bh;
t_1 = par.t1;
t_2 = par.t2;
theta = par.theta;
s=par.s;


Q = obj.values.Q;
H = obj.values.H;
num_stim = size(Q,2);
num_choice = size(Q,1);

%%% Extract trial data
STIM = obj.data.Stimulus(T);
CHOSE = obj.data.Choice(T);
OUTCOME = obj.data.Outcome(T);

if step == "learn"
    %%% Calculate prediction error and update variables

    Q(:,:,T+1)=Q(:,:,T);
    H(:,:,T+1)=H(:,:,T);
    
    dq =OUTCOME-Q(CHOSE,STIM,T);
    Q(CHOSE,STIM,T+1)=Q(CHOSE,STIM,T)+aq.*dq;
    
    choices = additional;


    for c = 1:length(choices)
        dh(c) = double(choices(c)==CHOSE)-H(c,STIM,T); % 1 if action was made, 0 otherwise
        H(c,STIM,T+1)= H(c,STIM,T)+ah.*dh(c);
    end
    
    obj.values.Q = Q;
    obj.values.H = H;
    obj.values.dq(CHOSE,STIM,T)=dq;
    obj.values.dh(:,STIM,T)=dh;
    output1 = obj;

    output2 = [];
end


if step == "mu"
    %%% Extract drift-rate
    output1 = bh.*H(:,STIM,T);
    output2 = bh.*H(:,STIM,T)+bq.*Q(:,STIM,T);
end



if step == "obs" 
    dt = 0.001;
   
    if ~isempty(additional{1})
        RT = additional{1};
        force = 1;
    else
        force = 0;
    end
    
    %%% Accumulation 
     t0 = round(t_1/dt);
     noise =  randn(num_stim,2/dt);
     noise(:,1:t0)=0;
    [mu1,mu2] = Habit1_Race(obj,"mu",T);
    th = round((t_1-dt)/dt:t_2/dt);
    tp = round((t_2+dt)/dt:2/dt);
     V = zeros(num_stim,2/dt);
    if th(1)==0
        th(1)=[];
    end
    V(:,th)=repmat(mu1,1,length(th));
    V(:,tp)=repmat(mu2,1,length(tp));
    X=cumsum(dt.*V+sqrt(dt).*noise,2);

    %%% Free or Time-Controlled, calculate choice and RT
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
            obj.data.Reaction_Time(T) = nan;
            obj.data.Choice(T)=nan;
        else
          [obj.data.Reaction_Time(T),obj.data.Choice(T)] = min(RT);
        end
    elseif force == 1
        if RT<=t_1
    	    obj.data.Choice(T) = randi(num_choice);
        else
            X_force = X(:,round((RT)/dt));
            [~,obj.data.Choice(T)] = max((X_force == max(X_force)));
        end
        obj.data.Reaction_Time(T)=RT;
    end

    output1=obj; output2 = [];

end

end