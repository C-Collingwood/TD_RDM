function [output1,output2] = Habit_race(obj,step,T,additional)
%% This function holds all the observation model and learning rules for the HABIT RACE model
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
ah = par.ah;
bq = par.bq;
bh = par.bh;
t_h = par.th;
t_q = par.tq;
theta = par.theta;
s=par.s;



Q = obj.values.Q;
H = obj.values.H;
num_stim = size(Q,2);
num_choice = size(Q,1);


STIM = obj.data.Stimulus(T);
CHOSE = obj.data.Choice(T);
OUTCOME = obj.data.Outcome(T);

if step == "learn"
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
    output1 = bh.*H(:,STIM,T);
    output2 = bh.*H(:,STIM,T)+bq.*Q(:,STIM,T);
end



if step == "obs" % observation model
    dt = 0.001;
   
    if ~isempty(additional{1})
        RT = additional{1};
        force = 1;
    else
        force = 0;
    end
    
    %%% Accumulation 
     t0 = round(t_h/dt);
     noise =  randn(num_stim,2/dt);
     noise(:,1:t0)=0;
     mu1 = bh.*H(:,STIM,T);
     mu2 =  bh.*H(:,STIM,T)+bq.*Q(:,STIM,T);
    t1 = round((t_h-dt)/dt:t_q/dt);
    t2 = round((t_q+dt)/dt:2/dt);
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
            obj.data.Reaction_Time(T) = nan;
            obj.data.Choice(T)=nan;
        else
          [obj.data.Reaction_Time(T),obj.data.Choice(T)] = min(RT);
        end
    elseif force == 1
        if RT<=t_h
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