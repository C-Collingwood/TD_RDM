function [output1,output2] = Qlearn2(par,variable,step,additional)
%% This function holds all the observation model and learning rules for Miller
%%% Input:
%%% par: Parameter values, in the order [a_q, beta, kappa, a_h]
%%% variable = all Q, H values at that timestep
%%% step = which subequation do you need
%%% options: "obs" calculates utility, "learn" updates Q values,
%%% "priorlog" calculates negloglike with additional prior
%%% additional provides extra inputs required for that step;
%%% if "obs", additional has [last_choice]
%%% if "learn", additional has [reward, choice taken]
%%% if "priorlog", additional has [population parameters]
%%% cont_val = max RT and min RT values for contaminant distribution
aq1 = par(1);
bq1 = par(2);
aq2 = par(3);
bq2 = par(4);
t_h = par(5);
t_q = par(6);
s=par(7);
theta = par(8);

if length(par)==11
p0=par(9);
uf=(par(11)-par(10))./2;
end

dt = 0.001;
Q1 = variable{1};
Q2 = variable{2};
num_stim = length(Q1);


if step == "obs" % observation model
    % %         Adapted to free v forced choice
    t0 = round(t_h/dt);
    noise =  randn(num_stim,2/dt);
    noise(:,1:t0)=0;
    t = additional{1};
    mu1 = bq2.*Q2;
    mu2 = bq2.*Q2+bq1.*Q1;
    t1 = round((t_h-dt)/dt:t_q/dt);
    t2 = round((t_q+dt)/dt:2/dt);
    v = zeros(num_stim,2/dt);
    if t1(1)==0
        t1(1)=[];
    end

    v(:,t1)=repmat(mu1,1,length(t1));
    v(:,t2)=repmat(mu2,1,length(t2));
    X=cumsum(dt.*v+sqrt(dt).*noise,2);


    if isempty(t) % free choice
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
            output2 = nan;
            output1 = nan;
        else
            
            if length(par)==11 & rand<p0 %% including contaminants
                a = par(10);b=par(11);
               % RT = rand(1,4).*(b-a)+a; % add time before choice
                %[output2,output1] = min(RT);
                output2=rand.*(b-a)+a;
                output1=randi(4);
               %%% output2 = output2+rand.*(b-a)+a; % time to add;
            else 
                [output2,output1] = min(RT);
            end
        end
    else %% t is predetermined
        if t<=t_h
            output1 = randi(4);
        else
            X_force = X(:,round((t)/dt));
            [~,output1] = max((X_force == max(X_force)));
        end
        output2=t;
    end
end

if step == "learn"
    rew = additional(1);
    choice = additional(2);
    dq1 =rew-Q1(choice);
    output1 = Q1(choice)+aq1.*dq1;
    dq2 =rew-Q2(choice);
    output2 = zeros(1,num_stim);
    output2(choice) = Q2(choice)+aq2.*dq2;
end




end