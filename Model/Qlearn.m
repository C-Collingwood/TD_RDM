function [output1,output2] = Qlearn(par,variable,step,additional)
       %% This function holds all the observation model and learning rules for Qlearn
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
    aq = par(1);
    bq = par(2);
    t_h = par(5);
    s = par(7);
    theta = par(8);
    Q = variable{1};
    num_stim = length(Q);
    dt = 0.001;
    t0 = (t_h/dt);
    num_stim = length(Q);
  if length(par)==11
p0=par(9);
uf=(par(11)-par(10))./2;
end


    if step == "obs" % observation model
        % %         Adapted to free v forced choice
        
        t0 = round(t_h/dt);
        noise = randn(num_stim,2/dt);
        noise(:,1:t0)=0;
        t = additional{1};
        mu = bq.*Q;
        t1 = round((t_h-dt)/dt:2/dt);
        v = zeros(num_stim,2/dt);
        if t1(1)==0
            t1(1)=[];
        end
        v(:,t1)=repmat(mu,1,length(t1));
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
                X_force = X(:,round(t/dt));
                [~,output1] = max((X_force == max(X_force)));
        	end
            output2=t;
        end
    end
    if step == "learn"
        rew = additional(1);
        choice = additional(2);
        dq =rew-Q(choice);
        output1 = Q(choice)+aq.*dq;
        output2 = zeros(num_stim,1);
    end
end