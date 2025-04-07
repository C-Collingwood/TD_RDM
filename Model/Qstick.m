function [output1,output2] = Qstick(par,variable,step,additional)
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

    aq = par(1);
    bq = par(2);
    bh = par(4);
    t_h = par(5);
    t_q = par(6); 
    s = par(7);
    theta = par(8);
    Q = variable{1};
    L = variable{2};
    num_stim = length(Q);
    dt = 0.001;
    noise = randn(num_stim,2/dt);
   
        noise(:,1:t_h/dt)=0;
    if step == "obs" % observation model
        % %         Adapted to free v forced choice
        t = additional{1};
        mu1 = bh.*L;
        mu2 = bh.*L+bq.*Q;
        noise(:,1:t_h/dt)=0;


        t1 = round((t_h+dt)/dt:t_q/dt);
        t2 = round((t_q+dt)/dt:2/dt);
        v = zeros(num_stim,2/dt);
        
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
            [output2,output1] = min(RT);
        else %% t is predetermined 
           X_force = X(:,round(t/dt));
           [~,output1] = max((X_force == max(X_force)));
           output2=t;
        end

    end
    if step == "learn"
        rew = additional(1);
        choice = additional(2);
        dq =rew-Q(choice);
        output1 = Q(choice)+aq.*dq;
        output2 = zeros(num_stim,1); % Replacing H with Last Choice
        if ~isnan(choice)
        output2(choice)=1;
        end
    end
%     if step == "priorlog"
%         aq_pop = additional(:,1);
%         beta_pop = additional(:,2);
%         kappa_pop = additional(:,3);
%         ah_pop = additional(:,4);
%         output1 = -log(normpdf(aq,aq_pop(1),aq_pop(2)))-log(normpdf(bq,beta_pop(1),beta_pop(2)))-log(normpdf(bh,kappa_pop(1),kappa_pop(2)))-log(normpdf(ah,ah_pop(1),ah_pop(2)));
%     end

    
    
end