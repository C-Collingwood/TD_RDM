function [output1,output2] = Simple_HabitsQ(par,variable,step,additional)
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
ah = par(3);
bh = par(4);
t_h = par(5);
t_q = par(6);
s=par(7);
theta = par(8);
dt = 0.001;
Q = variable{1};
H = variable{2};
num_stim = length(Q);

%    noise = [zeros(num_stim,t0),randn(num_stim,2/dt-t0)];
if step == "obs" % observation model
    % %         Adapted to free v forced choice
     t0 = round(t_h/dt);
     noise =  randn(num_stim,2/dt);
     noise(:,1:t0)=0;
    t = additional{1};
    mu1 = bh.*H;
    mu2 =  bq.*Q; %  bh.*H+bq.*Q;% 
%    t05 = 1:t0;
    t1 = round((t_h-dt)/dt:t_q/dt);
    t2 = round((t_q+dt)/dt:2/dt);
    v = zeros(num_stim,2/dt);
 %   v(:,t05)=zeros(num_stim,length(t05));
    v(:,t1)=repmat(mu1,1,length(t1));
    v(:,t2)=repmat(mu2,1,length(t2));
    X=cumsum(dt.*v+sqrt(dt).*noise,2);
    % 
    % 
    % % figure
    % % hold on
    % % plot(X')
    % % xline(t./dt,':')
    % % xline(t_q./dt,'-')
    % % legend({'1','2','3','4','RT','t_q'})

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
        [output2,output1] = min(RT);
        end
    else %% t is predetermined
        if t<=t_h+5*dt
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
    dq =rew-Q(choice);
    output1 = Q(choice)+aq.*dq;
    for stim = 1:num_stim
        dh(stim) = double(stim==choice)-H(stim); % 1 if stim was chosen, 0 otherwise
        output2(stim) = H(stim)+ah.*dh(stim);
    end
end




end