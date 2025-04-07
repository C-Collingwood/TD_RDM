function [chose,r,bl,RT,stim_seen, first,trial_num,Q,H] = Simulate(model,par,trials,err)
%% This function creates data for the Hardwick experiment
%%% Input:
%%% model: string name for the model used to determine choices
%%% par: the parameters to be used for the model
%%% par_name: optional input to label the parameters, if included, will
%%% plot data
%%% trials: number of free and forced trials
%%% start: determines start values for Q and H

%% set up parameters:
fun = str2func(model);
%%% Create data structure
trials_free = trials(1);
trials_forced = trials(2);
stim = randi(4,1,trials_free);
stim_forced = randi(4,1,trials_forced);
time_forced = randi(1600,1,trials_forced)./1000; %milliseconds

rew = 1; 
% switch model
%     case 'Qlearn' % Provide [aq, bq,tr, s]
%         par(1:2) = par_in(1:2);
%        par(5) = par_in(3);
%         par(7) = par_in(4);
%     case 'Qstick' % Provide [aq bq bh tr th s]
%         par(1:2) = par_in(1:2);
%         par(4:7) = par_in(3:6);
%     case 'Simple_Habits' %Provide [aq bq ah bh tr th s]
%         par = par_in;
% end
% par(8)=theta(1);


switch model
    case 'Qlearn'
        par([3 4 6])=0;
    case 'Qstick'
        par(3)=0;
    case 'Simple_Habits'
    case 'Qlearn2'
end
par(5) = round(par(5),4);
par(6) = round(par(6),4);
%%%%par = par.*[0.1,10,0.01,1,0.1,0.1,1,1];
% Parameters
chose = zeros(1,trials_free);
RT = zeros(1,trials_free);
r = zeros(1,trials_free);
Q=0.5*ones(4,4,trials_free+trials_forced); %Action, Stimulus, Time
H = zeros(4,4,trials_free+trials_forced);

% Sq = ones(4,trials);
% Sh = ones(4,trials);
% P=zeros(:,trials);
% U = zeros(2,1);
bl = ones(1,trials_free+trials_forced);
stim_seen = zeros(1,trials_free+trials_forced);

%% Run free trials
i=1;
for t = 1:trials_free
    corr = 0;
    j = 1;
    while corr ==0 & j<10 %stop a trial going on forever
        if j ==1
            first(i)=1;
        else
            first(i)=0;
        end
        trial_num(i)=t;
        j = j+1;
        bl(i)=1;
        % New value
        Q(:,:,i+1)=Q(:,:,i);
        H(:,:,i+1)=H(:,:,i);
        %         Sh(:,i+1)=Sh(:,i);
        %         Sq(:,i+1)=Sq(:,i);

        % Observation Model %%% Need to change to match Sam's race model
        variable{1} = Q(:,stim(t),i);
        variable{2} = H(:,stim(t),i);
       %%% stim(t)
        %     variable{3} = Sq(:,i);
        %     variable{4} = Sh(:,i);
        [chose(i),RT(i)] = fun(par,variable,"obs",{[]});

        if chose(i)==stim(t)
            r(i) =rew;
            corr = 1;
        else
            r(i)=err;
        end
        % Learning Model
        if ~isnan(chose(i))
            [Q(chose(i),stim(t),i+1),H(:,stim(t),i+1)] = fun(par,variable,"learn",[r(i),chose(i)]);
        end

        stim_seen(i)=stim(t);
        i=i+1;
    end
end

%% Run Criterion 1
n_cor = zeros(1,4);
j = 1;
t_end = trial_num(end);
t = 1;
k=1;
while any(n_cor ~=[5,5,5,5]) && j<100
    if r(i-1)==1 | k ==5
        first(i)=1;
        if stim_seen(i-1)<4 
            stim_a=stim_seen(i-1)+1;
        else
            stim_a=1;
        end
        trial_num(i)=t_end+t;
        t = t+1;
        k = 1;
    else
      first(i)=0;
      trial_num(i) = trial_num(i-1);
       stim_a = stim_seen(i-1);
      k = k+1;
    end
    bl(i)=2;
   
    

    % New value
    Q(:,:,i+1)=Q(:,:,i);
    H(:,:,i+1)=H(:,:,i);
    % Observation Model %%% Need to change to match Sam's race model
    variable{1} = Q(:,stim_a,i);
    variable{2} = H(:,stim_a,i);
    [chose(i),RT(i)] = fun(par,variable,"obs",{[]});
    if chose(i)==stim_a
        r(i) =rew;
       
    else
        r(i)=err;
    end
    if r(i)==rew && n_cor(stim_a)<5
        n_cor(stim_a) = n_cor(stim_a)+1;
    elseif r(i)~=rew
        n_cor(stim_a)=0;
    end
    % Learning Model#
    if ~isnan(chose(i))
        [Q(chose(i),stim_a,i+1),H(:,stim_a,i+1)] = fun(par,variable,"learn",[r(i),chose(i)]);
    end
    stim_seen(i)=stim_a;
    i=i+1;
    j = j+1;
end


%% Run Criterion 2
n_cor = zeros(1,4);
j = 1;
t_end = trial_num(end);
t = 1;
k = 1;
while any(n_cor ~=[5,5,5,5]) && j<100
    if r(i-1)==1 | k ==5
        first(i)=1; %%% Add new stimulus
        if stim_seen(i-1)<4
            stim_b=stim_seen(i-1)+1;
        else
            stim_b=1;
        end

        trial_num(i)=t_end+t;
        t=t+1;
        k = 1;
    else
      first(i)=0;
      trial_num(i) = trial_num(i-1);
       stim_b = stim_seen(i-1);
       k = k+1;
    end
    bl(i)=3;
   
   
    % New value
    Q(:,:,i+1)=Q(:,:,i);
    H(:,:,i+1)=H(:,:,i);
    % Observation Model 
    variable{1} = Q(:,stim_b,i);
    variable{2} = H(:,stim_b,i);
    [chose(i),RT(i)] = fun(par,variable,"obs",{[]});
    
    
    if (chose(i)==4 & stim_b == 3) | (chose(i) == 3 & stim_b == 4) | (chose(i)==2 & stim_b == 2) | (chose(i)==1 & stim_b == 1)
        r(i) =rew;
      
    else
        r(i)=err;
    end

    % stim_b
    % chose(i)
    % r(i)

    if r(i)==rew && n_cor(stim_b)<5
            n_cor(stim_b) = n_cor(stim_b)+1;
        
    elseif r(i) ~=rew
        n_cor(stim_b)=0;
    end
    % Learning Model
    if ~isnan(chose(i))
        [Q(chose(i),stim_b,i+1),H(:,stim_b,i+1)] = fun(par,variable,"learn",[r(i),chose(i)]);
    end
    stim_seen(i)=stim_b;
    i=i+1;
    j = j+1;
end


%% Run Forced Response
t_end = trial_num(end);
for t = 1:trials_forced
    first(end+1)=1;
    trial_num(end+1)=t+t_end;
    bl(i) = 4;
    % New value
    Q(:,:,i+1)=Q(:,:,i);
    H(:,:,i+1)=H(:,:,i);
    % Observation Model %%% Need to change to match Sam's race model
    variable{1} = Q(:,stim_forced(t),i);
    variable{2} = H(:,stim_forced(t),i);
    [chose(i),RT(i)] = fun(par,variable,"obs",{time_forced(t)});


    if (chose(i)==4 && stim_forced(t) == 3) || (chose(i) == 3 && stim_forced(t) == 4) || (chose(i)==2 && stim_forced(t) == 2) || (chose(i)==1 && stim_forced(t) == 1)
        r(i) =rew;
    else
        r(i)=err;
    end
    % Learning Model
    if ~isnan(chose(i))
        [Q(chose(i),stim_forced(t),i+1),H(:,stim_forced(t),i+1)] = fun(par,variable,"learn",[r(i),chose(i)]);
    end
    stim_seen(i)=stim_forced(t);
    i = i+1;
end


end