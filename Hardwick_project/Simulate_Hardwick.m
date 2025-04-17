function [MODEL, data] = Simulate_Hardwick(model,trials,parameter,map,options)
%% Simulate Paradigm from Hardwick et al. 2019
%%% This function creates a 'MODEL' object with a simulated set of choice
%%% and RT data.
%%% Input:
% model: Habit_race (default), RL2_race, RL_race (STRING)
% trials: number of trials [free minimal forced_minimal free extended forced extended]
% parameter: a vector of parameters in the order 
%           [alpha, beta, t_start, t_switch, theta]
% map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
%      (e.g., [1,1,1;2,2,3;3,3,2])

%%% Options:
% error: Value for output on error trials (= 0, default).
% reward: Value for output on correct trials (=1, default).

%%% Returns:
% MODEL: object containing data and model information.


%% Input checks

arguments
    model string {mustBeMember(model,["Habit_race","RL2_race","RL_race"])}
    trials (1,4) double
    parameter (1,:) double
    map (:,3) double
    options.error (1,1) double = 0
    options.reward (1,1) double = 1
end
switch model
    case {"Habit_race"}
        if length(parameter)~=7
            error('Expecting parameter vector in order [aq, ah, bq, bh, th, tq, theta]')
        end
    case {"RL2_race"}
        if length(parameter)~=7
            error('Expecting parameter vector in order [aq2, aq1, bq2, bq1, tq1 tq2, theta]')
        end
    case {"RL_race"}
        if length(parameter)~=4
            error('Expecting parameter vector in order [aq, bq, tq1, theta]')
        end
end

ERR = options.error; REW = options.reward;


%% CODE:

%%% Create Model Object

MODEL = Model(model,parameter,sum(trials)*2,map);
data_vec = nan(1,sum(trials)*2);
MODEL = initialise(MODEL,data_vec,data_vec,data_vec,data_vec,data_vec,[],sum(trials));

%%% Add trial information
MODEL.data.TrialNum=data_vec';
MODEL.data.Block = data_vec';
MODEL.data.First = data_vec';
MODEL.data.Exp = data_vec';


%%% Set up simulation:
n_stim = length(unique(map(:,1)));
choice_opt = unique(map(:,2:3));

fun = str2func(model);



%%% Run trials
i = 1;
for exp = 1:2
    switch exp
        case 1
            stim_free = randi(n_stim,1,trials(1));
            stim_force = randi(n_stim,1,trials(2));
            exp_trials = trials(1:2);

            time_force = randi(1600,1,trials(2))./1000; %milliseconds

        case 2
            stim_free = randi(n_stim,1,trials(3));
            stim_force = randi(n_stim,1,trials(4));
            exp_trials = trials(3:4);
            time_force = randi(1600,1,trials(4))./1000; %milliseconds

            f = fieldnames(MODEL.values);
            for val = 1:length(f)
                if contains(f{val},"Q")
                    MODEL.values.(f{val})(:,:,i:end)=0.5;
                elseif contains(f{val},"d") % don't need to reset delta
                else
                    MODEL.values.(f{val})(:,:,i:end)=0;
                end
            end
    end


    %% Block 1

    for t = 1:exp_trials(1)
        corr = 0;
        j = 1;
        while corr == 0 & j<10 % Can't try more than 10 times / trial
            %%% Update trial data
            if j ==1
                MODEL.data.First(i)=1;
            else
                MODEL.data.First(i)=0;
            end
            MODEL.data.TrialNum(i)=t;
            MODEL.data.Block(i)=1;
            MODEL.data.Exp(i)=exp;
            %%% Update model
            MODEL.data.Stimulus(i) = stim_free(t);
            MODEL.data.Forced(i) = 0;
            MODEL.data.Remap(i) = 0;

            %%% Observation Function
            [MODEL] = fun(MODEL,"obs",i,{[]});
            CHOSE = MODEL.data.Choice(i);
            if CHOSE==map(map(:,1)==stim_free(t),2)
                MODEL.data.Outcome(i)=REW;
                corr = 1;
            else
                MODEL.data.Outcome(i)=ERR;
            end
            
            %%% Learning Model
            if ~isnan(CHOSE)
                [MODEL] = fun(MODEL,"learn",i,choice_opt);
            end
            i = i+1;
            j = j+1;
        end
    end
    
    %% Criterion A + B
    for bl = 2:3

        n_corr = zeros(1,4);
        j = 1;
        t_end = MODEL.data.TrialNum(i-1);
        t = 1;
        k = 1;

        while any(n_corr ~=[5,5,5,5]) & j<100

            %%% Update Trial Data
            MODEL.data.Block(i)=bl;
            MODEL.data.Exp(i)=exp;
            MODEL.data.Forced(i) = 0;
            MODEL.data.TrialNum(i)=t_end+t;

            if MODEL.data.Outcome(i-1)==REW | k ==5
                MODEL.data.First(i)=1;




                if MODEL.data.Stimulus(i-1)<4
                    MODEL.data.Stimulus(i)=MODEL.data.Stimulus(i-1)+1;
                else
                    MODEL.data.Stimulus(i)=1;
                end

                switch bl
                    case 2 %%% Criterion A
                        MODEL.data.Remap(i) = 0;
                    case 3 %%% Criterion B
                        MODEL.data.Remap(i) = 1;
                end
                t = t+1; k = 1;
            else

                MODEL.data.First(i)=0;

                MODEL.data.Stimulus(i)=MODEL.data.Stimulus(i-1);
            end






            %%% Observation Function
            [MODEL] = fun(MODEL,"obs",i,{[]});
            STIM = MODEL.data.Stimulus(i);
            CHOSE = MODEL.data.Choice(i);
            if CHOSE==map(map(:,1)==STIM,bl)
                MODEL.data.Outcome(i)=REW;
                if  n_corr(STIM) <5
                    n_corr(STIM)=n_corr(STIM)+1;
                end
            else
                MODEL.data.Outcome(i)=ERR;
                n_corr(STIM)=0;
            end
            %%% Learning Model
            if ~isnan(CHOSE)
                [MODEL] = fun(MODEL,"learn",i,choice_opt);
            end

            i = i+1; j = j+1;

        end
    end

    %% Forced Trials
    t_end = MODEL.data.TrialNum(i-1);
    for t = 1:exp_trials(2)
        %%% Update Trial Data
        MODEL.data.First(i)=1;
        MODEL.data.Block(i)=4;
        MODEL.data.Exp(i)=exp;
        MODEL.data.TrialNum(i)=t_end+t;
        MODEL.data.Forced(i)=1;
        MODEL.data.Remap(i) = 1;
        MODEL.data.Stimulus(i) = stim_force(t);

        %%% Observation Function
        [MODEL] = fun(MODEL,"obs",i,{time_force(t)});
        CHOSE = MODEL.data.Choice(i);
        if CHOSE==map(map(:,1)==stim_force(t),3)
            MODEL.data.Outcome(i)=REW;

        else
            MODEL.data.Outcome(i)=ERR;
        end
        %%% Learning Model
        if ~isnan(CHOSE)
            [MODEL] = fun(MODEL,"learn",i,choice_opt);
        end
        i = i+1;
    end
end

%%% Remove Excess Rows
idx = isnan(MODEL.data.Stimulus);
MODEL.data(isnan(MODEL.data.Stimulus),:)=[];
f = fieldnames(MODEL.values);
for i = 1:length(f)
    MODEL.values.(f{i})(:,:,idx)=[];
end

data = MODEL.data;
end

