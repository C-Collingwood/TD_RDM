function [MODEL] = Simulate(model,trials,parameter,map,non_decision,options)
%% Simulate free and time-controlled trials
%%% This function creates a 'MODEL' object with a simulated set of choice
%%% and RT data, for the four models provided.

%%% Input:
% model: Habit1_Race, Habit2_Race, RL2_Race, RL_Race (STRING)
% trials: total number of trials
% parameter: a vector of parameters in the order 
%           [alpha(s), beta(s), t_2, theta]
% map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
%      (e.g., [1,1,1;2,2,3;3,3,2])
% non_decision: t1 parameter, (= 0.1s, default)

%%% Options:
% stimulus: A vector of the stimuli seen in each trial, default = random.
% time_cont: logical vector of flags for forced (=1) or free trials (=0, default)
% forced_rt: a vector (1,trial number) with set reaction time (seconds) at
%            the elements corresponding to forced trials, , default = uniform[0,2].
% remap_trial: index of trials where map switches from A to B, (=1, default)
% new_cond: Trial number where a new condition starts, (return to map A,
%           all variables reset to initial), (= [], default).
% error: Value for output on error trials (= 0, default).
% reward: Value for output on correct trials (=1, default).

%%% Returns:
% MODEL: object containing data and model information.


%% Input checks

arguments
    model string {mustBeMember(model,["Habit1_Race","Habit2_Race","RL2_Race","RL_Race"])}
    trials (1,1) double
    parameter (1,:) double
    map (:,3) double
    non_decision (1,1) double = 0.1
    options.stimulus (1,:) double = nan(1,trials)
    options.time_cont (1,:) double = zeros(1,trials)
    options.force_rt (1,:) double = nan(1,trials)
    options.remap_trial double = 1
    options.new_cond (1,:) double = []
    options.error (1,1) double = 0
    options.reward (1,1) double = 1
end
if length(options.stimulus)~=trials
    error('Stimulus vector must be the length of trial number')
end
if length(options.time_cont)~=trials
    error('Time-Controlled flag vector must be the length of trial number')
end
if length(options.force_rt)~=trials
    error(['Force-RT vector must be the length of trial number, ' ...
        'elements corresponding to free trials should be set to 0 or nan'])
end
switch model
    case {"Habit1_Race"}
        if length(parameter)~=6
            error('Expecting parameter vector in order [aq, ah, bq, bh, t2, theta]')
        end
    case {"Habit2_Race"}
        if length(parameter)~=7
            error('Expecting parameter vector in order [aq, ah, bq, bh1, bh2, t2, theta]')
        end
    case {"RL2_Race"}
        if length(parameter)~=6
            error('Expecting parameter vector in order [aq2, aq1, bq2, bq1, t2, theta]')
        end
    case {"RL_Race"}
        if length(parameter)~=3
            error('Expecting parameter vector in order [aq, bq, theta]')
        end
end

if ~all(isnan(options.force_rt)) & any(options.force_rt(options.time_cont==1)<=0|isnan(options.force_rt(options.time_cont==1)<=0))
    error('Reaction times missing for some forced trials')
end


STIM = options.stimulus; FORCE = options.time_cont; FORCE_RT = options.force_rt;
REMAP = options.remap_trial; NEW_COND = options.new_cond;
ERR = options.error; REW = options.reward;


%% CODE:

%%% Create Model Object

MODEL = Model_obj(model,parameter,non_decision,trials,map);
data_vec = nan(1,trials);
MODEL = initialise(MODEL,data_vec,STIM,data_vec,data_vec,FORCE,[],REMAP);

%%% Set up simulation:
n_stim = length(unique(map(:,1)));
choice_opt = unique(map(:,2:3));
n_forced = sum(FORCE);

fun = str2func(model);
if all(isnan(STIM))
    STIM = randi(n_stim,1,trials);
end
if ~all(isnan(FORCE_RT))
    FORCE_RT(FORCE==1)=randi(2000,1,n_forced)./1000;
end


%%% Run trials
for t = 1:trials
    if any(NEW_COND==t)
        f = fieldnames(MODEL.values);
        for i = 1:length(f)
            if contains(f{i},"Q")
                MODEL.values.(f{i})(:,:,t:end)=0.5;
            elseif contains(f{i},"d") % don't need to reset delta
            else
                MODEL.values.(f{i})(:,:,t:end)=0;
            end
        end
    end

    %%% Observation function
    if FORCE(t) ==0
        [MODEL] = fun(MODEL,"obs",t,{[]});
    else
        [MODEL] = fun(MODEL,"obs",t,{FORCE_RT(t)});
    end
    
    if MODEL.data.Remap_Trial(t)==0
        corr_choice = map(map(:,1)==STIM(t),2);
    else
        corr_choice = map(map(:,1)==STIM(t),3);
    end


    if MODEL.data.Choice(t)==corr_choice
        MODEL.data.Outcome(t)=REW;
    else
        MODEL.data.Outcome(t)=ERR;
    end


    %%% Learning Model
    if ~isnan(MODEL.data.Choice(t))
        [MODEL] = fun(MODEL,"learn",t,choice_opt);
    end

end
end
