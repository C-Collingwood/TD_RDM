function [MODEL] = Simulate(model,trials,parameter,map,options)
%% Simulate free and forced trials

%%% This function creates a 'MODEL' object with a simulated set of choice
%%% and RT data.
%%% Input:
% model: Habit_race (default), RL2_race, RL_race (STRING)
% trials: number of trials
% parameter: a vector of parameters in the order 
%           [alpha, beta, t_start, t_switch, theta]
% map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
%      (e.g., [1,1,1;2,2,3;3,3,2])

%%% Options:
% stimulus: A vector of the stimuli seen in each trial, default = random.
% force: logical vector of flags for forced (=1) or free trials (=0, default)
% forced_rt: a vector (1,trial number) with set reaction time (seconds) at
%            the elements corresponding to forced trials, , default = uniform[0,2].
% remap_trial:index of trials where map switches from A to B, (=1, default)
% new_cond: Trial number where a new condition starts, (return to map A,
%           all variables reset to initial), (= [], default).
% error: Value for output on error trials (= 0, default).
% reward: Value for output on correct trials (=1, default).

%%% Returns:
% MODEL: object containing data and model information.


%% Input checks

arguments
    model string {mustBeMember(model,["Habit_race","RL2_race","RL_race"])}
    trials (1,1) double
    parameter (1,:) double
    map (:,3) double
    options.stimulus (1,:) double = nan(1,trials)
    options.force (1,:) double = zeros(1,trials)
    options.force_rt (1,:) double = nan(1,trials)
    options.remap_trial double = 1
    options.new_cond (1,:) double = []
    options.error (1,1) double = 0
    options.reward (1,1) double = 1
end
if length(options.stimulus)~=trials
    error('Stimulus vector must be the length of trial number')
end
if length(options.force)~=trials
    error('Force flag vector must be the length of trial number')
end
if length(options.force_rt)~=trials
    error(['Force RT vector must be the length of trial number, ' ...
        'elements corresponding to free trials should be set to 0 or nan'])
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

if ~all(isnan(options.force_rt)) & any(options.force_rt(options.force==1)<=0|isnan(options.force_rt(options.force==1)<=0))
    error('Reaction times missing for some forced trials')
end


STIM = options.stimulus; FORCE = options.force; FORCE_RT = options.force_rt;
REMAP = options.remap; NEW_COND = options.new_cond;
ERR = options.error; REW = options.reward;


%% CODE:

%%% Create Model Object

MODEL = Model(model,parameter,trials,map);
data_vec = nan(1,trials);
MODEL = initialise(MODEL,data_vec,STIM,data_vec,data_vec,FORCE,[],REMAP);

%%% Set up simulation:
n_stim = sum(unique(map(:,1)));
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

    if MODEL.data.Choice(t)==STIM(t)
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
