function [PARAM, NLL, LL, MODEL] = fitTwoDriftRM(choice,stimulus,rt,outcome,map,force,options)
%% Fit forced trials
%%% This function takes choice and reaction time data for FORCED CHOICE
%%% trials, and returns the best fitting parameters according to a
%%% two-drift race model.
%%% Input:
% choice: A vector of choices made across the trials
% stimulus: A vector of the stimuli seen in each trial
% rt: A vector of reaction times in each trial  (in seconds)
% outcome: A vector of outcomes (0/1) of each trial
% map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
%      (e.g., [1,1,1;2,2,3;3,3,2])
% force: logical vector of flags for forced (=1, default) or free trials (=0)
%%% Options:
% model: Habit_race (default), RL2_race, RL_race (STRING)
% remap_trial:index of trials where map switches from A to B, (=1, default)
% non_decision: non-decision time, default uses heuristic described by
%               approx_th.m. If non_decision <=0, heuristic used.
% include_fit: logical vector for trials to be included in fitting
%              procedure. (=ones(1,number trials), default).
% weight: Weighting towards forced trials, (0-1, =0.5,default)
% new_cond: Trial number where a new condition starts, (return to map A,
%           all variables reset to initial), = [], default.

%%% Returns:
% PARAM: The best fitting parameters, [aq, ah, bq, bh, th, tq, theta] 
% NLL: Negative loglikelihood for best fitting parameters
% LL: a vector of loglikelihoods for each trial.
% MODEL: object containing data about model iteration

%% Input checks
arguments
    choice (1,:) double
    stimulus (1,:) double
    rt (1,:) double
    outcome (1,:) double
    map (:,3) double
    force (1,:) 
    options.model string {mustBeMember(options.model,["Habit_race","RL2_race","RL_race"])} ="Habit_race"
   
    options.remap_trial (1,:) double = 1
    options.non_decision (1,1) double = NaN;
    options.include_fit (1,:) double = ones(size(choice))
    options.weight (1,1) double =0.5;
    options.new_cond(1,:) double = [];

end

if ~isequal(size(choice),size(stimulus)) | ~isequal(size(choice),size(rt)) | ~isequal(size(choice),size(outcome)) 
    error('Data vectors must be the same size')
end
if ~isequal(size(choice),size(options.include_fit))
    error('include_fit vector must be the same size as data vectors')
end
if ~isequal(size(choice),size(force))
    error('logical vector for forced must be the same size as data vectors')
end
if ~all(options.include_fit==1 | options.include_fit==0)
    error('include_fit vector must only contain 1 (forced) and 0 (free) values')
end

model = options.model; remap_trial = options.remap_trial; fit = options.include_fit; new_cond = options.new_cond;


%% CODE:

%%% Set up parameters
fun = str2func(model);
switch model
    case {"Habit_race","RL2_race"}
        PARAM =  [0.3,0.0015,10,1,0.15,0.4,3]; %[0.999999432387102		0.00499999894096654	4.30097989383138 2.40442082198110	0.237167319044520	0.418849992901298	1.14124285051639];%
        PARAM_boundary = [1,0.005,100,100,0.5,0.6,100; 0, 0, 0, 0, 0.1 ,0.1, 0.1];
    case "RL_race"
        PARAM = [0.3,10,0.15,3]; 
        PARAM_boundary = [1,100,0.5,100; 0, 0, 0.1 , 0.1];
end

LL = nan(1,length(choice));

%%% Create Model Object
if contains(model,"Habit") || contains(model,"RL2")
    free_idx = [1 2 3 4 6 7];
else
    free_idx = [1 2 4];
end
n_trials = length(choice);
MODEL = Model(model,PARAM,n_trials,map);
MODEL = initialise(MODEL,choice,stimulus,rt,outcome,force,PARAM_boundary(:,free_idx),remap_trial);
MODEL.weight = options.weight;

if ~isempty(new_cond) & ~isempty(remap_trial)
    remap_vec = ones(1,length(fit));
    remap_vec(1:remap_trial(1))=0;
    for i = 1:length(new_cond)
        r = find((remap_trial-new_cond(i))>0,1,"first");
        if~isempty(r)
            remap_vec(new_cond(i):remap_trial(r))=0;
        end
    end
end
MODEL.data.Remap_trial=remap_vec';

%%% Approximate the time parameters
MODEL = update_time(MODEL,options.non_decision);

%%% Don't fit trials with RT < non-decision time;
fit = remove_fast(MODEL,fit);

%%% Run fitting
PARAM = table2array(MODEL.par);
init = PARAM(free_idx);

fun = @(x)runfitting(x,MODEL.free_par,fit,MODEL,new_cond);

LB = table2array(MODEL.par_bound("LB",:));
UB = table2array(MODEL.par_bound("UB",:));


[P_fit,NLL]=fmincon(fun,init,[],[],[],[],LB,UB);

PARAM(free_idx)=P_fit;

%%% Update MODEL obj
MODEL = update_param(MODEL,P_fit);


fun = @(init)runfitting(init,MODEL.free_par,fit,MODEL,new_cond);
[~,LL,MODEL]=fun(P_fit);

MODEL.fit.LL = LL;
MODEL.fit.NLL = NLL;






end