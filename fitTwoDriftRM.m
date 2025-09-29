function [PARAM, NLL, LL, MODEL] = fitTwoDriftRM(choice,stimulus,rt,outcome,map,time_cont,non_decision,options)
%% Fit time-controlled and free-RT trials
%%% This function takes choice and reaction time data, and returns the 
%%% best-fitting parameters according to a two-drift race model

%%% Input:
% choice: A vector of choices made across the trials
% stimulus: A vector of the stimuli seen in each trial
% rt: A vector of reaction times in each trial  (in seconds)
% outcome: A vector of outcomes (0/1) of each trial
% map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
%      (e.g., [1,1,1;2,2,3;3,3,2])
% time_cont: logical vector of flags for time_controlled (=1, default) or free trials (=0)
% non_decision: non-decision time, default = 0.1s
%%% Options:
% model: Habit_race (default), RL2_race, RL_race (STRING)
% remap_trial: index of trials where map switches from A to B, (=1, default)
% include_fit: logical vector for trials to be included in fitting
%              procedure. (=ones(1,number trials), default).
% weight: Weighting towards time_controlled trials, (0-1, =0.5,default)
% new_cond: Trial number where a new condition starts, (return to map A,
%           all variables reset to initial), = [], default.
% mult_start: How many start-points to run, = 0, default (non-parallel).

%%% Returns:
% PARAM: The best fitting parameters;  [aq, bq, t1, theta]   [aq2, aq1, bq2, bq1, t1, t2, theta]   [aq, ah, bq, bh, t1, t2, theta]  [aq, ah, bq, bh, t1, t2, theta] 
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
    time_cont (1,:) 
    non_decision (1,1) double = 0.1;
    options.model string {mustBeMember(options.model,["Habit1_Race","Habit2_Race","RL2_Race","RL_Race"])} ="Habit1_Race"
    options.remap_trial (1,:) double = 1
    options.include_fit (1,:) double = ones(size(choice))
    options.weight (1,1) double =0.5;
    options.new_cond(1,:) double = [];
    options.mult_start (1,1) double =0;
end

if ~isequal(size(choice),size(stimulus)) | ~isequal(size(choice),size(rt)) | ~isequal(size(choice),size(outcome)) 
    error('Data vectors must be the same size')
end
if ~isequal(size(choice),size(options.include_fit))
    error('include_fit vector must be the same size as data vectors')
end
if ~isequal(size(choice),size(time_cont))
    error('logical vector for time_cont must be the same size as data vectors')
end
if ~all(options.include_fit==1 | options.include_fit==0)
    error('include_fit vector must only contain 1 (time_cont) and 0 (free) values')
end

model = options.model; remap_trial = options.remap_trial; fit = options.include_fit; new_cond = options.new_cond;
weight = options.weight;

%% CODE:

%%% Set up parameters
t1 = non_decision;
if t1>=0.4
    t2 = t1+0.01;
else
    t2 = 0.4;
end

t2_min = t1+0.001;
switch model
    case "RL2_Race"
        PARAM =  [0.3, 0.0015, 10, 1, t2, 3];
        PARAM_boundary = [1, 1, 100, 100,  0.7,   100; ...
                          0, 0,  0,   0,  t2_min, 0.1];
    case "Habit1_Race"  %[aq, ah, bq, bh, t2, theta] 
        PARAM =  [0.3, 0.0015, 10, 1, t2, 3];
        PARAM_boundary = [1, 0.005, 100, 100, 0.7,   100; ...
                          0,   0,    0,   0, t2_min, 0.1];
    case {"Habit2_Race"} %[aq, ah, bq, bh1, bh2, t1, t2, theta] 
        PARAM =  [0.3, 0.0015, 10, 1, 1, t2,  3, 1];
        PARAM_boundary = [1, 0.005, 100, 100, 100, 0.7,   100; ...
                          0,   0,    0,   0,   0, t2_min, 0.1];
    case "RL_Race"
        PARAM = [0.3, 10, 3];
        PARAM_boundary = [1, 100, 100;...
                          0,  0, 0.1];
end


%%% Create Model Object
if weight ~=1
    free_idx = 1:length(PARAM_boundary);
else
    free_idx = 1:(length(PARAM_boundary)-1);
end

n_trials = length(choice);
MODEL = Model_obj(model,PARAM,t1,n_trials,map);
if weight==1
    MODEL.free_par(MODEL.free_par=="theta")=[];
end
MODEL = initialise(MODEL,choice,stimulus,rt,outcome,time_cont,PARAM_boundary(:,free_idx),remap_trial);
MODEL.weight = weight;

%%% Adding remapped trial to data structure if included
if ~isempty(new_cond) & ~isempty(remap_trial)
    remap_vec = ones(1,length(fit));
    remap_vec(1:remap_trial(1))=0;
    for i = 1:length(new_cond)
        r = find((remap_trial-new_cond(i))>=0,1,"first");
        if~isempty(r)
            remap_vec(new_cond(i):remap_trial(r))=0;
        end
    end
elseif ~isempty(remap_trial)
    remap_vec = zeros(1,length(fit)); % Set all to Map A
    for i = 1:2:length(remap_trial)
        if i ==length(remap_trial)
            remap_vec(remap_trial(i):end)=1; % Switch to Map B
        else
            remap_vec(remap_trial(i):remap_trial(i+1)-1)=1; % Switch to Map B
        end
    end
else
    remap_vec = ones(1,length(fit));
end
MODEL.data.Remap_Trial=remap_vec';

%%% Don't fit trials with RT < non-decision time;
fit = remove_fast(MODEL,fit);


%%% Set upper and lower boundaries
LB = table2array(MODEL.par_bound("LB",:));
UB = table2array(MODEL.par_bound("UB",:));
A = zeros(1,length(free_idx)); b = 0;
if contains(model,"RL2_race")
    A(1)=-1; A(2)=1;
end

    
%%% Run fitting
PARAM = table2array(MODEL.par);
init = PARAM(free_idx);

fun = @(x)runfitting(x,MODEL.free_par,fit,MODEL,new_cond);

if options.mult_start==0
    [P_fit,NLL]=fmincon(fun,init,A,b,[],[],LB,UB);
else
    problem = createOptimProblem('fmincon','objective',fun,'x0',init,'lb',LB,'ub',UB,'Aineq',A,'bineq',b);

ms = MultiStart("StartPointsToRun","bounds",Display="final",MaxTime=20000,UseParallel=true);
[P_fit,NLL] = run(ms,problem,options.mult_start);
end

PARAM(free_idx)=P_fit;

%%% Update MODEL obj
MODEL = update_param(MODEL,P_fit);


fun = @(x)runfitting(x,MODEL.free_par,fit,MODEL,new_cond);
[~,LL,MODEL]=fun(P_fit);

MODEL.fit.LL = LL;
MODEL.fit.NLL = NLL;

end