%% HARDWICK EXAMPLE:
% This section of the script fits the Habit1-Race model to the behavioural data of S1, 
% from Hardwick et al. (2019).
load('Hardwick_Data_2019')
addpath('Cost_Functions')
addpath('Model')

%%% Fitting to Participant, s1
s_data = data.s1;

%%% Clean data
s_data(s_data.Correct==0,:)=[];
s_data.Correct(s_data.Correct<0)=0;
s_data.RT(s_data.RT<0)=0.001;%[];
s_data(s_data.RT>5,:)=[];
s_data(s_data.Response==0,:)=[];
s_data(isnan(s_data.Response),:)=[];

CHOICE = s_data.Response;
STIMULUS = s_data.Stimulus;
RT = s_data.RT;
OUTCOME = s_data.Correct;

MAP = [1,1,1;
       2,2,2;
       3,3,4;
       4,4,3];


%%% Using optional arguments
MODEL = "Habit1_Race"; % (this is default value)
FORCE = zeros(1,height(s_data));
    FORCE(s_data.Block==4)=1; % Time-controlled trials were included in Block 4
REMAP_TRIAL = [find(s_data.Block==3&s_data.Exp==1,1,"first"),find(s_data.Block==3&s_data.Exp==2,1,"first")]; % Purely to add to the data structure - unused during fitting
NON_DECISION = approx_t1(s_data,[]); % (Using our heuristic), 
INCLUDE_FIT = ones(1,height(s_data));
     INCLUDE_FIT(s_data.First~=1) =0 ; % only include first responses in fitting
     INCLUDE_FIT((s_data.Block==1|s_data.Block==2) & s_data.Trial>50) = 0; % Only first 50 trials of initial learning included in fitting
WEIGHT = 0.95;
NEW_COND = find(s_data.Exp==2,1,"first");

%%% Run fitting
[PARAM, NLL, LL, MODEL_s1] = fitTwoDriftRM(CHOICE,STIMULUS,RT,OUTCOME,MAP,FORCE,NON_DECISION,...
                            model=MODEL,...
                            remap_trial = REMAP_TRIAL, ...
                            include_fit=INCLUDE_FIT, ...
                            weight=WEIGHT,...
                            new_cond = NEW_COND, ...
                            mult_start = 10);


%% SIMULATE EXAMPLE:
% This section of the script produces a MODEL object containing simulated
% surrogate data.
addpath('Cost_Functions')
addpath('Model')

%%% Required Input
MODEL = "Habit1_Race"; % (this is default value)
TRIALS = 150; % Total number of trials, both free-RT and time-controlled
PARAMETER =  [0.3, 0.0015, 10, 1, 0.45, 3]; %[aq, ah, bq, bh, t2, theta] 
MAP = [1,1,1;
       2,2,2;
       3,3,4;
       4,4,3];
NON_DECISION = 0.15;

%%% Using optional arguments
STIMULUS = randi(4,1,TRIALS);
TIME_CONT = [zeros(1,100), ones(1,50)]; %100 free-RT trials followed by 50 time-controlled trials.
FORCED_RT = [nan(1,100),randi(2000,1,50)./1000]; % the externally enforced RT for the 50 time-controlled trials, uniformly-distributed between 0-2 seconds.
REMAP_TRIAL = 90; % Reversal of mapping occurs on the 90th trial.
NEW_COND = []; % No reset to new conditions are included.
ERROR = 0;
REWARD = 1;

%%% Create simulated data
[MODEL_sim] = Simulate(MODEL,TRIALS,PARAMETER,MAP,NON_DECISION,...
              stimulus=STIMULUS,...
              time_cont=TIME_CONT, ...
              force_rt = FORCED_RT,...
              remap_trial=REMAP_TRIAL, ...
              new_cond = NEW_COND, ...
              error = ERROR, ...
              reward = REWARD);

%%% Extract simulated data from model object
DATA = MODEL_sim.data;
%%

%%% Run fitting on simulated data (using default optional arguments and
%%% data contained in the MODEL_sim object)
[PARAM_sim, NLL_sim, LL_sim, MODEL_sim_fit] = fitTwoDriftRM(DATA.Choice,DATA.Stimulus,DATA.Reaction_Time,DATA.Outcome,MODEL_sim.map,DATA.Time_Controlled,0.15, ...
    "remap_trial",REMAP_TRIAL);


%% Approximate T_1
function [th]=approx_t1(data,mmn)
    if isempty(mmn)
        mmn=200;
    end
    if isstruct(data)
        d = struct2table(flipstruct(data));
    else
        d = data;
    end
    d(isnan(d.RT),:)=[];
    rt=data.RT(data.Block==1&data.Trial>3000&data.First==1&d.Correct==1);
    th=min(rt(rt>0.1));
end