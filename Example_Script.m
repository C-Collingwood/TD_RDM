%%% HARDWICK EXAMPLE:

load('Hardwick_Data_2019')

%%% Fitting to subject 1

s_data = data.s1;

%%% Clean data
s_data(s_data.Correct==0,:)=[];
s_data.Correct(s_data.Correct<0)=0;
s_data.RT(s_data.RT<0)=0.001;%[];
s_data(s_data.RT>5,:)=[];


CHOICE = s_data.Response;
STIMULUS = s_data.Stimulus;
RT = s_data.RT;
OUTCOME = s_data.Correct;

MAP = [1,1,1;
       2,2,2;
       3,3,4;
       4,4,3];


%%% Using optional arguments
MODEL = "Habit_race"; % (this is default value)
FORCE = zeros(1,height(s_data));
    FORCE(s_data.Block==4)=1; % Forced choice trials were included in Block 4
REMAP_TRIAL = [find(s_data.Block==3&s_data.Exp==1,1,"first"),find(s_data.Block==3&s_data.Exp==2,1,"first")];
NON_DECISION = t_h(1);% 0; % ms, (Use inbuilt heuristic, this is default value), 
INCLUDE_FIT = ones(1,height(s_data));
     INCLUDE_FIT(s_data.First~=1) =0 ; % only include first responses in fitting
     %INCLUDE_FIT(s_data.Block==) = 0; % ignore Criterion A in fitting
     INCLUDE_FIT((s_data.Block==1|s_data.Block==2) & s_data.Trial>50) = 0; % Only first 50 trials of initial learning included in fitting
WEIGHT = 0.75;
NEW_COND = find(s_data.Exp==2,1,"first");

%%% Run fitting
 
[PARAM, NLL, LL, MODEL_obj] = fitTwoDriftRM(CHOICE,STIMULUS,RT,OUTCOME,MAP,FORCE,...
                            model=MODEL,...
                            remap_trial = REMAP_TRIAL, ...
                            non_decision=NON_DECISION,...
                            include_fit=INCLUDE_FIT, ...
                            weight=WEIGHT,...
                            new_cond = NEW_COND);