%% Running Model Recovery for Hardwick experiment using TOOLBOX
clear
%% Set up
%%% Set up functions
load('Hardwick_Data_2019')
addpath("Cost_Functions\","Model\")

%%% Set up flags
model_list = {"RL_race","RL2_race","Habit_race"};
min_flag=1; 
err = 0;
tests = 30;

%%% Create true parameters
par_num = [3 6 6];

 max_i = [0.4, 0.002,15,2.5, 0.2,0.4, 3];
 min_i = [0.3,0.001, 8, 0.5,0.1,0.2, 1];

pidx{1}=[1 3 5 7];
pidx{2}=1:7;
pidx{3}=1:7;

par_true = cell(1,3);

 for i = 1:length(model_list)
     par_true{i}=repmat(max_i(pidx{i})-min_i(pidx{i}),tests,1).*rand(tests,length((pidx{i})))+repmat(min_i(pidx{i}),tests,1);
     if  model_list{i}~="RL_race"
         idx = find(par_true{i}(:,6)<=par_true{i}(:,5));
         par_true{i}(idx,6) = par_true{i}(idx,5)+rand(length(idx),1);
     end

 end

 %% Create Data:
 %%% Set up variables
 MAP = [1,1,1;
     2,2,2;
     3,3,4;
     4,4,3];
 Trials = [1 500 4000 500];

 data = cell(1,tests*3);
 all_models_true = cell(1,tests*3);
 parfor i = 1:tests*3
     if i <=tests
         m = 1;
         j = i;
     elseif i <= tests*2
         m=2;
         j = i-tests;
     else
         m = 3;
         j = i-tests*2;
     end
     [all_models_true{i},data{i}]=Simulate_Hardwick(model_list{m},Trials,par_true{m}(j,:),MAP);
     true_model(i)=m;
 end

%%
%%% Fitting to all subjects

 parfor i = 1:tests%*3
   
     s_data = data{i};

     %%% Clean data
     % s_data(s_data.Outcome==0,:)=[];
     % s_data.Outcome(s_data.Outcome<0)=0;
     s_data.Reaction_Time(s_data.Reaction_Time<0)=0.001;%[];
     s_data(s_data.Reaction_Time>5,:)=[];
     s_data(s_data.Choice==0,:)=[];
     s_data(isnan(s_data.Choice),:)=[];
     TrialNum(i)=height(s_data);

     
     CHOICE = s_data.Choice;
     STIMULUS = s_data.Stimulus;
     RT = s_data.Reaction_Time;
     OUTCOME = s_data.Outcome;
     FORCE = s_data.Forced;
     MAP = [1,1,1;
         2,2,2;
         3,3,4;
         4,4,3];

     for m = 1%:3
         try
         %%% Using optional arguments
         MODEL = model_list{m}; %
         FORCE(s_data.Block==4)=1; % Forced choice trials were included in Block 4
         REMAP_TRIAL = [find(s_data.Block==3&s_data.Exp==1,1,"first"),find(s_data.Block==3&s_data.Exp==2,1,"first")];
         INCLUDE_FIT = ones(1,height(s_data));
         INCLUDE_FIT(s_data.First~=1) =0 ; % only include first responses in fitting
         INCLUDE_FIT((s_data.Block==1|s_data.Block==2) & s_data.TrialNum>50) = 0;
         WEIGHT = 0.75;
         NEW_COND = find(s_data.Exp==2,1,"first");

         %%% Run fitting

         [PARAM, NLL_m, LL_m, MODEL_obj] = fitTwoDriftRM(CHOICE,STIMULUS,RT,OUTCOME,MAP,FORCE,...
             model=MODEL,...
             remap_trial = REMAP_TRIAL, ...
             include_fit=INCLUDE_FIT, ...
             weight=WEIGHT,...
             new_cond = NEW_COND);

         par_out{m,i}=PARAM;
         NLL(m,i)=NLL_m;
         LL{m,i}=LL_m;
         all_models_fit{m,i}=MODEL_obj;
         catch
         end
     end
  
     
 end
 save("Hardwick_project\MR_075.mat")
 
 %%


 BIC = 2.*(NLL)+(repmat(par_num',1,((length(model_list)).*tests))).*log(repmat(TrialNum,length(model_list),1));
[~,pref_model] = min(BIC,[],1);

% for i = 1:length(true_model)
%     true_model_string(i)=(true_model(i));
%     pref_model_string(i)=(pref_model(i));
% end
model_table = table(true_model',pref_model','VariableNames',{'True','Recovered'});
figure
%%%%% Confusion and Inversion matrix 
model_table.con_val = repmat(1/tests,height(model_table),1);
h = heatmap(model_table,'True', 'Recovered',"ColorVariable","con_val",'ColorMethod','sum');


model_table.in_val = zeros(height(model_table),1);
G = groupsummary(model_table,"Recovered");
for i = 1:height(model_table)
    model_table.in_val(i)=1/G.GroupCount(G.Recovered==model_table.Recovered(i));
end
h = heatmap(model_table,'True', 'Recovered',"ColorVariable","in_val",'ColorMethod','sum');
