function [chose,r,bl,RT,stim_seen, first,trial_num,Q,H] = SimulateB4(model,par,data,Q,H,err)
%% This function creates data for the Hardwick experiment
%%% Input:
%%% model: string name for the model used to determine choices
%%% par: the parameters to be used for the model
%%% par_name: optional input to label the parameters, if included, will
%%% plot data
%%% trials: number of free and forced trials
%%% start: determines start values for Q and H

%% set up parameters:
 rew = 1; 
fun = str2func(model);
%%% Create data structure
data(data.Block ~=4,:)=[];
if size(Q,3)>1
    Q = Q(:,:,end);
end
if size(H,3)>1
    H = H(:,:,end);
end

stim_forced = data.Stimulus;
time_forced = data.RT;


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

% Parameters
chose = zeros(1,length(stim_forced));
bl = zeros(1,length(stim_forced));
RT = zeros(1,length(stim_forced));
r = zeros(1,length(stim_forced));
first = []; t_end = data.Trial(1);


%% Run Forced Response
i = 1;
for t = 1:height(data)
    first(end+1)=1;
    trial_num(t)=t_end+t;
    bl(i) = 4;
    % New value
    Q(:,:,i+1)=Q(:,:,i);
    H(:,:,i+1)=H(:,:,i);
    % Observation Model %%% Need to change to match Sam's race model
    variable{1} = Q(:,stim_forced(t),i);
    variable{2} = H(:,stim_forced(t),i);
    % 
  %   stim_forced(t)
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