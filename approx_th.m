function [th_es,tq_es]=approx_th(obj)
%% Estimating t_h from data
%%% This function provides an estimate of t_h based on the last time the
%%% forced trial habit curve was under 0.33, before it reaches its peak
%%% value.
%%% Input:
% obj = Model object containing all required data
%%% Output:
% th_es = estimate of t_h
% tq_es = time of peak habit value
    
    mmn =200;
    d = obj.data(obj.data.Forced==1,:);
    r_trial = find(d.Remap_trial==1);
    map = obj.map;
    %%% Finding idx of trials with a remapped stimulus
    hab_tr = zeros(2,height(d)); % Trials where a remapped stimulus was seen
    


    for i = 1:size(map,1)
        idx_h = find(d.Stimulus==map(i,1) & d.Choice ==map(i,2));
        idx_c = find(d.Stimulus==map(i,1) & d.Choice ==map(i,3));
        if ~isempty(r_trial)
        idx_h(r_trial(idx_h)==0)=[];
        idx_c(r_trial(idx_h)==0)=[];
        end
        hab_tr(1,idx_h)=1;
        hab_tr(2,idx_c)=1;
    end
    hab_tr = logical(hab_tr);
    hab = nan(1,round(max(d.Reaction_Time))*1000);
    hab(uint16(d.Reaction_Time(hab_tr(1,:))*1000))= 1;
    hab(uint16(d.Reaction_Time(hab_tr(2,:))*1000))= 0;
    hab_av=movmean(hab,mmn,'omitnan');
    [tq_es,tq_idx] = max(hab_av(150:end));
    th_es = find(hab_av(1:tq_idx)<=0.33,1,"last")./1000;
    
    if isempty(th_es)
       th_es = 0.15;
    end
    
end