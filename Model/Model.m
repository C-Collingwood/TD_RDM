classdef Model
    %% Class that contains all model data and information
    properties
        lr
        free_par
        par
        type
        values
        par_bound
        data
        map
        weight
        fit
    end
    methods
        function [obj] = Model(name,par,n_trial,map)
            %%% Creates object
            obj.type = name;
            n_choice = length(unique(map(:,2:3)));
            n_stim = length(unique(map(:,1)));
            out_vec = zeros(n_choice,n_stim,n_trial); 

            switch name
                case {"Habit_race"}
                    obj.lr=2;
                    obj.free_par = {'aq','ah','bq','bh','tq','theta'};
                    obj.par = array2table([par,1],'VariableNames',{'aq','ah','bq','bh','th','tq','theta','s'}); 
                    obj.values = struct("Q",0.5*ones(size(out_vec)),"H",out_vec,'dq',out_vec,'dh',out_vec);
                case {"RL2_race"}
                    obj.lr=2;
                    obj.free_par = {'aq2','aq1','bq2','bq1','tq2','theta'}; % All except t_q1
                    obj.par =  array2table([par,1],'VariableNames',{'aq2','aq1','bq2','bq1','tq1','tq2','theta','s'});  
                    obj.values = struct("Q2",0.5*ones(size(out_vec)),"Q1",0.5*ones(size(out_vec)),'dq2',out_vec,'dq1',out_vec);
                case {"RL_race"}
                    obj.lr=1;
                    obj.free_par = {'aq','bq','theta'};
                    obj.par = array2table([par,1],'VariableNames',{'aq','bq','tq','theta','s'}); 
                    obj.values = struct("Q",0.5*ones(size(out_vec)),'dq',out_vec);

            end
           % f = obj.par.Properties.VariableNames;
            obj.weight = 0.5;
            obj.map = map;
            obj.data=[];
            obj.par_bound = [];
        end


        function [obj] = initialise(obj,choice,stim,rt,outcome,force,par_bound,r_trial)
            n_trial = length(choice);
            %%% Creates object
            if ~isempty(par_bound)
            obj.par_bound=array2table(par_bound,"RowNames",{'UB','LB'},"VariableNames",obj.free_par);
            end
            remap_trial = zeros(1,n_trial); remap_trial(r_trial:end)=1;
            obj.data = table(choice',stim',rt',outcome',force',remap_trial',VariableNames= {'Choice','Stimulus','Reaction_Time','Outcome','Forced','Remap_trial'});
        
        end
        
        %%% Approximate T_start
        function [th_es,tq_es]=approx_th(obj)
            %%% This function provides an estimate of t_h based on the last time the
            %%% forced trial habit curve was under 0.25, before it reaches its peak
            %%% value.
            %%% Input:
            % obj = Model object containing all required data
            %%% Output:
            % th_es = estimate of t_h
            % tq_es = time of peak habit value

            mmn =200;
            d = obj.data(obj.data.Forced==1,:);
            r_trial = find(d.Remap_trial==1);
            Map = obj.map;
            %%% Finding idx of trials with a remapped stimulus
            hab_tr = zeros(2,height(d)); % Trials where a remapped stimulus was seen



            for i = 1:size(Map,1)
                idx_h = find(d.Stimulus==Map(i,1) & d.Choice ==Map(i,2));
                idx_c = find(d.Stimulus==Map(i,1) & d.Choice ==Map(i,3));
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
            [~,tq_idx] = max(hab_av(150:end));
            tq_idx=tq_idx+150;
            tq_es = tq_idx(end)./1000;
             th_es = find(0.25-hab_av(1:tq_idx)<=0,1,'last')./1000;
                
            
            if isempty(th_es)
                th_es = 0.15;
            end

        end



        %%% update model object to account for time parameter heuristics
        function [obj] = update_time(obj,nd)
            [t1,t2]=approx_th(obj);
            if ~isempty(nd) & nd>0
                t1 = nd;
            end
            switch obj.type
                case {'Habit_race'}
                    obj.par.th=t1;
                    obj.par.tq=t2;
                    obj.par_bound("LB","tq")={t1+0.001};
                case {"RL2_race"}
                    obj.par.tq1=t1;
                    obj.par.tq2=t2;
                    obj.par_bound("LB","tq2")={t1+0.001};
                case {"RL_race"}
                    obj.par.tq = t1;
            end
        end

        
        %%% Trials with RT during non-decision time are removed
        function [fit]=remove_fast(obj,fit)
                switch obj.type
                    case {"Habit_race"}
                        t_start = 'th';
                    case {"RL2_race"}
                        t_start = 'tq1';
                    case {"RL_race"}
                        t_start = 'tq';
                end
                
                fit(obj.data.Reaction_Time <= obj.par.(t_start)+0.001) = 0;
              
        end



        %%%% update model with fitted free parameters
        function [obj] = update_param(obj,p)
            f = obj.free_par;
            for i = 1:length(p)
                obj.par.(f{i})=p(i);
            end
        end

    end
end