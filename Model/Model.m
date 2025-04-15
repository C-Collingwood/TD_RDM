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