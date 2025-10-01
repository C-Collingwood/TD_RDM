classdef Model_obj
    %% Class that contains all model data and information
    %%% Class creation function is Model_obj() for the provided RL-Race,
    %%% RL2-Race, Habit1-Race and Habit2-Race functions.
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
        function [obj] = Model_obj(name,par,t1,n_trial,map)
            %%% Creates object
            obj.type = name;
            n_choice = length(unique(map(:,2:3)));
            n_stim = length(unique(map(:,1)));
            out_vec = zeros(n_choice,n_stim,n_trial); 

            switch name
                case {"Habit1_Race"}
                    obj.free_par = {'aq','ah','bq','bh','t2','theta'};
                    obj.par = array2table([par,t1,1],'VariableNames',{'aq','ah','bq','bh','t2','theta','t1','s'}); 
                    obj.values = struct("Q",0.5*ones(size(out_vec)),"H",out_vec,'dq',out_vec,'dh',out_vec);
                case {"Habit2_Race"}
                    obj.free_par = {'aq','ah','bq','bh1','bh2','t2','theta'};
                    obj.par = array2table([par,t1,1],'VariableNames',{'aq','ah','bq','bh1','bh2','t2','theta','t1','s'}); 
                    obj.values = struct("Q",0.5*ones(size(out_vec)),"H",out_vec,'dq',out_vec,'dh',out_vec);
                case {"RL2_Race"}
                    obj.free_par = {'aq2','aq1','bq2','bq1','t2','theta'}; % All except t_1
                    obj.par =  array2table([par,t1,1],'VariableNames',{'aq2','aq1','bq2','bq1','t2','theta','t1','s'});  
                    obj.values = struct("Q2",0.5*ones(size(out_vec)),"Q1",0.5*ones(size(out_vec)),'dq2',out_vec,'dq1',out_vec);
                case {"RL_Race"}
                    obj.free_par = {'aq','bq','theta'};
                    obj.par = array2table([par,t1,1],'VariableNames',{'aq','bq','theta','t1','s'}); 
                    obj.values = struct("Q",0.5*ones(size(out_vec)),'dq',out_vec);
            end
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
            obj.data = table(choice',stim',rt',outcome',force',remap_trial',VariableNames= {'Choice','Stimulus','Reaction_Time','Outcome','Time_Controlled','Remap_Trial'});
        end
        
     

        
        %%% Trials with RT during non-decision time are removed
        function [fit]=remove_fast(obj,fit)
                fit(obj.data.Reaction_Time <= obj.par.t1+0.001) = 0;
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