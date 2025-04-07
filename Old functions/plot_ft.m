function [ft,data_out,remap_prop,out_lines] = plot_ft(param,model,reps,data,l_flag,err,np)
mmn = 200;
x = 1:1600;
col = {'#E200FF','#0A00FF','#FF7B00'};
if isempty(param)
    t_h = 0.25;
else
    t_h = param(5);
end

cons_mm = []; remp_mm = []; hab_mm = [];mean_all=[];std_all=[];

for rep = 1:reps
    cons_tr= []; remap_tr=[]; hab_tr = []; nohab_tr=[];cons = nan(1,1600);remp = nan(1,1600); hab = nan(1,1600);
    if ~isempty(data)
        if isstruct(data)
            d = struct2table(flipstruct(data));
        else
            d = data;

        end

    else
        d = [];
        if length(param)==11
            [d.Response,d.Correct,d.Block,d.RT,d.Stimulus, d.First,d.Trial] = Simulatecont(model,param,[4000 500],err);
        else
            [d.Response,d.Correct,d.Block,d.RT,d.Stimulus, d.First,d.Trial] = Simulate(model,param,[4000 500],err);
            % end
        end
        remap_prop(rep) = sum(d.Block==3 & d.First==1);
    end

    if  any(d.Correct<0) & err == 0 %& %% we always want to be working between
    %0 and 1 here
        d(d.Correct==0,:)=[]; % Remove trials during punishment break
        d.Correct(d.Correct<0,:)=0;
    elseif any(d.Correct<0) 
        d.Correct(d.Correct<0)=0;
    end
    
    % d.Response(d.Response ==0)=nan;
    % d.RT(d.RT<=0)=0.001;
    % d.RT(d.RT>5)=5;
    data_out{rep}=d;
    hab = nan(1,1600);
    cons_tr = find(d.Block==4& (d.Stimulus ==1 |d.Stimulus == 2));
    remap_tr = find(d.Block==4& (d.Stimulus ==3 |d.Stimulus == 4));
    hab_tr = find(d.Block==4& ((d.Stimulus ==3 & d.Response==3 )|(d.Stimulus == 4&d.Response==4)));
    nohab = find(d.Block==4& ((d.Stimulus ==3 & d.Response~=3 )|(d.Stimulus == 4&d.Response~=4)));
    cons(round(d.RT(cons_tr)*1000))= d.Correct(cons_tr);
    remp(round(d.RT(remap_tr)*1000))=d.Correct(remap_tr);
    hab(round(d.RT(hab_tr)*1000))= 1;
    hab(round(d.RT(nohab)*1000))= 0;
    hab_av=movmean(hab,mmn,'omitnan');

    
    if ~isempty(param)
        t_q(rep)=param(6);%approx_tq(d,[]);
    elseif isempty(data)
        t_q(rep)=approx_tq(d,[]);
        t_h(rep)=approx_th(d,[],t_q(rep));
    else
        t_q(rep)=0;
        t_h(rep)=0;
    end
    %  t_q(rep) = find(hab_av(1:800)>remp_av(1:800),1,'last')./1000;
    % t_q(rep)=(find(hab_av(150:700) == max(hab_av(150:700)),1,'last')+149)./1000;
    %t_q(rep)=(find(hab_av(150:700) == max(hab_av(150:700)),1,'last')+149)./1000;
    cons_mm(rep,:) = movmean(cons,mmn,"omitnan");
    remp_mm(rep,:) = movmean(remp,mmn,"omitnan");
    hab_mm(rep,:) = movmean(hab,mmn,"omitnan");
end
cons_mm = cons_mm(:,1:1600);
remp_mm = remp_mm(:,1:1600);
hab_mm = hab_mm(:,1:1600);


if reps>1
    mean_all(1,:) =  mean(cons_mm);  mean_all(2,:) =  mean(remp_mm);  mean_all(3,:) =  mean(hab_mm);
    std_all(1,:) =  std(cons_mm);  std_all(2,:) =  std(remp_mm);  std_all(3,:) =  std(hab_mm);

else
    mean_all(1,:) =  (cons_mm);  mean_all(2,:) =  (remp_mm);  mean_all(3,:) =  (hab_mm);

end
if isempty(np)   
    hold on
    for j = 1:3
        if reps>1
            xconf = [x x(end:-1:1)] ;
            yconf = [mean_all(j,:)+std_all(j,:) mean_all(j,end:-1:1)-std_all(j,end:-1:1)];
            p = fill(xconf,yconf,'r');
            p.FaceColor=col{j};
            p.FaceAlpha = 0.4;
            p.EdgeColor = 'none';
        end
    plot(x,mean_all(j,:),'Color',col{j},'LineWidth',1)
    end
    yline(0.25,'--')
    xline(t_h*1000,'--')
    switch model
        case "Simple_Habits"
            xline(mean(t_q)*1000,'--')
        case "Qlearn2"
            xline(mean(t_q)*1000,'--')
        case "Qlearn"
    end
    xlim([0,1200])
 %   ylim([0 1.1])
    if l_flag ==1
        legend({'Consistent','Remapped','Habitual'});
    end

    ft = gcf;
else
    ft = [];
end
%if rep ==1
%    out_lines = [cons_mm;remp_mm;hab_mm];
%else
    out_lines = mean_all;
%end
end