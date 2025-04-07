function [th]=approx_th(data,mmn,tq)
    load trainedEstimator.mat trainedEstimator
    if isempty(mmn)
        mmn=200;
    end
    if ismember('Exp', fieldnames(data)) 
        data=data(data.Exp==2,:);
    end
    if isstruct(data)
        d = struct2table(flipstruct(data));
    else
        d = data;
    end
    hab = nan(1,1600);
    cons = nan(1,1600);
    remp = nan(1,1600);
    d.Correct(d.Correct<0)=0; %% So that our plots look right
    cons_tr = find(d.Block==4& (d.Stimulus ==1 |d.Stimulus == 2));
    remap_tr = find(d.Block==4& (d.Stimulus ==3 |d.Stimulus == 4));
    hab_tr = find(d.Block==4& ((d.Stimulus ==3 & d.Response==3 )|(d.Stimulus == 4&d.Response==4)));
    nohab = find(d.Block==4& ((d.Stimulus ==3 & d.Response~=3 )|(d.Stimulus == 4&d.Response~=4)));
    cons(round(d.RT(cons_tr)*1000))= d.Correct(cons_tr);
    remp(round(d.RT(remap_tr)*1000))=d.Correct(remap_tr);
    hab(round(d.RT(hab_tr)*1000))= 1;
    hab(round(d.RT(nohab)*1000))= 0;
    cons_av=movmean(cons,mmn,'omitnan');
    remp_av=movmean(remp,mmn,'omitnan');
    
    hab_av=movmean(hab,mmn,'omitnan');
    std_all = (std([cons_av(1:500);remp_av(1:500);hab_av(1:500)]));
    mi = ((movmean(MI(d(d.Block==4,:),50),100)));
    th_es = find(hab_av(1:tq*1000)<=0.25,1,"last")./1000;
    if isempty(th_es)
       th_es =  find(hab_av(1:tq*1000)<=0.33,1,"last")./1000;
    end
    if isempty(th_es)
       th_es = 0.001;
    end
    
    T = table(th_es,mi(1:500),std_all(1:500),cons_av(1:500),remp_av(1:500),hab_av(1:500),'VariableNames',{'Th_estimate','MI', 'std','cons','remap','habitual'});
    th = trainedEstimator.predictFcn(T);

end