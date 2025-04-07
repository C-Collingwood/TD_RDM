function [tq,mh]=approx_tq(data,mmn)
    if isempty(mmn)
        mmn=200;
    end
    if ismember('Exp', fieldnames(data)) 
        data=data(data.Exp==2,:);
    end
    if istable(data)
        d = data;
    elseif isstruct(data)
        d = struct2table(flipstruct(data));
    
    end
    d.RT(d.RT<=0)=0.001;
    d(d.RT>5,:)=[];
    d.Response(d.Response==0)=nan;
    d.Correct(d.Correct<0)=0;
    hab = nan(1,1600);
 
    hab_tr = find(d.Block==4& ((d.Stimulus ==3 & d.Response==3 )|(d.Stimulus == 4&d.Response==4)));
    nohab = find(d.Block==4& ((d.Stimulus ==3 & d.Response~=3 )|(d.Stimulus == 4&d.Response~=4)));

    hab(round(d.RT(hab_tr)*1000))= 1;
    hab(round(d.RT(nohab)*1000))= 0;
    hab_av=movmean(hab,mmn,'omitnan');
    tq =(find(hab_av(150:700) == max(hab_av(150:700)),1,'last')+149)./1000;
    mh = max(hab_av(150:700));


end
