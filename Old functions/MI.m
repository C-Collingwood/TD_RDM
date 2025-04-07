function [mi]=MI(data,window)

%%% This function calculates the mutual information between stimulus and
%%% response across RT during forced choice trials:
%%% Takes an (n,3) matrix, where 1 = RT, 2 = Stimulus seen and 3 = Response
%%% made
dt = 0.001;
rt_all = 0:dt:1.6;
IN = zeros(4,4,length(rt_all)); 
[rt,idx]=sort(data.RT);
stimulus = data.Stimulus(idx);
response = data.Response(idx);
RT = round(rt*1000);
RT(RT>1600)=1600;

for i = 1:length(rt)
    idx = round(RT(i));
    if ~isnan(response(i))
        IN(stimulus(i),response(i),idx)=IN(stimulus(i),response(i),idx)+1;
    end

end


for i = 1:length(rt)-window
    mm = RT(i):RT(i+window);
    ps = zeros(4,1);
    pr = zeros(1,4);
    psr=zeros(4,4);
    for stim = 1:4
        ps(stim,1)=sum(IN(stim,:,mm)>0,'all','omitnan')./window;
        for resp = 1:4
            pr(1,resp)=sum(IN(:,resp,mm)>0,'all','omitnan')./window;
            psr(stim,resp)=sum(IN(stim,resp,mm)>0,'all','omitnan')./window;
        end
    end
    mi(RT(i+round(window/2)))=sum(psr.*log(psr./(ps*pr)),'all','omitnan');


end




end