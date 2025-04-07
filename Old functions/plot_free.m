function [ft,out_lines,len] = plot_free(data,mmn,l_flag,np,rt)

if exist("rt","var")==0
    rt=0;
end
block_names = {'EL','CA','CB','FT'};
col = {'#E200FF','#0A00FF','#FF7B00'};
if isempty(mmn)
    mmn = 10;
end

if isstruct(data)
    d = struct2table(flipstruct(data));
else
    d = data;
end

if any(d.Correct<0)
    d(d.Correct==0,:)=[]; % Remove trials during punishment break
    d.Correct(d.Correct<0,:)=0;
end
d = d(d.Block~=4,:);

if rt==0
    block_change(1,:) = find(diff(d.Block)~=0);
    block_change(:,end+1)=height(d);
    block_change = [0, block_change];
    x = 0;
    hold on
    for i = 1:length(block_change)-1
        bl = d.Block(block_change(i)+1);
        d_bl = d(block_change(i)+1:block_change(i+1),:);

        cons{i}=nan(1,height(d_bl));
        remap{i}=nan(1,height(d_bl));
        cons{i}(d_bl.Stimulus ==1 |d_bl.Stimulus == 2)=d_bl.Correct(d_bl.Stimulus ==1 |d_bl.Stimulus == 2);
        remap{i}(d_bl.Stimulus ==3 |d_bl.Stimulus == 4)=d_bl.Correct(d_bl.Stimulus ==3 |d_bl.Stimulus == 4);


        x = x(end)+1:x(end)+length(cons{i});
        if np==0
        plot(x,movmean(cons{i},mmn,'omitnan','Endpoints','fill'),'Color',col{1},'LineWidth',1)
        plot(x,movmean(remap{i},mmn,'omitnan','Endpoints','fill'),'Color',col{2},'LineWidth',1)
        xline(x(end),'--')%,block_names{d_bl.Block(1)+1})
        ylim([0 1.1])
        end
        len(i)=(height(d_bl));
    end
    if np==0 &  l_flag ==1
        legend({'Consistent','Remapped'});

    
    end
    ft = gcf;
    out_lines={(cons),(remap)};
else
    free_RT = d.RT(d.Block~=4 & d.First==1);
    histogram(free_RT(isoutlier(free_RT)==0))
    %xlim([0.3 1])
    out_lines= {};
    ft = gcf;
    len = [];
end
end