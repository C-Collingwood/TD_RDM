function [LL] = forced_cost(mu_chose,mu_other,s)
%% Cost function for time-controlled trials


fun = @(x)fcost(x,mu_chose,mu_other,s);
P = integral(fun,-Inf,Inf);
if P<=0
    P=0.0000000001;
end
LL = log(P);
    %%% Function to be integrated
    function [cost]=fcost(x,mu_chose,mu_other,s)
        n_stim = length(mu_other);
        non_chose = zeros(length(x),n_stim);
        for i = 1:n_stim
            non_chose(:,i) = normcdf(x,mu_other(i),s);
        end
        cost = normpdf(x,mu_chose,s).*prod(non_chose,2)';
    end


end

