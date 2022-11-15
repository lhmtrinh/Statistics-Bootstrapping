function [CIs] = non_parametric(data,B,alpha,CI_level)
    %ES vec to store ES after each bootstrap
    ESvec=zeros(B,1);
    for j=1:B
        % sample with replacement
        sample = datasample(data,length(data));
        % VaR and ES
        VaR=quantile(sample,alpha); 
        temp=sample(sample<=VaR); 
        ESvec(j)=mean(temp);
    end
    % confidence interval
    CIs = quantile(ESvec,[(1-CI_level)/2,1-(1-CI_level)/2]);
end
