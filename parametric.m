function [CIs] = parametric(data,B,alpha,CI_level)
    % ES vec to store ES values
    ESvec=zeros(B,1);
    for j=1:B
        % sample with replacement
        sample = datasample(data,length(data));
        % estimate t-dist parameters
        estimate_param = MLE(sample,[3,0,1]);
        ESvec(j)=ES_t(estimate_param,alpha);
    end
    % confidence interval
    CIs = quantile(ESvec,[(1-CI_level)/2,1-(1-CI_level)/2]);
end
