function [CIs] = parametric(data,B,alpha,CI_level)
    ESvec=zeros(B,1);
    for j=1:B
        sample = datasample(data,length(data));
        estimate_param = MLE(sample,[3,0,1]);
        ESvec(j)=ES(estimate_param,alpha);
    end
    CIs = quantile(ESvec,[(1-CI_level)/2,1-(1-CI_level)/2]);
end
