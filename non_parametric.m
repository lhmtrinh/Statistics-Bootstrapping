function [CIs] = non_parametric(data,B,alpha,CI_level)
    ESvec=zeros(B,1);
    %CIvec=zeros(rep,1);
    for j=1:B
        sample = datasample(data,length(data));
        VaR=quantile(sample,alpha); 
        temp=sample(sample<=VaR); 
        ESvec(j)=mean(temp);
    end
    % left quantile of ES
    CIs = quantile(ESvec,[(1-CI_level)/2,1-(1-CI_level)/2]);
end
