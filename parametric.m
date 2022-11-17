function [CIs] = parametric(data,B,alpha,CI_level,parametric_model)
    % ES vec to store ES values
    ESvec=zeros(B,1);
    for j=1:B
        % sample with replacement
        sample = datasample(data,length(data));

        % estimate parameters based on the dist name
        estimate_param=0;
        if parametric_model == "T"
            estimate_param = MLE_t(sample,[3,0,1]);
            ESvec(j)=ES_t(estimate_param,alpha);
        end
        if parametric_model == "NCT"
            estimate_param = MLE_nct(sample,[6, -2, 0, 1]);
            ESvec(j)=ES_nct(estimate_param,alpha);
        end

        
    end
    % confidence interval
    CIs = quantile(ESvec,[(1-CI_level)/2,1-(1-CI_level)/2]);
end
