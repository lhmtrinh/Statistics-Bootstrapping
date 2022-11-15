function [bool_nonparam,length_nonparam, bool_param,length_param] = repetition(dist, T, rep, B, true_ES,alpha, confidence_level)
   
    bool_nonparam = zeros(1,rep);
    length_nonparam = zeros(1,rep);
    bool_param = zeros(1,rep);
    length_param = zeros(1,rep);
    
    f = waitbar(0, 'Starting');
    for i=1:rep
        waitbar(i/rep, f, sprintf('Progress for T=%d : %d %%',T,floor(i/rep*100)));
        % data generation process
        data= random(dist,T,1);
    
        % calculate confidence intervals for ES with non parametric bootstrap
        ci=non_parametric(data,B,alpha,confidence_level);
        % calculate the lenght of interval
        length_nonparam(i) = ci(2)-ci(1);
        % calculate coverage array
        bool_nonparam(i) = (true_ES>=ci(1))&(true_ES<=ci(2));
        
        % calculate confidence intervals for ES with non parametric bootstrap
        ci=parametric(data,B,alpha,confidence_level);
        % calculate the lenght of interval
        length_param(i) = ci(2)-ci(1);
        % calculate coverage array
        bool_param(i) = (true_ES>=ci(1)) & (true_ES<=ci(2));
    end
    close(f)
end

