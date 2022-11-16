function [bool_nonparam,length_nonparam, bool_param,length_param] = repetition(dist,param, T, rep, B, true_ES,alpha, confidence_level,parametric_model)
    % Arrays to store actual coverage array and confidence interval lengths
    % Size 1 x number of reps
    bool_nonparam = zeros(1,rep);
    length_nonparam = zeros(1,rep);
    bool_param = zeros(1,rep);
    length_param = zeros(1,rep);
    
    % waitbar so that we can know the progress
    f = waitbar(0, 'Starting');
    for i=1:rep
        waitbar(i/rep, f, sprintf('Progress for T=%d : %d %%',T,floor(i/rep*100)));

        % data generation process
        data = zeros(T,1);
        if dist=="Stable"
            pdist= makedist("Stable",param(1),param(2),param(3),param(4));
            data= random(pdist,T,1);
        end
        if dist=="T"
            % Create distribution object for random data generation in repetition()
            pdist = makedist("tLocationScale","nu",param(1),"mu",param(2),"sigma",param(3));
            data= random(pdist,T,1);
        end
        if dist=="NCT"
            % Generate T samples of noncentral t distribution, with
            % parameters df = param(1), noncentrality parameter mu =
            % param(2), location = param(3), scale = param(4)
            data = asymtrnd(T, param(1), param(2), param(3), param(4), i);
        end
    
        % calculate confidence intervals for ES with non parametric bootstrap
        ci=non_parametric(data,B,alpha,confidence_level);
        % calculate the lenght of interval
        length_nonparam(i) = ci(2)-ci(1);
        % calculate coverage array
        bool_nonparam(i) = (true_ES>=ci(1))&(true_ES<=ci(2));
        
        % calculate confidence intervals for ES with parametric bootstrap
        ci=parametric(data,B,alpha,confidence_level,parametric_model);
        % calculate the lenght of interval
        length_param(i) = ci(2)-ci(1);
        % calculate coverage array
        bool_param(i) = (true_ES>=ci(1)) & (true_ES<=ci(2));
    end
    close(f)
end

