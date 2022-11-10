function mle = MLE(x,initvec)
    tol=1e-5;
    opts=optimset('Disp','none','LargeScale','Off','TolFun',tol, ...
        'TolX ',tol,'Maxiter',200);
    % minimum of negative log likelihood
    % maximize loglikelihood
    mle =  fminunc(@(param) tloglik(param,x),initvec,opts);
end

function ll = tloglik(param,x)
    % param=[df, location, scale]
    v= param(1); 
    mu=param(2); 
    c= param(3);
    % prevent negative values
    if v<0.01 
        v=rand;
    end
    if c<0.01 
        c=rand;
    end
    K=beta(v/2,0.5)*sqrt(v); 
    z=(x-mu)/c ;
    ll = -log(c)-log(K)-((v+1)/2)*log(1+(z.^2)/v); 
    ll = -sum(ll);
end
