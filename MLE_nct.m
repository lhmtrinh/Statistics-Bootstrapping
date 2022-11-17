function mle = MLE_nct(x,initvec)
    tol=1e-5;
    opts=optimset('Disp','none','LargeScale','Off','TolFun',tol, ...
        'TolX ',tol,'Maxiter',200);
    % minimum of negative log likelihood
    % maximize loglikelihood
    mle = fminunc(@(param) nctloglik_dda(param,x), initvec, opts);
end

function ll = nctloglik_dda(param,x)
    % param=[df, noncentrality parameter, location, scale]
    df = param(1); 
    mu = param(2);
    loc = param(3); 
    c = param(4);
    % prevent negative values
    if df<0 
        df=rand;
    end
    if c<0 
        c=rand;
    end
    % prevent positive values
    if mu>0 
        mu=rand;
    end
    z = (x-loc)./c;
    pdfln = -log(c)+stdnctpdfln_j(z, df, mu);
    ll = -sum(pdfln);
end

