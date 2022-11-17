function ES = ES_nct(param,alpha)
    df = param(1);
    mu = param(2);
    loc = param(3);
    c = param(4);
    c01 = nctinv(alpha,df,mu); % left tail quantile, for loc-0 scale-1
    fun = @(x) x.*nctpdf(x,df,mu);
    ES = integral(fun,-inf,c01,"AbsTol",1e-5)/alpha;
    ES = ES*c+loc;
end

