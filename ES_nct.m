function ES = ES_nct(param,alpha)
    df = param(1);
    mu = param(2);
    loc = param(3);
    c = param(4);
    c01 = nctinv(alpha,df,mu); % left tail quantile, for loc-0 scale-1
    ES01 = -nctpdf(c01,df,mu)/nctcdf(c01,df,mu)*(df+c01^2)/(df-1);
    ES = loc+c*ES01; % true theoretical ES
end
