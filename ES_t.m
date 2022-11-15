function es = ES_t(param,alpha)
    df = param(1);
    location = param(2);
    scale = param(3);
    c01 = tinv(alpha , df); % left tail quantile, for loc-0 scale-1
    es = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
    es = location+scale*es;
end

