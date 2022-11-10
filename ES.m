function es = ES(param,alpha)
    c01 = tinv(alpha , param(1)); % left tail quantile, for loc-0 scale-1
    es = -tpdf(c01,param(1))/tcdf(c01,param(1)) * (param(1)+c01^2)/(param(1)-1);
    es = param(2)+param(3)*es;
end

