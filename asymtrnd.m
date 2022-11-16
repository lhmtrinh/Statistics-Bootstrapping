function sample = asymtrnd(n_samp, df, mu, loc, c, seed)
% function to sample from the non-central t distribution
% returns an [n_samp 1] vector

% params
% % n_samp:         sample size
% % df:             degrees of freedom
% % mu:             noncentrality parameter
% % loc:            location
% % c:              scale
% % seed:           random seed for reproducability

% Return NaN for elements corresponding to illegal parameter values.
% Non-integer degrees of freedom are allowed.
    df(df <= 0) = NaN;

% Infinite degrees of freedom is equivalent to a normal, and valid.
% Prevent Inf/Inf==NaN for the standardized chi-square in the denom.
    df(isinf(df)) = realmax;

% for reproducable results
    if nargin<4 , seed=rand*1000; end
    rng(seed, 'twister');

% sample from the normal and non-central chisq distribution
    norm_rv = normrnd(mu, 1, [n_samp 1]);
    %chisq_rv = chi2rnd(df, [n_samp 1]); % with 'normal' chisq dist
    noncent_chisq_rv = ncx2rnd(df, 0, [n_samp 1]); % with noncentral chisq dist

% obtain a sample of the non-central t dist by applying the theory
    sample = norm_rv ./ sqrt(noncent_chisq_rv/df);
    sample = loc+c*sample;

end % function