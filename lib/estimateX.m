function x = estimateX(y, cf, eqtype)
% estimate x form a y value
% use dose response curve- hill type and exponential
% eqtype:   'hill'
%           'exponential'

T = cf.T; B = cf.B; C = cf.C; n = cf.n;

switch eqtype
    case 'hill'
        x = C * (((y-B)/(T-y))^(1/n));
    case 'exponential'
        x = (1/C * log((T-B)/(T-y)))^(1/n);
end
end