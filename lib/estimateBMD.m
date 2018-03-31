function [BMD, yBMD, BMD_1, yBMD_1] = estimateBMD(cf, dose, yM, BMDref, eqtype)
% estimate BMD form the dose response equation
% BMDref: reference BMD 

% BMD with respect to control teli at dose 1.
y0 = cf(dose(1));
%yM = cf(dose(end));
yBMD = y0 + ((yM-y0)*BMDref/(100*yM));
BMD = estimateX(yBMD, cf, eqtype);

% BMD with respect to control teli = 1; so y0 = 1
yBMD_1 = 1 + ((yM-1)*BMDref/(100*yM));
BMD_1 = estimateX(yBMD_1, cf, eqtype); 
end