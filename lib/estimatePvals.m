function [modelSig, coeffSig]   = estimatePvals(coeffs, y, yPredict, dfe, jacbn)
% function that calculates the model p vals 
[stdErr,~,~,~,vcv]=regdata(coeffs,yPredict,y,jacbn);
[modelSig.pVals, modelSig.Fstat] = linhyptest(coeffs,vcv,[],[],dfe);
coeffSig.stdErr = stdErr;
coeffSig.tstat = coeffs./stdErr;
%coeffSig.pVals = 2*(1 - tcdf(real(coeffSig.tstat),dfe));
end