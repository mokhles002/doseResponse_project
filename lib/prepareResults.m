function[resultVls, hdrTexts] = prepareResults(cf, G,  modelSig)

%function[resultVls, hdrTexts] = prepareResults(cf, G, s, modelSig, coeffSig)
% This function prepare the results in a vector format
%{
if s == 1 || 3
    n = 1;
else
    n = cf.n;
end
%}
T = cf.T; B = cf.B; EC50 = cf.C; n = cf.n;
temp = struct2table(G);
tmpHdr3 = temp.Properties.VariableNames;
vls = struct2array (G);
tmpHdr1 = {'Bottom', 'Slope', 'EC50', 'Top'};
tmpHdr2 = {'F-statistics', 'p-value'};
hdrTexts = [tmpHdr1, tmpHdr2, tmpHdr3];
resultVls = [B; n; EC50; T; modelSig.Fstat; modelSig.pVals; vls'];
end