function [allSignificantResults, geneName_significant_model, pathName_significant_model, ytest_sig, idAllSig, xtest] = conductDRM_geneLevel (dose, geneVls, eqtype, sffx,folderName, geneName, pathName, pathUni)

% geneVls is 2D: gene*conc
%% genes
numGene = size(geneVls,1);
result_summary = nan(numGene,18);
xtest = (logspace(log10(dose(1,1)), log10(dose(end,1)),100))';
for i = 1:numGene
    try
    response = (geneVls(i,:));
    [cf, G, output] = fourParameterDRM(dose, response, eqtype);
    [modelSig, ~] = estimatePvals((coeffvalues(cf))', response, cf(dose), G.dfe, output.Jacobian);
    yPredict = cf(xtest);
    yM = max(yPredict);    
    [BMD, yBMD, BMD_1, yBMD_1] = estimateBMD(cf,dose,yM, 10, eqtype);
    [BMD5, ~] = estimateBMD(cf,dose,yM, 5, eqtype);
    BMDratio = BMD/BMD5;
    if i == 1
        [resultVls, hdrTexts] = prepareResults(cf, G, modelSig);
        rowNames = [hdrTexts'; {'BMD'; 'yBMD'; 'BMD_1'; 'yBMD_1';'TELI1.5';'TELImax';'BMDratio'}];
    else
        [resultVls] = prepareResults(cf, G, modelSig);
    end
    TELI_15 = estimateX (1.5, cf, eqtype);
    result_summary (i, :) = [resultVls; BMD; yBMD; BMD_1; yBMD_1; TELI_15; yM;BMDratio]; 
    yTest(:,i) = yPredict(:);
   
    catch ME
        fileID = fopen(strcat(folderName,'/errorException_',sffx, '.txt'),'a');
        %fprintf(fileID,'%d GeneName: %s',i, char(geneName(i)));
        fprintf(fileID,'%d\n',i);
    end
end
save (strcat(folderName,'/gene_modelFitResult_',sffx, '.mat'), 'rowNames');
ind = result_summary(:,6) <= 0.05; % get index for p-values <= 0.05
significantResults = result_summary(ind,:);
gName_significant_model = geneName(ind);
pName_significant_model = pathName(ind);
ytest_sig_tmp = yTest(:,ind);
idAllSig = 1:numGene;
idAllSig(~ind) = [];
ind = significantResults(:,12) ~= 0.0; % get index for BMD_10 == 0
allSignificantResults = significantResults(ind,:);
geneName_significant_model = gName_significant_model(ind);
pathName_significant_model = pName_significant_model(ind);
ytest_sig = ytest_sig_tmp(:,ind);
idAllSig(~ind) = [];
save (strcat(folderName,'/gene_modelFitResult_',sffx, '.mat'), 'allSignificantResults',...
    'geneName_significant_model', 'pathName_significant_model', 'idAllSig','-append');
[numGene,pathGrp,minVal, maxVal, meanVal, medianVal, ciVal] = grpstats(allSignificantResults,pathName_significant_model,{'numel','gname','min','max','mean','median','meanci'});
[~,~,ia] = intersect(pathUni,pathGrp,'stable');
pathUnique = pathGrp(ia); numGene = numGene(ia,:); minVal = minVal(ia,:);
maxVal = maxVal(ia,:); medianVal = medianVal(ia,:); meanVal = meanVal(ia,:); ciVal = ciVal(ia,:);
save (strcat(folderName,'/gene_modelFitResult_',sffx, '.mat'), 'numGene', 'pathUnique',...
        'minVal','maxVal', 'meanVal', 'ciVal', 'medianVal', '-append');
end

