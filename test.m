
%check for outlier
geneTELI_repli_n(:,:,1) = geneTELI_repli(:,1:6);
geneTELI_repli_n(:,:,2) = geneTELI_repli(:,7:12);
geneTELI_repli_n(:,:,3) = geneTELI_repli(:,13:18);

aa(:,:,1) = pathTELI_repli(:,1:6);aa(:,:,2) = pathTELI_repli(:,7:12);aa(:,:,3) = pathTELI_repli(:,13:18);
pathTELI_repli_n = aa;


for i = 1:3
    tmp_repli = geneTELI_repli_n(:,:,i);
    tmp_repli(tmp_repli > geneTELI(:,:,1)+1.38*geneTELI_std) = NaN;
    tmp_repli(tmp_repli < geneTELI(:,:,1)-1.38*geneTELI_std) = NaN;
    gneOtlier(i,:) = sum(isnan(tmp_repli));
    tmp_repli = pathTELI_repli_n(:,:,i);
    tmp_repli(tmp_repli > pathTELI(:,:,1)+1.38*pathTELI_std) = NaN;
    tmp_repli(tmp_repli < pathTELI(:,:,1)-1.38*pathTELI_std) = NaN;
    pthOtlier(i,:) = sum(isnan(tmp_repli));
end

pathTELI_min = min(aa,[],3);

%% 
addpath lib;
load data/data_vF.mat;
load data/timepoint_rawData.mat;

[result_summary_teli, geneNameSig_teli, pathNameSig_teli, yTestSig_teli, idSig_teli, xTest] = ...
        conductDRM_geneLevel (conc, geneTELI(:,:,1), 'hill', 'gene_teli','test', ...
        geneName, pathName, uniquePathName_repli(1:5));

[result_summary_t0, geneNameSig_t0, pathNameSig_t0, yTestSig_t0, idSig_t0] = ...
        conductDRM_geneLevel (conc, data_t0(:,:,1), 'hill', 'time_0','test', ...
        geneName, pathName, uniquePathName_repli(1:5));
[result_summary_t20, geneNameSig_t20, pathNameSig_t20, yTestSig_t20, idSig_t20] = ...
        conductDRM_geneLevel (conc, data_t20(:,:,1), 'hill', 'time_20','test', ...
        geneName, pathName, uniquePathName_repli(1:5));
[result_summary_t40, geneNameSig_t40, pathNameSig_t40, yTestSig_t40, idSig_t40] = ...
        conductDRM_geneLevel (conc, data_t40(:,:,1), 'hill', 'time_40','test', ...
        geneName, pathName, uniquePathName_repli(1:5));
[result_summary_t60, geneNameSig_t60, pathNameSig_t60, yTestSig_t60, idSig_t60] = ...
        conductDRM_geneLevel (conc, data_t60(:,:,1), 'hill', 'time_60','test', ...
        geneName, pathName, uniquePathName_repli(1:5));
[result_summary_t120, geneNameSig_t120, pathNameSig_t120, yTestSig_t120, idSig_t120] = ...
        conductDRM_geneLevel (conc, data_t120(:,:,1), 'hill', 'time_120','test', ...
        geneName, pathName, uniquePathName_repli(1:5));


clear chemName chemNameAll clr data_t0 data_t120 data_t20 data_t40 data_t60 dataInd_mean dataInd_repli
clear geneTELI geneTELI_repli geneTELI_std path_t0 path_t120 path_t20 path_t40 path_t60
clear pathTELI pathTELI_repli pathTELI_std
clear pathUnique_t0 pathUnique_t120 pathUnique_t20 pathUnique_t40 pathUnique_t60
save test/allBiomarkerResultSummary.mat


