function boxPlotGeneDRMSummary(allBiomarker_model,nCol,dose,pltPath)

 id_t0 = strcmp(allBiomarker_model.pathName_significant_t0,pltPath);
 id_t20 = strcmp(allBiomarker_model.pathName_significant_t20,pltPath);
 id_t40 = strcmp(allBiomarker_model.pathName_significant_t40,pltPath);
 id_t60 = strcmp(allBiomarker_model.pathName_significant_t60,pltPath);
 id_t120 = strcmp(allBiomarker_model.pathName_significant_t120,pltPath);
 id_TELI = strcmp(allBiomarker_model.pathName_significant_TELI,pltPath);
 hold on;
bplot(allBiomarker_model.result_t0(id_t0,nCol),1,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
bplot(allBiomarker_model.result_t20(id_t20,nCol),2,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
bplot(allBiomarker_model.result_t40(id_t40,nCol),3,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
bplot(allBiomarker_model.result_t60(id_t60,nCol),4,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
bplot(allBiomarker_model.result_t120(id_t120,nCol),5,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
bplot(allBiomarker_model.result_TELI(id_TELI,nCol),6,'horiz','nooutliers','whiskers',2.5,'width',0.5,'linewidth',1.2);
xL = [floor(log10(min(dose))), ceil(log10(max(dose)))];
xlim(10.^(xL)); 
ylim([0 7]);
set(gca,'xscale','log','xtick',10.^(xL(1)-1:2:xL(2)),'ytick',1:1:6,...%set(gca,'ytick',1:1:6,...
    'fontname','arial','fontsize',10,'ylim',[0.5 6.5],...
    'yticklabel',{'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min','T = 120 min','TELI'});
%ylabel('Data Descriptor','fontname','Arial','fontsize',12);
%xlabel('Concentration (\muM)','fontname','Arial','fontsize',12);
title(strcat(pltPath),'fontname','Arial','fontsize',12,'fontweight', 'normal');
box on; hold off;

%

end




