function plotPathDRMSummary(allpath_model,nCol,dose,pltPath)

id_t0 = strcmp(allpath_model.uniquePathName,pltPath);
hold on;
scatter(allpath_model.result_t0(id_t0,nCol),1,'or','filled','linewidth',1.2);
scatter(allpath_model.result_t20(id_t0,nCol),2,'or','filled','linewidth',1.2);
scatter(allpath_model.result_t40(id_t0,nCol),3,'or','filled','linewidth',1.2);
scatter(allpath_model.result_t60(id_t0,nCol),4,'or','filled','linewidth',1.2);
scatter(allpath_model.result_t120(id_t0,nCol),5,'or','filled','linewidth',1.2);
scatter(allpath_model.result_TELI(id_t0,nCol),6,'or','filled','linewidth',1.2);
xL = [floor(log10(min(dose))), ceil(log10(max(dose)))];
xlim(10.^(xL)); 
ylim([0 7]);
set(gca,'xscale','log','xtick',10.^(xL(1)-1:2:xL(2)),'ytick',1:1:6,...%set(gca,'ytick',1:1:6,...
    'fontname','arial','fontsize',10,'ylim',[0.5 6.5],...
    'yticklabel',{'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min','T = 120 min','TELI'});
%ylabel('Data Descriptor','fontname','Arial','fontsize',12);
xlabel('Concentration (\muM)','fontname','Arial','fontsize',12);
title(strcat(pltPath),'fontname','Arial','fontsize',12,'fontweight', 'normal');
box on; hold off;

%

end




