%% plot the summary of all level of dose response models
addpath lib/
load results/allPath_modelSummary.mat;
load results/allBiomarker_modelSummary.mat;
load results/allCellularLevelResults.mat;
lgnd = {'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min','T = 120 min'};
txtTtle = {'a) Biomarker Level','b) Pathway Level', 'c) Cellular Level'};
h2 = figure;
set(h2, 'Units','inches', 'Position',[0 0 7 8],...
    'PaperPosition',[0 0 7 8],'PaperUnits', 'Inches', 'PaperSize', [7 8],'color','w');
for i = 1:2
     annotation('textbox',[0.00 0.733-(i-1)*.29+0.238 1 0.029],'str', txtTtle(i),...
       'backgroundcolor',[0.9 0.9 0.9],'linestyle','none','margin',0,...
        'fontname','Arial','fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle')
    
    for j = 1:5
        subplot('position',[0.12+(j-1)*.177,  0.775-(i-1)*.29, 0.165, 0.168])
        if i == 1
            boxPlotGeneDRMSummary(allBiomarker_model,12,conc,uniquePathName_repli(j));
            if j==3, xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold'); end            
        end
        if i == 2
            plotPathDRMSummary(allpath_model,12,conc,uniquePathName_repli(j));
            xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold');
        end
        if ~(j==1), ylabel(''); set(gca,'yticklabel',''); end
        if ~(j==3), xlabel(''); end
        set(gca,'ticklength',[0.03, 0.03]);
    end
end
i = 3; j = 1;
annotation('textbox',[0.00 0.731-(i-1)*.29+0.238 1 0.029],'str', txtTtle(i),...
       'backgroundcolor',[0.9 0.9 0.9],'linestyle','none','margin',0,...
        'fontname','Arial','fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle');
subplot('position',[0.12+(j-1)*.177, 0.77-(i-1)*.29, 0.165, 0.168])
hold on; nCol = 12;
scatter(allpath_model.result_t0(6,nCol),1,'bv','filled','linewidth',1.2);
scatter(allpath_model.result_t20(6,nCol),2,'bv','filled','linewidth',1.2);
scatter(allpath_model.result_t40(6,nCol),3,'bv','filled','linewidth',1.2);
scatter(allpath_model.result_t60(6,nCol),4,'bv','filled','linewidth',1.2);
scatter(allpath_model.result_t120(6,nCol),5,'bv','filled','linewidth',1.2);
scatter(allpath_model.result_TELI(6,nCol),6,'bv','filled','linewidth',1.2);
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
ylim([0 7]);
set(gca,'xscale','log','xtick',10.^(xL(1)-1:2:xL(2)),'ytick',1:1:6,...%set(gca,'ytick',1:1:6,...
    'fontname','arial','fontsize',10,'ylim',[0.5 6.5],'ticklength',[0.03, 0.03],...
    'yticklabel',{'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min','T = 120 min','TELI'});
title('Toxicogenomic quantifier','fontname','Arial','fontsize',14,'fontweight','normal');
xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold');
box on; hold off;

i = 3; j = 3;id = [1:6];
subplot('position',[0.22+(j-1)*.17, 0.77-(i-1)*.29, 0.165, 0.168])
hold on; nCol = 12;
scatter(result_summary_network(id,nCol),1:6,'s',...
            'markeredgecolor','none','markerfacecolor',[1 0.5 0.0]);
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
ylim([0 7]);
set(gca,'xscale','log','xtick',10.^(xL(1)-1:2:xL(2)),'ytick',1:1:6,...%set(gca,'ytick',1:1:6,...
    'fontname','arial','fontsize',10,'ylim',[0.5 6.5],'ticklength',[0.03, 0.03],...
    'yticklabel',propertyName_network(id));
xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold');
title({'Network Properties'},'fontname','Arial','fontsize',14,'fontweight','normal');
box on; hold off;

i = 3; j = 5; id = [1, 3:7];
subplot('position',[0.15+(j-1)*.17, 0.77-(i-1)*.29, 0.165, 0.168])
hold on; nCol = 12;
scatter(result_summary_PARAFAC(:,nCol),1:2,'<m','filled','linewidth',1.2);
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
ylim([0 3]);
set(gca,'xscale','log','xtick',10.^(xL(1)-1:2:xL(2)),'ytick',1:1:2,...
    'fontname','arial','fontsize',10,'ticklength',[0.03, 0.03],...
    'yticklabel',{'Factor 1','Factor 2'});
httl = title('PARAFAC Model','fontname','Arial','fontsize',14,'fontweight','normal');
set(httl,'units','normalized');
set(httl,'Position', httl.Position - [0.15 0 0]);
xlblPos = xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold');
set(xlblPos,'units','normalized');
set(xlblPos,'Position', xlblPos.Position - [0.15 0 0]);
box on; hold off;

%plot legend
annotation('textbox',[0.12 0.097 0.3 0.029],'str', 'Biomarker Level',...%    'backgroundcolor',[0.9 0.9 0.9],
        'linestyle','none','margin',0,...
        'fontname','Arial','fontsize', 13,'HorizontalAlignment','center','VerticalAlignment','middle');
    
ha = axes('position',[0.12, 0.003, 0.8 0.12]);
hold on;
hs = plot(conc,conc,'d','color',[0 0.6 0],'linewidth',1.2,'markersize',6); hs.Visible = 'off';
hs = plot(conc, conc,'-','linewidth',1.2,'color','m'); hs.Visible = 'off';
hfill=bar(conc,'edgecolor',[0.0005    0.3593    0.7380],...
    'linewidth',1.2,'facecolor','none'); hfill.Visible = 'off';
hs = plot(conc, conc,'--','linewidth',1.2,'color',[0.3 0.3 0.3]); hs.Visible = 'off';

hold off;
legend({'Mean','Median', '25-75 percentile', '95% Confidence Interval'},...
    'fontname', 'arial','FontSize',12,'box','off','location','NorthWest');
ha.Visible = 'off';

hb = axes('position',[0.45, 0.003, 0.5 0.12]);
hold on;
hs = scatter(conc,conc,'or','filled','linewidth',1.2); hs.Visible = 'off';
hs = scatter(conc, conc,'bv','filled','linewidth',1.2); hs.Visible = 'off';
hs = scatter(conc,conc,'s','markerfacecolor',[1 0.5 0.0],'markeredgecolor','none'); hs.Visible = 'off';
hs = scatter(conc, conc,'<m','filled','linewidth',1.2); hs.Visible = 'off';

hold off;
legend({'Pathway level endpoint','Cellular level endpoint', 'Gene co-expression network properties', 'PARAFAC model properties'},...
    'fontname', 'arial','FontSize',12,'box','off','location','NorthWest');
hb.Visible = 'off';

print (h2,'-dpdf','-r300','results/figure5_allBMD.pdf');
