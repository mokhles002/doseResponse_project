%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This script prepares the dose-response plot at cellular level
%   The responses are plotted at five time-points and time-aggregated endpoint
%   Created by: Sheikh M. Rahman
%   Date: March, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load required data
load data/data_vF.mat;
load data/timepoint_rawData.mat;
load data/network_data_dose_response.mat;
load data/parafacLoads.mat;
addpath lib;            % add required script library
lnstts = 1;             % 1 to plot log(I) and 0 to plot I
ylbl = 'ln({\itI} )';   % Set the y-axis label
yLim = [-0.5 1.0];      % set y-axis limit
ytcklbl = -1:0.5:1; % set the y-axis label
%lnstts = 0; ylbl ='Induction Factor ({\itI} )'; ytcklbl = 0:1:4; yLim = [0 4];
shp = {'v','*','+','^','x','o',};   % assign shape for the dose-responses of various time-points
clr = [1, 0.1, 0; 0, 0.4, 0; 1, 0.5, 0;
    1, 0, 1; 0.52, 0.52, 1; 0, 0, 1;];    % assign color for various time-points

%% Plot figure 4
h2 = figure;
set(h2, 'Units','inches', 'Position',[0 0 6.5 7],'color','w',...
    'PaperPosition',[0 0 6.5 7],'PaperUnits', 'Inches', 'PaperSize', [6.5 7]);

% TOXICOGENOMIC QUANTIFIERS
k = 6;  % index for the cellular level response
subplot('position',[0.09, 0.75, 0.44, 0.205]);
hold on;
% find the dose-response fit and plot the observation as point and the model fit as line
[hAx,~,~] = plotyy(conc, nan(numel(conc),1),conc, nan(numel(conc),1),'scatter','scatter');
plotPathDRM_iInd_TELI  (hAx, conc, path_t0(k,:,1),clr(1,:),shp(1), 'hill', 'induction', lnstts);    % t = 0 min
plotPathDRM_iInd_TELI  (hAx, conc, path_t20(k,:,1),clr(2,:),shp(2),'hill', 'induction', lnstts);    % t = 20 min
plotPathDRM_iInd_TELI  (hAx, conc, path_t40(k,:,1),clr(3,:),shp(3), 'hill', 'induction', lnstts);   % t = 40 min
plotPathDRM_iInd_TELI  (hAx, conc, path_t60(k,:,1),clr(4,:),shp(4),  'hill', 'induction', lnstts);  % t = 60 min
plotPathDRM_iInd_TELI  (hAx, conc, path_t120(k,:,1),clr(5,:),shp(5),'hill', 'induction', lnstts);   % t = 120 min
plotPathDRM_iInd_TELI (hAx, conc, pathTELI(k,:,1),clr(6,:),shp(6),'hill', 'teli');                  % aggregated time-point

% set graph axis entries
set(hAx(1), 'ylim',yLim, 'xscale', 'log', 'fontname','Arial', 'fontsize', 11, 'box','off');
set(hAx(2), 'ylim',[0 3.5], 'xscale', 'log', 'fontname','Arial', 'fontsize', 11);

% left vertical axis for the static time-points
ylabel(hAx(1), ylbl,'fontname','Arial','fontsize',12, 'fontweight','bold');
xlabel(hAx(1), 'Concentration (\muM)','fontname','Arial','fontsize',12, 'fontweight','bold');
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
set(hAx(1),'xtick',10.^(xL(1):2:xL(2)),'xlim',10.^(xL),'ytick',ytcklbl,'TickLength',[0.02, 0.02]);
% right vertical axis for the aggregated time endpoint
ylabel(hAx(2),'TELI','fontname','Arial','fontsize',12, 'fontweight','bold');
set (hAx(2), 'xlim',10.^(xL),'xaxislocation','top','xticklabel',[],'ytick',1:1:4,'TickLength',[0.02, 0.02]);
hAx(2).XAxis.Visible = 'on';
set(hAx(2),'ycolor','b')
hold off;

% plot legend
ha = axes('position',[0.74, 0.745, 0.25 0.32]);
hold on;
for i = 1:6
    hs = scatter(conc, conc, 30,  'Marker',shp{i},'markeredgecolor',clr(i,:),...
        'markerfacecolor',clr(i,:)); hs.Visible = 'off';
end
hold off;
lgnd = {'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min',...
    'T = 120 min', 'Aggregated over 2-hr exposure'};
legend(lgnd,...
    'fontname', 'arial','FontSize',10,'box','off','location','Southeast');
ha.Visible = 'off';
% title of the plot
annotation('textbox',[0.00 0.97 1 0.03],'str', 'a) Toxicogenomic Quantifier',...
    'backgroundcolor',[0.9 0.9 0.9],'faceAlpha',0.5,'linestyle','none','margin',0,...
    'fontname','Arial','fontsize', 13,'HorizontalAlignment','center','VerticalAlignment','middle')

% NETWORK PROPERTIES
subplot('position',[0.107, 0.385, 0.295, 0.22]);
[hAx,~,~] = plotyy(conc, nan(numel(conc),1),conc, nan(numel(conc),1),'scatter','scatter');
plotDRM_networkProperties (hAx(2), conc, allNetworkProperties(:,1),clr(6,:),shp(6), 'hill')
plotDRM_networkProperties (hAx(1), conc, allNetworkProperties(:,3),clr(3,:),shp(3), 'hill')
plotDRM_networkProperties (hAx(1), conc, allNetworkProperties(:,4),clr(1,:),shp(1), 'hill')
set(hAx(1), 'ylim',[0 0.3], 'xscale', 'log', 'fontname','Arial', 'fontsize', 11, 'box','off');
set(hAx(2), 'ylim',[0 1500], 'xscale', 'log', 'fontname','Arial', 'fontsize', 11);

ylabel(hAx(1), [propertyName(3); propertyName(4)],'fontname','Arial','fontsize',12, 'fontweight','bold');
xlabel(hAx(1), 'Concentration (\muM)','fontname','Arial','fontsize',12, 'fontweight','bold');
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
set(hAx(1),'xtick',10.^(xL(1):2:xL(2)),'xlim',10.^(xL),'ytick',0:0.1:0.3,'TickLength',[0.024, 0.02]);

ylabel(hAx(2),propertyName(1),'fontname','Arial','fontsize',12, 'fontweight','bold');
set (hAx(2), 'xlim',10.^(xL),'xaxislocation','top','xticklabel',[],'ytick',0:500:1500,'TickLength',[0.024, 0.02]);
hAx(2).XAxis.Visible = 'on';
set(hAx(2),'ycolor',clr(6,:))
hold off;

ha = axes('position',[0.12, 0.365, 0.3, 0.25]);
hold on;
for i = [6, 1, 3]
    hs = scatter(conc, conc, 30,  'Marker',shp{i},'markeredgecolor',clr(i,:),...
        'markerfacecolor',clr(i,:)); hs.Visible = 'off';
end
hold off;
legend(propertyName([1,4,3]),...
    'fontname', 'arial','FontSize',10,'box','off','location','Southeast');
ha.Visible = 'off';

% Second plot of network properties
subplot('position',[0.62, 0.385, 0.3, 0.22]);
[hAx,~,~] = plotyy(conc, nan(numel(conc),1),conc, nan(numel(conc),1),'scatter','scatter');
plotDRM_networkProperties (hAx(2), conc, allNetworkProperties(:,5),clr(2,:),shp(2), 'hill')
plotDRM_networkProperties (hAx(1), conc, allNetworkProperties(:,2),clr(4,:),shp(4), 'hill')
plotDRM_networkProperties (hAx(1), conc, allNetworkProperties(:,6),clr(5,:),shp(5), 'hill')
set(hAx(1), 'ylim',[0.12 0.85], 'xscale', 'log', 'fontname','Arial', 'fontsize', 11, 'box','off');
set(hAx(2), 'ylim',[1.95 2.55], 'xscale', 'log', 'fontname','Arial', 'fontsize', 11);

ylabel(hAx(1), [propertyName(2); propertyName(6)],'fontname','Arial','fontsize',12, 'fontweight','bold');
xlabel(hAx(1), 'Concentration (\muM)','fontname','Arial','fontsize',12, 'fontweight','bold');
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
set(hAx(1),'xtick',10.^(xL(1):2:xL(2)),'xlim',10.^(xL),'ytick',0.2:0.2:0.8,'TickLength',[0.024, 0.02]);

ylAx = ylabel(hAx(2),propertyName(5),'fontname','Arial','fontsize',12, 'fontweight','bold');
ylAx.Units = 'normalized'; 
ylAx.Position = [1.16 .435 0];
set (hAx(2), 'xlim',10.^(xL),'xaxislocation','top','xticklabel',[],'ytick',2:0.1:2.5,'TickLength',[0.024, 0.02]);
hAx(2).XAxis.Visible = 'on';
set(hAx(2),'ycolor',clr(2,:))
hold off;

ha = axes('position',[0.74, 0.365, 0.2, 0.25]);
hold on;
for i = [2, 4, 5]
    hs = scatter(conc, conc, 30,  'Marker',shp{i},'markeredgecolor',clr(i,:),...
        'markerfacecolor',clr(i,:)); hs.Visible = 'off';
end
hold off;
legend(['Char. Path Length';propertyName([2, 6])],...
    'fontname', 'arial','FontSize',10,'box','off','location','SouthEast');
ha.Visible = 'off';
% title of the plot
annotation('textbox',[0.00 0.645 1 0.03],'str', 'b) Gene Co-expression Network Properties',...
    'backgroundcolor',[0.9 0.9 0.9],'faceAlpha',0.5,'linestyle','none','margin',0,...
    'fontname','Arial','fontsize', 13,'HorizontalAlignment','center','VerticalAlignment','middle')

% PARAFAC MODEL LOADINGS
hAx = axes ('position',[0.075, 0.068, 0.5, 0.2]);
hold on;
plotDRM_networkProperties (hAx, conc, prfcLoadsAve(:,1),clr(1,:),shp(1), 'hill');
plotDRM_networkProperties (hAx, conc, prfcLoadsAve(:,2),clr(2,:),shp(2), 'hill');

ylabel('Factor Loading','fontname','Arial','fontsize',12, 'fontweight','bold');
xlabel( 'Concentration (\muM)','fontname','Arial','fontsize',12, 'fontweight','bold');
xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
xlim(10.^(xL));
set(gca,'xtick',10.^(xL(1):2:xL(2)),'xlim',10.^(xL),'TickLength',[0.016, 0.02],...
    'ytick',0:0.2:0.8, 'xscale','log','ylim',[0 0.8], 'box','on');
% title of the plot
annotation('textbox',[0.00 0.28 1 0.03],'str', 'c) PARAFAC Model',...
    'backgroundcolor',[0.9 0.9 0.9],'faceAlpha',0.5,'linestyle','none','margin',0,...
    'fontname','Arial','fontsize', 13,'HorizontalAlignment','center','VerticalAlignment','middle')

% Legend
ha = axes('position',[0.68, 0.15, 0.1, 0.2]);
hold on;
for i = 1:2
    hs = scatter(conc, conc, 30,  'Marker',shp{i},'markeredgecolor',clr(i,:),...
        'markerfacecolor',clr(i,:)); hs.Visible = 'off';
end
hold off;
hl=legend({'Factor 1', 'Factor 2'},...
    'fontname','Arial','fontsize',10,'box', 'off','location','southeast');
ha.Visible = 'off';

% save the figure as PDF
print (h2,  '-dpdf', '-r300', 'results/figure4_molecular.pdf');
