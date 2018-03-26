%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This script prepares the dose-response plot at pathway level
%   The responses are plotted at five time-points and time-aggregated endpoint
%   Created by: Sheikh M. Rahman
%   Date: March, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load data/data_vF.mat; % load the data
load data/timepoint_rawData.mat; % load individual time-point data
addpath lib;            % add required script library
lnstts = 1;             % 1 to plot log(I) and 0 to plot I
ylbl = 'ln({\itI} )';   % Set the y-axis label
yLim = [-0.5 1.5];      % set y-axis limit
ytcklbl = -0.5:0.5:1.5; % set the y-axis label
shp = {'v','*','+','^','x','o',}; % assign shape for the dose-responses of various time-points
clr = [1, 0.1, 0; 0, 0.4, 0; 1, 0.5, 0;
    1, 0, 1; 0.52, 0.52, 1; 0, 0, 1;]; % assign color for various time-points
ttl_prfx = {'a','b','c','d','e'};  % Figure sub-id
lgnd = {'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min',...
    'T = 120 min', 'Aggregated over 2-hr exposure'}; % the legend entries

%% Plot figure 3
h2 = figure;
set(h2, 'Units','inches', 'Position',[0 0 6.5 7],'color','w',...
    'PaperPosition',[0 0 6.5 7],'PaperUnits', 'Inches', 'PaperSize', [6.5 7]);
k = 1;
for i = 1:3
    for j = 1:2
        if k<6
            subplot('position',[0.085+(j-1)*.5, 0.74-(i-1)*.335, 0.35, 0.225]);
            hold on;
            % find the dose-response fit and plot the observation as point
            % and the model fit as line
            [hAx,~,~] = plotyy(conc, nan(numel(conc),1),conc, nan(numel(conc),1),'scatter','scatter');          % set two vertical axes
            plotPathDRM_iInd_TELI  (hAx, conc, path_t0(k,:,1),clr(1,:),shp(1), 'hill', 'induction', lnstts);    % T = 0 min 
            plotPathDRM_iInd_TELI  (hAx, conc, path_t20(k,:,1),clr(2,:),shp(2),'hill', 'induction', lnstts);    % T = 20 min 
            plotPathDRM_iInd_TELI  (hAx, conc, path_t40(k,:,1),clr(3,:),shp(3), 'hill', 'induction', lnstts);   % T = 40 min 
            plotPathDRM_iInd_TELI  (hAx, conc, path_t60(k,:,1),clr(4,:),shp(4),  'hill', 'induction', lnstts);  % T = 60 min 
            plotPathDRM_iInd_TELI  (hAx, conc, path_t120(k,:,1),clr(5,:),shp(5),'hill', 'induction', lnstts);   % T = 120 min 
            plotPathDRM_iInd_TELI (hAx, conc, pathTELI(k,:,1),clr(6,:),shp(6),'hill', 'teli');                  % aggregated time-point
            
            % set graph axis entries
            set(hAx(1), 'ylim',yLim,'ytick',ytcklbl, 'xscale', 'log', 'fontname','Arial', 'fontsize', 11, 'box','off');
            set(hAx(2),  'xscale', 'log', 'fontname','Arial', 'fontsize', 11);
            % left vertical axis for the static time-points
            ylabel(hAx(1), ylbl,'fontname','Arial','fontsize',12, 'fontweight','bold');
            xlabel(hAx(1), 'Concentration (\muM)','fontname','Arial','fontsize',12, 'fontweight','bold');
            xL = [floor(log10(min(conc))), ceil(log10(max(conc)))];
            xlim(10.^(xL));
            set(hAx(1),'xtick',10.^(xL(1):2:xL(2)),'xlim',10.^(xL), 'box','off');
            % right vertical axis for the aggregated time endpoint
            ylabel(hAx(2),'TELI','fontname','Arial','fontsize',12, 'fontweight','bold');
            set (hAx(2), 'xlim',10.^(xL),'xaxislocation','top','xticklabel',[],'ylim',[0 4],'ytick',1:1:4);
            hAx(2).XAxis.Visible = 'on';
            set(hAx(2),'ycolor','b')
            
            % sub-plot title
            annotation('textbox',[0.085+(j-1)*.5, 0.74-(i-1)*.335+0.227, 0.35, 0.032],...
                'str', strcat(ttl_prfx{k},{') '}, pathUnique_t0{k},' stress'),...
                'backgroundcolor',[0.9 0.9 0.9],'faceAlpha',1,'linestyle','none','margin',0,...
                'fontname','Arial','fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle');
          %  box on;
            k = k+1;
        end
    end
end
hold off;

% set the legend
ha = axes('position',[0.73, 0.08, 0.26 0.32]);
hold on;
for i = 1:6
    hs = scatter(conc, conc, 30,  'Marker',shp{i},'markeredgecolor',clr(i,:),...
        'markerfacecolor',clr(i,:)); hs.Visible = 'off';
end
hold off;
legend(lgnd,...
    'fontname', 'arial','FontSize',12,'box','off','location','South');
ha.Visible = 'off';

% save the figure as pdf
print (h2,  '-dpdf', '-r300', 'results/figure3_pathway.pdf');
