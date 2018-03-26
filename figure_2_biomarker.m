%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This script prepares the dose-response plot at biomarker level
%   The responses are plotted at five time-points and time-aggregated endpoint
%   Created by: Sheikh M. Rahman
%   Date: March, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load data/data_vF.mat; % load the data
lgnd = {'T = 0 min','T = 20 min', 'T = 40 min', 'T = 60 min',...
    'T = 120 min', 'Aggregated over 2-hr exposure'}; % the legend entries
id = [1, 5, 9, 13, 25];     %   relevant id of the time points ? 0, 20, 40, 60, and 120 minutes
ttl_prfx = {'a','b','c','d','e', 'f'};  % Figure sub-id
%% Plot figure 2, dose-dependent responses at time points and time-aggregated index
h2 = figure;
set(h2, 'Units','inches', 'Position',[0 0 7 9],'color','w',...
    'PaperPosition',[0 0 7 9],'PaperUnits', 'Inches', 'PaperSize', [7 9]);

for i = 1:6 % 6 for 6 panels, 5 static time-points and 1 for time-aggregated endpoint
    for j = 1:5 % 5 column; 1 for each pathways
        subplot('position',[0.065+(j-1)*.187, 0.817-(i-1)*.148, 0.18, 0.118])
        if i == 6
            %plot TELI of all genes with 95% confidence interval
            plotRawData_with_CI(conc,geneTELI(:,:,1),pathName, uniquePathName_repli{j}, 1);
        else
            % plot ln(I) of all genes with 95% confidence interval at 5 static time points
            plotInductionData_with_CI(conc,(permute(dataInd_mean(:,id(i),:),[1,3,2])),pathName, uniquePathName_repli{j}, 1);
        end
         % remove labels, ticklabels etc.
        if ~(j==1), ylabel('');  end
        if i ~= 6, set(gca,'xticklabel',''); end
        if j ~= 1, set(gca,'yticklabel','');  end
        if ~and(i == 6, j==3), xlabel(''); end
        if j ~= 3, title('');end
        % add the column title
        if i == 1
            annotation('textbox',[0.065+(j-1)*.187, 0.983, 0.18, 0.01],'str', [uniquePathName_repli{j}],...
                'linestyle','none','margin',0,'fontweight','bold',...
                'fontname','Arial','fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle')
            annotation('textbox',[0.065+(j-1)*.187, 0.964, 0.18, 0.01],'str', ['stress'],...
                'linestyle','none','margin',0,'fontweight','bold',...
                'fontname','Arial','fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle')

        end
    end
    % add the panel title
    annotation('textbox',[0 0.817-(i-1)*.148+0.119 1 0.023],...
        'str', strcat(ttl_prfx{i},{') '}, lgnd{i}),...
        'backgroundcolor',[0.9 0.9 0.9],'faceAlpha',0.5,'linestyle','none','margin',0,...
        'fontname','Arial','fontsize', 13,'HorizontalAlignment','center','VerticalAlignment','middle')

end
% add legend
ha = axes('position',[0.25, 0.00, 0.26 0.03]);
plot (conc,geneTELI(1,:,1),'k-','linewidth',1.2,'visible','off');
legend({'Individual Genes'},...
    'fontname', 'arial','FontSize',11,'box','off','location','southoutside');
hb = axes('position',[0.6, 0.00, 0.26 0.03]);
hfill=ciplot(geneTELI(1,:,2),geneTELI(1,:,3),conc,[0.7 0.7 0.7],'none',0.5);
hfill.Visible = 'off';
hl = legend({'95% Confidence Interval'},...
    'fontname', 'arial','FontSize',11,'box','off','location','southoutside');
ha.Visible = 'off';
hb.Visible = 'off';
% save the figure as pdf
print (h2, '-painters', '-dpdf', '-r300', 'figure2_biomarker.pdf');
