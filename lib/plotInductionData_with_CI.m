function[yP, LL, UL] =  plotInductionData_with_CI(dose,response,pathName, pltPath,status, clr)
% status defines the color of the plot.
ind = ismember(pathName, pltPath);
yP = response(ind,:);
%yP(yP<0.1) = 0.1; yP(yP>10) = 10;
yP(yP < log(0.1)) = log(0.1); yP(yP>log(10)) = log(10);
%ave_data = mean(yP,1);
% find the lower and upper limit of the 95% confidence interval from the
LL = prctile(yP,2.5,1);
UL = prctile(yP,97.5,1);

hold on; box on;
if status == 0;
    ciplot(LL,UL,dose,clr,'none',0.5);
   % plot (dose,ave_data,'linewidth',2,'color',clr);
    plot (dose,yP,'linestyle',':','linewidth',2 ,'color',clr);
else
    ciplot(LL,UL,dose,[0.7 0.7 0.7],'none',0.5);
    % plot (dose,ave_data,'linewidth',2,'color','k');
    plot (dose,yP,'linestyle','-','linewidth',0.8);
end
    set(gca,'xscale','log','fontname','Arial','fontsize',10,'ticklength',[0.03 0.035]);
    ylabel('ln({\itI} )','fontname','Arial','fontsize',12,'fontweight','Bold');
    xlabel('Concentration (\muM)','fontname','Arial','fontsize',12,'fontweight','Bold');
    xl = [floor(log10(min(dose))), ceil(log10(max(dose)))];
    %yLmax = ceil(max(max(response(:,:))));
    %yLmin = floor(min(min(response(:,:))));
    xlim(10.^(xl));% ylim([-0.5, 10.5]); 
    ylim([-2.5, 2.5]);
     set(gca,'xtick',10.^((xl(1)-1):2:(xl(2)-1)),'ytick',-2:1:2);

      % set(gca,'xtick',10.^((xl(1)-1):2:(xl(2)-1)),'ytick',0:5:10);
% title(strcat (pltPath,' stress'),'fontname','Arial','fontsize',13,'fontweight','normal');
end

%{
Induction
 




for i = 1:5
    h2 = figure;
    set(h2, 'Units','inches', 'Position',[0 0 6 8],...
    'PaperPosition',[0 0 6 8],'PaperUnits', 'Inches', 'PaperSize', [6 8]);
    plotPathDoseResponse (geneTELI, uniquePathName_repli{i}, pathName, conc, geneName);
    xlim([1E-2, 1E4]);
    print (h2, '-dtiff', '-r300', strcat('time/geneDose_',uniquePathName_repli{i},'.tiff'));
end

%}
