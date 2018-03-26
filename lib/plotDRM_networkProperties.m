function plotDRM_networkProperties (hAx, dose, response,clr,shp, eqtype)
cf = fourParameterDRM(dose, response, eqtype);
xtest = (logspace(log10(dose(1,1)), log10(dose(end,1)),100))';
yPredict = cf(xtest);

hold(hAx, 'on');
scatter(hAx,dose, response(:,1),30,  'Marker',shp{:},'markeredgecolor',clr,...
    'markerfacecolor',clr);
plot(hAx,xtest,yPredict, 'color',clr, 'linewidth',1.4);

end