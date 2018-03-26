function plotPathDRM_iInd_TELI (hAx, dose, pathVls,clr,shp, eqtype, varType, lnstts)
if nargin < 8, lnstts = 0; end
response =  permute(pathVls,[2 3 1]);
cf = fourParameterDRM(dose, response(:,1), eqtype);
xtest = (logspace(log10(dose(1,1)), log10(dose(end,1)),100))';
yPredict = cf(xtest);
if strcmp(varType, 'teli')
    plotfigures_teli (hAx, dose, response, xtest, yPredict, clr, shp);
else
    plotfigures_induction (hAx,dose, response, xtest, yPredict, clr, shp, lnstts);
end
end

function plotfigures_induction (hAx,dose, response, xtest, ytest, clr, shp, lnstts)
% a function to plot the dose response function with the 95% confidence
% interval;
% Created by: Sheikh M Rahman
% Date created: 9/23/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lnstts ~= 0    
    y = log(response);
    ytest = log(ytest);
else y = response; 
end

hold(hAx(1), 'on');
scatter(hAx(1),dose, y(:,1),30,  'Marker',shp{:},'markeredgecolor',clr,...
    'markerfacecolor',clr);
plot(hAx(1),xtest,ytest, 'color',clr, 'linewidth',1.4);

end

function plotfigures_teli (hAx, dose, response, xtest, ytest, clr, shp)
% a function to plot the dose response function with the 95% confidence
% interval;
% Created by: Sheikh M Rahman
% Date created: 9/23/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = response; 

hold(hAx(2), 'on');
scatter(hAx(2),dose, y(:,1),30,  'Marker',shp{:},'markeredgecolor',clr,...
    'markerfacecolor',clr);
plot(hAx(2),xtest,ytest, 'color',clr, 'linewidth',1.4);

end