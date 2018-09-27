function [times, plotprod] = getPlotData(age, sf, productionRate, nuclide, scaling_model)
  %calculating the production rates/scaling factors through time
  %create a time vector
  numberpoints=300;
  times=linspace(0,-age,numberpoints); %find 300 points for the plot
  
  scaleFactorField = '';
  switch nuclide
    case 'al'
        scaleFactorField = 'Sel26';
    case 'be'
        scaleFactorField = 'Sel10';
    case 'c'
        scaleFactorField = 'Sel14';
    case 'he'
        scaleFactorField = 'Sel3';
    case 'ne'
        scaleFactorField = 'Sel21';
  end
  
  plotscalefactors=getcurrentsf(sf,times(1),scaling_model,nuclide);
  prod=productionRate/plotscalefactors.(scaleFactorField)(1);
  
  plotprod = zeros(1, numberpoints);
  for k=2:numberpoints;
      plotscalefactors=getcurrentsf(sf,times(k),scaling_model,nuclide);
      plotprod(k)=prod*plotscalefactors.(scaleFactorField);
  end
end