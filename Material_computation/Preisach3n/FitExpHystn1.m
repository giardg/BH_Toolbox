function [fitresult, gof, coefsFit, fitAimantation] = FitExpHystn1(interpScaledHext, interpScaledAimant)
%FitExpHyst(interpScaledHext,interpScaledAimant)
%  Create a fit.
%
%  Data for 'FitExpHyst' fit:
%      X Input : interpScaledHext
%      Y Output: interpScaledAimant
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%


%% Fit: 'FitExpHyst'.
[xData, yData] = prepareCurveData( interpScaledHext, interpScaledAimant );

% Set up fittype and options.
ft = fittype( 'a*atan((x+aa)/aaa)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [0.5 0.5 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

coefsFit = coeffvalues(fitresult);
fitAimantation = coefsFit(1)*atan((interpScaledHext+coefsFit(2))/coefsFit(3));



