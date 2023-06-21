function [fitresult, gof, coefsFit, fitAimantation] = FitExpHystn4(interpScaledHext, interpScaledAimant)
%FitExpHyst(interpScaledHext,interpScaledAimant)
%  Create a fit.
%
%  Data for 'FitExpHystn3' fit:
%      X Input : interpScaledHext
%      Y Output: interpScaledAimant
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%


%% Fit: 'FitExpHystn4'.
[xData, yData] = prepareCurveData( interpScaledHext, interpScaledAimant );

% Set up fittype and options.
ft = fittype( 'a*atan((x+aa)/aaa) + b*atan((x+bb)/bbb) + c*atan((x+cc)/ccc) + d*atan((x+dd)/ddd)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0 0 0 0 0 0 0];
opts.StartPoint = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

coefsFit = coeffvalues(fitresult);
fitAimantation = coefsFit(1)*atan((interpScaledHext+coefsFit(2))/coefsFit(3)) + coefsFit(4)*atan((interpScaledHext+coefsFit(5))/coefsFit(6)) + coefsFit(7)*atan((interpScaledHext+coefsFit(8))/coefsFit(9)) + coefsFit(10)*atan((interpScaledHext+coefsFit(11))/coefsFit(12));



