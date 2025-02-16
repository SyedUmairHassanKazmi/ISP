function [fitresult, gof] = createFit1(Months, Po, Avg_Max_Fixed_Power, Avg_Max_Monthly_Power, Avg_Desired_Power)
%CREATEFIT1(MONTHS,PO)
%  Create a fit.
%
%  Data for 'Optimum Power' fit:
%      X Input : Months
%      Y Output: Po
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Jan-2021 14:16:32


%% Fit: 'Optimum Power'.
[xData, yData] = prepareCurveData( Months, Po );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999996940403878;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Optimum Power' );
h = plot( fitresult, xData, yData );
legend( h, 'Power vs. Months', 'Optimum Monthly Power', 'Location', 'NorthEast' );

% Label axes
xlabel Months
ylabel('Monthly Power Yield (KW-hr / m^2 / Month)')
xlim([1 12])
title('Monthly Optimum Power Yield')
grid on
hold on 
plot(Months,Avg_Max_Monthly_Power*ones(size(Months)), 'LineStyle','--' , 'DisplayName','Average Monthly Power', 'LineWidth' , 1.5)
plot(Months,Avg_Max_Fixed_Power*ones(size(Months)), 'LineStyle','-.' , 'DisplayName','Fixed Optimal Power','LineWidth' , 1.5)
plot(Months,Avg_Desired_Power*ones(size(Months)), 'LineStyle',':' , 'DisplayName','Fixed Desired Power', 'LineWidth' , 1.5)

end


