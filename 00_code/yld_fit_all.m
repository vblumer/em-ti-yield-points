%% ----------------------------- DESCRIPTION ----------------------------- &&
% Fit yield functions to datapoints computed in DAMASK
%  ----------------------------- END DESCRIPTION -------------------------  &
format compact


%% ----------------------------- READ DATA POINTS ------------------------ &&
pth               = '../99_results/results_2023-09-22_16-23-27/';
yield_surf_points = 1e-6*readmatrix([pth,'yld_results.txt']);
% -------------------------------------------------------------------------
yield_surf_points = [yield_surf_points;-yield_surf_points(2:end-1,:)];
n                 = length(yield_surf_points);
target_values     = zeros(n, 1);
%  ----------------------------- END READ DATA POINTS -------------------- &&

%% ----------------------------- MISES YIELD FUNCTION -------------------- &&
% define the implicit curve function
implicit_curve_mises = @(param, xy) param*(1/2*((xy(:,1)-xy(:,2)).^2+xy(:,1).^2+xy(:,2).^2))-1;

% initial guess for parameters
initial_param_mises = 10;

% fit the curve using lsqcurvefit
fit_param_mises     = lsqcurvefit(implicit_curve_mises,initial_param_mises, yield_surf_points, target_values);
% -------------------------------------------------------------------------
% error computations - root mean square error
model_predictions = implicit_curve_mises(fit_param_mises, yield_surf_points);
residuals         = target_values - model_predictions;
rmse_mises        = sqrt(mean(residuals.^2));
% -------------------------------------------------------------------------
% display the fitted parameters
% disp('Fitted Parameters, root mean square error:');
% fit_param_mises
% rmse_mises
%  ----------------------------- END MISES YIELD FUNCTION ---------------- &&

%% ----------------------------- HILL48 YIELD FUNCTION ------------------- &&
% define the implicit curve function
implicit_curve_hill48 = @(params, xy) (xy(:,1).^2*(params(1)+params(2))-2*params(2)*xy(:,1).*xy(:,2)+xy(:,2).^2*(params(2)+params(3)))-1;
% params(1) - G
% params(2) - H
% params(3) - F
% initial guess for parameters
initial_param_hill48 = [0.1,0.2,0.3];

% fit the curve using lsqcurvefit
fit_param_hill48 = lsqcurvefit(implicit_curve_hill48, initial_param_hill48, yield_surf_points, target_values);
% -------------------------------------------------------------------------
% error computations - root mean square error
model_predictions = implicit_curve_hill48(fit_param_hill48, yield_surf_points);
residuals = target_values - model_predictions;
rmse_hill48 = sqrt(mean(residuals.^2));
% -------------------------------------------------------------------------
% % display the fitted parameters 
% fprintf('Fitted Parameters, root mean square error:\n');
% fit_param_hill48
% rmse_hill48
% -----------------------------  END HILL48 YIELD FUNCTION --------------- &&

%% ----------------------------- PLOT------------------------------------- &&
figure;
% -------------------------------------------------------------------------
% Generate a grid of x and y values
x = linspace(-1750, 1750, 1000); % Adjust the range as needed
y = linspace(-1750, 1750, 1000); % Adjust the range as needed
[X, Y] = meshgrid(x, y);
% -------------------------------------------------------------------------
% MISES - Evaluate the implicit curve function for Mises
Z = implicit_curve_mises(fit_param_mises, [X(:), Y(:)]);
Z = reshape(Z, size(X));
% Create a contour plot of the implicit curve
contour(X, Y, Z, [0, 0],"Edgecolor","#0073ec", 'LineWidth', 4);
hold on;
% -------------------------------------------------------------------------
% HILL48 - Evaluate the implicit curve function for Hill48
Z = implicit_curve_hill48(fit_param_hill48, [X(:), Y(:)]);
Z = reshape(Z, size(X));
% Create a contour plot of the implicit curve
contour(X, Y, Z, [0, 0],"Edgecolor","#EC7900", 'LineWidth', 4);
hold on;
% -------------------------------------------------------------------------
scatter(yield_surf_points(:,1), yield_surf_points(:,2), 70,'k',  'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
% -------------------------------------------------------------------------
xlabel('$\sigma_{11}\,$[MPa]', 'Interpreter', 'latex', 'FontSize', 16); % Enable LaTeX rendering and set font size
ylabel('$\sigma_{22}\,$[MPa]', 'Interpreter', 'latex', 'FontSize', 16);
axis equal;
% title('Implicit Curve: param - sqrt(1/2 * ((x-y)^2 + x^2 + y^2)) = 0');
% -------------------------------------------------------------------------
% Set font size for axis numbers (tick labels)
set(gca, 'FontSize', 30,'FontName','Times');
% -------------------------------------------------------------------------
xlim([-1600, 1600]); % Adjust the limits as needed for your specific case
ylim([-1600, 1600]);
% -------------------------------------------------------------------------
grid on;
% -------------------------------------------------------------------------
legend1=legend('von Mises', 'Hill48');
set(legend1,     'Position',[0.226 0.665 0.239 0.264]);
legend('boxoff')
hold off;
% -------------------------------------------------------------------------
% Save the plot as a PNG image
% saveas(gcf, [pth,'yld_surf_plot.png']);
% -------------------------------------------------------------------------
%  ----------------------------- END PLOT------------------------------------- &&
