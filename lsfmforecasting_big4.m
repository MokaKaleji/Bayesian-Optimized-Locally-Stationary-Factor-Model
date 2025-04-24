%% lsfm_bo_forecasting.m
% Author: Moka Kaleji
% Affiliation: Master Thesis in Econometrics, University of Bologna
% Description:
%   Executes forecasting for macroeconomic indicators (GDP, Unemployment,
%   Inflation, Interest Rate) using LSFM outputs optimized via Bayesian
%   estimation. Computes performance metrics (MSFE, RMSE) and visualizes
%   forecast accuracy and error distributions.

clear; close all; clc;

%% Dataset Selection
options = {'Monthly (MD1959.xlsx)', 'Quarterly (QD1959.xlsx)'};
[choiceIndex, ok] = listdlg('PromptString','Select the dataset to load:',...
                             'SelectionMode','single',...
                             'ListString',options,...
                             'Name','Dataset Selection',...
                             'ListSize',[400,200]);
if ~ok
    error('Dataset selection cancelled. Exiting.');
end

%% Load Data and Key Variables
switch choiceIndex
    case 1  % Monthly frequency
        filepath = '/Users/moka/Research/Thesis/Live Project/Processed_Data/MD1959.xlsx';
        raw = readtable(filepath); x = table2array(raw);
        T = size(x,1);
        key_vars = [1,24,105,77];
    case 2  % Quarterly frequency
        filepath = '/Users/moka/Research/Thesis/Live Project/Processed_Data/QD1959.xlsx';
        raw = readtable(filepath); x = table2array(raw(:,2:end));
        T = size(x,1);
        key_vars = [1,58,116,136];
    otherwise
        error('Invalid dataset choice.');
end
var_names = {'GDP','Unemployment','Inflation','Interest Rate'};

%% Load Optimized Estimation Results
% Contains: q_opt, h_opt, p_opt, Fhat_train, Lhat_train, mean_train, std_train, T_train, N, key_vars
load('lsfm_estimation_results.mat');

%% Forecast Settings
H = 8;  % Forecast horizon
fprintf('Using optimized parameters: q=%d, h=%.3f, p=%d\n', q_opt, h_opt, p_opt);

%% Prepare Test Sample
x_test = x(T_train+1:T_train+H, :);
x_test_norm = (x_test - mean_train) ./ std_train;

%% Forecast via LSFM + VAR
yhat_norm = forecast_LSFM(Fhat_opt, Lhat_opt, T_train, N_var, H, p_opt);
yhat = yhat_norm .* std_train + mean_train;

%% Compute Forecast Error Metrics
sq_err = (x_test_norm(:,key_vars) - yhat_norm(:,key_vars)).^2;
MSFE_horizon = mean(sq_err,2);         % MSFE by horizon across key vars
MSFE_vars    = mean(sq_err,1);         % MSFE by variable across horizons
MSFE_overall = mean(MSFE_vars);        % Overall MSFE

fprintf('Overall MSFE (key vars): %.6f\n', MSFE_overall);
for idx = 1:length(key_vars)
    fprintf('MSFE %s: %.6f\n', var_names{idx}, MSFE_vars(idx));
end

%% Visualization of Forecasting Results
figure('Position', [100, 100, 1200, 800]);

% MSFE over forecast horizon
subplot(3,3,1);
plot(1:H, MSFE_horizon, '-o', 'LineWidth', 1.5);
title('MSFE by Horizon'); xlabel('Horizon'); ylabel('MSFE'); grid on;

% Actual vs Forecast for each indicator
for idx=1:4
    subplot(3,3,1+idx);
    plot(1:T_train, (x(1:T_train,key_vars(idx))), 'b-', 'DisplayName','Train Actual'); hold on;
    plot(T_train+1:T_train+H, x_test(:,key_vars(idx)), 'k-', 'DisplayName','Test Actual');
    plot(T_train+1:T_train+H, yhat(:,key_vars(idx)), 'r--','DisplayName','Forecast');
    hold off; title(var_names{idx}); xlabel('Time'); ylabel('Value'); grid on;
    if idx==1, legend('Location','Best'); end
end
% Error distribution histogram
subplot(3,3,6);
histogram(MSFE_vars,10);
title('MSFE Distribution Across Indicators'); xlabel('MSFE'); ylabel('Frequency'); grid on;

% Time series of factors
subplot(3,3,7);
plot(Fhat_opt, 'LineWidth', 1.5);
title('Estimated Factors over Time'); xlabel('Time'); ylabel('Value');
legend(arrayfun(@(k)['Factor ' num2str(k)], 1:size(Fhat_opt,2), 'UniformOutput', false), 'Location', 'Best'); grid on;

sgtitle(sprintf('LSFM Bayesian-Optimized Forecasting: q=%d, h=%.3f, p=%d', q_opt,h_opt,p_opt));

fprintf('Forecasting complete. Results saved to %s\n', outname);

%% forecast_LSFM (Helper Function)
function yhat = forecast_LSFM(Fhat, Lhat, T, N, H, p)
    q = size(Fhat,2);
    if p >= T, p = max(1,floor(T/2)); end
    mdl = varm(q,p); Est = estimate(mdl,Fhat);
    Ff = forecast(Est,H,Fhat);
    Lambda = squeeze(Lhat(:,:,T));
    yhat = zeros(H,N);
    for h1=1:H, yhat(h1,:) = (Lambda*Ff(h1,:)')'; end
end
