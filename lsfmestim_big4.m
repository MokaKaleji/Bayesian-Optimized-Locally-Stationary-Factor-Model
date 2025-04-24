%% lsfm_bo_estimation.m
% Author: Moka Kaleji • Contact: mohammadkaleji1998@gmail.com
% Affiliation: Master Thesis in Econometrics, University of Bologna
% Description:
%   Extends LSFM estimation by selecting optimal model hyperparameters
%   (number of factors q, kernel bandwidth h, VAR lag p) via Bayesian
%   optimization. The objective minimizes mean squared forecast error (MSFE)
%   on key macroeconomic indicators (GDP, Unemployment, Inflation,
%   Interest Rate) using a cross-validation (CV) framework.

clear; close all; clc;

%% Dataset Selection
% Present data periodicity options and obtain user selection
options = {'Monthly (MD1959.xlsx)', 'Quarterly (QD1959.xlsx)'};
[choiceIndex, ok] = listdlg('PromptString', 'Select the dataset to load:', ...
                             'SelectionMode', 'single', ...
                             'ListString', options, ...
                             'Name', 'Dataset Selection', ...
                             'ListSize', [400 200]);
if ~ok
    error('Dataset selection cancelled. Exiting script.');
end

%% Load Data and Specify Key Variables
% Depending on frequency, read processed data and set time dimension T.
switch choiceIndex
    case 1  % Monthly data
        filepath = '/Users/moka/Research/Thesis/Live Project/Processed_Data/MD1959.xlsx';
        raw = readtable(filepath);
        x = table2array(raw);
        T = size(x,1);
        key_vars = [1, 24, 105, 77];  % Indices for: GDP, Unemp, Infl, IntRate
    case 2  % Quarterly data
        filepath = '/Users/moka/Research/Thesis/Live Project/Processed_Data/QD1959.xlsx';
        raw = readtable(filepath);
        x = table2array(raw(:,2:end));  % Exclude date column
        T = size(x,1);
        key_vars = [1, 58, 116, 136];
    otherwise
        error('Invalid choiceIndex.');
end
[N_obs, N_var] = size(x);

%% Define Training Horizon via Dialog
defaultTrain = '640';
prompt = sprintf('Dataset has %d observations. Enter training size (T_train):', T);
userInput = inputdlg(prompt, 'Training Size', [3 100], {defaultTrain});
if isempty(userInput)
    error('Training size not provided. Exiting.');
end
T_train = str2double(userInput{1});
assert(T_train>0 && T_train<T, 'T_train must lie in (0, %d)', T);

%% Standardize Training Sample
x_train = x(1:T_train, :);
mean_train = mean(x_train);
std_train  = std(x_train);
x_train_norm = (x_train - mean_train) ./ std_train;

%% Diagnostic: Variance Explained (Pre-Optimization)
q_max = 12;    % Maximum latent factors for diagnostic
h_temp = 0.10; % Fixed bandwidth
[CChat_d, ~, ~, Sigmahat] = lsfm(x_train_norm, q_max, h_temp);
% Compute eigenvalues over time for variance decomposition
eigvals = zeros(N_var, T_train);
for t = 1:T_train
    Sigma_t = squeeze(Sigmahat(:,:,t));
    eigvals(:,t) = sort(real(eig(Sigma_t)), 'descend');
end
cum_var = cumsum(eigvals,1) ./ sum(eigvals,1);
mean_var_explained = mean(cum_var,2);

%% Bayesian Optimization Setup
% Split last T_val observations as validation set
T_val = 24;
x_val_norm = x_train_norm(end-T_val+1:end, :);
x_cv_norm  = x_train_norm(1:end-T_val, :);

% Objective: minimize MSFE on key_vars using CV splits
dim_q = optimizableVariable('q',[3 q_max],'Type','integer');
dim_h = optimizableVariable('h',[0.05 0.25],'Type','real');
dim_p = optimizableVariable('p',[3 8],'Type','integer');
vars = [dim_q, dim_h, dim_p];
objfun = @(par) compute_cv_msfe(par.q, par.h, par.p, x_cv_norm, x_val_norm, key_vars);

% Execute Bayesian optimization
tmax = 75; % Max iterations
results = bayesopt(objfun, vars, ...
    'MaxObjectiveEvaluations', tmax, ...
    'AcquisitionFunctionName','expected-improvement-plus', ...
    'Verbose',1, 'UseParallel',false);

% Extract optimal hyperparameters
best = results.XAtMinObjective;
q_opt = best.q;
h_opt = best.h;
p_opt = best.p;

fprintf('Optimal: q=%d, h=%.3f, p=%d\n', q_opt, h_opt, p_opt);
fprintf('Min CV-MSFE: %.6f\n', results.MinObjective);

%% Re-estimate LSFM with Optimal q & h
[CChat_opt, Fhat_opt, Lhat_opt] = lsfm(x_train_norm, q_opt, h_opt);

%% Pre-Plot Diagnostic: Variance Explained by q_opt
fprintf('Variance explained by q_opt=%d: %.2f%%\n', q_opt, mean_var_explained(q_opt)*100);

%% Save Outputs
save('lsfm_estimation_results.mat','q_opt','h_opt','p_opt',... 
     'Fhat_opt','Lhat_opt','mean_train','std_train','T_train','N_var','key_vars');

fprintf('Estimation with optimized parameters complete. Results saved.\n');

%% compute_cv_msfe (Nested Helper)
function msfe=compute_cv_msfe(q,h,p,x_train_cv,x_val_norm,key_vars)
    % Fit LSFM on training CV segment
    [~,F_cv,L_cv] = lsfm(x_train_cv,q,h);
    % Forecast normalized series on validation horizon
    H = size(x_val_norm,1);
    yhat = forecast_LSFM(F_cv,L_cv,size(x_train_cv,1),size(x_val_norm,2),H,p);
    % Compute MSFE across key variables & periods
    se=(x_val_norm(:,key_vars)-yhat(:,key_vars)).^2;
    msfe = mean(se(:));
end

%% lsfm.m (Resolver)
% Same as previous implementation ensuring local covariance, eigen-decomp,
% and factor extraction.
function [CChat,Fhat,Lhat,Sigmahat]=lsfm(X,R,h)
    [T,N]=size(X);
    u=(1:T)'/T;
    CChat=zeros(T,N);Fhat=zeros(T,R);Lhat=zeros(N,R,T);
    Sigmahat=zeros(N,N,T);
    for t=1:T
        z=(u-u(t))/h;
        w=(1/sqrt(2*pi))*exp(-0.5*z.^2);
        w=w/sum(w);
        for i=1:N
            for j=1:i
                Sigmahat(i,j,t)=sum(w.*(X(:,i).*X(:,j)));
                Sigmahat(j,i,t)=Sigmahat(i,j,t);
            end
        end
    end
    for t=1:T
        S=squeeze(Sigmahat(:,:,t));S=(S+S')/2;
        [V,D]=eig(S);D=max(real(diag(D)),0);S=V*diag(D)*V';
        [A,~]=eigs(S,R,'largestabs');Lhat(:,:,t)=A;
        Fhat(t,:)=(A'*A)\(A'*X(t,:)');
        CChat(t,:)=(A*Fhat(t,:)')';
    end
end

%% forecast_LSFM.m (Resolver)
function yhat=forecast_LSFM(Fhat,Lhat,T,N,H,p)
    q=size(Fhat,2);
    if p>=T, p=floor(T/2); end
    mdl=varm(q,p);Est=estimate(mdl,Fhat);
    Ff=forecast(Est,H,Fhat);
    Lambda=squeeze(Lhat(:,:,T));
    yhat=zeros(H,N);
    for h1=1:H, yhat(h1,:)=(Lambda*Ff(h1,:)')'; end
end
