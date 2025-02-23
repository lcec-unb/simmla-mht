%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Multi-regressors comparison by hiperparameters bayesian optimzation
%
%                             AUTHOR: João Luiz Junho Pereira
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
warning off

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

% Parameters
Runs = 30;  % Number of runs per regressor
kfolds = 5; % Number of cross validation folds

% Dataset
A = load('DS.txt'); % (fi, H [A/m], omega [Hz], a [m]) and  (Tc(°C), t(s))
X = A(:,1:4);
Y = A(:,5:6);

% Cross-validation folds
kIdx = crossvalind('Kfold', size(X,1), kfolds);

%% SECTION 2 - MODELS SELECTION - Run for each response i
i = 1; % Desired response
    clc
    clear L MM1 MM2 MM3 MM4 MM5 MM6 MM7 RrmseK RrmseKsd M1 M2 M3 M4 M5 M6 M7

    % MODELS HYPERPARAMETER OPTIMIZATION 
    % Regression Gaussian Process
    MM1 = fitrgp(X,Y(:,i),'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
    L(i,1) = MM1.HyperparameterOptimizationResults.MinObjective;

    % Support Vector Machines 
    MM2 = fitrsvm(X,Y(:,i),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'));
    L(i,2) = MM2.HyperparameterOptimizationResults.MinObjective;

    % Decision Trees
    MM3 = fitrtree(X,Y(:,i),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'));
    L(i,3) = MM3.HyperparameterOptimizationResults.MinObjective;

    % Linear Regression
    hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
    [MM4,FitInfo,HyperparameterOptimizationResults] = fitrlinear(X,Y(:,i),...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',hyperopts);
    L(i,4) = HyperparameterOptimizationResults.MinObjective;

    % Ensemble of learners
    t = templateTree('Reproducible',true);
    MM5 = fitrensemble(X,Y(:,i),'OptimizeHyperparameters','auto','Learners',t, ...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));
    L(i,5) = MM5.HyperparameterOptimizationResults.MinObjective;

    % Kernels
    [MM6,FitInfo,HyperparameterOptimizationResults] = fitrkernel(X,Y(:,i),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));
    L(i,6) = HyperparameterOptimizationResults.MinObjective;

    % Artificial Neural Networks   
    MM7 = fitrnet(X,Y(:,i),"OptimizeHyperparameters","auto", ...
    "HyperparameterOptimizationOptions",struct("AcquisitionFunctionName","expected-improvement-plus"));
    L(i,7) = MM7.HyperparameterOptimizationResults.MinObjective;

     
    % MODELS COMPARISON THROUGH CROSS VALIDATION
    for k = 1 : kfolds
    % Splitting into trainning and test data
    % Trainning
    Predictors = X(kIdx~=k,:);
    Response = Y(kIdx~=k,i);
    % Testing
    TEX = X(kIdx==k,:);
    TEY = Y(kIdx==k,i);
    % Normalizing
        % extracting normalization factors from training data
         minimums = min(Predictors, [], 1);
         ranges = max(Predictors, [], 1) - minimums;
        % normalizing training data	
         Predictors = (Predictors - repmat(minimums, size(Predictors, 1), 1)) ./ repmat(ranges, size(Predictors, 1), 1);
        % normalizing testing data based on factors extracted from training data
         TEX = (TEX - repmat(minimums, size(TEX, 1), 1)) ./ repmat(ranges, size(TEX, 1), 1);

    for j = 1 : Runs
    % RGP 
    M1 = fitrgp(...
    Predictors, ...
    Response, ...
    'BasisFunction', MM1.BasisFunction, ...
    'KernelFunction', MM1.KernelFunction, ...
    'Sigma', MM1.Sigma, ...
    'Standardize', false);
    RMSEteste(j,1) = rmse(TEY,predict(M1,TEX));

    % SVM
    responseScale = iqr(Response);
    if ~isfinite(responseScale) || responseScale == 0.0
    responseScale = 1.0;
    end
    boxConstraint = responseScale/1.349;
    epsilon = responseScale/13.49;
    M2 = fitrsvm(...
        Predictors, ...
        Response, ...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', MM2.KernelParameters.Scale, ...
        'BoxConstraint', boxConstraint, ...
        'Epsilon', epsilon, ...
        'Standardize', false);
    RMSEteste(j,2) = rmse(TEY,predict(M2,TEX));

    % DT
    M3 = fitrtree(...
    Predictors, ...
    Response, ...
    'MinLeafSize', MM3.HyperparameterOptimizationResults.XAtMinObjective.MinLeafSize, ...
    'Surrogate', 'off');
    RMSEteste(j,3) = rmse(TEY,predict(M3,TEX));

    % LR 
    M4 = fitrlinear(Predictors,Response,'Lambda',MM4.Lambda,...
    'Learner',MM4.Learner,'Solver','sparsa','Regularization','lasso');
    RMSEteste(j,4) = rmse(TEY,predict(M4,TEX));

    % ENS
    template = templateTree(...
    'MinLeafSize', MM5.HyperparameterOptimizationResults.XAtMinObjective.MinLeafSize, ...
    'NumVariablesToSample', 'all');
    Method = char(MM5.HyperparameterOptimizationResults.XAtMinObjective.Method);
    M5 = fitrensemble(...
    Predictors, ...
    Response, ...
    'Method', Method, ...
    'NumLearningCycles', MM5.HyperparameterOptimizationResults.XAtMinObjective.NumLearningCycles, ...
    'Learners', template);
    RMSEteste(j,5) = rmse(TEY,predict(M5,TEX));

    % KERNELS
    M6 = fitrkernel(...
    Predictors, ...
    Response, ...
    'Learner', MM6.Learner, ...
    'NumExpansionDimensions', MM6.NumExpansionDimensions, ...
    'Lambda', MM6.Lambda, ...
    'KernelScale', MM6.KernelScale, ...
    'IterationLimit', 1000);
    RMSEteste(j,6) = rmse(TEY,predict(M6,TEX));

    % ANN
    M7 = fitrnet(...
    Predictors, ...
    Response, ...
    'LayerSizes', MM7.LayerSizes, ...
    'Activations', MM7.Activations, ...
    'Lambda', MM7.HyperparameterOptimizationResults.XAtMinObjective.Lambda, ...
    'IterationLimit', 1000, ...
    'Standardize', false);
    RMSEteste(j,7) = rmse(TEY,predict(M7,TEX));
    end
    RrmseK(k,:) = mean(RMSEteste);
    RrmseKsd(k,:) = std(RMSEteste);
    end
Rrmse = mean(RrmseK);
Rrmsesd = mean(RrmseKsd);
[a,b]=min(Rrmse);
% SELECTION RESULTS
disp('Best Model')
if b == 1; disp('RGP'); end
if b == 2; disp('SVM'); end
if b == 3; disp('DT'); end
if b == 4; disp('LR'); end
if b == 5; disp('ENSEMBLE'); end
if b == 6; disp('KERNEL'); end
if b == 7; disp('ANN'); end
disp('Loss in the BO')
disp(L)
disp('Mean Rmse in test')
disp(Rrmse);
disp('Std Rmse in test')
disp(Rrmsesd);