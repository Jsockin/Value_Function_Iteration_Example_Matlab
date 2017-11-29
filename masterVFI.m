%{ 
    This program is the master file for running value function iteration 
    under various grid specifications including:
            - exogenous grid
            - "accelerator"
            - multigrid
            - stochastic grid

        Author: Jason Sockin    Date: November 2017                                    
%}

%--------------------------------------------
% Section 0: Housekeeping
%--------------------------------------------

clear;
clc;

%--------------------------------------------
% Section 1: Parametrization
%--------------------------------------------

% Model Parameters
alpha   = 0.33;        % Production function C-D share
psi     = 0.5;         % Consumption utility share exponent
beta    = 0.96;        % Discount Factor
delta   = 0.1;         % Depreciation
epsilon = 10^-6;     % Error Threshold

% Values and transition matrix for technology shock z 
zVals = [-0.0673 -0.0336 0 0.0336 0.0673];
piZ   = [0.9727 0.0273 0      0      0;
          0.0041 0.9806 0.0153 0      0;
          0      0.0082 0.9836 0.0082 0;
          0      0      0.0153 0.9806 0.0041;
          0      0      0      0.0273 0.9727];

% Values and transition matrix for production shock A
aVals = [0.9  1   1.1];
piA   = [0.9  0.1 0;
          0.05 0.9 0.05;
          0    0.1 0.9];
      
% Additional parameters for looping
fminconOptions  = optimoptions('fmincon','Algorithm','sqp','Display','none');
imagThreshold   = 10^-4;
l1_bar          = 1;

% Suppress warning output from fmincon during parallelization
warningID = 'optimlib:fwdFinDiffInsideBnds:StepReduced';
warning('off',warningID);

%--------------------------------------------
% Section 2: Steady State
%--------------------------------------------

z_ss    = 0;
A_ss    = 1;
l1_ss   = fzero(@(l1) steadyStateLabor(alpha,beta,delta,psi,l1),0.5);
k_ss    = (((1/beta)+delta-1) / (alpha * l1_ss^(1-alpha))) ^ (1/(alpha-1));
l2_ss   = ((1-psi)/psi) * ((k_ss^alpha*l1_ss^(1-alpha) - delta*k_ss) / ( (1-alpha)*k_ss^alpha*l1_ss^(-alpha)));  
c2_ss   = A_ss * l2_ss;
c1_ss   = exp(z_ss)*k_ss^alpha*l1_ss^(1-alpha) + 0.9*k_ss - k_ss;

%--------------------------------------------
% Section 3: Capital Grid
%--------------------------------------------

percentBand = 30;
kmin        = k_ss - (percentBand/100)*k_ss;
kmax        = k_ss + (percentBand/100)*k_ss;

kPoints     = 1:250;
kVals       = (kmin + ((kPoints-1)/(length(kPoints)-1))*(kmax-kmin))';

%----------------------------------------------------
% Section 4: Value Function Iterate w/ Fixed Grid
%----------------------------------------------------
%{
result_exogrid = runExogenousGrid(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,imagThreshold,fminconOptions);
plotFunctions('exogrid','exogrid',kVals,zVals,aVals);
%}
%----------------------------------------------------
% Section 5: Accelerator
%----------------------------------------------------
%{
result_exogrid = runExogenousGrid(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,imagThreshold,fminconOptions);
plotFunctions('accelerator','acc',kVals,zVals,aVals);
plotComparison('exogrid','exogrid','accelerator','acc',kVals,kVals,zVals,aVals);
%}
%----------------------------------------------------
% Section 6: Multigrid
%----------------------------------------------------
%{
kVals_exogrid   = kVals;
gridPoints      = [100 500 5000];
kPoints         = 1:gridPoints(1);
kVals           = (kmin + ((kPoints-1)/(length(kPoints)-1))*(kmax-kmin))';

result_multi = runMultigrid(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,gridPoints,imagThreshold,fminconOptionsm);
plotFunctions('multigrid','multi',kVals,zVals,aVals);
plotComparison('exogrid','exogrid','multigrid','multi',kVals_exogrid,kVals,zVals,aVals);
%}
%----------------------------------------------------
% Section 7: Stochastic Grid
%----------------------------------------------------
%{
percentBand = 25;
kmin        = k_ss - (percentBand/100)*k_ss;
kmax        = k_ss + (percentBand/100)*k_ss;

kPoints = 1:250;
kVals = (kmin + ((kPoints-1)/(length(kPoints)-1))*(kmax-kmin))';

stochDraws = 500;
rng(1);

result_stoch = runStochastic(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,stochDraws,imagThreshold,fminconOptions);
plotFunctions('stochastic','stoch',kVals,zVals,aVals);
plotComparison('exogrid','exogrid','stochastic','stoch',kVals_exogrid,kVals,zVals,aVals);
%}