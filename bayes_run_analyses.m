%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   bayes_run_analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 09 August 2025
% This is the soup-to-nuts code. Click go, wait a couple hours,
% and you'll have the figures from Piecuch (2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

nameYouWant='myAnalysis'; % whatever you want to name your files
numYouWant=6; % number of separate chains you want to run
% 1. perform Bayesian data analysis
% for kk=1:numYouWant % <-- if you can't/don't want to run in parallel
%parfor kk=1:numYouWant
%    bayes_main_code(nameYouWant,kk); 
%end
% 2. perform posterior prediction to get solutions on grid
%bayes_posterior_prediction(nameYouWant,1:numYouWant)
% 3. make figures
bayes_make_figures(nameYouWant,1:numYouWant)
% 4. plot residuals
bayes_plot_residuals(nameYouWant,1:numYouWant)