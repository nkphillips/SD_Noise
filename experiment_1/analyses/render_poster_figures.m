%%% render_poster_figures
% One-shot driver to render the six Phase 2 poster figures using the
% fresh 04.24.2026 estimates produced by Phase 1.

clear; clc; close all;

addpath('functions');
addpath('plotting');

analysis_date = '04.24.2026';
n_back_list   = [1 2 3];

disp('Rendering Poster Figure 1 (sigma grids)...');
plotResponseBiasSigma(analysis_date, n_back_list);

disp('Rendering Poster Figure 2 (sigma summary)...');
plotSigmaSummary(analysis_date, n_back_list);

disp('Rendering Poster Figure 3 (amplitude summary)...');
plotAmplitudeSummary(analysis_date, n_back_list);

disp('Rendering Poster Figure 4 (FWHM summary)...');
plotWidthSummary(analysis_date, n_back_list);

disp('Rendering Poster Figure 5 (A vs FWHM scatter)...');
plotAmplitudeWidthScatter(analysis_date, n_back_list);

disp('Rendering Poster Figure 6 (full 3x3 bias+DoG)...');
plotBiasDogGrid(analysis_date, n_back_list);

disp('Done.');
