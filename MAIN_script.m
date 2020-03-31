% Main script to call other code in this folder
% Created by Jack Baker 12/17/2019
%
% Results from these calculations are in:
%
% Baker and Chen (2020), "Ground motion spatial correlation fitting methods 
% and estimation uncertainty" (in review).

clear; close all; clc;

homeDir = [pwd]; % reference directory for outputs (the current directory by default


%% process NGA-West2 data to get residuals
if 1==0 % do new analysis (otherwise, if analysis is complete, we can load that data instead)
    cd('data')
    getResiduals % select data and run mixed effects regression to get residuals 
    cd(homeDir) % switch back to original directory
end


%% compute event-specific variograms for each earthquake (this script sets most of the user-defined parameters)
compute_variograms

%% get illustrative results for two earthquakes
study_two_events

%% generate synthetic data using station configurations, and fit variograms
if 1==1 % do new analysis
    load main_data
    [h, synthetic] = fn_synthetic_range (station_lat, station_long, recIdx, recsPerEQ, eventIdx, EQ_name_string, options);
    save synthetic_data h synthetic
end

%% plots of results from the synthetic data 
synthetic_data_plots

%% analyze estimation and true variance
variance_decomposition

%% Posterior analysis of empirical data
posterior_distributions

%% study role of WLS coefficient value
wls_coeff_study 

