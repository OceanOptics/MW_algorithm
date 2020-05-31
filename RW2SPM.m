
%% script example for the MW algorithm

close all
clear all; 
clc

%-------------------------------------------------------------------------%
%% INPUTS:

L         = [0.0949,0.0794]; %Gordon et al., 1988 constants
Q_filter  = 0.5; % saturation treshold 

% shape parameters
S         = 0.006:0.001:0.014;
Y         = 0:0.15:1.8;

% mass-specific coefficients
a_nap443     = 0.01:0.01:0.06;
a_nap750     = 0.013:0.001:0.015; 
bbp700       = 0.002:0.001:0.021;

temp = 29;

%-------------------------------------------------------------------------%
% LOADS sample data of 10 measurements from Nechad et al., 2015 (Indonesian waters)
load ('sample_data.mat') % loads Rw its std, SPM, and wavelength

%% Water-leaving reflectance to below water remote sensing reflectance

[U, rrs, nm] = RW2rrs(RW, nm, L);

[SPM_mw, err_mw, IOP_table] = MW_algorithm(nm, std, rrs, U, L, S, Y, a_nap443, bbp700, a_nap750, temp, Q_filter);

%% OUTPUTS:

%SPM_modeled:
SPM_mw;

%SPM_uncertainty:
err_mw;

% constrained shape paremeters and mass-specific coeffs: 
IOP_S       = IOP_table(:,:,6);
IOP_Y       = IOP_table(:,:,5);
IOP_ap443   = IOP_table(:,:,3);
IOP_ap750   = IOP_table(:,:,4);
IOP_bbp700  = IOP_table(:,:,2);




