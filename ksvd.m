% this script is KVSVD for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././DICTOL-master'))

% files are in the folder called utils
run ksvd_train

run ksvd_test