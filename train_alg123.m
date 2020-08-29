% this script is KVSVD for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././DICTOL-master'))

% files are in the folder called utils
run alg1_tr
run alg1_cv

run alg2_tr
run alg2_cv

run alg3_tr
run alg3_cv

% These files only contains nn, logistic regression, MF, ZF
% SVM results are in seperate files

% How to train nn, logistic regression, svm, are in trainLR_NN_SVM.m