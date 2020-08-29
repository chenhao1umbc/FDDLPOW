% this script is KVSVD for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././DICTOL-master'))

% files are in the folder called utils
run lrsdl_train

run lrsdl_test

% These files only contains nn, logistic regression, MF, ZF
% SVM results are in seperate files

% How to train nn, logistic regression, svm, are in trainLR_NN_SVM.m