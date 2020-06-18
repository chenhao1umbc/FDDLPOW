% this script is used to test LRSDL+SVM for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('/extra/chenhao1/DICTOL-master/LRSDL_FDDL'))

k = 10;
k0 = 5;
lambda1 = 0.001;
lambda2 = 0.01;
lambda3 = 0.02;

[acc_lrsdl, rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda1, lambda2, lambda3);