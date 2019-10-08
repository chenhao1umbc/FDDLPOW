% all tuning with 20db Dynamic Ratio, 3 signal mixture using cv data

beta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10 100];
%% overall whitening performance plot
result_beta = [0.841777777777778,0.841944444444445,0.841833333333333,0.838388888888889,...
    0.839277777777778,0.835555555555556,0.825888888888889,0.825666666666667,0.807722222222222];
result_betaWEEK = [0.525333333333333,0.525833333333333,0.525500000000000,0.515166666666667,...
    0.517833333333333,0.506666666666667,0.477666666666667,0.477000000000000,0.423166666666667];
semilogx(beta, result_beta, '-x')
hold on
semilogx(beta, result_betaWEEK, '-x')

% Beta = 0
Projected_WithinClass_covariance = 0.5767
Frobenius_norm_W = 8.3

% Beta = 0.1 --working
Projected_WithinClass_covariance = 0.4117
whitening_term_without_beta = 0.7038
Frobenius_norm_W = 8.2252

% Beta = 1
Projected_WithinClass_covariance = 0.1978
whitening_term_without_beta = 0.7751
Frobenius_norm_W = 4.1527

% Beta = 100 -- really whitened but not visually whitened
Projected_WithinClass_covariance = 0.7807
whitening_term_without_beta = 0.0124
Frobenius_norm_W = 2.7767

% Beta = 1e-5 -- gives the best performance
Projected_WithinClass_covariance = 0.9221
Frobenius_norm_W = 9.7503
whitening_term_without_beta = 1.1306


%% perclass whitening

% Beta = 0
Projected_WithinClass_covariance = 0.5767
Frobenius_norm_W = 8.3

% Beta = 0.1  -- working
Projected_WithinClass_covariance = 0.1865
whitening_term_without_beta = 0.0474
Frobenius_norm_W = 6.8290

% Beta = 1
Projected_WithinClass_covariance = 1.0201
whitening_term_without_beta = 0.0048
Frobenius_norm_W = 4.2480

% Beta = 100 -- some whitened some just shrink
Projected_WithinClass_covariance = 1.3769
whitening_term_without_beta = 0.0022
Frobenius_norm_W = 2.0892

% Beta = 1e-4 -- gives the best performance
Projected_WithinClass_covariance = 0.8784
Frobenius_norm_W = 9.3017
whitening_term_without_beta = 1.2376