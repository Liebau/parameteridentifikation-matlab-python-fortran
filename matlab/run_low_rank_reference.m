% RUN_LOW_RANK_REFERENCE Kleine Genauigkeitsstudie fuer Rang und beta.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

sizes = [32, 128];
ranks = [4, 8, 12, 16, 20, 24];
betas = [1.0e-3, 1.0e-6];

rows = numel(sizes)*numel(ranks)*numel(betas);
implementation = repmat("MATLAB", rows, 1);
nColumn = zeros(rows, 1);
rankColumn = zeros(rows, 1);
betaColumn = zeros(rows, 1);
remainderBound = zeros(rows, 1);
kernelMaximumError = zeros(rows, 1);
kernelRelativeFrobeniusError = zeros(rows, 1);
operatorMaximumError = zeros(rows, 1);
operatorRelativeFrobeniusError = zeros(rows, 1);
reconstructionMaximumError = zeros(rows, 1);
reconstructionRelativeError = zeros(rows, 1);

row = 0;
for n = sizes
    for rank = ranks
        diagnostics = gwrparameter.low_rank_diagnostics(n, rank);
        for beta = betas
            row = row+1;
            result = gwrparameter.solve_tikhonov_low_rank(n, rank, beta);
            nColumn(row) = n;
            rankColumn(row) = rank;
            betaColumn(row) = beta;
            remainderBound(row) = diagnostics.remainderBound;
            kernelMaximumError(row) = diagnostics.kernelMaximumError;
            kernelRelativeFrobeniusError(row) = ...
                diagnostics.kernelRelativeFrobeniusError;
            operatorMaximumError(row) = diagnostics.operatorMaximumError;
            operatorRelativeFrobeniusError(row) = ...
                diagnostics.operatorRelativeFrobeniusError;
            reconstructionMaximumError(row) = result.maximumError;
            reconstructionRelativeError(row) = result.relativeError;
        end
    end
end

results = table(implementation, nColumn, rankColumn, betaColumn, ...
    remainderBound, kernelMaximumError, ...
    kernelRelativeFrobeniusError, operatorMaximumError, ...
    operatorRelativeFrobeniusError, reconstructionMaximumError, ...
    reconstructionRelativeError, ...
    'VariableNames', {'implementation', 'n', 'rank', 'beta', ...
    'remainder_bound', 'kernel_max_error', ...
    'kernel_relative_frobenius_error', 'operator_max_error', ...
    'operator_relative_frobenius_error', ...
    'reconstruction_max_error', 'reconstruction_relative_error'});

disp(results(results.n == 128 & results.beta == 1.0e-6, :));
outputFile = fullfile(matlabRoot, 'low_rank_matlab.csv');
writetable(results, outputFile);
fprintf('Niedrigrangstudie gespeichert: %s\n', outputFile);
