% RUN_REFERENCE Reproduziert die Tabellen des GWR-Sage-Beispiels.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

ns = [1, 2, 4, 8, 16, 32];
nsExtended = [ns, 128];

fprintf('Direkte Loesung (ohne Regularisierung)\n');
fprintf('%5s %14s %14s %14s\n', ...
    'n', 'cond(A)', 'Fehler', 'rel. Fehler');
for n = ns
    result = gwrparameter.analyze(n);
    fprintf('%5d %14.6e %14.6e %13.4f%%\n', n, ...
        result.condition, result.directError, ...
        100*result.directRelativeError);
end

print_regularized(ns, 1.0e-3, "weighted");
print_regularized(ns, 1.0e-3, "euclidean");
print_regularized(nsExtended, 1.0e-6, "weighted");
print_regularized(nsExtended, 1.0e-6, "euclidean");

function print_regularized(ns, beta, mode)
fprintf('\nTikhonov: beta=%.0e, Adjungierter=%s\n', beta, mode);
fprintf('%5s %14s %14s\n', 'n', 'Fehler', 'rel. Fehler');
for n = ns
    result = gwrparameter.analyze(n, beta, mode);
    fprintf('%5d %14.6e %13.4f%%\n', n, ...
        result.regularizedError, ...
        100*result.regularizedRelativeError);
end
end
