function results = measure_performance(sizes, beta)
%MEASURE_PERFORMANCE Misst getrennte Rechenphasen mit MATLAB timeit.

if nargin < 1
    sizes = [64, 128, 256, 512, 1024];
end
if nargin < 2
    beta = 1.0e-6;
end

sizes = sizes(:);
nCases = numel(sizes);
assemblySeconds = zeros(nCases, 1);
directSolveSeconds = zeros(nCases, 1);
tikhonovSetupSeconds = zeros(nCases, 1);
tikhonovSolveSeconds = zeros(nCases, 1);

% Warnungen werden einmal fuer die gesamte Messung deaktiviert. Eine
% Warnungsabfrage innerhalb der zeitgemessenen Funktion wuerde besonders bei
% kleinen Matrizen mehr kosten als die eigentliche lineare Algebra.
warningState = warning;
restoreWarnings = onCleanup(@() warning(warningState)); %#ok<NASGU>
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

for index = 1:nCases
    n = sizes(index);
    problem = gwrparameter.build_problem(n);
    [lhs, rhs] = make_tikhonov_system(problem.A, problem.y, beta);

    assemblySeconds(index) = timeit( ...
        @() gwrparameter.build_problem(n));
    directSolveSeconds(index) = timeit( ...
        @() solve_system(problem.A, problem.y));
    tikhonovSetupSeconds(index) = timeit( ...
        @() make_tikhonov_system(problem.A, problem.y, beta));
    tikhonovSolveSeconds(index) = timeit( ...
        @() solve_system(lhs, rhs));
end

implementation = repmat("MATLAB", nCases, 1);
threadModeValue = string(getenv('GWR_BENCHMARK_THREAD_MODE'));
if strlength(threadModeValue) == 0
    threadModeValue = "unspecified";
end
threadMode = repmat(threadModeValue, nCases, 1);
betaColumn = repmat(beta, nCases, 1);
results = table(implementation, threadMode, sizes, betaColumn, assemblySeconds, ...
    directSolveSeconds, tikhonovSetupSeconds, tikhonovSolveSeconds, ...
    'VariableNames', {'implementation', 'thread_mode', 'n', 'beta', ...
    'assembly_seconds', 'direct_solve_seconds', ...
    'tikhonov_setup_seconds', 'tikhonov_solve_seconds'});

end

function solution = solve_system(A, rhs)
solution = A\rhs;
end

function [lhs, rhs] = make_tikhonov_system(A, y, beta)
lhs = beta*eye(size(A, 1))+A*A;
rhs = A*y;
end
