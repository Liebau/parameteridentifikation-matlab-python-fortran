% RUN_BAUMANN_COMPARISON Vergleicht den direkten G-Loeser mit der duennen SVD.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

sizes = read_integer_list( ...
    'GWR_BAUMANN_SIZES', [64, 128, 256, 512, 1024, 8192]);
rank = read_positive_integer('GWR_BAUMANN_RANK', 16);
beta = read_positive_scalar('GWR_BAUMANN_BETA', 1.0e-6);
outputFile = getenv('GWR_BENCHMARK_OUTPUT');

nCases = numel(sizes);
maximumDifferences = zeros(nCases, 1);
reconstructionMaximumErrors = zeros(nCases, 1);
reducedResiduals = zeros(nCases, 1);
conditionGrams = zeros(nCases, 1);
conditionReduced = zeros(nCases, 1);
baumannSecondsColumn = zeros(nCases, 1);
svdSecondsColumn = zeros(nCases, 1);

fprintf('MATLAB-Version: %s\n', version);
fprintf('rank=%d, beta=%.1e\n', rank, beta);
fprintf(['n,N,max_abs_difference_to_svd,reconstruction_max_error,' ...
    'reduced_residual_inf,condition_gram,condition_reduced,' ...
    'baumann_workflow_ms,svd_workflow_ms,speedup_over_svd\n']);

for index = 1:nCases
    n = sizes(index);
    baumann = solve_baumann(n, rank, beta);
    svdResult = gwrparameter.solve_tikhonov_low_rank(n, rank, beta);
    exactLambda = exp(baumann.x);

    maximumDifference = max(abs(baumann.alpha-svdResult.alpha));
    reconstructionMaximumError = max(abs(baumann.alpha-exactLambda));

    % Beide Wege vor der Zeitmessung einmal ausfuehren.
    solve_baumann(n, rank, beta);
    gwrparameter.solve_tikhonov_low_rank(n, rank, beta);
    baumannSeconds = timeit(@() solve_baumann(n, rank, beta));
    svdSeconds = timeit( ...
        @() gwrparameter.solve_tikhonov_low_rank(n, rank, beta));

    fprintf(['%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,' ...
        '%.6f,%.6f,%.6f\n'], ...
        n, n+1, maximumDifference, reconstructionMaximumError, ...
        baumann.reducedResidual, cond(baumann.gramMatrix), ...
        cond(baumann.reducedMatrix), 1.0e3*baumannSeconds, ...
        1.0e3*svdSeconds, svdSeconds/baumannSeconds);

    maximumDifferences(index) = maximumDifference;
    reconstructionMaximumErrors(index) = reconstructionMaximumError;
    reducedResiduals(index) = baumann.reducedResidual;
    conditionGrams(index) = cond(baumann.gramMatrix);
    conditionReduced(index) = cond(baumann.reducedMatrix);
    baumannSecondsColumn(index) = baumannSeconds;
    svdSecondsColumn(index) = svdSeconds;

    assert( ...
        maximumDifference <= 2.0e-11+2.0e-11*max(abs(svdResult.alpha)), ...
        'Baumann- und SVD-Loesung weichen fuer n=%d zu stark ab.', n);
end

implementation = repmat("MATLAB", nCases, 1);
threadModeValue = string(getenv('GWR_BENCHMARK_THREAD_MODE'));
if strlength(threadModeValue) == 0
    threadModeValue = "unspecified";
end
threadMode = repmat(threadModeValue, nCases, 1);
threadCount = repmat(maxNumCompThreads, nCases, 1);
nColumn = sizes(:);
NColumn = nColumn+1;
rankColumn = repmat(rank, nCases, 1);
betaColumn = repmat(beta, nCases, 1);
speedupOverSvd = svdSecondsColumn./baumannSecondsColumn;

results = table(implementation, threadMode, threadCount, nColumn, NColumn, ...
    rankColumn, betaColumn, maximumDifferences, ...
    reconstructionMaximumErrors, reducedResiduals, conditionGrams, ...
    conditionReduced, baumannSecondsColumn, svdSecondsColumn, ...
    speedupOverSvd, ...
    'VariableNames', {'implementation', 'thread_mode', 'thread_count', ...
    'n', 'N', 'rank', 'beta', 'max_abs_difference_to_svd', ...
    'reconstruction_max_error', 'reduced_residual_inf', 'condition_gram', ...
    'condition_reduced', 'baumann_workflow_seconds', ...
    'svd_workflow_seconds', 'speedup_over_svd'});

if ~isempty(outputFile)
    outputDirectory = fileparts(outputFile);
    if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end
    writetable(results, outputFile);
    fprintf('Vergleich gespeichert: %s\n', outputFile);
end

function result = solve_baumann(n, rank, beta)
pointCount = n+1;
x = (0:n).'/n;

weights = ones(pointCount, 1)/n;
weights([1, end]) = 0.5/n;
data = expm1(x+1)./(x+1);

factor = zeros(pointCount, rank);
factor(:, 1) = 1.0;
for k = 2:rank
    factor(:, k) = factor(:, k-1).*x/sqrt(k-1);
end

gramMatrix = factor.'*(weights.*factor);
reducedRhs = factor.'*(weights.*data);
reducedMatrix = beta*eye(rank)+gramMatrix*gramMatrix;
coefficients = reducedMatrix\reducedRhs;
alpha = factor*coefficients;

result = struct();
result.x = x;
result.alpha = alpha;
result.gramMatrix = gramMatrix;
result.reducedMatrix = reducedMatrix;
result.coefficients = coefficients;
result.reducedResidual = norm( ...
    reducedMatrix*coefficients-reducedRhs, inf);
end

function values = read_integer_list(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    values = defaultValue;
    return;
end
values = sscanf(textValue, '%f').';
if isempty(values) || any(~isfinite(values)) || ...
        any(values < 1) || any(values ~= floor(values))
    error('GWR:InvalidEnvironment', ...
        '%s muss eine Liste positiver ganzer Zahlen sein.', variableName);
end
end

function value = read_positive_integer(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    value = defaultValue;
    return;
end
value = str2double(textValue);
if ~isfinite(value) || value < 1 || value ~= floor(value)
    error('GWR:InvalidEnvironment', ...
        '%s muss eine positive ganze Zahl sein.', variableName);
end
end

function value = read_positive_scalar(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    value = defaultValue;
    return;
end
value = str2double(textValue);
if ~isfinite(value) || value <= 0
    error('GWR:InvalidEnvironment', ...
        '%s muss ein positiver endlicher Skalar sein.', variableName);
end
end
