function results = measure_low_rank_performance(sizes, rank, beta)
%MEASURE_LOW_RANK_PERFORMANCE Misst getrennte Phasen der Niedrigrangrechnung.

if nargin < 1
    sizes = 128;
end
if nargin < 2
    rank = 16;
end
if nargin < 3
    beta = 1.0e-6;
end

if ~isscalar(rank) || ~isfinite(rank) || rank < 1 || rank ~= floor(rank)
    error('gwrparameter:InvalidRank', ...
        'rank muss eine positive ganze Zahl sein.');
end
if ~isscalar(beta) || ~isfinite(beta) || beta <= 0
    error('gwrparameter:InvalidBeta', ...
        'beta muss ein positiver endlicher Skalar sein.');
end

sizes = sizes(:);
nCases = numel(sizes);
factorAssemblySeconds = zeros(nCases, 1);
operatorApplySeconds = zeros(nCases, 1);
thinSvdSeconds = zeros(nCases, 1);
reducedSolveSeconds = zeros(nCases, 1);
additionalRhsSeconds = zeros(nCases, 1);
regularizedWorkflowSeconds = zeros(nCases, 1);
reconstructionMaximumError = zeros(nCases, 1);
reconstructionRelativeError = zeros(nCases, 1);

for index = 1:nCases
    n = sizes(index);
    factors = gwrparameter.build_taylor_factors(n, rank);
    data = expm1(factors.x+1)./(factors.x+1);
    transformedData = sqrt(factors.weights).*data;
    additionalData = sqrt(factors.weights).*cos(factors.x);
    values = cos(factors.x);

    [leftVectors, singularMatrix, ~] = svd( ...
        factors.weightedFactor, 'econ');
    eigenvalues = diag(singularMatrix).^2;
    reference = gwrparameter.solve_tikhonov_low_rank(n, rank, beta);

    factorAssemblySeconds(index) = timeit( ...
        @() gwrparameter.build_taylor_factors(n, rank));
    operatorApplySeconds(index) = timeit( ...
        @() gwrparameter.apply_low_rank(factors, values));
    thinSvdSeconds(index) = timeit( ...
        @() compute_thin_svd(factors.weightedFactor));
    reducedSolveSeconds(index) = timeit( ...
        @() apply_reduced_solution(leftVectors, eigenvalues, ...
        transformedData, factors.weights, beta));
    additionalRhsSeconds(index) = timeit( ...
        @() apply_reduced_solution(leftVectors, eigenvalues, ...
        additionalData, factors.weights, beta));
    regularizedWorkflowSeconds(index) = timeit( ...
        @() gwrparameter.solve_tikhonov_low_rank(n, rank, beta));
    reconstructionMaximumError(index) = reference.maximumError;
    reconstructionRelativeError(index) = reference.relativeError;
end

implementation = repmat("MATLAB", nCases, 1);
threadModeValue = string(getenv('GWR_BENCHMARK_THREAD_MODE'));
if strlength(threadModeValue) == 0
    threadModeValue = "unspecified";
end
threadMode = repmat(threadModeValue, nCases, 1);
threadCount = repmat(maxNumCompThreads, nCases, 1);
rankColumn = repmat(rank, nCases, 1);
betaColumn = repmat(beta, nCases, 1);

results = table(implementation, threadMode, threadCount, sizes, rankColumn, ...
    betaColumn, factorAssemblySeconds, operatorApplySeconds, thinSvdSeconds, ...
    reducedSolveSeconds, additionalRhsSeconds, regularizedWorkflowSeconds, ...
    reconstructionMaximumError, reconstructionRelativeError, ...
    'VariableNames', {'implementation', 'thread_mode', 'thread_count', 'n', ...
    'rank', 'beta', 'factor_assembly_seconds', 'operator_apply_seconds', ...
    'thin_svd_seconds', 'reduced_solve_seconds', 'additional_rhs_seconds', ...
    'regularized_workflow_seconds', 'reconstruction_max_error', ...
    'reconstruction_relative_error'});

end

function decomposition = compute_thin_svd(weightedFactor)
[leftVectors, singularMatrix, rightVectors] = svd(weightedFactor, 'econ');
decomposition = {leftVectors, singularMatrix, rightVectors};
end

function alpha = apply_reduced_solution( ...
        leftVectors, eigenvalues, transformedData, weights, beta)
filterFactors = eigenvalues./(beta+eigenvalues.^2);
transformedAlpha = leftVectors*( ...
    filterFactors.*(leftVectors.'*transformedData));
alpha = transformedAlpha./sqrt(weights);
end
