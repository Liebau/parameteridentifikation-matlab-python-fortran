function result = solve_tikhonov_low_rank(n, rank, beta)
%SOLVE_TIKHONOV_LOW_RANK Loest Tikhonov ueber die duenne SVD von sqrt(W)*P.

if ~isscalar(beta) || ~isfinite(beta) || beta <= 0
    error('gwrparameter:InvalidBeta', ...
        'beta muss ein positiver endlicher Skalar sein.');
end

factors = gwrparameter.build_taylor_factors(n, rank);
data = expm1(factors.x+1)./(factors.x+1);
exactLambda = exp(factors.x);
transformedData = sqrt(factors.weights).*data;

[leftVectors, singularMatrix, ~] = svd( ...
    factors.weightedFactor, 'econ');
singularValues = diag(singularMatrix);
eigenvalues = singularValues.^2;
filterFactors = eigenvalues./(beta+eigenvalues.^2);
transformedAlpha = leftVectors*( ...
    filterFactors.*(leftVectors.'*transformedData));
alpha = transformedAlpha./sqrt(factors.weights);

result = struct();
result.n = n;
result.rank = rank;
result.beta = beta;
result.alpha = alpha;
result.singularValues = singularValues;
result.maximumError = max(abs(alpha-exactLambda));
result.relativeError = result.maximumError/max(abs(exactLambda));

end
