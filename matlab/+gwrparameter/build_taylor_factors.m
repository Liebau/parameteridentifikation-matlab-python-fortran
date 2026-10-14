function factors = build_taylor_factors(n, rank)
%BUILD_TAYLOR_FACTORS Erzeugt P und sqrt(W)*P ohne dichte Operatormatrix.

if ~isscalar(n) || ~isfinite(n) || n < 1 || n ~= floor(n)
    error('gwrparameter:InvalidSubdivision', ...
        'n muss eine positive ganze Zahl sein.');
end
if ~isscalar(rank) || ~isfinite(rank) || rank < 1 || rank ~= floor(rank)
    error('gwrparameter:InvalidRank', ...
        'rank muss eine positive ganze Zahl sein.');
end

x = (0:n).'/n;
weights = ones(n+1, 1)/n;
weights([1, end]) = 1/(2*n);

factor = zeros(n+1, rank);
factor(:, 1) = 1.0;
for k = 2:rank
    factor(:, k) = factor(:, k-1).*x/sqrt(k-1);
end

factors = struct();
factors.n = n;
factors.rank = rank;
factors.x = x;
factors.weights = weights;
factors.factor = factor;
factors.weightedFactor = sqrt(weights).*factor;

end
