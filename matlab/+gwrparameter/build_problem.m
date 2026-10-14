function problem = build_problem(n)
%BUILD_PROBLEM Diskretisiert den Integraloperator mit der Trapezregel.

if ~isscalar(n) || ~isfinite(n) || n < 1 || n ~= floor(n)
    error('gwrparameter:InvalidSubdivision', ...
        'n muss eine positive ganze Zahl sein.');
end

x = (0:n).'/n;
weights = ones(n+1, 1)/n;
weights([1, end]) = 1/(2*n);

kernel = exp(x*x.');
A = kernel.*weights.';

problem = struct();
problem.n = n;
problem.x = x;
problem.weights = weights;
problem.A = A;
problem.y = expm1(x+1)./(x+1);
problem.exactLambda = exp(x);

end
