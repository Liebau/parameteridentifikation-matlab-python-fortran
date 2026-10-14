function result = analyze(n, beta, adjointMode)
%ANALYZE Berechnet direkte und optional regularisierte Rekonstruktionen.
%
% adjointMode = "weighted"  : A_star = A, wie auf der GWR-Seite
% adjointMode = "euclidean" : A_star = A.', zum dokumentierten Vergleich

if nargin < 2
    beta = 0.0;
end
if nargin < 3
    adjointMode = "weighted";
end

if ~isscalar(beta) || ~isfinite(beta) || beta < 0
    error('gwrparameter:InvalidBeta', ...
        'beta muss ein endlicher nichtnegativer Skalar sein.');
end

problem = gwrparameter.build_problem(n);

% Die direkte Lösung ist für große n erwartungsgemäß nahezu singulär. Die
% Diagnose wird über cond und rcond zurückgegeben; Warnungen würden die
% Referenztabellen unnötig überdecken.
warningState = warning;
restoreWarnings = onCleanup(@() warning(warningState)); %#ok<NASGU>
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

alphaDirect = problem.A\problem.y;

result = struct();
result.n = n;
result.condition = cond(problem.A, 2);
result.reciprocalCondition = rcond(problem.A);
result.alphaDirect = alphaDirect;
result.directError = max(abs(alphaDirect-problem.exactLambda));
result.directRelativeError = ...
    result.directError/max(abs(problem.exactLambda));
result.beta = beta;
result.adjointMode = string(adjointMode);

if beta > 0
    switch string(adjointMode)
        case "weighted"
            adjoint = problem.A;
        case "euclidean"
            adjoint = problem.A.';
        otherwise
            error('gwrparameter:InvalidAdjointMode', ...
                'adjointMode muss "weighted" oder "euclidean" sein.');
    end

    lhs = beta*eye(n+1)+adjoint*problem.A;
    rhs = adjoint*problem.y;
    alphaRegularized = lhs\rhs;

    result.alphaRegularized = alphaRegularized;
    result.regularizedError = ...
        max(abs(alphaRegularized-problem.exactLambda));
    result.regularizedRelativeError = ...
        result.regularizedError/max(abs(problem.exactLambda));
end

end
