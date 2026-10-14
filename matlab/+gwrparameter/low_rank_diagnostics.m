function diagnostics = low_rank_diagnostics(n, rank)
%LOW_RANK_DIAGNOSTICS Vergleicht dichte und separierbare kleine Testmatrizen.

factors = gwrparameter.build_taylor_factors(n, rank);
exactKernel = exp(factors.x*factors.x.');
approximateKernel = factors.factor*factors.factor.';
exactOperator = exactKernel.*factors.weights.';
approximateOperator = approximateKernel.*factors.weights.';

kernelDifference = approximateKernel-exactKernel;
operatorDifference = approximateOperator-exactOperator;

diagnostics = struct();
diagnostics.n = n;
diagnostics.rank = rank;
diagnostics.remainderBound = exp(1)/factorial(rank);
diagnostics.kernelMaximumError = max(abs(kernelDifference), [], 'all');
diagnostics.kernelRelativeFrobeniusError = ...
    norm(kernelDifference, 'fro')/norm(exactKernel, 'fro');
diagnostics.operatorMaximumError = max(abs(operatorDifference), [], 'all');
diagnostics.operatorRelativeFrobeniusError = ...
    norm(operatorDifference, 'fro')/norm(exactOperator, 'fro');

end
