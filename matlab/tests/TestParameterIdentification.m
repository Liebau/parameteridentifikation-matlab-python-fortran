classdef TestParameterIdentification < matlab.unittest.TestCase

    methods (Test)
        function buildsTrapezoidalMatrix(testCase)
            problem = gwrparameter.build_problem(2);

            expectedX = [0; 0.5; 1.0];
            expectedWeights = [0.25; 0.5; 0.25];
            expectedA = exp(expectedX*expectedX.').*expectedWeights.';

            testCase.verifyEqual(problem.x, expectedX, 'AbsTol', 1e-15);
            testCase.verifyEqual(problem.weights, expectedWeights, ...
                'AbsTol', 1e-15);
            testCase.verifyEqual(problem.A, expectedA, 'AbsTol', 1e-15);
        end

        function reproducesStableDirectCase(testCase)
            result = gwrparameter.analyze(4);

            testCase.verifyEqual(result.condition, 9.1002107e5, ...
                'RelTol', 2e-7);
            testCase.verifyEqual(result.directError, 1.08964, ...
                'AbsTol', 1e-4);
        end

        function weightedRegularizationImprovesFineGrid(testCase)
            direct = gwrparameter.analyze(32);
            regularized = gwrparameter.analyze(32, 1.0e-6, "weighted");

            testCase.verifyLessThan(regularized.regularizedError, 0.02);
            testCase.verifyLessThan( ...
                regularized.regularizedError, direct.directError);
        end

        function adjointConventionsRemainExplicit(testCase)
            weighted = gwrparameter.analyze(32, 1.0e-3, "weighted");
            euclidean = gwrparameter.analyze(32, 1.0e-3, "euclidean");

            testCase.verifyLessThan(weighted.regularizedError, 0.06);
            testCase.verifyGreaterThan(euclidean.regularizedError, 1.0);
        end

        function rejectsInvalidSubdivision(testCase)
            testCase.verifyError( ...
                @() gwrparameter.build_problem(0), ...
                'gwrparameter:InvalidSubdivision');
        end

        function taylorErrorDecreasesWithRank(testCase)
            rank8 = gwrparameter.low_rank_diagnostics(32, 8);
            rank16 = gwrparameter.low_rank_diagnostics(32, 16);

            testCase.verifyLessThan(rank16.kernelMaximumError, 2.0e-13);
            testCase.verifyLessThan( ...
                rank16.kernelMaximumError, rank8.kernelMaximumError);
        end

        function matrixFreeApplicationMatchesDenseOperator(testCase)
            problem = gwrparameter.build_problem(32);
            factors = gwrparameter.build_taylor_factors(32, 18);
            values = cos(problem.x);

            dense = problem.A*values;
            matrixFree = gwrparameter.apply_low_rank(factors, values);
            testCase.verifyEqual(matrixFree, dense, ...
                'RelTol', 2e-15, 'AbsTol', 2e-15);
        end

        function lowRankTikhonovMatchesDenseSolution(testCase)
            dense = gwrparameter.analyze(32, 1.0e-6, "weighted");
            lowRank = gwrparameter.solve_tikhonov_low_rank( ...
                32, 20, 1.0e-6);

            testCase.verifyEqual(lowRank.alpha, dense.alphaRegularized, ...
                'RelTol', 2e-9, 'AbsTol', 2e-9);
        end
    end
end
