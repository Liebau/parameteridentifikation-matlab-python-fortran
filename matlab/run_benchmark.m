% RUN_BENCHMARK Misst Matrixaufbau und Loeser getrennt.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

results = gwrparameter.measure_performance( ...
    [64, 128, 256, 512, 1024], 1.0e-6);
disp(results);

outputFile = fullfile(matlabRoot, 'benchmark_matlab.csv');
writetable(results, outputFile);
fprintf('Benchmark gespeichert: %s\n', outputFile);
fprintf('MATLAB-Version: %s\n', version);
fprintf('Thread-Modus: %s\n', getenv('GWR_BENCHMARK_THREAD_MODE'));
