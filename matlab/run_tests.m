% RUN_TESTS Startet die MATLAB-Unittests der Parameteridentifikation.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

suite = testsuite(fullfile(matlabRoot, 'tests'), ...
    'IncludeSubfolders', true);
results = run(suite);
disp(results);
assert(all([results.Passed]), ...
    'Mindestens ein MATLAB-Test ist fehlgeschlagen.');
