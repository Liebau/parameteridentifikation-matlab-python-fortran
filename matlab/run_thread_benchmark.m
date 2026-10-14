% RUN_THREAD_BENCHMARK Parameterisierter Lauf fuer Threadskalierung.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

n = read_positive_integer('GWR_BENCHMARK_SIZE');
threadCount = read_positive_integer('GWR_BENCHMARK_THREADS');
outputFile = getenv('GWR_BENCHMARK_OUTPUT');
if isempty(outputFile)
    error('GWR:MissingEnvironment', ...
        'GWR_BENCHMARK_OUTPUT muss gesetzt sein.');
end

previousThreadCount = maxNumCompThreads(threadCount);
restoreThreads = onCleanup( ...
    @() maxNumCompThreads(previousThreadCount)); %#ok<NASGU>
activeThreadLimit = maxNumCompThreads;
if activeThreadLimit ~= threadCount
    error('GWR:ThreadLimitMismatch', ...
        'Angefordert: %d Threads, aktiv: %d Threads.', ...
        threadCount, activeThreadLimit);
end

results = gwrparameter.measure_performance(n, 1.0e-6);
results = addvars(results, repmat(threadCount, height(results), 1), ...
    'After', 'thread_mode', 'NewVariableNames', 'thread_count');

outputDirectory = fileparts(outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writetable(results, outputFile);

disp(results);
fprintf('Benchmark gespeichert: %s\n', outputFile);
fprintf('MATLAB-Version: %s\n', version);
fprintf('Thread-Limit: %d\n', activeThreadLimit);
fprintf('Thread-Modus: %s\n', getenv('GWR_BENCHMARK_THREAD_MODE'));

function value = read_positive_integer(variableName)
textValue = getenv(variableName);
value = str2double(textValue);
if isempty(textValue) || ~isfinite(value) || value < 1 || value ~= fix(value)
    error('GWR:InvalidEnvironment', ...
        '%s muss eine positive ganze Zahl sein.', variableName);
end
end
