% RUN_LOW_RANK_BENCHMARK Parameterisierter Niedrigrang-Zeitvergleich.

clearvars;
clc;

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);

sizes = read_integer_list('GWR_BENCHMARK_SIZES', 128);
rank = read_positive_integer('GWR_LOW_RANK_RANK', 16);
beta = read_positive_scalar('GWR_BENCHMARK_BETA', 1.0e-6);
requestedThreadCount = read_positive_integer( ...
    'GWR_BENCHMARK_THREADS', maxNumCompThreads);
previousThreadCount = maxNumCompThreads(requestedThreadCount);
restoreThreads = onCleanup( ...
    @() maxNumCompThreads(previousThreadCount)); %#ok<NASGU>
activeThreadCount = maxNumCompThreads;
if activeThreadCount ~= requestedThreadCount
    error('GWR:ThreadLimitMismatch', ...
        'Angefordert: %d Threads, aktiv: %d Threads.', ...
        requestedThreadCount, activeThreadCount);
end
outputFile = getenv('GWR_BENCHMARK_OUTPUT');
if isempty(outputFile)
    outputFile = fullfile(matlabRoot, 'low_rank_benchmark_matlab.csv');
end

results = gwrparameter.measure_low_rank_performance(sizes, rank, beta);
outputDirectory = fileparts(outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writetable(results, outputFile);

disp(results);
fprintf('Niedrigrangbenchmark gespeichert: %s\n', outputFile);
fprintf('MATLAB-Version: %s\n', version);
fprintf('Thread-Limit: %d\n', activeThreadCount);
fprintf('Thread-Modus: %s\n', getenv('GWR_BENCHMARK_THREAD_MODE'));

function values = read_integer_list(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    values = defaultValue;
    return;
end
values = sscanf(textValue, '%f').';
if isempty(values) || any(~isfinite(values)) || ...
        any(values < 1) || any(values ~= floor(values))
    error('GWR:InvalidEnvironment', ...
        '%s muss eine Liste positiver ganzer Zahlen sein.', variableName);
end
end

function value = read_positive_integer(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    value = defaultValue;
    return;
end
value = str2double(textValue);
if ~isfinite(value) || value < 1 || value ~= floor(value)
    error('GWR:InvalidEnvironment', ...
        '%s muss eine positive ganze Zahl sein.', variableName);
end
end

function value = read_positive_scalar(variableName, defaultValue)
textValue = getenv(variableName);
if isempty(textValue)
    value = defaultValue;
    return;
end
value = str2double(textValue);
if ~isfinite(value) || value <= 0
    error('GWR:InvalidEnvironment', ...
        '%s muss ein positiver endlicher Skalar sein.', variableName);
end
end
