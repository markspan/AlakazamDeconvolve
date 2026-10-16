function [pool, workers] = workerPool(jobs, workload, threshold, mode, bytesPerWorker)
%WORKERPOOL  Parallel workers for a step done channel by channel, or none.
%   [POOL, WORKERS] = Unfold.workerPool(JOBS, WORKLOAD, THRESHOLD, MODE,
%   BYTESPERWORKER) is a running process pool to hand JOBS independent
%   per-channel jobs to, with WORKERS its size, or [] and 1 when the jobs are
%   better done here, one after another.
%
%   MODE 'off' never uses workers. 'on' uses them whenever the Parallel
%   Computing Toolbox can give two or more (the tests use it). 'auto', the
%   default, uses them only when the work is large, WORKLOAD (a cost the
%   caller measures) at least THRESHOLD, and Alakazam's settings allow
%   parallel processing (Settings > Processing, Process recordings in
%   parallel): starting workers takes a while, which pays only on a long job.
%   On a worker itself (Apply to All Raw Files runs each recording on one)
%   no workers are started: the work stays on that worker.
%
%   HOW MANY: the smallest of JOBS, the physical cores, and what the free
%   memory allows at BYTESPERWORKER each, with 2 GB left for the app. A
%   process pool already running is used as it is. A thread pool is not:
%   its workers cannot take this MATLAB's path, which EEGLAB and the toolbox
%   need. The path is handed to the workers every time, since a pool started
%   earlier need not have the toolbox on it.
%
%   See also UNFOLD.FITBINS, PARPOOL.
    pool = [];
    workers = 1;
    mode = lower(char(string(mode)));
    if strcmp(mode, 'off') || jobs < 2
        return;
    end
    if ~(license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')))
        return;
    end
    if onAWorker()
        return;
    end
    if strcmp(mode, 'auto') && (workload < threshold || ~settingsAllow())
        return;
    end
    try
        running = gcp('nocreate');
    catch
        return;
    end
    if ~isempty(running) && isa(running, 'parallel.ThreadPool')
        return;
    end
    if isempty(running)
        wanted = min([jobs, feature('numcores'), byMemory(bytesPerWorker)]);
        if wanted < 2
            return;
        end
        TransTools.BusyGate('message', sprintf('Starting %d parallel workers...', wanted));
        try
            running = parpool('Processes', wanted);
        catch err
            warning('Alakazam:Unfold:noPool', ['I could not start parallel workers (%s), so ' ...
                'the channels are done one at a time.'], err.message);
            return;
        end
    end
    if running.NumWorkers < 2
        return;
    end
    pool = running;
    workers = pool.NumWorkers;
    wait(parfevalOnAll(pool, @path, 0, path));
end

function tf = settingsAllow()
%SETTINGSALLOW  Alakazam's own switch for parallel processing; on where the
%   settings cannot be read (a script without the app).
    try
        tf = logical(AlakazamSettings.get('processing', 'applyToAll', 'parallel'));
    catch
        tf = true;
    end
end

function tf = onAWorker()
%ONAWORKER  True inside a parallel worker, which cannot start a pool of its
%   own.
    try
        tf = ~isempty(getCurrentTask());
    catch
        tf = false;
    end
end

function n = byMemory(bytesPerWorker)
%BYMEMORY  How many workers the free memory holds; any number where it
%   cannot be read.
    n = inf;
    try
        free = availableMemory();
        if isfinite(free)
            n = floor((free - 2e9) / bytesPerWorker);
        end
    catch
    end
end
