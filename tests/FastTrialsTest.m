classdef FastTrialsTest < matlab.unittest.TestCase
%FASTTRIALSTEST  The faster ways to Deconvolve's overlap-corrected trials
%   give the trials the toolbox's way gives, and say how far they are.
%
%   Alakazam's own arithmetic (trialsMethod 'fast', the default) is held to the toolbox's
%   uf_erpimage on a recording where every rule that decides a trial is at
%   work: a cut, artefact stretches left out of the model, a stimulus on
%   the sample of a response, and formulas with a covariate and a spline.
%   Both report their progress in the app's busy indicator (a stand-in for
%   it records what they report). The fit and the toolbox's trials on
%   parallel workers (Slow: starting workers takes a while) give what they
%   give one channel after another, to rounding error.
%
%   Run with: runtests('tests/FastTrialsTest.m').
%
%   See also UNFOLD.FITBINS, UNFOLD.WORKERPOOL, OVERLAPCORRECTEDTRIALSTEST.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            plugin = fileparts(fileparts(mfilename('fullpath')));
            root = getenv('ALAKAZAM_ROOT');   % Alakazam's source, which the plugin runs in
            if isempty(root)
                root = fullfile(fileparts(plugin), 'Alakazam');
            end
            testCase.assumeTrue(isfolder(fullfile(root, 'src')), ...
                'Alakazam''s source is not found (set ALAKAZAM_ROOT).');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(plugin, 'tests')));
            for p = {fullfile(root, 'src'), fullfile(root, 'src', 'Transformations'), ...
                     fullfile(plugin, 'Deconvolve'), ...
                     fullfile(root, 'src', 'Transformations', 'DefineBins'), ...
                     fullfile(root, 'src', 'Transformations', 'Average'), ...
                     fullfile(root, 'src', 'Support'), fullfile(root, 'src', 'IO'), ...
                     fullfile(root, 'src', 'Reports'), fullfile(root, 'src', 'Views'), ...
                     fullfile(plugin, 'tests')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end

    methods (Test)
        function theFastTrialsAreTheToolboxs(testCase)
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = FastTrialsTest.everyRuleAtWork();
            options = DeconvolveTest.options('output', 'trials', 'artifactThresholdUv', 150, ...
                'artifactChannels', 'all', 'parallel', 'off');

            fast = Deconvolve(EEG, options);   % no trialsMethod: the default
            options.trialsMethod = 'toolbox';
            toolbox = Deconvolve(EEG, options);

            info = toolbox.etc.alz.unfold;
            testCase.assertNotEmpty(info.excludedIntervals, 'The eye channel''s swings are left out.');
            testCase.assertGreaterThan(info.trialsDropped, 0, 'Some trials touch what was left out.');
            testCase.assertEqual(info.trialsSharingASample, 1, 'One stimulus shares a sample.');
            FastTrialsTest.verifySameTrials(testCase, fast, toolbox);
            testCase.verifyEqual(fast.etc.alz.unfold.trialsMethod, 'fast', ...
                'Options without a trialsMethod (templates, older nodes) take Alakazam''s way.');
            testCase.verifyEqual(info.trialsMethod, 'toolbox');
        end

        function aCovariateAndASplineToo(testCase)
        %ACOVARIATEANDASPLINETOO  An event's own fitted response is its row of
        %   the design times the betas, columns of a covariate and of a
        %   spline included, as uf_erpimage computes it.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            rng(8);
            for k = 1:numel(EEG.event)
                EEG.event(k).rt = 400 + 50 * randn();
            end
            formulas = struct('bin', {'Frequent', 'Rare'}, ...
                'formula', {'y ~ 1 + rt', 'y ~ 1 + spl(rt, 4)'});
            options = DeconvolveTest.options('output', 'trials', 'formulas', formulas, ...
                'parallel', 'off', 'trialsMethod', 'toolbox');

            toolbox = Deconvolve(EEG, options);
            options.trialsMethod = 'fast';
            fast = Deconvolve(EEG, options);

            FastTrialsTest.verifySameTrials(testCase, fast, toolbox);
        end

        function bothSayHowFarTheyAre(testCase)
        %BOTHSAYHOWFARTHEYARE  The toolbox's trials report each channel done,
        %   and the fast ones their trials, in the app's busy indicator, each
        %   ending at 1, where the indicator turns back into a spinner.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            writeText(fullfile(folder, 'busyGate.m'), strjoin({ ...
                'function busyGate(action, varargin)', ...
                'calls = getappdata(groot, ''FastTrialsTestGate'');', ...
                'value = [];', ...
                'if ~isempty(varargin), value = varargin{1}; end', ...
                'calls(end + 1, :) = {action, value};', ...
                'setappdata(groot, ''FastTrialsTestGate'', calls);', ...
                'end'}, newline));
            testCase.addTeardown(@rehash);   % after the path is restored: no stand-in stays cached
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
            rehash;
            testCase.addTeardown(@() rmappdata(groot, 'FastTrialsTestGate'));
            EEG = DeconvolveTest.untaggedRecording();
            options = DeconvolveTest.options('output', 'trials', 'parallel', 'off', ...
                'trialsMethod', 'toolbox');

            setappdata(groot, 'FastTrialsTestGate', cell(0, 2));
            Deconvolve(EEG, options);
            calls = getappdata(groot, 'FastTrialsTestGate');
            messages = calls(strcmp(calls(:, 1), 'message'), 2);
            testCase.verifyTrue(any(contains(messages, 'channel 1 of 2')), strjoin(messages, ' | '));
            testCase.verifyTrue(any(contains(messages, 'channel 2 of 2')));
            progress = [calls{strcmp(calls(:, 1), 'progress'), 2}];
            testCase.verifyEqual(progress(end), 1);
            testCase.verifyTrue(any(progress > 0 & progress < 1), 'Partway is shown, not only the end.');

            setappdata(groot, 'FastTrialsTestGate', cell(0, 2));
            options.trialsMethod = 'fast';
            Deconvolve(EEG, options);
            calls = getappdata(groot, 'FastTrialsTestGate');
            messages = calls(strcmp(calls(:, 1), 'message'), 2);
            testCase.verifyTrue(any(contains(messages, 'fast')), strjoin(messages, ' | '));
            progress = [calls{strcmp(calls(:, 1), 'progress'), 2}];
            testCase.verifyEqual(progress(end), 1);
        end

        function fitBinsMakesThemAlakazamsWayUnlessAsked(testCase)
        %FITBINSMAKESTHEMALAKAZAMSWAYUNLESSASKED  Unfold.fitBins called
        %   without a TrialsMethod, as a script would, takes the default too.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            [~, info] = Unfold.fitBins(UnfoldBinsTest.recording(), 'WindowMs', UnfoldBinsTest.WindowMs, ...
                'ArtifactThresholdUv', 0, 'Output', 'trials', 'Parallel', 'off');
            testCase.verifyEqual(info.trialsMethod, 'fast');
        end

        function aSmallJobStaysOnThisMatlab(testCase)
        %ASMALLJOBSTAYSONTHISMATLAB  Starting workers costs more than a short
        %   fit takes, so 'auto' leaves a small one here.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            trials = Deconvolve(DeconvolveTest.untaggedRecording(), ...
                DeconvolveTest.options('output', 'trials', 'trialsMethod', 'toolbox'));
            testCase.verifyEqual(trials.etc.alz.unfold.fitWorkers, 1);
            testCase.verifyEqual(trials.etc.alz.unfold.trialWorkers, 1);
        end
    end

    methods (Test, TestTags = {'Slow'})
        function workersGiveTheSameFitAndTrials(testCase)
        %WORKERSGIVETHESAMEFITANDTRIALS  The fit and the toolbox's trials on
        %   parallel workers ('parallel' 'on'), against one channel after
        %   another: the same numbers, the same trials. And none started
        %   from a worker.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')), ...
                'The Parallel Computing Toolbox is not installed.');
            hadPool = ~isempty(gcp('nocreate'));
            if ~hadPool
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
            EEG = FastTrialsTest.everyRuleAtWork();
            options = DeconvolveTest.options('output', 'trials', 'artifactThresholdUv', 150, ...
                'artifactChannels', 'all', 'trialsMethod', 'toolbox');

            options.parallel = 'off';
            here = Deconvolve(EEG, options);
            options.parallel = 'on';
            spread = Deconvolve(EEG, options);

            testCase.verifyGreaterThan(spread.etc.alz.unfold.fitWorkers, 1);
            testCase.verifyGreaterThan(spread.etc.alz.unfold.trialWorkers, 1);
            % To rounding error: a worker computes on one thread, and sums
            % some terms in another order than this MATLAB's threads do.
            testCase.verifyEqual(spread.data, here.data, 'AbsTol', 1e-9 * max(abs(here.data(:))));
            testCase.verifyEqual(spread.trials, here.trials);
            testCase.verifyEqual(spread.etc.alz.unfold.notes, here.etc.alz.unfold.notes);

            % On a worker (Apply to All Raw Files runs each recording on
            % one) the channels stay on that worker, without a word: a
            % worker cannot start workers of its own.
            said = fetchOutputs(parfeval(gcp('nocreate'), @evalc, 1, ...
                '[~, workers] = Unfold.workerPool(4, inf, 0, ''on'', 1); disp(workers)'));
            testCase.verifyEqual(strtrim(said), '1', said);
        end
    end

    methods (Static, Access = private)
        function EEG = everyRuleAtWork()
        %EVERYRULEATWORK  The test recording with an eye channel whose swings
        %   the artefact scan leaves out, a cut, and a stimulus moved onto
        %   the sample of a response.
            EEG = DeconvolveTest.withEyeChannel(DeconvolveTest.untaggedRecording());
            stimuli = find(strcmp({EEG.event.type}, 'S1'));
            responses = find(strcmp({EEG.event.type}, 'response'));
            EEG.event(responses(10)).latency = EEG.event(stimuli(10)).latency;
            EEG.event(end + 1) = struct('type', 'boundary', 'latency', 11500);
            [~, order] = sort([EEG.event.latency]);
            EEG.event = EEG.event(order);
        end

        function verifySameTrials(testCase, fast, toolbox)
        %VERIFYSAMETRIALS  The same trials of the same events, the same
        %   numbers to rounding error, and the same layout.
            testCase.verifyEqual(size(fast.data), size(toolbox.data));
            scale = max(abs(toolbox.data(:)));
            testCase.verifyEqual(fast.data, toolbox.data, 'AbsTol', 1e-9 * scale);
            testCase.verifyEqual(fast.trials, toolbox.trials);
            testCase.verifyEqual(fast.bindesc, toolbox.bindesc);
            testCase.verifyEqual(fast.epoch, toolbox.epoch);
            testCase.verifyEqual(fast.times, toolbox.times);
            a = fast.etc.alz.unfold;
            b = toolbox.etc.alz.unfold;
            testCase.verifyEqual([a.trials a.trialsDropped a.trialsSharingASample], ...
                [b.trials b.trialsDropped b.trialsSharingASample]);
        end
    end
end

% ======================================================================= %
function writeText(file, text)
    fid = fopen(file, 'w');
    fprintf(fid, '%s', text);
    fclose(fid);
end
