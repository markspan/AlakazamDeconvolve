classdef TemplateReplayTest < matlab.unittest.TestCase
%TEMPLATEREPLAYTEST  The plugin's two templates: every bin script in them
%   parses, and each, replayed on its own data the way Apply Template
%   replays it, gives what it was checked to give.
%
%   FaceSaccadesDeconvolution was replayed on the authors' face recording
%   when it was written, and the result written down in the manual (its
%   section "Overlap correction by deconvolution"): 3323 of the 3829 events
%   keep an overlap-corrected trial, and the template's Average of them
%   peaks at 17.9 uV between 300 and 500 ms at Oz. ReadingDeconvolution has
%   no recorded result; it is checked to run end to end on EYE-EEG's
%   reading data and to produce what it is for. It runs Eye tracking, a
%   plugin of its own, whose repository is looked for beside this one.
%
%   THE DATA IS NOT IN THE REPOSITORY. The recordings are looked for under
%   the data folder (ALAKAZAM_DATA, Alakazam's Data folder, or D:\data), as
%   opendata/face_saccades_opendata_fig10.set and reading/reading_eeg.set.
%   A case skips rather than fails where its recording, or a toolbox it
%   needs (EEGLAB, Unfold, EYE-EEG), is missing. The replays are tagged
%   Slow: a deconvolution of a whole recording takes minutes.
%
%   Run with: runtests('tests/TemplateReplayTest.m').
%
%   See also DECONVOLVE, FACESACCADESTEMPLATETEST.

    properties
        Alakazam = ''   % Alakazam's source folder, which the plugin runs in
        Plugin = ''     % this repository
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            plugin = fileparts(fileparts(mfilename('fullpath')));
            root = getenv('ALAKAZAM_ROOT');   % Alakazam's source, which the plugin runs in
            if isempty(root)
                root = fullfile(fileparts(plugin), 'Alakazam');
            end
            testCase.assumeTrue(isfolder(fullfile(root, 'src')), ...
                'Alakazam''s source is not found (set ALAKAZAM_ROOT).');
            testCase.Alakazam = root;
            testCase.Plugin = plugin;
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'src', 'Transformations'), 'IncludeSubfolders', true));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src', 'Support')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(plugin, 'Deconvolve')));
            try
                EEGLabEnvironment.ensure();
            catch
                % The assumptions in each case say what is missing.
            end
        end
    end

    methods (Test)
        function everyTemplateBinScriptParses(testCase)
        %EVERYTEMPLATEBINSCRIPTPARSES  The bin language refuses text before
        %   its first statement and knows an epoch statement, so a script
        %   that parsed when a template was written could stop parsing. The
        %   templates are the scripts users start from.
            files = dir(fullfile(testCase.Plugin, 'Deconvolve', 'templates', '*.alztemplate'));
            testCase.assertNotEmpty(files, 'No templates found to check.');
            checked = 0;
            for f = 1:numel(files)
                template = jsondecode(fileread(fullfile(files(f).folder, files(f).name)));
                steps = templateSteps(template);
                for k = 1:numel(steps)
                    if ~strcmp(steps(k).transformId, 'Deconvolve') || ~isstruct(steps(k).params) ...
                            || ~isfield(steps(k).params, 'binScript')
                        continue;
                    end
                    checked = checked + 1;
                    try
                        DefineBinsEngine.parseSpec(char(steps(k).params.binScript));
                    catch err
                        testCase.verifyFail(sprintf('%s, step %d: the bin script does not parse: %s', ...
                            files(f).name, k, err.message));
                    end
                end
            end
            testCase.verifyEqual(checked, numel(files), 'Each template has one Deconvolve step.');
        end
    end

    methods (Test, TestTags = {'Slow'})
        function theFaceSaccadeTemplateGivesFigureElevensNumbers(testCase)
        %THEFACESACCADETEMPLATEGIVESFIGUREELEVENSNUMBERS  3323 of the 3829
        %   events keep an overlap-corrected trial, and the template's
        %   Average of them peaks at 17.9 uV between 300 and 500 ms at Oz,
        %   within 0.3 uV of the fitted waveform's 17.7.
            file = testCase.dataFile('opendata', 'face_saccades_opendata_fig10.set');
            nodes = testCase.templateNodes('FaceSaccadesDeconvolution.alztemplate');
            results = testCase.replay(nodes, file);

            trials = results{stepIndex(nodes, 'Deconvolve')};
            average = results{stepIndex(nodes, 'Average')};
            testCase.verifyEqual(trials.trials, 3323, 'Overlap-corrected trials kept.');
            testCase.verifyEqual(peakIn(average, 'Stimulus', 'Oz', [300 500]), 17.9, 'AbsTol', 0.06);
        end

        function theReadingTemplateRunsEndToEnd(testCase)
        %THEREADINGTEMPLATERUNSENDTOEND  No recorded result: the eye track
        %   joins onto the EEG, with its synchronisation recorded, and the
        %   deconvolution returns a waveform per bin with numbers in it.
            eye = fullfile(fileparts(testCase.Plugin), 'AlakazamEyetracking', 'EyeTracking');
            testCase.assumeTrue(isfolder(eye), ['The reading template runs Eye tracking, ' ...
                'whose repository (AlakazamEyetracking) is not beside this one.']);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(eye));
            file = testCase.dataFile('reading', 'reading_eeg.set');
            nodes = testCase.templateNodes('ReadingDeconvolution.alztemplate');
            results = testCase.replay(nodes, file);

            joined = results{stepIndex(nodes, 'EyeTracking')};
            testCase.verifyTrue(isfield(joined.etc, 'alz') && isfield(joined.etc.alz, 'eyeTracking'), ...
                'EyeTracking should record how good the synchronisation was.');
            testCase.verifyTrue(any(strcmpi({joined.chanlocs.type}, 'EYE')), ...
                'The gaze and pupil channels should be typed EYE.');

            fitted = results{stepIndex(nodes, 'Deconvolve')};
            testCase.verifyEqual(char(fitted.DataFormat), 'Averaged');   % a string, as Average sets it
            testCase.verifyNotEmpty(fitted.bindesc);
            testCase.verifyTrue(any(isfinite(fitted.data(:))), 'The fit returned no numbers.');
        end
    end

    methods (Access = private)
        function nodes = templateNodes(testCase, name)
        %TEMPLATENODES  A template as Apply Template reads it
        %   (Alakazam.readTemplate): a node list, or an older flat list of
        %   steps, each step then the child of the one before.
            file = fullfile(testCase.Plugin, 'Deconvolve', 'templates', name);
            testCase.assertTrue(isfile(file), sprintf('%s is not in the plugin''s templates.', name));
            raw = jsondecode(fileread(file));
            if isfield(raw, 'nodes')
                items = asItems(raw.nodes);
                parents = cellfun(@(it) double(it.parent), items);
            else
                items = asItems(raw.steps);
                parents = (1:numel(items)) - 1;
                parents(1) = -1;
            end
            nodes = struct('transformId', cellfun(@(it) char(it.transformId), items, 'UniformOutput', false), ...
                'params', cellfun(@(it) it.params, items, 'UniformOutput', false), ...
                'parent', num2cell(parents));
        end

        function results = replay(testCase, nodes, file)
        %REPLAY  Every step on its parent's result, as onApplyTemplate walks
        %   the node list, after the toolboxes the steps need are confirmed.
        %   Each dataset carries a cache path in .File, as in the app.
            testCase.requireToolboxesFor(nodes);
            cache = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            raw = loadRecording(file);
            [~, name] = fileparts(file);
            raw.File = fullfile(cache, [name '.mat']);
            results = cell(1, numel(nodes));
            for k = 1:numel(nodes)
                if nodes(k).parent < 1
                    input = raw;
                else
                    input = results{nodes(k).parent};
                end
                results{k} = TransTools.invoke(nodes(k).transformId, input, nodes(k).params);
                results{k}.File = resultCacheFile(input.File, nodes(k).transformId);
            end
        end

        function requireToolboxesFor(testCase, nodes)
            ids = {nodes.transformId};
            testCase.assumeTrue(~isempty(which('pop_loadset')), 'EEGLAB is not on the path.');
            if ismember('Deconvolve', ids)
                testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            end
            if ismember('EyeTracking', ids)
                testCase.assumeTrue(EyeEeg.isAvailable(), 'EYE-EEG is not installed.');
            end
        end

        function file = dataFile(testCase, folder, name)
        %DATAFILE  A recording under the data folder, or the case skips.
            for c = {getenv('ALAKAZAM_DATA'), fullfile(testCase.Alakazam, 'Data'), 'D:\data'}
                file = fullfile(c{1}, folder, name);
                if ~isempty(c{1}) && isfile(file)
                    return;
                end
            end
            testCase.assumeFail(sprintf(['%s is not under %s in ALAKAZAM_DATA, Alakazam''s ' ...
                'Data folder or D:\\data.'], name, folder));
        end
    end
end

% ======================================================================= %
function EEG = loadRecording(file)
%LOADRECORDING  A .set file as WorkSpace.loadSETFile prepares it.
    [folder, name, ext] = fileparts(file);
    EEG = pop_loadset([name ext], folder);
    EEG = eeg_checkset(EEG);
    EEG = recordRawFile(EEG, file);
    EEG.DataType = 'TIMEDOMAIN';
    EEG.DataFormat = inferDataFormat(EEG);
    if strcmpi(EEG.DataFormat, 'CONTINUOUS')
        EEG.times = ((1:EEG.pnts) - 1) / EEG.srate;
    elseif strcmpi(EEG.DataFormat, 'EPOCHED') && ~isempty(EEG.times) && max(abs(EEG.times(:))) < 10
        EEG.times = EEG.times * 1000;
    end
    EEG = deriveBinsFromEpochs(EEG);
end

function items = asItems(list)
%ASITEMS  A decoded JSON list as a cell array, whichever shape it came in.
    if iscell(list)
        items = reshape(list, 1, []);
    else
        items = reshape(num2cell(list), 1, []);
    end
end

function steps = templateSteps(template)
%TEMPLATESTEPS  A template's steps, which older files call nodes.
    if isfield(template, 'steps')
        steps = template.steps;
    else
        steps = template.nodes;
    end
    if iscell(steps)
        steps = [steps{:}];
    end
end

function k = stepIndex(nodes, transformId)
%STEPINDEX  The template's first step of this transformation.
    k = find(strcmp({nodes.transformId}, transformId), 1);
    assert(~isempty(k), 'The template has no %s step.', transformId);
end

function value = peakIn(EEG, binLabel, channel, windowMs)
%PEAKIN  The largest value of one bin's waveform at a channel in a window.
    c = strcmpi({EEG.chanlocs.labels}, channel);
    b = strcmp({EEG.bindesc.label}, binLabel);
    t = EEG.times >= windowMs(1) & EEG.times <= windowMs(2);
    value = max(EEG.data(c, t, b));
end
