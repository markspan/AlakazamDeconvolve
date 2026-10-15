classdef UnfoldTimeBasisTest < matlab.unittest.TestCase
%UNFOLDTIMEBASISTEST  Deconvolve's time basis: uf_timeexpandDesignmat's
%   'method', one beta per sample (the toolbox's default), cubic splines or a
%   Fourier set, with every result read from uf_condense.
%
%   The waveforms are checked against the same model run as an Unfold script
%   runs it; the rest is what Deconvolve does around it: no trials (the
%   toolbox's uf_erpimage cannot make them for such a fit), no basis without
%   the deconvolution, and one beta per sample for options stored before
%   the setting existed.
%
%   Run with: runtests('tests/UnfoldTimeBasisTest.m').
%
%   See also UNFOLD.FITBINS, DECONVOLVE, UNFOLDMODELFITTEST.

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
                     fullfile(root, 'src', 'Support'), fullfile(root, 'src', 'IO'), ...
                     fullfile(plugin, 'tests')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end

    methods (Test)
        function aBasisNeedsTheDeconvolution(testCase)
            EEG = DeconvolveTest.untaggedRecording();
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('timeBasis', 'splines', ...
                'overlapCorrection', false)), 'Alakazam:Unfold:TimeBasisNeedsDeconvolution');
        end

        function trialsAreRefusedWithABasisBeforeFitting(testCase)
            EEG = DeconvolveTest.untaggedRecording();
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('timeBasis', 'fourier', ...
                'output', 'trials')), 'Alakazam:Unfold:BasisTrials');
            testCase.verifySubstring(Unfold.basisTrialsMessage(), 'uf_erpimage');
        end
    end

    methods (Test, TestTags = {'External'})
        function eachBasisGivesTheToolboxsWaveforms(testCase)
        %EACHBASISGIVESTHETOOLBOXSWAVEFORMS  With y ~ 1 a bin's waveform is
        %   its intercept as uf_condense returns it, the betas of the splines
        %   or frequencies turned back into samples: the same model run as a
        %   script gives the same.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            for basis = {{'splines', 20}, {'fourier', 10}}
                [method, count] = basis{1}{:};
                fitted = Deconvolve(EEG, DeconvolveTest.options('timeBasis', method, ...
                    'timeBasisCount', count, 'excludeAroundCuts', false, 'baselineMs', []));
                condensed = UnfoldTimeBasisTest.script(EEG, method, count);

                testCase.verifyEqual(fitted.etc.alz.unfold.timeBasis, method);
                testCase.verifyEqual(fitted.data(:, :, 1), condensed.beta(:, :, 1), 'AbsTol', 1e-6, method);
                testCase.verifyEqual(fitted.data(:, :, 2), condensed.beta(:, :, 2), 'AbsTol', 1e-6, method);
            end
        end

        function storedOptionsFitOneBetaPerSample(testCase)
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            fitted = Deconvolve(DeconvolveTest.untaggedRecording(), DeconvolveTest.options());
            testCase.verifyEqual(fitted.etc.alz.unfold.timeBasis, 'stick');
        end
    end

    methods (Static)
        function condensed = script(EEG, method, count)
        %SCRIPT  EEG as an EEGLAB dataset, with the model Deconvolve fits on
        %   it (DeconvolveTest's bins and the response modelled), expanded by
        %   METHOD with COUNT functions, fitted and condensed as a script does.
            script = eeg_emptyset();
            script.data = EEG.data;
            script.srate = EEG.srate;
            script.nbchan = size(EEG.data, 1);
            script.pnts = size(EEG.data, 2);
            script.trials = 1;
            script.xmax = (script.pnts - 1) / script.srate;
            script.chanlocs = EEG.chanlocs;
            script.event = rmfield(EEG.event, setdiff(fieldnames(EEG.event), {'type', 'latency'}));
            evalc(['script = uf_designmat(script, ''eventtypes'', {''S1'', ''S2'', ''response''}, ' ...
                '''formula'', {''y ~ 1'', ''y ~ 1'', ''y ~ 1''});']);
            evalc(['script = uf_timeexpandDesignmat(script, ''timelimits'', [-0.2 0.8], ' ...
                '''method'', method, ''timeexpandparam'', count);']);
            evalc('script = uf_glmfit(script);');
            condensed = uf_condense(script);
        end
    end
end
