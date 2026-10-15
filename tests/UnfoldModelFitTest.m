classdef UnfoldModelFitTest < matlab.unittest.TestCase
%UNFOLDMODELFITTEST  Deconvolve's model fit: the toolbox's own
%   uf_checkmodelfit, off unless asked, as R2, partial R2 or their
%   cross-validated versions.
%
%   The first check is that the number is the toolbox's, against the same
%   model run as an Unfold script runs it. The rest checks what Deconvolve
%   adds: the terms named by their bins, the folds cut at the code chosen,
%   and a refusal, with the reason, where the toolbox cannot measure a fit.
%
%   Run with: runtests('tests/UnfoldModelFitTest.m').
%
%   See also UNFOLD.FITBINS, DECONVOLVE, OVERLAPCORRECTEDTRIALSTEST.

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
        function aFitNeedsTheDeconvolution(testCase)
            EEG = DeconvolveTest.untaggedRecording();
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('modelFit', 'R2', ...
                'overlapCorrection', false)), 'Alakazam:Unfold:FitNeedsDeconvolution');
        end

        function crossValidationNeedsAFoldEvent(testCase)
            EEG = DeconvolveTest.untaggedRecording();
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('modelFit', 'crossValR2')), ...
                'Alakazam:Unfold:NoFoldEvent');
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('modelFit', 'crossValR2', ...
                'foldEvent', {'pause'})), 'Alakazam:Unfold:NoFoldEvent', 'A code the recording lacks.');
        end
    end

    methods (Test, TestTags = {'External'})
        function theFitIsTheToolboxsOwn(testCase)
        %THEFITISTHETOOLBOXSOWN  The same model as an Unfold script writes
        %   it, and uf_checkmodelfit's R2 of it: Deconvolve's is that.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            [fitted, info] = Deconvolve(EEG, DeconvolveTest.options('modelFit', 'R2', ...
                'excludeAroundCuts', false)); %#ok<ASGLU>
            fit = fitted.etc.alz.unfold.modelFit;

            script = UnfoldModelFitTest.script(EEG);
            r2 = [];
            evalc('r2 = uf_checkmodelfit(script, ''method'', ''R2'');');

            testCase.verifyEqual(fit.method, 'R2');
            testCase.verifyEqual(fit.r2, r2, 'AbsTol', 1e-6);
            testCase.verifyEqual(fit.channels, {'Cz', 'Pz'});
            testCase.verifyTrue(any(contains(fitted.etc.alz.unfold.notes, 'Model fit (uf_checkmodelfit, R2)')));
        end

        function eachTermsPartialFitIsNamedByItsBin(testCase)
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            fitted = Deconvolve(EEG, DeconvolveTest.options('modelFit', 'partialR2'));
            fit = fitted.etc.alz.unfold.modelFit;

            terms = unique({fit.partial.term});
            testCase.verifyEqual(sort(terms), sort({'Frequent (Intercept)', 'Rare (Intercept)', ...
                'response (Intercept)'}));
            testCase.verifyNumElements(fit.partial, 2 * 3, 'Every term on every channel.');
            testCase.verifyTrue(all([fit.partial.partialR2] > 0 & [fit.partial.partialR2] < 1));
            testCase.verifyTrue(any(contains(fitted.etc.alz.unfold.notes, 'partial R2 of each term')));
        end

        function theCollinearityIsTheToolboxsOwn(testCase)
        %THECOLLINEARITYISTHETOOLBOXSOWN  uf_vif of the time-expanded design,
        %   as a script computes it, every column; and each term's largest
        %   named by its bin.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            script = UnfoldModelFitTest.script(EEG);
            fitted = Deconvolve(EEG, DeconvolveTest.options('collinearity', 'expanded', ...
                'excludeAroundCuts', false));
            vif = fitted.etc.alz.unfold.collinearity;
            expected = [];
            evalc('expected = uf_vif(script.unfold.Xdc);');

            testCase.verifyEqual(vif.perColumn, reshape(expected, 1, []), 'RelTol', 1e-9);
            testCase.verifyEqual(sort(vif.terms), sort({'Frequent (Intercept)', ...
                'Rare (Intercept)', 'response (Intercept)'}));
            testCase.verifyEqual(vif.largest, cellfun(@(t) max(vif.perColumn(strcmp(vif.columns, t))), ...
                vif.terms));
            testCase.verifyTrue(any(contains(fitted.etc.alz.unfold.notes, 'Collinearity (uf_vif')));
            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options('collinearity', 'expanded', ...
                'overlapCorrection', false)), 'Alakazam:Unfold:VifNeedsDeconvolution');
        end

        function crossValidationCutsTheFoldsAtTheCode(testCase)
        %CROSSVALIDATIONCUTSTHEFOLDSATTHECODE  Two cuts, with a window left
        %   out around each, make three folds; fitted across the cuts, the
        %   model has responses there and the toolbox refuses, as the message
        %   says.
            testCase.assumeTrue(Unfold.isAvailable(), 'The Unfold toolbox is not installed.');
            EEG = DeconvolveTest.untaggedRecording();
            for cut = [7000 14000]           % the recording is 20000 samples long
                EEG.event(end + 1) = struct('type', 'boundary', 'latency', cut);
            end
            [~, order] = sort([EEG.event.latency]);
            EEG.event = EEG.event(order);
            args = {'modelFit', 'crossValR2', 'foldEvent', {'boundary'}};

            fitted = Deconvolve(EEG, DeconvolveTest.options(args{:}, 'excludeAroundCuts', true));
            fit = fitted.etc.alz.unfold.modelFit;
            testCase.verifyEqual(size(fit.r2Folds), [2 3]);
            testCase.verifyEqual(fit.r2, mean(fit.r2Folds, 2), 'AbsTol', 1e-12);

            testCase.verifyError(@() Deconvolve(EEG, DeconvolveTest.options(args{:}, ...
                'excludeAroundCuts', false)), 'Alakazam:Unfold:FoldOverlap');
        end
    end

    methods (Static)
        function script = script(EEG)
        %SCRIPT  EEG as an EEGLAB dataset, with the model Deconvolve fits on it
        %   (DeconvolveTest's bins and the response modelled), built and fitted
        %   as an Unfold script does.
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
            evalc('script = uf_timeexpandDesignmat(script, ''timelimits'', [-0.2 0.8]);');
            evalc('script = uf_glmfit(script);');
        end
    end
end
