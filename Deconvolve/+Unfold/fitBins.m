function [EEG, info] = fitBins(input, varargin)
%FITBINS  One overlap-corrected waveform per bin, from continuous data.
%   [EEG, INFO] = Unfold.fitBins(INPUT) fits the model Unfold.binModel builds
%   from DefineBins' bins and returns it in the shape Average returns:
%   DataFormat "Averaged", channels x samples x bins, the same EEG.bindesc,
%   so Measure, ScalpDistribution, GrandAverage and the reports read it
%   without knowing it came from a regression. INFO carries the model, the
%   notes and what was excluded.
%
%   THE ALTERNATIVE TO EPOCH-AND-AVERAGE, not a replacement for it. Averaging
%   assumes each epoch holds the response to its own event and nothing else.
%   Where events follow each other faster than the response decays, that
%   assumption is wrong and the average of a bin carries a smear of whatever
%   came before and after it. This fits all the bins at once against the
%   whole continuous recording, so a sample explained by two events is
%   attributed to both, and what comes out is the response to each bin with
%   the others' overlap removed (Ehinger & Dimigen, 2019; Dimigen & Ehinger,
%   2021, J Vis 21(1):3).
%
%   IT NEEDS DATA WITHOUT A DC OFFSET, and refuses data that has one. A
%   time-expanded design has no constant term, so whatever mean voltage a
%   channel carries has to be explained by the event responses, and the
%   betas come back in the thousands of microvolts. Epoch-and-average never
%   shows this because Baseline subtracts the offset epoch by epoch, so a
%   DC-coupled recording (BioSemi and the like, tens of millivolts of
%   electrode offset) averages perfectly and deconvolves into nonsense. The
%   remedy belongs upstream, in DCDetrend or a high-pass Filter, not here:
%   silently de-meaning someone's data would hide the one fact that decides
%   whether the numbers mean anything. Without overlap correction (below)
%   the offset is no obstacle: each epoch's sample has its own intercept.
%
%   WITHOUT OVERLAP CORRECTION (OverlapCorrection false) the same design is
%   fitted the toolbox's other way, a mass-univariate regression on epochs:
%   uf_epoch cuts the events' windows, leaving out every one that touches a
%   stretch the artefact scan marked, spans a cut or runs off the recording,
%   and uf_glmfit_nodc fits the design at every sample of them. With y ~ 1
%   a bin's waveform is then the mean of its epochs, as Average gives it,
%   but with Unfold's rounding of each latency to the nearest sample; with a
%   formula it is the same regression as the deconvolution, without the
%   neighbours' overlap taken out. It is the toolbox's own comparison of the
%   two (its tutorial on deconvolved and not deconvolved results).
%
%   Options; where the toolbox has a default, it is the default here:
%     WindowMs           [-200 800], the response window per event. It should
%                        cover the whole response, including anything that
%                        precedes the event (a saccade's own motor activity
%                        precedes the fixation it produces).
%     BaselineMs         [] (the default) leaves the betas exactly as the
%                        solver returned them, as uf_condense does. A window
%                        [start stop] in ms, or 'pre-event' (the part of
%                        WindowMs before the event), baseline-corrects the
%                        fitted waveforms over it. A beta is a regression
%                        coefficient, so its zero is wherever the model put
%                        it, and reading one beside an average (which
%                        Baseline has corrected) means correcting it too. The
%                        window is taken as Unfold's uf_plotParam takes its
%                        baseline: the samples from start up to, but not
%                        including, stop (baselineSamples).
%     OtherEvents        which event codes in no bin to fit as nuisance and
%                        drop from the result: 'all' (default), a cellstr of
%                        codes, or {} for none. See Unfold.binModel for why
%                        this is a choice per code.
%     Formulas           each bin's formula in Unfold's notation, as a struct
%                        array of .bin (label) and .formula; 'y ~ 1' for a
%                        bin without one (see Unfold.binModel).
%     Covariates         the older option: event fields added as linear
%                        terms to every bin without a formula of its own.
%                        Unfold.eventCovariates lists what a dataset offers.
%     EvaluateAt         for Output 'terms': where each continuous or spline
%                        term is evaluated, as text, "sac_amplitude = 0.5 1 2;
%                        rating = 1 5". A term not named is evaluated at ten
%                        quantiles of its own values (uf_predictContinuous's
%                        own default), which differ from one recording to
%                        the next, so name the values when the terms are to
%                        be combined across subjects.
%     ArtifactThresholdUv, ArtifactWindowMs, ArtifactStepMs
%                        150 uV in a 2000 ms window stepped by 100 ms, which
%                        are the toolbox's own defaults and the paper's
%                        larger dataset. Bad intervals are EXCLUDED BY
%                        ZEROING the design matrix there, not by cutting
%                        epochs: the surrounding data still contributes, and
%                        an event near a bad stretch does not lose its whole
%                        epoch. Set the threshold to 0 to skip detection.
%                        What happens around the recording's own boundaries
%                        (cuts) is ExcludeAroundCuts.
%                        The threshold is PEAK-TO-PEAK within each moving
%                        window: uf_continuousArtifactDetect hands it to
%                        ERPLAB's basicrap without a threshold type, whose
%                        default is 'peak-to-peak'. Its own header says
%                        "[-lim +lim] is marked", which reads like an
%                        absolute limit and is not what the code does
%                        (checked in basicrap.m, unfold 1.3.1). So a standing
%                        offset does not trip it; drifts, blinks and, in
%                        free viewing, eye movements within a window do.
%                        THE TOOLBOX IS CALLED AS A SCRIPT CALLS IT, and its
%                        results are kept as they are (1.3.1, unchanged in
%                        its current version). The scan is handed the
%                        dataset's own events, boundaries included, so
%                        ERPLAB's scan looks at each stretch between cuts on
%                        its own (forArtefactScan). Its windows start one
%                        sample in and stop when the next would run off the
%                        end, so the first sample and less than one step at
%                        the end of a stretch are not looked at. With
%                        ExcludeAroundCuts, the scan's marks and the zones
%                        around cuts are joined by
%                        uf_combineWinrej, which ends a joined stretch where
%                        the later-starting one ends; where one lies wholly
%                        inside another, the joined stretch ends with the
%                        inner one. With the scan between cuts a marked
%                        window never spans a cut, so the one case left is a
%                        marked stretch lying wholly inside a cut's zone,
%                        on one side of it: the zone then ends where that
%                        stretch does. It needs an artefact window shorter
%                        than the longer side of the response window (the
%                        zone reaches that far either side of a cut), so
%                        with the defaults, 2000 ms against 800 ms, it
%                        cannot happen.
%     ExcludeAroundCuts  false (the default, as the toolbox does nothing about
%                        them): the data either side of a cut in the
%                        recording (a boundary event) is fitted as it is,
%                        across the join. true leaves the longer side of
%                        WindowMs out on each side of every cut, through
%                        uf_continuousArtifactExclude with the scan's marks,
%                        so no event's fitted window spans two moments that
%                        were never adjacent. Either way a trial whose window
%                        spans a cut is left out, as EEGLAB's pop_epoch
%                        leaves out an epoch that holds a boundary.
%     SolverIterations   400, the toolbox's own limit on the lsmr solver's
%                        iterations (uf_glmfit's lsmriterations). A fit that
%                        reaches it says so in the notes; raising it lets a
%                        slow but sound fit finish, while a nearly collinear
%                        design needs a different model instead.
%     Solver             'default' (the default): each fitting function's
%                        own, lsmr for uf_glmfit and pinv for uf_glmfit_nodc;
%                        'matlab', MATLAB's own exact solver, which needs no
%                        iteration limit but, as the toolbox warns, can need a
%                        great deal of memory for a long recording; or
%                        'glmnet', a regularised fit (lasso, ridge or elastic
%                        net) whose strength is chosen by cross-validation,
%                        so its betas are shrunk towards zero and are not the
%                        least-squares waveforms; with overlap correction it
%                        gives no trials (Unfold.glmnetTrialsMessage says
%                        why), and the two together are refused. The toolbox's other two,
%                        'par-lsmr' and 'pinv' for a deconvolution, are not
%                        offered: its own help calls the first no faster and
%                        not recommended, and the second unstable.
%     GlmnetAlpha        for Solver 'glmnet': 1 (the default, the toolbox's)
%                        is lasso, 0 ridge, between them elastic net.
%     OverlapCorrection  true (the default): deconvolution. false: a
%                        regression on epochs (see above).
%     Channels           which channels the artefact scan looks at: [] (the
%                        default) for every channel, as
%                        uf_continuousArtifactDetect scans by default;
%                        'scalp' for the scalp EEG only (eegChannelMask); or
%                        channel indices. An EOG channel's range is several
%                        times the EEG's, so scanning it can mark much of the
%                        recording bad and take the bins' data with it, which
%                        is what 'scalp' is for.
%     CodingSchema       how a factor is coded, uf_designmat's
%                        'codingschema': 'reference' (the default, the
%                        toolbox's: treatment coding, the intercept is the
%                        first level in sorted order) or 'effects' (the
%                        intercept is the mean over the levels; uf_designmat
%                        then also centres each continuous term on its event
%                        type's mean). A bin's waveform is the same under
%                        either, since it is the prediction at the bin's own
%                        mix of levels; the terms differ. Effects coding with
%                        an interaction of a factor and a continuous term is
%                        refused (requireRebuildableInteractions).
%     EvaluateSpacing    for Output 'terms': where a continuous or spline term
%                        not named in EvaluateAt is evaluated,
%                        uf_predictContinuous's 'auto_method': 'quantile' (the
%                        default, the toolbox's: ten quantiles of its values),
%                        'linear' (ten values evenly spaced from its minimum
%                        to its maximum) or 'average' (its mean alone).
%     TimeBasis          how each response is expanded over the window,
%                        uf_timeexpandDesignmat's 'method': 'stick' (the
%                        default, the toolbox's: one beta per sample),
%                        'splines' (TimeBasisCount cubic splines across the
%                        window, so neighbouring samples share their betas
%                        and the waveforms come out smooth) or 'fourier' (the
%                        first TimeBasisCount frequencies). uf_condense turns
%                        the betas back into samples, and every result is
%                        read from it. With overlap correction only, and
%                        without trials: uf_erpimage cannot make them for
%                        such a fit (Unfold.basisTrialsMessage says why).
%     TimeBasisCount     how many splines or frequencies: 30, the toolbox's
%                        default ('timeexpandparam').
%     Collinearity       the variance inflation factor of every column of
%                        the time-expanded design as fitted (Xdc), by the
%                        toolbox's own uf_vif: 'none' (the default: the
%                        toolbox computes none unless asked) or 'expanded':
%                        what the overlap does to telling the bins apart.
%                        With overlap correction only; its memory grows with
%                        the square of the design's columns. The largest of
%                        each term is in the notes, and every column's in
%                        INFO.collinearity (collinearity). The toolbox notes
%                        that splines give high values that are usually no
%                        problem. uf_vif of the design X alone, its default
%                        when handed a dataset, is not offered: every
%                        modelled event belongs to one event type, so the
%                        types' intercepts add up to a constant and their
%                        factors are infinite (in a script the event list's
%                        unmodelled events are rows of zeros, which keep
%                        the matrix just short of singular).
%     ModelFit           how well the model explains the recording, by the
%                        toolbox's own uf_checkmodelfit: 'none' (the default:
%                        the toolbox computes no fit unless asked), 'R2' (the
%                        coefficient of determination per channel, over the
%                        samples the model covers), 'partialR2' (what each
%                        term adds to it, refitting the model once without
%                        each), 'crossValR2' or 'crossValpartialR2' (the same,
%                        each fold fitted on the rest and tested on itself).
%                        Refits use uf_glmfit's own default solver. With
%                        overlap correction only: it needs the time-expanded
%                        design. The result is in the notes and, in full, in
%                        INFO.modelFit (modelFit).
%     FoldEvent          for the cross-validated fits: the event code(s) of
%                        the recording the folds are cut at (uf_cv_getFolds),
%                        such as a break between blocks; the toolbox refuses a
%                        fold cut where the model has a response.
%     Marginal           for Output 'terms': 'MEM' (the default, as
%                        uf_addmarginal's) or 'AME' (see below).
%     MissingValues      what happens to an event whose formula names a
%                        field it has no number for: one of the toolbox's
%                        uf_imputeMissing methods, 'median' (the default, as
%                        the toolbox's), 'mean', 'marginal' (a random draw
%                        from the others) or 'drop' (the event left out of
%                        the model, so not overlap-corrected either), applied
%                        between the design and its time expansion, as
%                        toolboxWorkflow.rst places it; or 'refuse'. A factor
%                        without a level is refused whatever the choice:
%                        uf_designmat cannot build it, and uf_imputeMissing
%                        fills in numbers. TOOLBOX BEHAVIOURS KEPT AS THEY
%                        ARE (1.3.1): it fills a spline's basis column by
%                        column, so the row it makes is not the spline at any
%                        one value; its warning that more than 5% are
%                        missing counts against every event in the model, not
%                        the type's own; and 'marginal' fails when two
%                        predictors miss different numbers of values, since
%                        it reuses one list of drawn values from predictor to
%                        predictor (refused here with that reason). A
%                        waveform per bin is held at the values the events
%                        have, not at filled-in ones (Unfold.binModel).
%     Output             'average' (the default): one waveform per bin, as
%                        described above. 'trials': one OVERLAP-CORRECTED
%                        TRIAL per binned event instead, in the epoched shape
%                        DefineBins cuts (DataFormat 'EPOCHED', EEG.epoch with
%                        each trial's .bini, bindesc.trials), so EpochView,
%                        Average and the data-quality report take it like any
%                        epoch node. Each trial is the recording around its
%                        event with every other event's fitted response
%                        subtracted: what the toolbox's own uf_erpimage
%                        makes, its modelled responses plus residuals (see
%                        TrialsMethod for how).
%                        At the least-squares fit a bin's residuals average
%                        to zero at every lag, so averaging its trials gives
%                        back its fitted waveform. A trial whose window
%                        touches a stretch the model left out, or runs off
%                        the recording, is dropped by uf_epoch, since nothing
%                        was subtracted there; the count is in the
%                        provenance. Without overlap
%                        correction the trials are the epochs as uf_epoch
%                        cut them, nothing subtracted. 'terms': one waveform
%                        per model term instead (see below).
%     TrialsMethod       how the overlap-corrected trials are made:
%                        'fast' (the default, the user's choice): by Alakazam
%                        (fastTrials), the trials uf_erpimage makes from the
%                        same model, computed only at the trials' own
%                        samples and for every channel at once: the
%                        recording less the model's prediction there, plus
%                        the event's own fitted response, kept and cut as
%                        uf_epoch and pop_epoch keep and cut them. It agrees
%                        with the toolbox's to rounding error
%                        (FastTrialsTest), but it is Alakazam's arithmetic
%                        rather than the toolbox's call, so a change to
%                        uf_erpimage would not reach it. 'toolbox': by
%                        uf_erpimage, channel by channel, as above. Each call
%                        epochs the whole recording around every modelled
%                        event, three times, so on a long recording with
%                        many events this takes long.
%     Parallel           'auto' (the default), 'off' or 'on': whether the
%                        steps done one channel at a time run on parallel
%                        workers (Unfold.workerPool): the fit with the
%                        default solver, each worker running the toolbox's
%                        uf_glmfit on its channels, and the toolbox's trials.
%                        'auto' uses workers only for a large job and only
%                        where Alakazam's settings allow parallel processing.
%                        The results are the same either way, to rounding
%                        error (a worker computes on one thread); not stored
%                        with the options, since it is no part of the
%                        analysis.
%
%   A BIN'S WAVEFORM is the model's prediction with every continuous and
%   spline term at its average over ALL the modelled events that carry it
%   (Unfold.binModel's plan.pooled), the same for every bin, and every factor
%   at the bin's own mix of levels. With y ~ 1 that is the intercept, as it
%   always was. With a covariate it is each bin's response at the same value
%   of it, so a covariate whose values differ between bins (a bigger saccade
%   in one condition, a slower response in another) is held constant rather
%   than showing up as a difference between them. That is what makes it a
%   control. For a spline the average is taken over its basis, which is
%   Unfold's own average marginal effect and is right for a circular one
%   too, where a mean angle would not be. The spline bases are evaluated by
%   the toolbox's own function (EEG.unfold.splines{s}.splinefunction with its
%   knots and removedSplineIdx). A spline that is not circular is only
%   averaged over the values every bin using it shares (plan.pooled's
%   .common), since outside a bin's own values its spline is extrapolating;
%   bins that share none are each averaged over their own, and the notes say
%   which of the two happened. A 2D spline (2dspl) is averaged the same way
%   over the pooled pairs of its two fields. An interaction with a
%   continuous term (cat(side) * rt) has that term at the pooled mean too,
%   like the term's own column. Where every bin carries the same values, as
%   with y ~ 1 or a covariate one bin alone uses, this is also exactly what
%   the bin's overlap-corrected trials average to.
%
%   BINS THAT SHARE EVENTS: each set of bins the events can be in is one
%   event type of the model (Unfold.binModel), and a bin's waveform is the
%   average of those types' waveforms weighted by how many of its events
%   each holds. That is Average's own rule, an event counting in every bin it
%   is in, so "Rare" nested in "All stimuli" is the rare events' response,
%   not its difference from the others.
%
%   THE TERMS OUTPUT is Unfold's own view of the model: uf_condense, then
%   uf_predictContinuous (each continuous and spline term at the EvaluateAt
%   values, or ten quantiles), then uf_addmarginal, which makes every term a
%   whole waveform with the model's other terms added in. With Marginal
%   'MEM', the toolbox's default, those other terms are at their mean value;
%   with 'AME' they are their average marginal effect, a spline averaged
%   over its events' own values, as a bin's waveform has it. The two differ
%   only for a spline, and most for a circular one, whose mean angle can be
%   a direction no event had. So the intercept is the response at each
%   factor's reference level, a factor level is the response at that level,
%   and a spline at a value is the response at that value. One "bin" per
%   term of every modelled set of
%   bins, labelled "<bin>: <term>" ("<bin>: x = 1, z = 2" for a 2D spline),
%   in Average's shape, so Measure, GrandAverage and the reports read them.
%   A TOOLBOX BEHAVIOUR KEPT AS IT IS, by design, and the toolbox's own
%   warning says so: uf_predictContinuous takes a continuous term's
%   quantiles, and uf_addmarginal its mean, over the events whose value is
%   not exactly zero, since zero is what the other event types' rows of X
%   hold. A covariate whose own values include zero (a rating of 0) is
%   summarised without them; name the values in EvaluateAt where that
%   matters.
%
%   WHAT THE RESULT DOES NOT CARRY is a standard error of its own. Average's
%   .stErr is the spread of the trials that went into a mean, and a
%   regression coefficient has no trials to spread; it has a standard error
%   from the model's residuals, which the lsmr solver used here does not
%   return. Rather than leave the field out, which AverageView indexes
%   unguarded, it is filled with zeros, exactly as that view does for any
%   series with no error available, and the provenance in EEG.etc.alz.unfold
%   records that there is none. Across-subject error bars in the reports are
%   unaffected: those come from the grand average, not from this.
%
%   See also UNFOLD.BINMODEL, UNFOLD.ENSURE, AVERAGE, DEFINEBINS.
    FIT_METHODS = {'none', 'R2', 'partialR2', 'crossValR2', 'crossValpartialR2'};
    parsed = inputParser();
    parsed.addParameter('WindowMs', [-200 800], @(v) isnumeric(v) && numel(v) == 2 && v(1) < v(2));
    parsed.addParameter('BaselineMs', [], ...
        @(v) isempty(v) || (ischar(v) || isstring(v)) || (isnumeric(v) && numel(v) == 2 && v(1) < v(2)));
    parsed.addParameter('OtherEvents', 'all', @(v) isempty(v) || ischar(v) || iscellstr(v) || isstring(v));
    parsed.addParameter('Covariates', {}, @(v) isempty(v) || iscellstr(v) || isstring(v));
    parsed.addParameter('Formulas', [], @(v) isempty(v) || isstruct(v));
    parsed.addParameter('EvaluateAt', '', @(v) ischar(v) || isstring(v));
    parsed.addParameter('ArtifactThresholdUv', 150, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    parsed.addParameter('ArtifactWindowMs', 2000, @(v) isnumeric(v) && isscalar(v) && v > 0);
    parsed.addParameter('ArtifactStepMs', 100, @(v) isnumeric(v) && isscalar(v) && v > 0);
    parsed.addParameter('Channels', [], @(v) isnumeric(v) || ((ischar(v) || isstring(v)) && strcmpi(v, 'scalp')));
    parsed.addParameter('Marginal', 'MEM', ...
        @(v) (ischar(v) || isstring(v)) && any(strcmpi(char(string(v)), {'MEM', 'AME'})));
    parsed.addParameter('MissingValues', 'median', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'refuse', 'median', 'mean', 'marginal', 'drop'})));
    parsed.addParameter('SolverIterations', 400, @(v) isnumeric(v) && isscalar(v) && v >= 1 && v == round(v));
    parsed.addParameter('Solver', 'default', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'default', 'matlab', 'glmnet'})));
    parsed.addParameter('GlmnetAlpha', 1, @(v) isnumeric(v) && isscalar(v) && v >= 0 && v <= 1);
    parsed.addParameter('OverlapCorrection', true, @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
    parsed.addParameter('ExcludeAroundCuts', false, @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
    parsed.addParameter('CodingSchema', 'reference', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'reference', 'effects'})));
    parsed.addParameter('EvaluateSpacing', 'quantile', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'quantile', 'linear', 'average'})));
    parsed.addParameter('ModelFit', 'none', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), FIT_METHODS)));
    parsed.addParameter('FoldEvent', {}, @(v) isempty(v) || ischar(v) || iscellstr(v) || isstring(v));
    parsed.addParameter('Collinearity', 'none', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'none', 'expanded'})));
    parsed.addParameter('TimeBasis', 'stick', @(v) (ischar(v) || isstring(v)) && ...
        any(strcmpi(char(string(v)), {'stick', 'splines', 'fourier'})));
    parsed.addParameter('TimeBasisCount', 30, @(v) isnumeric(v) && isscalar(v) && v >= 1 && v == round(v));
    parsed.addParameter('Output', 'average', ...
        @(v) (ischar(v) || isstring(v)) && any(strcmpi(char(string(v)), {'average', 'trials', 'terms'})));
    parsed.addParameter('TrialsMethod', 'fast', ...
        @(v) (ischar(v) || isstring(v)) && any(strcmpi(char(string(v)), {'toolbox', 'fast'})));
    parsed.addParameter('Parallel', 'auto', ...
        @(v) (ischar(v) || isstring(v)) && any(strcmpi(char(string(v)), {'auto', 'off', 'on'})));
    parsed.parse(varargin{:});
    opts = parsed.Results;
    opts.TrialsMethod = lower(char(string(opts.TrialsMethod)));
    opts.Parallel = lower(char(string(opts.Parallel)));
    opts.BaselineMs = resolveBaseline(opts.BaselineMs, opts.WindowMs);
    opts.Marginal = upper(char(string(opts.Marginal)));
    opts.MissingValues = lower(char(string(opts.MissingValues)));
    opts.Solver = lower(char(string(opts.Solver)));
    opts.OverlapCorrection = logical(opts.OverlapCorrection);
    opts.ExcludeAroundCuts = logical(opts.ExcludeAroundCuts);
    opts.CodingSchema = lower(char(string(opts.CodingSchema)));
    opts.EvaluateSpacing = lower(char(string(opts.EvaluateSpacing)));
    opts.ModelFit = FIT_METHODS{strcmpi(FIT_METHODS, char(string(opts.ModelFit)))};
    opts.FoldEvent = reshape(cellstr(string(opts.FoldEvent)), 1, []);
    requireFitPossible(opts, input);
    opts.TimeBasis = lower(char(string(opts.TimeBasis)));
    if ~strcmp(opts.TimeBasis, 'stick') && ~opts.OverlapCorrection
        throw(MException('Alakazam:Unfold:TimeBasisNeedsDeconvolution', '%s', sprintf([ ...
            'A time basis (%s) shapes the deconvolution''s time expansion, and a regression on ' ...
            'epochs has none: it fits every sample on its own. Would you switch Overlap ' ...
            'correction on, or keep one beta per sample?'], opts.TimeBasis)));
    end
    opts.Collinearity = lower(char(string(opts.Collinearity)));
    if strcmp(opts.Collinearity, 'expanded') && ~opts.OverlapCorrection
        throw(MException('Alakazam:Unfold:VifNeedsDeconvolution', '%s', ['The collinearity ' ...
            'of the time-expanded design needs that design, which only the deconvolution ' ...
            'builds. Would you switch Overlap correction on, or leave the check out?']));
    end
    if ~strcmp(opts.TimeBasis, 'stick') && strcmpi(char(string(opts.Output)), 'trials')
        throw(MException('Alakazam:Unfold:BasisTrials', '%s', Unfold.basisTrialsMessage()));
    end
    if opts.OverlapCorrection && strcmp(opts.Solver, 'glmnet') ...
            && strcmpi(char(string(opts.Output)), 'trials')
        throw(MException('Alakazam:Unfold:GlmnetTrials', '%s', Unfold.glmnetTrialsMessage()));
    end

    % Only the time-expanded design lacks a constant term; a regression on
    % epochs has an intercept at every sample, which takes up an offset as
    % an average does.
    if opts.OverlapCorrection
        requireCentredData(input);
    end
    plan = Unfold.binModel(input, 'OtherEvents', opts.OtherEvents, 'Covariates', opts.Covariates, ...
        'Formulas', opts.Formulas, 'MissingValues', opts.MissingValues, ...
        'OverlapCorrection', opts.OverlapCorrection, 'CodingSchema', opts.CodingSchema);
    if isempty(plan.eventTypes)
        throw(MException('Alakazam:Unfold:NothingToFit', ...
            ['None of this dataset''s bins hold any events, so there is no model to fit. ' ...
             'Check the bin definitions against the events the recording actually has.']));
    end
    Unfold.ensure('Deconvolution (rERP)');

    srate = double(input.srate);

    % Events in no bin that keep a near-constant lag to another modelled type
    % make the design nearly collinear, the usual reason the solver does not
    % converge; they are named in the notes before the fit, not only after it.
    % Only a deconvolution has to tell them apart.
    locked = struct('code', {}, 'other', {}, 'lagMs', {}, 'sdMs', {}, 'note', {});
    if opts.OverlapCorrection
        locked = Unfold.timeLockedEvents(plan, srate, opts.WindowMs);
        plan.notes = [plan.notes, {locked.note}];
    end

    % 1. The design: one event type per bin, each with its own formula (see
    %    Unfold.binModel).
    work = Unfold.designMatrix(input, plan);
    requireRebuildableInteractions(work.unfold, plan);

    % Missing values, filled in or their events left out by the toolbox's
    % own uf_imputeMissing, before the time expansion, as the workflow has it.
    [work, plan] = imputeMissing(work, plan);

    % Every channel unless asked otherwise, as the toolbox scans by default;
    % 'scalp' leaves the peripheral channels out (see this file's header).
    channels = opts.Channels;
    if ischar(channels) || isstring(channels)
        channels = scalpChannels(input);
    elseif isempty(channels)
        channels = 1:size(input.data, 1);
    end

    % 2 to 4: the time expansion, what not to model, and the fit; or, without
    % overlap correction, the epochs and a regression on them.
    % 5. Condensing: every result reads the betas uf_condense returns, with
    %    the time basis applied, rather than the solver's raw ones.
    if opts.OverlapCorrection
        [work, plan, excluded] = deconvolve(work, plan, input, opts, channels, srate, locked);
        if ~strcmp(opts.ModelFit, 'none')
            [fit, note] = modelFit(work, input, plan, opts);
            plan.notes{end + 1} = note;
        end
        condensed = uf_condense(work);
        model = condensed.unfold;
        model.beta = byColumn(condensed, 'beta');
    else
        [work, plan, excluded, epoched] = regressOnEpochs(work, plan, input, opts, channels);
        condensed = uf_condense(work);
        model = condensed.unfold;
        model.beta = byColumn(condensed, 'beta_nodc');
        X = zeros(numel(plan.events), size(model.X, 2));   % X's rows back on the plan's events
        X(epoched, :) = model.X;
        model.X = X;
    end

    % The waveforms are checked in every case: a fit that produced no
    % numbers produces no trials or terms either, and says so the same way.
    [waveforms, times, evaluationNotes] = fittedWaveforms(input, plan, model, condensed.times, opts);
    plan.notes = [plan.notes, evaluationNotes];
    if ~strcmp(opts.Collinearity, 'none')
        [vif, note] = collinearity(work, plan);
        plan.notes{end + 1} = note;
    end
    info = modelInfo(input, plan, opts, excluded, srate);
    info.fitWorkers = TransTools.FieldOr(plan, 'fitWorkers', 1);
    if exist('fit', 'var')
        info.modelFit = fit;
    end
    if exist('vif', 'var')
        info.collinearity = vif;
    end
    switch lower(char(string(opts.Output)))
        case 'trials'
            if opts.OverlapCorrection
                [EEG, info] = packageTrials(input, plan, work, info, opts, excluded, srate, times);
            else
                [EEG, info] = packageEpochs(input, plan, work, epoched, info, opts, srate, times);
            end
        case 'terms'
            [EEG, info] = packageTerms(input, plan, condensed, info, opts, times);
        otherwise
            EEG = package(input, plan, waveforms, times);
    end
    EEG = recordInfo(EEG, info, plan);
end

% ======================================================================= %
function [work, plan, excluded] = deconvolve(work, plan, input, opts, channels, srate, locked)
%DECONVOLVE  Steps 2 to 4 with overlap correction: the time expansion, the
%   stretches left out of the model, and the fit by uf_glmfit.
    % 2. Time expansion: the design matrix gains one column per predictor per
    %    time point in the window, which is what makes the fit a
    %    deconvolution rather than a regression on epochs.
    work = uf_timeexpandDesignmat(work, 'timelimits', opts.WindowMs / 1000, ...
        'method', opts.TimeBasis, 'timeexpandparam', opts.TimeBasisCount);

    % 3. What not to model: artifacts by the paper's method, and, when asked
    %    (ExcludeAroundCuts), the recording's own cuts. A boundary is a join,
    %    so the samples around it are two different moments stitched together
    %    and an event window that spans one is fitted across a discontinuity.
    %    The toolbox does nothing about them, which is the default; the zone
    %    is joined to the scan's marks through the toolbox's own combiner.
    excluded = zeros(0, 2);
    if opts.ExcludeAroundCuts
        excluded = boundaryIntervals(input, opts.WindowMs, srate);
    end
    detected = artefactScan(work, input, opts, channels);
    if ~isempty(detected)
        excluded = combineIntervals(excluded, detected);
    end
    if ~isempty(excluded)
        % Checked before the fit, not after: a design with most of its rows
        % zeroed takes minutes to not converge, and the answer it then
        % gives is unusable anyway. Refusing here costs the user the wait
        % and tells them which knob to turn.
        requireEnoughDataLeft(excluded, size(input.data, 2), opts, numel(channels));
        work = uf_continuousArtifactExclude(work, 'winrej', excluded);
    end

    % 4. The fit, by the toolbox's default solver (lsmr) with its iteration
    %    limit as a setting, or by MATLAB's own exact solver, or regularised
    %    by glmnet. lsmr warns rather than fails when it runs out of
    %    iterations, and an under-converged fit looks like a result, so the
    %    warning is caught and carried into the notes instead of scrolling
    %    past in the log.
    method = opts.Solver;
    if strcmp(method, 'default')
        method = 'lsmr';
    end
    % lsmr solves one channel at a time on a design every channel shares, so
    % a long fit can be spread over parallel workers (fitOnWorkers); the
    % exact solver and glmnet solve all channels together and are not.
    pool = [];
    plan.fitWorkers = 1;
    if strcmp(method, 'lsmr')
        nch = size(work.data, 1);
        bytes = bytesOf(rmfield(work, 'data')) + 6 * 8 * size(work.data, 2) + 1.5e9;
        [pool, plan.fitWorkers] = Unfold.workerPool(nch, nnz(work.unfold.Xdc) * nch, 5e7, ...
            opts.Parallel, bytes);
    end
    if isempty(pool)
        lastwarn('');
        work = uf_glmfit(work, 'method', method, 'lsmriterations', opts.SolverIterations, ...
            'glmnetalpha', opts.GlmnetAlpha);
        [solverWarning, ~] = lastwarn();
        solverWarnings = {solverWarning};
    else
        [work, solverWarnings] = fitOnWorkers(work, pool, plan.fitWorkers, method, opts);
    end
    if strcmp(method, 'lsmr') && any(contains(lower(solverWarnings), 'did not converge'))
        plan.notes{end + 1} = sprintf(['The solver ran out of iterations (Solver iterations, %d) ' ...
            'before it converged, for at least one channel, so these waveforms are an unfinished ' ...
            'estimate. That usually means the design is close to collinear (bins whose events ' ...
            'keep a near-constant lag) or that a lot of the data was excluded as artefact; ' ...
            'where it is neither, raising Solver iterations lets it finish.'], opts.SolverIterations);
        if ~isempty(locked)
            plan.notes{end} = sprintf('%s The likely cause here: %s.', plan.notes{end}, ...
                strjoin(arrayfun(@(p) sprintf('"%s" locked to "%s"', p.code, p.other), ...
                locked, 'UniformOutput', false), ', '));
        end
    end
end

function [work, warned] = fitOnWorkers(work, pool, workers, method, opts)
%FITONWORKERS  The fit, channel by channel on parallel workers: each worker
%   runs the toolbox's own uf_glmfit, with the same solver and limits, on
%   the dataset holding one channel. The solver takes one channel at a time
%   anyway, on the design every channel shares (the rows it drops and the
%   scaling it applies depend on the design alone), so the betas are the
%   ones a single call for every channel gives, to rounding error: a worker
%   computes on one thread, and sums some terms in another order (about
%   1e-12 of the values, FastTrialsTest). What each worker warned is
%   returned in WARNED, so a channel that ran out of iterations still says
%   so. Progress is the channels done.
    nch = size(work.data, 1);
    shared = parallel.pool.Constant(rmfield(work, 'data'));
    futures = parallel.FevalFuture.empty;
    for ch = 1:nch
        futures(ch) = parfeval(pool, @fitChannel, 3, shared, work.data(ch, :), method, ...
            opts.SolverIterations, opts.GlmnetAlpha);
    end
    stopTheRest = onCleanup(@() cancel(futures));
    betas = cell(nch, 1);
    customs = cell(nch, 1);
    warned = cell(1, nch);
    TransTools.BusyGate('progress', 0);
    for k = 1:nch
        [ch, beta, custom, said] = fetchNextOrCause(futures);
        betas{ch} = beta;
        customs{ch} = custom;
        warned{ch} = said;
        TransTools.BusyGate('message', sprintf('Fitting the model on %d workers: %d of %d channels...', ...
            workers, k, nch));
        TransTools.BusyGate('progress', k / nch);
    end
    TransTools.BusyGate('progress', 1);
    if all(~cellfun(@isempty, betas))
        work.unfold.beta_dc = cat(1, betas{:});
    else
        work.unfold.beta_dc = [];
    end
    if any(~cellfun(@isempty, customs))
        work.unfold.beta_dcCustomrow = cat(1, customs{:});
    end
    work.unfold.channel = 1:nch;
end

function [beta, custom, warned] = fitChannel(shared, row, method, iterations, alpha) %#ok<INUSD> read in evalc
%FITCHANNEL  On a worker: uf_glmfit of one channel, ROW, on the fitted
%   dataset SHARED holds without its data; its betas, any custom-row betas,
%   and the last warning it gave.
    one = shared.Value;
    one.data = row;
    one.nbchan = 1;
    lastwarn('');
    evalc(['one = uf_glmfit(one, ''method'', method, ''lsmriterations'', iterations, ' ...
        '''glmnetalpha'', alpha);']);
    beta = one.unfold.beta_dc;
    custom = [];
    if isfield(one.unfold, 'beta_dcCustomrow')
        custom = one.unfold.beta_dcCustomrow;
    end
    [warned, ~] = lastwarn();
end

function varargout = fetchNextOrCause(futures)
%FETCHNEXTORCAUSE  fetchNext, raising a worker's own error rather than the
%   pool's wrapper around it, so the message is the one the toolbox gave.
    try
        [varargout{1:nargout}] = fetchNext(futures);
    catch err
        if ~isempty(err.cause)
            rethrow(err.cause{1});
        end
        rethrow(err);
    end
end

function bytes = bytesOf(value) %#ok<INUSD> read by whos
%BYTESOF  The memory VALUE takes, sparse matrices counted as stored.
    info = whos('value');
    bytes = info.bytes;
end

function [fitted, plan, excluded, epoched] = regressOnEpochs(work, plan, input, opts, channels)
%REGRESSONEPOCHS  Steps 2 to 4 without overlap correction: the toolbox's
%   own mass-univariate regression on epochs, uf_epoch then
%   uf_glmfit_nodc, as toolboxWorkflow.rst sets it out.
%
%   uf_epoch leaves out every event whose window touches a stretch handed to
%   it (WINREJ) and, through EEGLAB's pop_epoch, every one whose window runs
%   off the recording. It is handed the scan's marks and the recording's
%   cuts, each a single sample, so an epoch is left out when it spans a
%   cut, as DefineBins leaves one out. EPOCHED lists the rows of plan.events
%   that became epochs, in the order of the epochs (read back from a field
%   the copy's events carry, alzRow), and the plan keeps only those
%   (Unfold.keepEvents). A regression on epochs has no neighbours to tell
%   apart, so the epochs left out are simply not in the fit.
%
%   The method is uf_glmfit_nodc's own default, pinv, unless MATLAB's exact
%   solver or glmnet is chosen.
    excluded = artefactScan(work, input, opts, channels);
    if ~isempty(excluded)
        requireEnoughDataLeft(excluded, size(input.data, 2), opts, numel(channels));
    end
    winrej = [excluded; cutSamples(input)];

    candidates = nnz(plan.kept);
    binned = numel(unique([plan.membership{:}]));
    copy = forEpoching(work);
    epochs = [];
    evalc('epochs = uf_epoch(copy, ''winrej'', winrej, ''timelimits'', opts.WindowMs / 1000);');
    epoched = [epochs.urevent.alzRow];
    keep = false(1, numel(plan.events));
    keep(epoched) = true;
    plan = Unfold.keepEvents(plan, keep);
    plan.binnedLeftOut = binned - numel(unique([plan.membership{:}]));
    left = candidates - nnz(plan.kept);
    if left > 0
        plan.notes{end + 1} = sprintf(['%d event(s) were left out (uf_epoch): their %g to %g ms ' ...
            'window touches a stretch left out as artefact or a cut, or runs off the recording.'], ...
            left, opts.WindowMs(1), opts.WindowMs(2));
    end
    if isempty(epoched)
        throw(MException('Alakazam:Unfold:NoEpochs', '%s', sprintf([ ...
            'None of the events has a whole %g to %g ms window of data clear of artefacts and ' ...
            'cuts, so there is nothing to regress on. Would you lower or switch off the artefact ' ...
            'threshold, or shorten the window?'], opts.WindowMs(1), opts.WindowMs(2))));
    end

    method = opts.Solver;
    if strcmp(method, 'default')
        method = 'pinv';
    end
    fitted = uf_glmfit_nodc(epochs, 'method', method, 'glmnetalpha', opts.GlmnetAlpha);
end

function detected = artefactScan(work, input, opts, channels)
%ARTEFACTSCAN  The stretches the toolbox's scan marks, by its own
%   uf_continuousArtifactDetect, as a winrej array; none with the threshold
%   at 0.
    detected = zeros(0, 2);
    if opts.ArtifactThresholdUv > 0
        detected = uf_continuousArtifactDetect(forArtefactScan(work, input.event), ...
            'amplitudeThreshold', opts.ArtifactThresholdUv, ...
            'windowsize', opts.ArtifactWindowMs, ...
            'stepsize', opts.ArtifactStepMs, ...
            'channels', channels);
    end
end

function EEG = forEpoching(EEG)
%FOREPOCHING  A copy shaped the way uf_epoch, through EEGLAB's pop_epoch,
%   expects a dataset: every field of EEGLAB's own empty dataset that the
%   Alakazam dataset lacks (pop_epoch reads .setname, among others), the
%   time axis left for EEGLAB to rebuild in its milliseconds (Alakazam keeps
%   a continuous recording's in seconds), and each event's row of the plan
%   in .alzRow, which the toolbox carries into the epochs' .urevent.
    empty = eeg_emptyset();
    for f = reshape(fieldnames(empty), 1, [])
        if ~isfield(EEG, f{1})
            EEG.(f{1}) = empty.(f{1});
        end
    end
    EEG.times = [];
    for k = 1:numel(EEG.event)
        EEG.event(k).alzRow = k;
    end
end

% ======================================================================= %
function [work, plan] = imputeMissing(work, plan)
%IMPUTEMISSING  The events without a number, as the toolbox's own
%   uf_imputeMissing handles them (plan.missingValues: 'median', 'mean',
%   'marginal' or 'drop'; under 'refuse' Unfold.binModel has refused them
%   already). Its printed lines are kept out of the command window, but its
%   warning that a predictor misses more than 5% of its values goes into the
%   notes. With 'drop' it zeroes each such event's row of the design, and
%   Unfold.binModel has left the same events out of the bins (plan.kept),
%   which is checked here: a difference would mean the two disagree about
%   which events the fit holds.
    missingRows = reshape(any(isnan(work.unfold.X), 2), 1, []);
    if ~any(missingRows)
        return;
    end
    method = plan.missingValues;
    try
        said = evalc('work = uf_imputeMissing(work, ''method'', method);');
    catch err
        if strcmp(method, 'marginal') && contains(err.message, 'Unable to perform assignment')
            throw(MException('Alakazam:Unfold:Marginal', '%s', sprintf([ ...
                'The Unfold toolbox could not draw the missing values ("%s"). Its ' ...
                'uf_imputeMissing (1.3.1) keeps the values it drew for one predictor and draws ' ...
                'into the same list for the next, so when two predictors miss different numbers ' ...
                'of values the second draw does not fit. Here: %s. Would you choose the median, ' ...
                'the mean or drop under Missing values instead?'], err.message, ...
                missingList(plan.missing))));
        end
        rethrow(err);
    end
    warned = regexp(said, '[^\n]*are missing! This could bias your analysis', 'match');
    for k = 1:numel(warned)
        plan.notes{end + 1} = ['The Unfold toolbox warns: ' ...
            regexprep(strtrim(warned{k}), '^\[?Warning:\s*', '')]; %#ok<AGROW>
    end
    if strcmp(method, 'drop') && ~isequal(missingRows, ~plan.kept)
        throw(MException('Alakazam:Unfold:Dropped', '%s', sprintf([ ...
            'uf_imputeMissing left out %d event(s) where %d were expected, so the bins would ' ...
            'count events the fit does not hold. This is a fault in Alakazam, not in the model.'], ...
            nnz(missingRows), nnz(~plan.kept))));
    end
end

function text = missingList(missing)
%MISSINGLIST  "rt, 2 of the 60 events of "A"; size, 1 of ..." for a message.
    text = strjoin(arrayfun(@(m) sprintf('%s, %d of the %d events of "%s"', m.field, m.n, ...
        m.of, m.label), missing, 'UniformOutput', false), '; ');
end

function [data, times, notes] = fittedWaveforms(input, plan, unfold, times, opts)
%FITTEDWAVEFORMS  One fitted waveform per bin, baseline-corrected, with the
%   combination bins computed from them. Each is the model's prediction at
%   the pooled values of its continuous and spline terms (see this file's
%   header), read from the documented EEG.unfold fields: X, colnames,
%   cols2eventtypes, cols2variablenames, variablenames, variabletypes,
%   splines, eventtypes, and the betas uf_condense returns (byColumn),
%   passed as UNFOLD.beta, whose third dimension runs over X's columns, with
%   X's rows on the plan's events; TIMES is uf_condense's, in seconds. NOTES
%   say where a spline could not be averaged over every pooled value (see
%   referencePrediction).
%
%   A BIN IS THE EVENT-WEIGHTED MEAN OF ITS SETS: the waveform of each event
%   type (a set of bins its events share, Unfold.binModel's plan.cellTypes)
%   weighted by how many of the bin's events it holds. For bins that share
%   no events each bin is one set and this is its own waveform; for nested
%   bins it is what Average gives, each event counted in every bin it is in.
    times = reshape(double(times) * 1000, 1, []);   % Unfold keeps seconds
    nchan = size(unfold.beta, 1);
    nbin = numel(input.bindesc);
    data = nan(nchan, numel(times), nbin);
    notes = {};

    perCell = cell(1, numel(plan.cellTypes));
    for c = 1:numel(plan.cellTypes)
        [perCell{c}, cellNotes] = referencePrediction(unfold, plan, c);
        if ~isempty(cellNotes)   % one note per spline, not one per type
            notes = reshape(unique([notes, cellNotes], 'stable'), 1, []);
        end
    end

    ordinary = find(~comboMask(input.bindesc));
    fitted = false(1, nbin);
    for k = 1:numel(plan.binLabels)
        % A set whose events were all left out (uf_imputeMissing's 'drop')
        % weighs nothing, and has no waveform to weigh.
        cells = plan.binCells{k};
        cells = cells(plan.cellCounts(cells) > 0);
        if isempty(cells) || any(cellfun(@isempty, perCell(cells)))
            continue;   % a bin with no events: left as NaN, and noted by binModel
        end
        waveform = zeros(nchan, numel(times));
        for c = reshape(cells, 1, [])
            waveform = waveform + plan.cellCounts(c) / plan.binCounts(k) * perCell{c};
        end
        data(:, :, ordinary(k)) = waveform;
        fitted(ordinary(k)) = true;
    end
    requireFiniteBetas(data, fitted, input.bindesc, opts);
    data = applyBaseline(data, fitted, times, opts.BaselineMs);
    data = resolveComboBins(data, input.bindesc);
end

function EEG = package(input, plan, data, times)
%PACKAGE  The fitted waveforms in Average's own output shape.
    nchan = size(data, 1);
    nbin = size(data, 3);
    EEG = input;
    EEG.data = data;
    EEG.times = times;
    EEG.pnts = numel(times);
    EEG.trials = 1;
    EEG.xmin = times(1) / 1000;
    EEG.xmax = times(end) / 1000;
    EEG.DataFormat = "Averaged";
    EEG.ntrials = sum(plan.cellCounts);  % each event once, as Average counts epochs
    EEG.stErr = zeros(size(data));      % see this file's header: there is none
    EEG.aSME = nan(nchan, nbin);
    ordinary = find(~comboMask(input.bindesc));
    for k = 1:numel(ordinary)
        EEG.bindesc(ordinary(k)).n = plan.binCounts(k);
    end
end

function [EEG, info] = packageTrials(input, plan, work, info, opts, excluded, srate, times)
%PACKAGETRIALS  Overlap-corrected trials, as the toolbox's own uf_erpimage
%   computes them (erpImageTrials), laid out by DefineBins' own cutEpochs,
%   so the result is a DefineBins epoch node in every field but its data:
%   one trial per binned event (an event in two bins is one trial carrying
%   both tags), EEG.epoch, the anchor events' .epoch and latencies,
%   bindesc.trials, and the neighbour table EpochView sorts by. Average,
%   EpochView and the data-quality report then read it without knowing it
%   came from a model, and there is one definition of an epoched dataset
%   rather than two to keep in step.
    events = unique([plan.membership{:}]);           % one trial per binned event
    % The stretches the model left out, and each cut as a stretch of one
    % sample, so a trial spanning a cut is left out, as pop_epoch leaves out
    % an epoch holding a boundary (the model's event list has none).
    if strcmp(opts.TrialsMethod, 'fast')
        [trials, rows] = fastTrials(work, plan, [excluded; cutSamples(input)]);
        workers = 1;
        madeBy = 'the fast route';
    else
        [trials, rows, workers] = erpImageTrials(work, plan, [excluded; cutSamples(input)], opts);
        madeBy = 'uf_erpimage';
    end
    % An event Missing values 'drop' left out keeps its place in the event
    % list, its rows of the design zeroed, so uf_erpimage cuts it a window
    % too; it is not in the fit, and not one of the bins' events.
    inFit = plan.kept(rows);
    % EVENTS ON ONE SAMPLE. uf_erpimage keeps one epoch for events that fall
    % on the same sample (uf_epoch rounds each latency to one): the epoch of
    % the one earlier in the event list, with that event's own fitted
    % response, under the label of the aligned event it finds first by type.
    % Those two need not be the same event, and what it returns does not say
    % which it is, so no event sharing its sample with another modelled
    % event has a trial here. They are counted and noted.
    lat = round([plan.events.latency]);
    [~, ~, sample] = unique(lat);
    perSample = accumarray(sample(:), 1);
    shared = reshape(perSample(sample) > 1, 1, []);
    sharedBinned = nnz(shared & ismember(plan.eventSource, events) & plan.kept);
    trials = trials(:, :, inFit & ~shared(rows));
    [kept, order] = sort(plan.eventSource(rows(inFit & ~shared(rows))));   % the input's events, in its order
    trials = trials(:, :, order);
    if size(trials, 2) ~= numel(times) || ~all(ismember(kept, events)) ...
            || numel(unique(kept)) ~= numel(kept)
        throw(MException('Alakazam:Unfold:Trials', '%s', sprintf([ ...
            '%s returned %d trials of %d samples that do not line up with the %d ' ...
            'binned events and the %d samples of the fit. This is a fault in Alakazam, not in ' ...
            'the model.'], madeBy, size(trials, 3), size(trials, 2), numel(events), numel(times))));
    end

    if isempty(kept)
        throw(MException('Alakazam:Unfold:NoTrials', '%s', sprintf([ ...
            'None of the %d binned events has a whole window of data the model covers, so ' ...
            'there are no overlap-corrected trials to return. A trial is dropped when its ' ...
            '%g to %g ms window runs off the recording or touches a stretch left out of the ' ...
            'model (artefact, or around a cut). Would you lower or switch off the artefact ' ...
            'threshold, or ask for one waveform per bin instead?'], ...
            numel(events), opts.WindowMs(1), opts.WindowMs(2))));
    end
    EEG = asEpochNode(input, plan, kept, trials, opts, srate, times);

    info.output = 'trials';
    info.trialsMethod = opts.TrialsMethod;
    info.trialWorkers = workers;
    info.trialCandidates = numel(events);
    info.trials = numel(kept);
    info.trialsDropped = numel(events) - numel(kept);
    info.trialsSharingASample = sharedBinned;
    if sharedBinned > 0
        note = sprintf(['%d binned event(s) fall on the same sample as another event of the ' ...
            'model, so they have no overlap-corrected trial: uf_erpimage keeps one epoch for ' ...
            'events at the same time, and which event''s fitted response it holds cannot be ' ...
            'read from what it returns. The fast route leaves them out too, so both give the ' ...
            'same trials.'], sharedBinned);
        info.notes{end + 1} = note;
        warning('Alakazam:Unfold:design', '%s', note);
    end
end

function [trials, rows, workers] = erpImageTrials(work, plan, excluded, opts)
%ERPIMAGETRIALS  One overlap-corrected trial per binned event, by the
%   toolbox's own uf_erpimage, called as its tutorial on ERP images calls it
%   for "modelled plus residuals": 'type' 'modelled' with 'addResiduals' 1
%   is each event's own fitted response (its row of X times the betas) plus
%   what the model leaves unexplained around it (the recording less every
%   event's fitted response, Xdc times the betas), which is the recording
%   with every OTHER event's response taken out. Its 'winrej' is the
%   stretches the model left out (EXCLUDED), as its help asks, so uf_epoch
%   leaves out a trial whose window touches one, and EEGLAB's pop_epoch one
%   whose window runs off the recording. It is aligned to the binned event
%   types only, so the nuisance events are subtracted but are not trials.
%
%   ONE CHANNEL AT A TIME, since uf_erpimage returns the mean of the
%   channels it is given, each from a copy of the fitted dataset holding that
%   channel alone, its data and its betas (erpImageOfChannel). Every step
%   uf_erpimage takes is per channel, so this is what it returns for that
%   channel of the whole dataset, without epoching every channel once for
%   each channel.
%
%   PROGRESS AND WORKERS. Each call epochs the whole recording around every
%   modelled event three times, which on a long recording with many events
%   takes a while, so the channels done are shown in the app's busy
%   indicator, and on a large job (Unfold.workerPool, opts.Parallel) the
%   calls go to parallel workers, the same call on the same copy. WORKERS
%   is how many there were (1: none).
%
%   WHICH EVENT EACH TRIAL IS comes from uf_erpimage's own sort vector: each
%   event carries its row of the model as .alzRow (forEpoching, which also
%   gives the copy the fields pop_epoch reads), and sorting by that field
%   at time zero ('sort_by' with 'sort_time' [0 0], the third use its help
%   gives) returns it per trial. ROWS are rows of plan.events.
%
%   TOOLBOX BEHAVIOUR KEPT AS IT IS (1.3.1): uf_epoch rounds each latency
%   to the nearest sample; and where events fall on the same sample,
%   uf_erpimage keeps one epoch for them, with a warning (packageTrials
%   says why such events have no trial).
    copy = forEpoching(work);   % shaped for EEGLAB's pop_epoch, rows in .alzRow
    copy.urevent = [];          % uf_erpimage numbers the events itself (.urevent)
    alignto = reshape(plan.cellTypes(plan.cellCounts > 0), 1, []);
    nchan = size(work.data, 1);
    base = rmfield(copy, 'data');
    base.unfold.beta_dc = [];
    custom = [];
    if isfield(work.unfold, 'beta_dcCustomrow') && ~isempty(work.unfold.beta_dcCustomrow)
        custom = work.unfold.beta_dcCustomrow;
        base.unfold.beta_dcCustomrow = [];
    end
    chanlocs = [];
    if isfield(work, 'chanlocs') && numel(work.chanlocs) >= nchan
        chanlocs = work.chanlocs;
    end
    part = @(ch) {work.data(ch, :), work.unfold.beta_dc(ch, :, :), rowOrEmpty(custom, ch), ...
        rowOrEmpty(chanlocs, ch), alignto, excluded};

    nEvents = numel(work.event);
    bytes = bytesOf(base) + 8 * size(work.data, 2) * 4 + ...
        3 * 8 * nEvents * size(work.unfold.beta_dc, 2) + 1.5e9;
    [pool, workers] = Unfold.workerPool(nchan, size(work.data, 2) * nEvents, 1e9, opts.Parallel, bytes);

    images = cell(1, nchan);
    orders = cell(1, nchan);
    TransTools.BusyGate('progress', 0);
    if isempty(pool)
        for ch = 1:nchan
            TransTools.BusyGate('message', sprintf('Overlap-corrected trials: channel %d of %d...', ...
                ch, nchan));
            args = part(ch);
            [images{ch}, orders{ch}] = erpImageOfChannel(base, args{:});
            TransTools.BusyGate('progress', ch / nchan);
        end
    else
        shared = parallel.pool.Constant(base);
        futures = parallel.FevalFuture.empty;
        for ch = 1:nchan
            args = part(ch);
            futures(ch) = parfeval(pool, @erpImageOfChannel, 2, shared, args{:});
        end
        stopTheRest = onCleanup(@() cancel(futures));
        for k = 1:nchan
            [ch, image, order] = fetchNextOrCause(futures);
            images{ch} = image;
            orders{ch} = order;
            TransTools.BusyGate('message', sprintf(['Overlap-corrected trials on %d workers: ' ...
                '%d of %d channels...'], workers, k, nchan));
            TransTools.BusyGate('progress', k / nchan);
        end
    end
    TransTools.BusyGate('progress', 1);

    trials = [];
    rows = [];
    for ch = 1:nchan
        if isempty(images{ch})
            % Every window was left out: packageTrials says so, and why.
            trials = zeros(nchan, 0, 0);
            rows = zeros(1, 0);
            return;
        end
        if ch == 1
            rows = reshape(orders{1}, 1, []);
            trials = nan(nchan, size(images{1}, 2), numel(rows));
        elseif ~isequal(reshape(orders{ch}, 1, []), rows)
            throw(MException('Alakazam:Unfold:Trials', '%s', sprintf([ ...
                'uf_erpimage returned other trials for channel %d than for channel 1. This ' ...
                'is a fault in Alakazam, not in the model.'], ch)));
        end
        trials(ch, :, :) = reshape(images{ch}, 1, size(images{ch}, 2), []);
    end
end

function [image, order] = erpImageOfChannel(base, row, beta, custom, chanloc, alignto, excluded) %#ok<INUSD> read in evalc
%ERPIMAGEOFCHANNEL  uf_erpimage's modelled responses plus residuals for
%   one channel: BASE, the fitted dataset without its data and betas (or a
%   parallel.pool.Constant holding it), given the channel's data ROW, its
%   betas, the betas of any column the solver added (glmnet's
%   'glmnet-DC-Correction') and its location. IMAGE and ORDER are empty when
%   every window was left out. Runs here or on a worker alike.
    if isa(base, 'parallel.pool.Constant')
        base = base.Value;
    end
    one = base;
    one.data = row;
    one.nbchan = 1;
    if ~isempty(chanloc)
        one.chanlocs = chanloc;
    end
    one.unfold.beta_dc = beta;
    if ~isempty(custom)
        one.unfold.beta_dcCustomrow = custom;
    end
    one.unfold.channel = 1;
    image = [];
    order = [];
    try
        % Quiet: it reports every step it takes.
        evalc(['[image, order] = uf_erpimage(one, ''channel'', 1, ''type'', ''modelled'', ' ...
            '''addResiduals'', 1, ''alignto'', alignto, ''winrej'', excluded, ' ...
            '''sort_by'', ''alzRow'', ''sort_time'', [0 0], ''figure'', 0);']);
    catch err
        if ~contains(err.message, 'Did not find any events to align to')
            rethrow(err);
        end
        image = [];
        order = [];
    end
end

function value = rowOrEmpty(values, k)
%ROWOREMPTY  Row K of VALUES, or [] when there are none.
    if isempty(values)
        value = [];
    elseif isstruct(values)
        value = values(k);
    else
        value = values(k, :);
    end
end

function [trials, rows] = fastTrials(work, plan, excluded)
%FASTTRIALS  The overlap-corrected trials uf_erpimage gives (erpImageTrials),
%   computed by Alakazam, for every channel at once and only at the trials'
%   own samples: TrialsMethod 'fast'.
%
%   WHAT A TRIAL IS, as uf_erpimage computes it with 'type' 'modelled' and
%   'addResiduals' 1: the event's own fitted response, its row of X times
%   the betas uf_condense returns, plus the recording less the whole model's
%   prediction (Xdc times the betas) over the event's window. The prediction
%   is needed only in the trials' windows, a few hundred thousand samples
%   rather than the whole recording three times per channel, which is where
%   the time goes.
%
%   WHICH EVENTS, AND WHICH SAMPLES, as the toolbox's own route has them:
%   the events of the binned types; each latency rounded, as uf_epoch rounds
%   it; an event left out when its window touches a stretch the model left
%   out (EXCLUDED), by uf_epoch's own test, over the window it passes to
%   pop_epoch (the fit's lags and one sample more); and the samples EEGLAB's
%   epoch cuts, from round(start x rate) to round(end x rate - 1) around the
%   rounded latency, an event whose window runs off the recording left out.
%   Events that share a sample are kept here and left out by packageTrials,
%   as from uf_erpimage's trials. ROWS are rows of plan.events, ascending,
%   as uf_erpimage's sort by .alzRow gives them.
    nchan = size(work.data, 1);
    srate = work.srate;
    times = work.unfold.times;                      % the fit's lags, in seconds
    limits = [times(1), times(end) + (times(2) - times(1))];   % uf_epoch's window
    first = round(limits(1) * srate);               % EEGLAB's epoch
    last = round(limits(2) * srate - 1);
    lags = (first:last)';
    if numel(lags) ~= numel(times)
        throw(MException('Alakazam:Unfold:Trials', '%s', sprintf([ ...
            'The trials'' windows hold %d samples but the fit has %d lags. This is a fault ' ...
            'in Alakazam, not in the model.'], numel(lags), numel(times))));
    end

    alignto = reshape(plan.cellTypes(plan.cellCounts > 0), 1, []);
    types = arrayfun(@(e) char(string(e.type)), work.event, 'UniformOutput', false);
    aligned = find(ismember(types, alignto));
    latency = round(double([work.event(aligned).latency]));
    % uf_epoch's own test, event by event against every left-out stretch.
    from = latency + limits(1) * srate;
    to = latency + limits(2) * srate;
    touches = false(size(latency));
    for w = 1:size(excluded, 1)
        win = excluded(w, :);
        touches = touches | (win(1) < from & win(2) >= from) | (win(1) >= from & to >= win(1));
    end
    inside = latency + first >= 1 & latency + last <= size(work.data, 2);
    keep = ~touches & inside;
    rows = reshape(aligned(keep), 1, []);
    latency = latency(keep);
    if isempty(rows)
        trials = zeros(nchan, 0, 0);
        rows = zeros(1, 0);
        return;
    end

    % The betas as uf_erpimage reads them: uf_condense's, channels x lags x
    % columns, and Xdc's columns lag within column, one block per column.
    condensed = [];
    evalc('condensed = uf_condense(work);');
    beta = condensed.beta;
    nLags = size(beta, 2);
    nCols = size(beta, 3);
    expanded = reshape(permute(beta, [2 3 1]), nLags * nCols, nchan);      % for Xdc
    own = reshape(permute(beta, [3 2 1]), nCols, nLags * nchan);           % for X
    X = work.unfold.X;
    XdcT = work.unfold.Xdc';    % columns of a sparse matrix are quick to take, rows are not

    n = numel(rows);
    trials = zeros(nchan, nLags, n);
    chunk = 64;
    TransTools.BusyGate('message', sprintf('Overlap-corrected trials (fast): %d trials...', n));
    TransTools.BusyGate('progress', 0);
    for s = 1:chunk:n
        k = s:min(s + chunk - 1, n);
        samples = latency(k) + lags;                                   % lags x trials
        predicted = full(XdcT(:, samples(:))' * expanded);             % (lags x trials) x channels
        recorded = double(work.data(:, samples(:)))';
        residual = permute(reshape(recorded - predicted, nLags, numel(k), nchan), [3 1 2]);
        fitted = permute(reshape(full(X(rows(k), :)) * own, numel(k), nLags, nchan), [3 2 1]);
        trials(:, :, k) = fitted + residual;
        TransTools.BusyGate('progress', k(end) / n);
    end
    TransTools.BusyGate('progress', 1);
end

function [EEG, info] = packageEpochs(input, plan, fitted, epoched, info, opts, srate, times)
%PACKAGEEPOCHS  Without overlap correction, the trials are the epochs as
%   uf_epoch cut them, one per binned event in the fit, laid out as
%   packageTrials lays out the corrected ones. Nothing is subtracted: there
%   is no model of the neighbours to subtract.
    kept = unique([plan.membership{:}]);              % binned events in the fit
    if isempty(kept)
        throw(MException('Alakazam:Unfold:NoTrials', ...
            'None of the binned events became an epoch, so there are no trials to return.'));
    end
    epochOf = zeros(1, numel(plan.events));
    epochOf(epoched) = 1:numel(epoched);
    rowOf = arrayfun(@(e) find(plan.eventSource == e, 1), kept);
    EEG = asEpochNode(input, plan, kept, fitted.data(:, :, epochOf(rowOf)), opts, srate, times);

    info.output = 'trials';
    info.trialCandidates = numel(kept) + plan.binnedLeftOut;
    info.trials = numel(kept);
    info.trialsDropped = plan.binnedLeftOut;
end

function EEG = asEpochNode(input, plan, kept, trials, opts, srate, times)
%ASEPOCHNODE  TRIALS, one per event in KEPT, baseline-corrected as the
%   waveforms are, and laid out by DefineBins' cutEpochs, time zero on the
%   sample Unfold placed each event on. cutEpochs is handed one channel
%   only: it would otherwise cut every channel of the raw recording just for
%   that to be replaced by the trials a line later.
    if ~isempty(opts.BaselineMs)
        inWindow = baselineSamples(times, opts.BaselineMs);
        trials = trials - mean(trials(:, inWindow, :), 2);
    end
    layout = input;
    layout.data = input.data(1, :);
    centre = zeros(1, numel(input.event));
    centre(kept) = round([input.event(kept).latency]);
    lags = round(times / 1000 * srate);
    [EEG, bindesc] = DefineBinsEngine.cutEpochs(layout, keptBins(input.bindesc, plan, kept), ...
        struct('lo', lags(1), 'hi', lags(end) + 1, 'unit', 'samples'), centre);
    EEG.bindesc = bindesc;
    EEG.data = trials;
    EEG.times = times;      % uf_condense's own, which the waveforms carry too
end

function bindesc = keptBins(bindesc, plan, kept)
%KEPTBINS  Each ordinary bin's events cut to those that became trials, which
%   is what cutEpochs makes trials of, with their reaction times kept
%   aligned: cutEpochs derives .trials from .events element by element, and a
%   sort by reaction time pairs .rt with .trials the same way.
    ordinary = find(~comboMask(bindesc));
    for k = 1:numel(ordinary)
        b = ordinary(k);
        members = plan.membership{k};
        members = members(ismember(members, kept));
        rts = nan(1, numel(members));
        if isfield(bindesc, 'events') && isfield(bindesc, 'rt') ...
                && numel(bindesc(b).rt) == numel(bindesc(b).events)
            [found, at] = ismember(members, bindesc(b).events);
            rts(found) = bindesc(b).rt(at(found));
        end
        bindesc(b).events = members;
        bindesc(b).rt = rts;
        bindesc(b).n = numel(members);
    end
end

function info = modelInfo(input, plan, opts, excluded, srate)
%MODELINFO  What was fitted and what was left out, for EEG.etc.alz.unfold.
    info = struct('output', 'average', 'window', opts.WindowMs, 'baseline', opts.BaselineMs, ...
        'binLabels', {plan.binLabels}, ...
        'binTypes', {plan.binTypes}, 'binCounts', plan.binCounts, ...
        'nuisanceTypes', {plan.nuisanceTypes}, ...
        'formulas', struct('type', plan.typeLabels, 'formula', plan.formulas), ...
        'notes', {plan.notes}, ...
        'excludedIntervals', excluded, ...
        'excludedSeconds', sum(diff(excluded, 1, 2)) / srate, ...
        'recordingSeconds', size(input.data, 2) / srate, ...
        'artifact', struct('thresholdUv', opts.ArtifactThresholdUv, ...
                           'windowMs', opts.ArtifactWindowMs, 'stepMs', opts.ArtifactStepMs, ...
                           'channels', {opts.Channels}), ...
        'solverIterations', opts.SolverIterations, ...
        'solver', opts.Solver, 'glmnetAlpha', opts.GlmnetAlpha, ...
        'overlapCorrection', opts.OverlapCorrection, ...
        'excludeAroundCuts', opts.ExcludeAroundCuts, ...
        'codingSchema', opts.CodingSchema, ...
        'modelFitMethod', opts.ModelFit, ...
        'collinearityOf', opts.Collinearity, ...
        'timeBasis', opts.TimeBasis, 'timeBasisCount', opts.TimeBasisCount, ...
        'missingValues', plan.missingValues, 'missing', {plan.missing}, ...
        'dropped', strcmp(plan.missingValues, 'drop') * nnz(plan.missingRows), ...
        'epochsLeftOut', TransTools.FieldOr(plan, 'binnedLeftOut', 0), ...
        'cellLabels', {plan.cellLabels}, 'binCells', {plan.binCells}, ...
        'hasStandardError', false);
end

function EEG = recordInfo(EEG, info, plan)
%RECORDINFO  The provenance on the dataset, and the model's notes in the log.
    if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc)
        EEG.etc = struct();
    end
    if ~isfield(EEG.etc, 'alz') || ~isstruct(EEG.etc.alz)
        EEG.etc.alz = struct();
    end
    EEG.etc.alz.unfold = info;

    for k = 1:numel(plan.notes)
        warning('Alakazam:Unfold:design', '%s', plan.notes{k});
    end
end

function window = resolveBaseline(requested, responseWindow)
%RESOLVEBASELINE  The baseline window: [] for none, the window as given, or
%   for 'pre-event' the part of the response window before the event. A
%   response window that starts at or after the event has no pre-event part,
%   so 'pre-event' then leaves the correction off rather than invent one.
    if isempty(requested)
        window = [];
        return;
    end
    if ischar(requested) || isstring(requested)
        if responseWindow(1) >= 0
            window = [];
        else
            window = [responseWindow(1) 0];
        end
        return;
    end
    window = reshape(double(requested), 1, 2);
    if window(1) < responseWindow(1) || window(2) > responseWindow(2)
        throw(MException('Alakazam:Unfold:Baseline', '%s', sprintf([ ...
            'The baseline window (%g to %g ms) reaches outside the response window (%g to ' ...
            '%g ms), so part of it was never fitted and there is nothing there to average. ' ...
            'Would you either move the baseline inside the response window or widen the ' ...
            'response window to cover it?'], window(1), window(2), ...
            responseWindow(1), responseWindow(2))));
    end
end

function data = applyBaseline(data, fitted, times, window)
%APPLYBASELINE  Each fitted waveform, shifted so its baseline window averages
%   zero. Combination bins are resolved afterwards, and a difference of two
%   corrected waveforms is the correction of their difference, so they inherit
%   it rather than needing their own.
    if isempty(window)
        return;
    end
    inWindow = baselineSamples(times, window);
    for b = find(fitted)
        slice = data(:, :, b);
        data(:, :, b) = slice - mean(slice(:, inWindow), 2);
    end
end

function inWindow = baselineSamples(times, window)
%BASELINESAMPLES  The samples a baseline window covers, as Unfold's
%   uf_plotParam takes its baseline: from the start up to, but not
%   including, the stop, so the default -200 to 0 ms leaves out the sample
%   at 0 ms. A window too short to hold a sample has no mean to subtract.
    inWindow = times >= window(1) & times < window(2);
    if ~any(inWindow)
        throw(MException('Alakazam:Unfold:Baseline', '%s', sprintf([ ...
            'The baseline window (%g to %g ms) holds no sample at this sampling rate, so ' ...
            'there is nothing to average over. Would you widen it?'], window(1), window(2))));
    end
end

function requireFitPossible(opts, input)
%REQUIREFITPOSSIBLE  Refuse a model fit that cannot be computed, before the
%   model itself is: uf_checkmodelfit works on the time-expanded design, and
%   its cross-validation cuts the folds at an event code of the recording.
    if strcmp(opts.ModelFit, 'none')
        return;
    end
    if ~opts.OverlapCorrection
        throw(MException('Alakazam:Unfold:FitNeedsDeconvolution', '%s', sprintf([ ...
            'A model fit (%s) is measured on the deconvolution''s time-expanded design, by the ' ...
            'toolbox''s uf_checkmodelfit, so it needs Overlap correction on. Would you switch ' ...
            'it on, or leave the model fit out?'], opts.ModelFit)));
    end
    if ~startsWith(opts.ModelFit, 'crossVal')
        return;
    end
    types = arrayfun(@(e) char(string(e.type)), input.event, 'UniformOutput', false);
    if isempty(opts.FoldEvent) || ~any(ismember(types, opts.FoldEvent))
        throw(MException('Alakazam:Unfold:NoFoldEvent', '%s', sprintf([ ...
            'A cross-validated fit (%s) cuts the recording into folds at an event code, and ' ...
            '%s. Would you choose a code that marks the breaks between blocks, say, under ' ...
            '"Folds cut at"?'], opts.ModelFit, foldEventText(opts.FoldEvent))));
    end
end

function text = foldEventText(codes)
    if isempty(codes)
        text = 'none was chosen';
    else
        text = sprintf('this recording has no event "%s"', strjoin(codes, '", "'));
    end
end

function [fit, note] = modelFit(work, input, plan, opts)
%MODELFIT  How well the model explains the recording, by the toolbox's own
%   uf_checkmodelfit, called as a script calls it, on the fitted dataset
%   with the recording's own events, which is where uf_cv_getFolds looks for
%   the fold events (the refits use the design alone). FIT holds the method,
%   the fold event, the channels' labels, .r2 per channel (the mean over the
%   folds when cross-validated, .r2Folds then channels x folds) and, for a
%   partial fit, .partial: one row per channel, term (and fold) with the
%   term's partial R2 and the full model's R2. NOTE says it in a line.
    copy = work;
    copy.event = input.event;
    args = {'method', opts.ModelFit};
    if startsWith(opts.ModelFit, 'crossVal')
        args = [args, {'fold_event', opts.FoldEvent}];
    end
    result = [];
    try
        evalc('result = uf_checkmodelfit(copy, args{:});');   % quiet: it reports every fold
    catch err
        if contains(err.message, 'folding events have modelled overlap')
            throw(MException('Alakazam:Unfold:FoldOverlap', '%s', sprintf([ ...
                'The toolbox''s cross-validation cuts the folds at "%s", and at least one of those ' ...
                'events falls where the model has an event''s response; a fold cut there would ' ...
                'split a response in two, so it refuses. Would you choose a code that marks a ' ...
                'pause with nothing modelled in it, such as the breaks between blocks (or the ' ...
                'cuts, with "Around cuts" leaving a window out)?'], strjoin(opts.FoldEvent, '", "'))));
        end
        rethrow(err);
    end
    labels = arrayfun(@(k) sprintf('%d', k), 1:size(work.data, 1), 'UniformOutput', false);
    if isfield(input, 'chanlocs') && numel(input.chanlocs) == size(work.data, 1)
        labels = cellfun(@(l) char(string(l)), {input.chanlocs.labels}, 'UniformOutput', false);
    end
    fit = struct('method', opts.ModelFit, 'foldEvent', {opts.FoldEvent}, 'channels', {labels}, ...
        'r2', [], 'r2Folds', [], 'partial', []);
    if istable(result)
        termOf = termLabels(work.unfold, plan);
        rows = table2struct(result);
        partial = struct('channel', {}, 'term', {}, 'partialR2', {}, 'r2', {}, 'fold', {});
        for k = 1:numel(rows)
            fold = NaN;
            if isfield(rows(k), 'fold')
                fold = double(rows(k).fold);
            end
            name = char(string(rows(k).variablenames));
            term = name;
            if isKey(termOf, name)
                term = termOf(name);
            end
            partial(end + 1) = struct('channel', labels{rows(k).channel}, 'term', term, ...
                'partialR2', double(rows(k).r2_ca), 'r2', double(rows(k).r2_total), 'fold', fold); %#ok<AGROW>
        end
        fit.partial = partial;
        terms = unique({partial.term}, 'stable');
        perTerm = cellfun(@(t) median([partial(strcmp({partial.term}, t)).partialR2]), terms);
        total = [partial.r2];
        fit.r2 = arrayfun(@(c) mean(total(strcmp({partial.channel}, labels{c}))), 1:numel(labels))';
        parts = arrayfun(@(k) sprintf('%s %.3g', terms{k}, perTerm(k)), 1:numel(terms), ...
            'UniformOutput', false);
        note = sprintf('%s; partial R2 of each term, the median over channels%s: %s.', ...
            fitLine(fit), foldsText(partial), strjoin(parts, ', '));
    else
        if size(result, 2) > 1
            fit.r2Folds = result;
        end
        fit.r2 = mean(result, 2, 'omitnan');
        note = [fitLine(fit) '.'];
    end
end

function [vif, note] = collinearity(work, plan)
%COLLINEARITY  Every column's variance inflation factor, by the toolbox's
%   own uf_vif, of the time-expanded design as fitted (handed Xdc). VIF
%   holds the matrix, each column's label and value, and each term's
%   largest; NOTE gives the terms.
    unfold = work.unfold;
    values = [];
    evalc('values = uf_vif(unfold.Xdc);');   %#ok<NASGU>  quiet: it says it is inverting
    termOfColumn = unfold.cols2variablenames(unfold.Xdc_terms2cols);
    matrix = 'Xdc';
    values = reshape(values, 1, []);
    termOf = termLabels(unfold, plan);
    names = cell(1, numel(termOfColumn));
    for c = 1:numel(termOfColumn)
        v = termOfColumn(c);
        names{c} = sprintf('variable %d', v);
        if v <= numel(unfold.variablenames) && isKey(termOf, unfold.variablenames{v})
            names{c} = termOf(unfold.variablenames{v});
        end
    end
    terms = unique(names, 'stable');
    largest = cellfun(@(t) max(values(strcmp(names, t))), terms);
    vif = struct('matrix', matrix, 'columns', {names}, 'perColumn', values, ...
        'terms', {terms}, 'largest', largest);
    which = 'the time-expanded design as fitted, Xdc';
    parts = arrayfun(@(k) sprintf('%s %.3g', terms{k}, largest(k)), 1:numel(terms), ...
        'UniformOutput', false);
    note = sprintf(['Collinearity (uf_vif, %s): the largest variance inflation factor of each ' ...
        'term: %s.'], which, strjoin(parts, ', '));
end

function text = fitLine(fit)
%FITLINE  "Model fit (uf_checkmodelfit, R2): R2 0.31, the median over 64
%   channels (0.12 to 0.45)".
    how = fit.method;
    if startsWith(fit.method, 'crossVal')
        how = sprintf('%s, folds cut at "%s"', fit.method, strjoin(fit.foldEvent, '", "'));
    end
    text = sprintf('Model fit (uf_checkmodelfit, %s): R2 %.3g, the median over %d channel(s) (%.3g to %.3g)', ...
        how, median(fit.r2, 'omitnan'), numel(fit.r2), min(fit.r2), max(fit.r2));
end

function text = foldsText(partial)
    text = '';
    folds = unique([partial.fold]);
    folds = folds(isfinite(folds));
    if ~isempty(folds)
        text = sprintf(' and %d folds', numel(folds));
    end
end

function termOf = termLabels(unfold, plan)
%TERMLABELS  Each of the toolbox's variable names, as uf_checkmodelfit
%   reports them, as "<bin> <term>": the bin or code its event type stands
%   for, and the name without the "2_" Unfold puts before a later type's.
    termOf = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for v = 1:numel(unfold.variablenames)
        t = unfold.cols2eventtypes(find(unfold.cols2variablenames == v, 1));
        if isempty(t) || t > numel(plan.typeLabels)
            continue;
        end
        termOf(unfold.variablenames{v}) = sprintf('%s %s', plan.typeLabels{t}, ...
            regexprep(unfold.variablenames{v}, '^\d+_', ''));
    end
end

function winrej = cutSamples(EEG)
%CUTSAMPLES  Each cut in the recording (a boundary event) as a winrej row of
%   one sample, for uf_epoch: an epoch whose window touches it is left out,
%   as EEGLAB's pop_epoch leaves out an epoch that holds a boundary.
    winrej = zeros(0, 2);
    if ~isfield(EEG, 'event') || isempty(EEG.event) || ~isfield(EEG.event, 'type')
        return;
    end
    types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
    cuts = round(double([EEG.event(strcmpi(types, 'boundary')).latency]));
    winrej = [cuts(:), cuts(:)];
end

function intervals = boundaryIntervals(EEG, windowMs, srate)
%BOUNDARYINTERVALS  The samples around each cut in the recording, as a winrej
%   array. A response window reaching across a boundary is fitted across the
%   join between two moments that were never adjacent, so the window's worth
%   of samples on each side of the cut is left out of the model.
    intervals = zeros(0, 2);
    if ~isfield(EEG, 'event') || isempty(EEG.event) || ~isfield(EEG.event, 'type')
        return;
    end
    types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
    cuts = double([EEG.event(strcmpi(types, 'boundary')).latency]);
    if isempty(cuts)
        return;
    end
    reach = ceil(max(abs(windowMs)) * srate / 1000);
    npnts = size(EEG.data, 2);
    intervals = [max(1, round(cuts(:)) - reach), min(npnts, round(cuts(:)) + reach)];
    intervals = intervals(intervals(:, 1) < intervals(:, 2), :);
end

function intervals = combineIntervals(a, b)
%COMBINEINTERVALS  Two winrej arrays as one, with overlaps merged, so the
%   share of the recording left out is counted once rather than twice. By
%   the toolbox's own uf_combineWinrej, as a script joins them, so the
%   exclusion is the one Unfold makes (this file's header says how it joins
%   a stretch lying inside another).
    if isempty(a); intervals = b; return; end
    if isempty(b); intervals = a; return; end
    intervals = uf_combineWinrej(a, b);
end

function EEG = forArtefactScan(EEG, events)
%FORARTEFACTSCAN  A copy shaped the way the ERPLAB scan expects a dataset.
%   The scan is ERPLAB's crap.m, and it reads two fields that EEGLAB always
%   writes but an Alakazam dataset need not carry: {chanlocs.type}, to warn
%   about EYE-EEG channels, and .epoch, to check the data is continuous. Both
%   are absent rather than empty on a dataset whose importer never wrote them,
%   and the scan then fails with a bare "Unrecognized field name" that says
%   nothing about what went wrong. Filled in on the copy handed to the scan, so
%   nothing is written into the dataset the user keeps.
%
%   THE SCAN GETS THE DATASET'S OWN EVENTS (EVENTS), boundaries included, as
%   it does in an Unfold script, where uf_continuousArtifactDetect is handed
%   the dataset itself: ERPLAB's scan then looks at each stretch between cuts
%   on its own. The model's event list, which the design was built from,
%   carries no boundaries, and with it the scan ran across every cut as if
%   the recording were one stretch.
    EEG.event = events;
    if ~isempty(EEG.chanlocs) && ~isfield(EEG.chanlocs, 'type')
        [EEG.chanlocs.type] = deal('');
    end
    if ~isfield(EEG, 'epoch')
        EEG.epoch = [];
    end
end

function channels = scalpChannels(EEG)
%SCALPCHANNELS  Which channels the artefact scan looks at with Channels
%   'scalp': the scalp EEG. An EOG channel swings several times as far as
%   the EEG, so including it means the scan reports the eyes rather than
%   the data.
    channels = 1:size(EEG.data, 1);
    if ~isfield(EEG, 'chanlocs') || numel(EEG.chanlocs) ~= numel(channels)
        return;
    end
    mask = eegChannelMask(EEG.chanlocs);
    if any(mask)
        channels = find(mask);
    end
end

function requireCentredData(EEG)
%REQUIRECENTREDDATA  Refuse a recording whose channels carry a DC offset.
%   The test is scale-free on purpose: a channel whose mean is further from
%   zero than its own variation is dominated by its offset, whatever the
%   units or the gain. Judged on the median channel, so one dead or saturated
%   electrode does not decide it.
    data = double(EEG.data);
    if isempty(data)
        return;
    end
    mu = mean(data, 2);
    sd = std(data, 0, 2);
    ratio = abs(mu) ./ max(sd, eps);
    if median(ratio) <= 1
        return;
    end
    [~, worst] = max(abs(mu));
    label = sprintf('channel %d', worst);
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= worst
        label = char(string(EEG.chanlocs(worst).labels));
    end
    throw(MException('Alakazam:Unfold:DcOffset', '%s', sprintf([ ...
        'This recording carries a DC offset that deconvolution cannot work with: the median ' ...
        'channel sits %.1f times its own standard deviation away from zero, and %s averages ' ...
        '%.0f uV. Would you run DCDetrend (order 0 removes the mean, order 1 a linear drift) ' ...
        'or a high-pass Filter on the continuous data first, and deconvolve that?\n\nThe ' ...
        'reason it matters here and not in an average: this fits one waveform per bin against ' ...
        'the whole recording, with no constant term to absorb a standing voltage, so that ' ...
        'offset has to come out of the event responses and the result is thousands of ' ...
        'microvolts of nothing. Epoch-and-average hides it because Baseline subtracts the ' ...
        'offset epoch by epoch, which is why the same file averages perfectly well.'], ...
        median(ratio), label, mu(worst))));
end

function requireEnoughDataLeft(excluded, npnts, opts, nchannels)
%REQUIREENOUGHDATALEFT  Refuse when the artefact scan has taken most of the
%   recording. The threshold is peak-to-peak within a window that is marked
%   whole, so frequent large swings (a drift, blinks, or the eye movements a
%   reading or free-viewing task is made of) can mark nearly everything, and
%   the fit that follows is a slow way of producing NaN.
    marked = sum(diff(excluded, 1, 2));
    fraction = marked / max(npnts, 1);
    if fraction <= 0.5
        return;
    end
    throw(MException('Alakazam:Unfold:TooMuchExcluded', '%s', sprintf([ ...
        'The artefact scan marked %.0f%% of this recording as bad (%d segments, at %g uV ' ...
        'peak-to-peak over %d channel(s)), which leaves too little for the model to be ' ...
        'fitted against. The limit is peak-to-peak within each moving window, and the ' ...
        'whole window is marked, so frequent large swings (a drift, blinks, or the eye ' ...
        'movements a reading or free-viewing task consists of) trip it almost everywhere. ' ...
        'Would you either high-pass or detrend the data first, raise the threshold, or set ' ...
        'it to 0 to leave artefact rejection out of this step?'], ...
        100 * fraction, size(excluded, 1), opts.ArtifactThresholdUv, nchannels)));
end

function requireFiniteBetas(data, fitted, bindesc, opts)
%REQUIREFINITEBETAS  Never hand back a waveform of NaNs. The solver returns
%   them without complaint when the system it was given has no answer, and a
%   plot of NaN looks like a bug in the app rather than a fact about the fit,
%   so the fit says so itself.
    bad = false(1, numel(fitted));
    for b = find(fitted)
        bad(b) = ~all(isfinite(data(:, :, b)), 'all');
    end
    if ~any(bad)
        return;
    end
    labels = arrayfun(@(d) char(string(d.label)), bindesc(bad), 'UniformOutput', false);
    throw(MException('Alakazam:Unfold:NotFinite', '%s', sprintf([ ...
        'The fit returned no usable numbers for %s, I''m afraid: the solver produced NaN ' ...
        'rather than a waveform. That happens when the system it was given has no answer ' ...
        'to give, which in practice means one of three things: too much of the data was ' ...
        'excluded as artefact (the threshold here was %g uV), the recording carries an ' ...
        'offset or drift the model cannot absorb (try DCDetrend or a high-pass Filter), or ' ...
        'two bins hold events whose timing never varies relative to each other, which no ' ...
        'amount of deconvolution can separate.'], strjoin(labels, ', '), ...
        opts.ArtifactThresholdUv)));
end

function mask = comboMask(bindesc)
    mask = false(1, numel(bindesc));
    if isfield(bindesc, 'combo')
        mask = ~cellfun(@isempty, {bindesc.combo});
    end
end

function [waveform, notes] = referencePrediction(unfold, plan, c)
%REFERENCEPREDICTION  The waveform of event type C of the plan (a set of
%   bins, plan.cellTypes): the betas of the type's own columns
%   (cols2eventtypes), weighted by its events' own average design row for
%   the intercept and factors, and with every continuous or spline term at
%   the pooled values (see this file's header). [] when the type is not in
%   the model.
%
%   A SPLINE IS ONLY AVERAGED WHERE EVERY BIN USING IT HAS DATA (the pooled
%   values inside plan.pooled's .common range), because a bin's spline is
%   fitted over its own values only, and outside them it is extrapolating,
%   which a spline does badly. Where the bins share no values at all there
%   is nothing a model can hold constant, and each bin is evaluated over its
%   own; both are said in NOTES. A circular spline covers the whole circle
%   and never extrapolates, so it takes every pooled value. A 2D spline is
%   averaged over the pooled PAIRS of its two fields, inside the box of
%   values every type using it shares.
%
%   AN INTERACTION WITH A CONTINUOUS TERM is held at the pooled values too
%   (pooledInteractions): its column is the product of a factor's dummy and
%   the term, and the term's part of it is set to the pooled mean, as the
%   term's own column is. Read from the type's own average row, it carried
%   the bin's own values into the waveform after all.
    waveform = [];
    notes = {};
    eventType = plan.cellTypes{c};
    label = plan.cellLabels{c};
    t = typeIndex(unfold, eventType);
    if isempty(t)
        return;
    end
    rows = strcmp({plan.events.type}, eventType) & plan.kept;   % not those 'drop' left out
    if ~any(rows)
        return;
    end
    cols = reshape(find(unfold.cols2eventtypes == t), 1, []);
    weights = mean(unfold.X(rows, cols), 1);
    variableOf = unfold.cols2variablenames(cols);
    single = cellfun(@isempty, {plan.pooled.pair});
    for v = reshape(unique(variableOf, 'stable'), 1, [])
        kind = unfold.variabletypes{v};
        if ~any(strcmp(kind, {'continuous', 'spline'}))
            continue;
        end
        name = regexprep(unfold.variablenames{v}, '^\d+_', '');
        at = variableOf == v;
        if strcmp(kind, 'continuous')
            % The column holds the value less the centre uf_designmat took
            % off it: 0 with reference coding, the type's own mean with
            % effects coding. Read off the design itself, row by row.
            centre = columnCentre(unfold.X(rows, cols(at)), plan.events(rows), name);
            weights(at) = mean(plan.pooled(single & strcmp({plan.pooled.name}, name)).values) - centre;
            continue;
        end
        spl = splineFor(unfold, unfold.colnames(cols(at)));
        if size(spl.knots, 1) == 2
            [basis, note] = surfaceBasis(spl, plan.pooled(~single & strcmp({plan.pooled.name}, name)), ...
                plan.events(rows), label);
        else
            [basis, note] = splineBasis(spl, plan.pooled(single & strcmp({plan.pooled.name}, name)), ...
                plan.events(rows), name, label);
        end
        notes = [notes, note]; %#ok<AGROW>
        weights(at) = mean(basis, 1);
    end
    interactions = strcmp(unfold.variabletypes(variableOf), 'interaction');
    vars = plan.variables{c};
    if any(interactions) && any(~[vars.categorical] & ~[vars.spline])
        weights(interactions) = pooledInteractions(plan, c, numel(cols), interactions);
    end
    waveform = sum(unfold.beta(:, :, cols) .* reshape(weights, 1, 1, []), 3);
end

function centre = columnCentre(column, events, name)
%COLUMNCENTRE  What uf_designmat subtracted from a continuous term's values
%   to make its column: the value less the column, on every event with
%   both; 0 when there is none to read it from.
    values = arrayfun(@(e) eventNumber(e.(name)), events);
    column = reshape(column, size(values));
    both = isfinite(values) & isfinite(column);
    centre = 0;
    if any(both)
        centre = median(values(both) - column(both));
    end
end

function value = eventNumber(value)
    if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
        value = NaN;
    end
    value = double(value);
end

function requireRebuildableInteractions(unfold, plan)
%REQUIREREBUILDABLEINTERACTIONS  With effects coding, refuse an interaction
%   of a factor with a continuous term. A bin's waveform holds such a term
%   at its pooled value, rebuilding the interaction's columns through
%   uf_designmat on the type's events with that value (pooledInteractions);
%   with effects coding uf_designmat centres a continuous term on the mean of
%   the events it is given, which are then all the pooled value, so the
%   rebuilt columns would be zero rather than the term at that value.
    if ~strcmp(plan.codingSchema, 'effects')
        return;
    end
    for v = find(strcmp(unfold.variabletypes, 'interaction'))
        parts = regexprep(strsplit(unfold.variablenames{v}, ':'), '^\d+_', '');
        t = unique(unfold.cols2eventtypes(unfold.cols2variablenames == v));
        vars = [plan.variables{t}];
        linear = {vars(~[vars.categorical] & ~[vars.spline]).name};
        hit = parts(ismember(parts, linear));
        if ~isempty(hit)
            throw(MException('Alakazam:Unfold:EffectsInteraction', '%s', sprintf([ ...
                'The formula of "%s" has an interaction of a factor with the continuous term ' ...
                '"%s", which cannot be held at its pooled value with effects coding: the ' ...
                'toolbox centres a continuous term on the mean of the events it is given, and ' ...
                'rebuilding the interaction at one value would centre it to zero. Would you use ' ...
                'reference coding for this model, or take the interaction out?'], ...
                plan.typeLabels{t(1)}, hit{1})));
        end
    end
end

function [basis, notes] = splineBasis(spl, pooled, events, name, label)
%SPLINEBASIS  A spline's basis at the pooled values of its field, cut to the
%   range every type using it shares (see referencePrediction).
    notes = {};
    values = pooled.values;
    if ~contains(func2str(spl.splinefunction), 'cyclical')
        if isempty(pooled.common)
            values = [events.(name)];
            values = values(isfinite(values));   % a missing value is not one of them
            notes{end + 1} = sprintf(['The bins using "%s" share no values of it, so it ' ...
                'cannot be held constant across them: "%s" is evaluated over its own.'], ...
                name, label);
        elseif any(values < pooled.common(1) | values > pooled.common(2))
            values = values(values >= pooled.common(1) & values <= pooled.common(2));
            notes{end + 1} = sprintf(['The spline of "%s" is averaged over its values ' ...
                'from %.4g to %.4g, the range every bin using it shares: outside a bin''s ' ...
                'own values its spline is extrapolating.'], name, ...
                pooled.common(1), pooled.common(2));
        end
    end
    basis = spl.splinefunction(values, spl.knots);
    basis(:, spl.removedSplineIdx) = [];
end

function [basis, notes] = surfaceBasis(spl, pooled, events, label)
%SURFACEBASIS  A 2D spline's basis at the pooled pairs of its two fields,
%   built the way uf_designmat_spline builds it (each pair's two bases
%   multiplied out), cut to the box every type using it shares.
    notes = {};
    values = pooled.values;
    fields = sprintf('"%s" and "%s"', pooled.pair{1}, pooled.pair{2});
    if isempty(pooled.common)
        values = [[events.(pooled.pair{1})]; [events.(pooled.pair{2})]];
        values = values(:, all(isfinite(values), 1));   % a pair missing a value is left out
        notes{end + 1} = sprintf(['The bins using the 2D spline of %s share no values of it, so ' ...
            'it cannot be held constant across them: "%s" is evaluated over its own.'], fields, label);
    else
        inside = all(values >= pooled.common(:, 1) & values <= pooled.common(:, 2), 1);
        if ~all(inside)
            values = values(:, inside);
            notes{end + 1} = sprintf(['The 2D spline of %s is averaged over the pairs of values ' ...
                'every bin using it shares (%.4g to %.4g, %.4g to %.4g): outside a bin''s own ' ...
                'values its spline is extrapolating.'], fields, pooled.common(1, 1), ...
                pooled.common(1, 2), pooled.common(2, 1), pooled.common(2, 2));
        end
    end
    first = spl.splinefunction(values(1, :), spl.knots(1, :));
    second = spl.splinefunction(values(2, :), spl.knots(2, :));
    basis = zeros(size(values, 2), size(first, 2) * size(second, 2));
    for k = 1:size(values, 2)
        product = first(k, :)' * second(k, :);
        basis(k, :) = product(:)';
    end
    basis(:, spl.removedSplineIdx) = [];
end

function weights = pooledInteractions(plan, c, ncols, interactions)
%POOLEDINTERACTIONS  The average design row of the interaction columns of
%   type C with every linear continuous field of its formula set to its
%   pooled mean: the type's own events, those fields replaced, through the
%   toolbox's own uf_designmat (Unfold.designMatrix), so the columns are
%   built exactly as the fitted ones were. A factor-by-factor interaction
%   has no such field and comes out as its own average row, as before.
    rows = strcmp({plan.events.type}, plan.cellTypes{c}) & plan.kept;
    events = plan.events(rows);
    vars = plan.variables{c};
    single = cellfun(@isempty, {plan.pooled.pair});
    for v = vars(~[vars.categorical] & ~[vars.spline])
        level = mean(plan.pooled(single & strcmp({plan.pooled.name}, v.name)).values);
        [events.(v.name)] = deal(level);
    end
    one = struct('events', events, 'eventTypes', {plan.cellTypes(c)}, ...
        'formulas', {plan.formulas(c)}, 'typeLabels', {plan.cellLabels(c)}, ...
        'variables', {plan.variables(c)}, 'codingSchema', plan.codingSchema); %#ok<NASGU>  read inside evalc
    designed = [];
    evalc('designed = Unfold.designMatrix(struct(''event'', events), one);');   % quiet: it prints
    if size(designed.unfold.X, 2) ~= ncols
        throw(MException('Alakazam:Unfold:Interaction', '%s', sprintf([ ...
            'Rebuilding the columns of "%s" at the pooled values gave %d columns where the fit ' ...
            'has %d, so its interaction cannot be evaluated. This is a fault in Alakazam, not in ' ...
            'the model.'], plan.cellLabels{c}, size(designed.unfold.X, 2), ncols)));
    end
    weights = mean(designed.unfold.X(:, interactions), 1);
end

function t = typeIndex(unfold, eventType)
%TYPEINDEX  Which entry of EEG.unfold.eventtypes holds EVENTTYPE, [] if
%   none. Not every entry is an event type: glmnet adds a column of its own,
%   'glmnet-DC-Correction', its intercept for the whole recording, whose
%   entry is {NaN} (uf_glmfit, through uf_designmat_addcol).
    t = [];
    for k = 1:numel(unfold.eventtypes)
        names = unfold.eventtypes{k};
        if ~iscell(names)
            names = {names};
        end
        if any(cellfun(@(x) (ischar(x) || isstring(x)) && strcmp(x, eventType), names))
            t = k;
            return;
        end
    end
end

function spl = splineFor(unfold, colnames)
%SPLINEFOR  The entry of EEG.unfold.splines whose columns are COLNAMES, the
%   columns of X it produced. Matched on the columns rather than the name,
%   since two event types can each have a spline of the same field.
    for s = 1:numel(unfold.splines)
        spl = unfold.splines{s};
        if isequal(reshape(cellstr(spl.colnames), 1, []), reshape(cellstr(colnames), 1, []))
            return;
        end
    end
    throw(MException('Alakazam:Unfold:SplineNotFound', '%s', sprintf([ ...
        'No spline in the fitted model produced the columns %s, so the waveform cannot be ' ...
        'evaluated. This is a fault in Alakazam, not in the model.'], strjoin(cellstr(colnames), ', '))));
end

function beta = byColumn(condensed, field)
%BYCOLUMN  uf_condense's betas (FIELD: 'beta' for a deconvolution,
%   'beta_nodc' for a regression on epochs; channels x samples x terms) put
%   back on the columns of EEG.unfold.X they belong to. uf_condense returns
%   one entry for each column whose variable has a type, in column order;
%   that leaves out a column a solver adds of its own (glmnet's intercept for
%   the whole recording, whose variable type is NaN), which is NaN here and
%   is no event type's column.
    unfold = condensed.unfold;
    typed = find(cellfun(@(type) ~all(isnan(type)), unfold.variabletypes));
    cols = find(ismember(unfold.cols2variablenames, typed));
    values = condensed.(field);
    if size(values, 3) ~= numel(cols)
        throw(MException('Alakazam:Unfold:Condensed', '%s', sprintf([ ...
            'uf_condense returned %d terms where the design has %d typed columns, so its ' ...
            'betas cannot be placed on the design. This is a fault in Alakazam, not in the ' ...
            'model.'], size(values, 3), numel(cols))));
    end
    beta = nan(size(values, 1), size(values, 2), size(unfold.X, 2));
    beta(:, :, cols) = values;
end

function [EEG, info] = packageTerms(input, plan, result, info, opts, times)
%PACKAGETERMS  One waveform per model term of every binned event type, by
%   Unfold's own route (see this file's header), in Average's shape, from
%   RESULT, the fit as uf_condense returns it.
    args = {'auto_method', opts.EvaluateSpacing, 'auto_n', 10};   % 'quantile' and 10, the toolbox's
    predictAt = Unfold.predictionValues(opts.EvaluateAt, result.unfold, plan.typeLabels);
    if ~isempty(predictAt)
        args = [args, {'predictAt', predictAt}];
    end
    lastwarn('');
    % 'MEM' (the toolbox's default) adds the other terms at their mean
    % value, 'AME' as their average marginal effect (see this file's header).
    marginal = uf_addmarginal(uf_predictContinuous(result, args{:}), 'type', opts.Marginal);
    [warningText, ~] = lastwarn();
    if contains(lower(warningText), 'interaction')
        plan.notes{end + 1} = ['The model has interactions, which uf_addmarginal does not ' ...
            'fold into the other terms: an interaction term is its own beta plus the ' ...
            'intercept, not the response to that combination of levels.'];
        info.notes = plan.notes;
    end

    events = arrayfun(@(p) char(string(p.event)), marginal.param, 'UniformOutput', false);
    labelOf = containers.Map(plan.eventTypes, plan.typeLabels);
    keep = find(ismember(events, plan.cellTypes));
    labels = cell(1, numel(keep));
    counts = zeros(1, numel(keep));
    surfaces = plan.pooled(~cellfun(@isempty, {plan.pooled.pair}));
    terms = struct('label', {}, 'bin', {}, 'name', {}, 'type', {}, 'value', {});
    for j = 1:numel(keep)
        p = marginal.param(keep(j));
        binLabel = labelOf(events{keep(j)});
        labels{j} = termLabel(binLabel, p, ...
            referenceLevels(result.unfold, plan.events, events{keep(j)}), surfaces);
        counts(j) = nnz(strcmp({plan.events.type}, events{keep(j)}) & plan.kept);
        terms(j) = struct('label', labels{j}, 'bin', binLabel, 'name', char(string(p.name)), ...
            'type', char(string(p.type)), 'value', double(p.value));
    end

    field = 'beta';          % uf_condense's name for a deconvolution's betas,
    if ~isfield(marginal, field)
        field = 'beta_nodc';     % and for a regression on epochs'
    end
    data = marginal.(field)(:, :, keep);
    requireFiniteBetas(data, true(1, numel(keep)), struct('label', labels), opts);
    data = applyBaseline(data, true(1, numel(keep)), times, opts.BaselineMs);

    EEG = input;
    EEG.data = data;
    EEG.times = times;
    EEG.pnts = numel(times);
    EEG.trials = 1;
    EEG.xmin = times(1) / 1000;
    EEG.xmax = times(end) / 1000;
    EEG.DataFormat = "Averaged";
    EEG.ntrials = sum(plan.cellCounts);
    EEG.stErr = zeros(size(data));      % see this file's header: there is none
    EEG.aSME = nan(size(data, 1), numel(keep));
    EEG.bindesc = struct('index', num2cell(1:numel(keep)), 'label', labels, ...
        'combo', repmat({[]}, 1, numel(keep)), 'n', num2cell(counts));

    info.output = 'terms';
    info.terms = terms;
    info.evaluateAt = char(string(opts.EvaluateAt));
    info.evaluateSpacing = opts.EvaluateSpacing;
    info.marginal = opts.Marginal;
end

function text = referenceLevels(unfold, events, eventType)
%REFERENCELEVELS  ", factor = level" for each factor of EVENTTYPE: the level
%   the intercept stands for, which is the one without a column of its own
%   (Unfold's reference coding makes the first level, in sorted order, the
%   reference); with effects coding ", factor at the mean of its levels".
%   '' for a type without factors.
    text = '';
    effects = strcmp(char(string(TransTools.FieldOr(unfold, 'codingschema', 'reference'))), 'effects');
    t = typeIndex(unfold, eventType);
    if isempty(t)
        return;
    end
    cols = reshape(find(unfold.cols2eventtypes == t), 1, []);
    variables = reshape(unique(unfold.cols2variablenames(cols), 'stable'), 1, []);
    rows = strcmp({events.type}, eventType);
    for v = variables
        if ~strcmp(unfold.variabletypes{v}, 'categorical')
            continue;
        end
        name = regexprep(unfold.variablenames{v}, '^\d+_', '');
        if ~isfield(events, name)
            continue;
        end
        if effects
            text = sprintf('%s, %s at the mean of its levels', text, name);
            continue;
        end
        levels = unique(cellfun(@(x) char(string(x)), {events(rows).(name)}, 'UniformOutput', false));
        own = regexprep(unfold.colnames(cols(unfold.cols2variablenames(cols) == v)), ['^(\d+_)?' name '_'], '');
        reference = setdiff(levels, own);
        if ~isempty(reference)
            text = sprintf('%s, %s = %s', text, name, reference{1});
        end
    end
end

function label = termLabel(binLabel, p, references, surfaces)
%TERMLABEL  "<bin>: <term>" for one of uf_condense's parameters, without the
%   "2_" Unfold puts before a second event type's names: the intercept with
%   the factor levels it stands for, a continuous or spline term with the
%   value it was evaluated at, a 2D spline with the value of each of its
%   fields (SURFACES, Unfold.binModel's pooled entries with a .pair, name
%   them: the toolbox's own name runs the two together), anything else by
%   Unfold's own column name.
    name = regexprep(char(string(p.name)), '^\d+_', '');
    surface = surfaces(strcmp({surfaces.name}, name));
    if strcmpi(name, '(Intercept)')
        label = sprintf('%s: (Intercept)%s', binLabel, references);
    elseif contains(char(string(p.type)), 'converted') && numel(p.value) == 2 && ~isempty(surface)
        label = sprintf('%s: %s = %.4g, %s = %.4g', binLabel, surface(1).pair{1}, p.value(1), ...
            surface(1).pair{2}, p.value(2));
    elseif contains(char(string(p.type)), 'converted') && isfinite(p.value(1))
        label = sprintf('%s: %s = %.4g', binLabel, name, p.value(1));
    else
        label = sprintf('%s: %s', binLabel, name);
    end
end

function data = resolveComboBins(data, bindesc)
%RESOLVECOMBOBINS  Difference bins, as Average computes them: the signed sum
%   of the bins they reference, in the dependency order TransTools.ComboOrder
%   works out for both, so a difference of differences works too. They are
%   not predictors (see Unfold.binModel); subtracting two fitted waveforms is
%   the same operation subtracting two averages is. A combination that
%   cannot be resolved is left NaN, which DefineBins' own parse-time check
%   makes unreachable from a script.
    for s = TransTools.ComboOrder(bindesc)
        acc = zeros(size(data, 1), size(data, 2));
        for t = 1:numel(s.parts)
            acc = acc + s.coeffs(t) * data(:, :, s.parts(t));
        end
        data(:, :, s.target) = acc;
    end
end
