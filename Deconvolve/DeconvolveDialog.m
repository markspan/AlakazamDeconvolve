function options = DeconvolveDialog(EEG, stored)
%DECONVOLVEDIALOG  Modal editor for Deconvolve's settings.
%   OPTIONS = DeconvolveDialog(EEG, STORED) shows the bins to fit, each
%   bin's formula, the response window, the artefact threshold, which event
%   codes outside every bin are modelled, and what comes out, seeded from
%   STORED (a previous run's options, or [] on first use), and returns the
%   options Deconvolve takes, or [] on Cancel.
%
%   THE BINS ARE A SETTING HERE, not a property of the dataset. Deconvolve
%   runs on the continuous recording (see Deconvolve for why it cannot be
%   otherwise), and a continuous recording does not normally carry bins, so
%   "Define bins..." opens DefineBins' own editor, in DefineBins' own
%   language, with its epoch fields hidden because nothing is being cut here.
%   A dataset that does already carry tags is used as it stands until a script
%   is given.
%
%   EACH BIN HAS A FORMULA, typed in Unfold's own notation in the table, one
%   row per bin: 'y ~ 1' by default, or with linear terms, factors (cat),
%   splines (spl, circspl, 2dspl) and interactions. The bins say which events are
%   fitted together; the formula says what explains their response. Beside
%   it is what a formula can use: this recording's event fields, by kind and
%   unit, and the notation itself. Options stored before formulas existed
%   carry a list of covariates instead; those are shown here as the formulas
%   they always meant, and OK stores the formulas.
%
%   IT SHOWS THE MODEL IT WOULD FIT, not just the numbers. The list names
%   every bin with its event count, the columns the toolbox builds for each
%   formula (Unfold.designMatrix on the events alone, so the preview is the
%   same call the fit makes), the nuisance event types, and anything
%   Unfold.binModel has to say about the design (a bin with no events, a
%   field a formula names that a bin's events lack). That is the one thing
%   a user cannot work out from the dialog's fields, and the one thing that
%   decides whether the answer will mean anything. It refreshes when the
%   bins, a formula or the chosen event codes change, because all three
%   change the model.
%
%   IT WARNS ABOUT EVENTS LOCKED TOGETHER. A code in no bin that is modelled
%   and keeps a near-constant lag to a bin or another modelled code (a
%   fixation 12 ms after its saccade) makes the design nearly collinear: the
%   solver may not converge, and the two responses are not told apart (see
%   Unfold.timeLockedEvents). Such a pair is listed in the model and raised in
%   an alert whenever the model gains one, since unticking the code is the
%   usual remedy and is made right here (manual issue M10). The alert does not
%   stop OK: a lock can be what the user means to model.
%
%   IT STARTS FROM THE TOOLBOX'S DEFAULTS where Unfold has one: the artefact
%   scan on every channel, the terms with the other terms at their mean
%   value ('MEM'), the betas as fitted, without a baseline, and a missing
%   number filled in with the median (uf_imputeMissing). Each is a choice
%   here: the scalp EEG only, the average marginal effect ('AME'), a
%   baseline window, and the mean, a random draw, leaving the event out or
%   refusing. So are the toolbox's other ways to fit: MATLAB's exact solver
%   or glmnet's regularised one instead of lsmr, and a regression on epochs
%   (uf_epoch, uf_glmfit_nodc) instead of a deconvolution.
%
%   OK checks the settings the way Deconvolve will apply them, so a design
%   the fit would refuse is reported here rather than after the dialog closes.
%
%   See also DECONVOLVE, DEFINEBINSDIALOG, UNFOLD.BINMODEL, UNFOLD.FITBINS,
%   UNFOLD.DESIGNMATRIX.
    options = [];

    seed = defaults();
    if isstruct(stored)
        for name = fieldnames(seed)'
            if isfield(stored, name{1}) && ~isempty(stored.(name{1}))
                seed.(name{1}) = stored.(name{1});
            end
        end
    end

    % An option set stored before the scan's channels, the marginal effect and
    % missing values were choices was run with the scalp EEG, 'AME' and a
    % refusal, and a replay still runs it so (see Deconvolve), so the dialog
    % shows it so too. Only a first use starts from the toolbox's defaults.
    if isstruct(stored) && ~isfield(stored, 'artifactChannels')
        seed.artifactChannels = 'scalp';
    end
    if isstruct(stored) && ~isfield(stored, 'marginal')
        seed.marginal = 'AME';
    end
    if isstruct(stored) && ~isfield(stored, 'missingValues')
        seed.missingValues = 'refuse';
    end
    if isstruct(stored) && ~isfield(stored, 'excludeAroundCuts')
        seed.excludeAroundCuts = true;      % what Deconvolve did before it was a choice
    end

    % Which unbinned codes to model, kept out of the generic seeding loop for
    % the same reason as the baseline: an empty list is an answer ("none"),
    % and the loop would read it as absent and restore the default. 'all' is
    % stored as the word, not as today's list of codes, so a replay on another
    % recording still means every code that recording has.
    otherSelection = Unfold.otherEventsChoice(stored);

    % The baseline is stored as the window itself, with [] meaning "leave the
    % betas as the solver returned them", so the checkbox and the two fields
    % are read back out of that one value. Off on first use, as Unfold
    % returns its betas; the fields then offer the pre-event window.
    baselineOn = false;
    baseline = [seed.windowMs(1) 0];
    if isstruct(stored) && ~isfield(stored, 'baselineMs')
        baselineOn = true;      % absent, Deconvolve corrects over the pre-event window
    elseif isstruct(stored)
        if ~isempty(stored.baselineMs)
            baselineOn = true;
            baseline = reshape(double(stored.baselineMs), 1, 2);
        end
    end
    if baseline(1) >= baseline(2)
        baseline = [seed.windowMs(1) seed.windowMs(1) + 100];
        baselineOn = false;     % a window starting at the event has no baseline
    end

    hasTags = datasetCarriesBins(EEG);
    binScript = strtrim(char(string(seed.binScript)));
    if isempty(binScript) && ~hasTags
        % Nothing stored and nothing on the dataset, so the bins have to come
        % from somewhere: the script last run in this workspace is almost
        % certainly the one meant here, since the usual way into this dialog
        % is having already binned and averaged the same recording.
        binScript = lastDefineBinsScript();
    end

    % Every formula typed in this session, by bin label, including those of
    % bins the script no longer has: rewriting the bins and then writing them
    % back should not cost the formulas that went with them.
    remembered = storedFormulas(seed.formulas);
    % The older covariates are turned into formulas once, the first time
    % there are bins to turn them into (see seedFormulaTable).
    legacyCovariates = cellstr(string(seed.covariates));
    if ~isempty(remembered)
        legacyCovariates = {};   % the formulas already say what they meant
    end
    selectedRow = 1;
    lockedWarned = {};      % the time-locked pairs last alerted about (see warnAboutLockedEvents)

    [accentColor, bgColor] = dialogChromeColors();
    fig = uifigure('Name', 'Deconvolve', 'Position', fitOnScreen([80 40 1200 920]), 'Color', bgColor);
    root = uigridlayout(fig, [2 1], 'RowHeight', {40, '1x'}, 'Padding', [0 0 0 0], 'RowSpacing', 0);
    titleBar = uilabel(root, 'Text', '  Deconvolve: overlapping responses separated', 'FontSize', 14, ...
        'FontWeight', 'bold', 'FontColor', [1 1 1], 'BackgroundColor', accentColor, ...
        'VerticalAlignment', 'center');
    if exist('manualButton', 'file') == 2   % the ? to the manual; not in an older Alakazam
        manualButton(titleBar, 'Deconvolve');
    end
    outer = uigridlayout(root, [6 1], 'RowHeight', {'fit', 'fit', 'fit', '1x', 'fit', 44});

    uilabel(outer, 'WordWrap', 'on', 'Text', [ ...
        'Fits every bin at once against the whole continuous recording, so where two events are ' ...
        'close enough for their responses to overlap, each bin keeps its own and gives up the ' ...
        'other''s. The result is one waveform per bin, in the shape Average produces, one ' ...
        'overlap-corrected trial per event, in the shape DefineBins cuts, or one waveform per ' ...
        'term of the model, so the rest of Alakazam reads any of them unchanged.']);

    binsRow = uigridlayout(outer, [1 2], 'ColumnWidth', {'1x', 130}, 'Padding', [0 4 0 4]);
    binsLabel = uilabel(binsRow, 'WordWrap', 'on', 'Text', '');
    uibutton(binsRow, 'Text', 'Define bins...', 'ButtonPushedFcn', @(~, ~) onDefineBins(), ...
        'Tooltip', 'Write the bins to fit, in DefineBins'' language');

    % Two columns: what the design and the result are on the left, how the
    % model is fitted and checked on the right.
    settingsRow = uigridlayout(outer, [1 2], 'ColumnWidth', {'fit', 'fit'}, ...
        'Padding', [0 0 0 0], 'ColumnSpacing', 24);
    settings = uigridlayout(settingsRow, [12 4], 'ColumnWidth', {190, 90, 210, 90}, ...
        'RowHeight', repmat({'fit'}, 1, 12), 'Padding', [0 0 0 0], 'RowSpacing', 4);
    fitting = uigridlayout(settingsRow, [8 4], 'ColumnWidth', {150, 80, 190, 80}, ...
        'RowHeight', repmat({'fit'}, 1, 8), 'Padding', [0 0 0 0], 'RowSpacing', 4);
    uilabel(settings, 'Text', 'Window start (ms):');
    % The window decides which events overlap, so the model is shown again.
    startField = uieditfield(settings, 'numeric', 'Value', seed.windowMs(1), ...
        'ValueChangedFcn', @(~, ~) showModel());
    uilabel(settings, 'Text', 'Artefact threshold (uV, 0 = off):');
    thresholdField = uieditfield(settings, 'numeric', 'Value', seed.artifactThresholdUv, ...
        'Limits', [0 Inf]);
    uilabel(settings, 'Text', 'Window stop (ms):');
    stopField = uieditfield(settings, 'numeric', 'Value', seed.windowMs(2), ...
        'ValueChangedFcn', @(~, ~) showModel());
    uilabel(settings, 'Text', 'Measured in a window of (ms):');
    artWindowField = uieditfield(settings, 'numeric', 'Value', seed.artifactWindowMs, ...
        'Limits', [0 Inf], 'LowerLimitInclusive', 'off');
    uilabel(settings, 'Text', 'Baseline start (ms):');
    baseStartField = uieditfield(settings, 'numeric', 'Value', baseline(1));
    uilabel(settings, 'Text', 'Stepped by (ms):');
    artStepField = uieditfield(settings, 'numeric', 'Value', seed.artifactStepMs, ...
        'Limits', [0 Inf], 'LowerLimitInclusive', 'off');
    uilabel(settings, 'Text', 'Baseline stop (ms):');
    baseStopField = uieditfield(settings, 'numeric', 'Value', baseline(2));

    % A beta's zero is wherever the model put it, so a fitted waveform has to
    % be baseline-corrected before it can be read beside an average that
    % Baseline has already corrected. Off until ticked, as Unfold leaves its
    % betas, with the pre-event part of the response window offered.
    baselineBox = uicheckbox(settings, 'Text', 'Baseline-correct the result', ...
        'Value', baselineOn, 'ValueChangedFcn', @(~, ~) onBaselineToggled());
    baselineBox.Layout.Row = 4;
    baselineBox.Layout.Column = [3 4];
    onBaselineToggled();

    % What comes out. The waveforms are the same either way (Average of the
    % trials gives them back); the trials are for looking at them one by one
    % and for the noise figures only trials can give; the terms are the
    % model's own view, one waveform per factor level and per value of a
    % continuous term.
    outputLabel = uilabel(settings, 'Text', 'Result:');
    outputLabel.Layout.Row = 5;
    outputLabel.Layout.Column = 1;
    outputDropdown = uidropdown(settings, 'Tag', 'output', ...
        'Items', {'One waveform per bin (as Average gives)', ...
                  'Overlap-corrected trials (as DefineBins cuts)', ...
                  'One waveform per model term'}, ...
        'ItemsData', {'average', 'trials', 'terms'}, 'Value', outputSeed(seed.output), ...
        'ValueChangedFcn', @(~, ~) onOutputChanged(), ...
        'Tooltip', ['Trials hold each event''s recording with every other event''s fitted ' ...
         'response subtracted: run Average on them for the waveforms, and EpochView shows ' ...
         'them as an ERP image without the overlap. Model terms are every factor level and ' ...
         'every continuous or spline term at chosen values, each a whole waveform with the ' ...
         'other terms added as set below.']);
    outputDropdown.Layout.Row = 5;
    outputDropdown.Layout.Column = [2 4];
    evaluateLabel = uilabel(settings, 'Text', 'Terms evaluated at:');
    evaluateLabel.Layout.Row = 6;
    evaluateLabel.Layout.Column = 1;
    evaluateField = uieditfield(settings, 'text', 'Tag', 'evaluateAt', ...
        'Value', char(string(seed.evaluateAt)), ...
        'Placeholder', 'e.g. sac_amplitude = 0.5 1 2 4; rt = 300 500', ...
        'Tooltip', ['Where each continuous or spline term is drawn, as "name = values", ' ...
         'separated by semicolons. A term not named here is drawn at ten quantiles ' ...
         'of its own values, the toolbox''s default.']);
    evaluateField.Layout.Row = 6;
    evaluateField.Layout.Column = [2 4];

    % Where a term not named above is evaluated: uf_predictContinuous's own
    % 'auto_method', ten quantiles of its values unless asked otherwise.
    spacingLabel = uilabel(settings, 'Text', 'A term not named, at:');
    spacingLabel.Layout.Row = 7;
    spacingLabel.Layout.Column = 1;
    spacingDropdown = uidropdown(settings, 'Tag', 'evaluateSpacing', ...
        'Items', {'Ten quantiles of its values (the toolbox''s default)', ...
                  'Ten values evenly from its minimum to its maximum', ...
                  'Its mean alone'}, ...
        'ItemsData', {'quantile', 'linear', 'average'}, ...
        'Value', choiceSeed(seed.evaluateSpacing, {'quantile', 'linear', 'average'}), ...
        'Tooltip', ['Where a continuous or spline term not named in "Terms evaluated at" is ' ...
         'drawn: the toolbox''s uf_predictContinuous takes ten quantiles of its values, ' ...
         'ten evenly spaced values, or the mean alone.']);
    spacingDropdown.Layout.Row = 7;
    spacingDropdown.Layout.Column = [2 4];

    % How each term's waveform carries the rest of the model: uf_addmarginal's
    % own default, or the average over the events' values. They differ only
    % for a spline, and most for a circular one.
    marginalLabel = uilabel(settings, 'Text', 'Other terms added at:');
    marginalLabel.Layout.Row = 8;
    marginalLabel.Layout.Column = 1;
    marginalDropdown = uidropdown(settings, 'Tag', 'marginal', ...
        'Items', {'Their mean value (MEM, the toolbox''s default)', ...
                  'Their average over the events (AME)'}, ...
        'ItemsData', {'MEM', 'AME'}, 'Value', choiceSeed(seed.marginal, {'MEM', 'AME'}), ...
        'Tooltip', ['Each term''s waveform has the model''s other terms added in. MEM takes ' ...
         'a spline at its mean value; AME averages it over the events'' own values. For an ' ...
         'angle (circspl) the mean can be a direction no event had.']);
    marginalDropdown.Layout.Row = 8;
    marginalDropdown.Layout.Column = [2 4];
    onOutputChanged();

    % How a factor is coded in the design: uf_designmat's 'codingschema'.
    % The design is shown again, since the columns are named by it.
    codingLabel = uilabel(settings, 'Text', 'Factors coded as:');
    codingLabel.Layout.Row = 9;
    codingLabel.Layout.Column = 1;
    codingDropdown = uidropdown(settings, 'Tag', 'codingSchema', ...
        'Items', {'Reference: the intercept is the first level (the toolbox''s default)', ...
                  'Effects: the intercept is the mean over the levels'}, ...
        'ItemsData', {'reference', 'effects'}, ...
        'Value', choiceSeed(seed.codingSchema, {'reference', 'effects'}), ...
        'ValueChangedFcn', @(~, ~) showModel(), ...
        'Tooltip', ['How the toolbox''s uf_designmat codes a factor (cat() in a formula). ' ...
         'A bin''s waveform is the same either way; the terms differ: with effects coding ' ...
         'the intercept is the mean over a factor''s levels, and each continuous term is ' ...
         'centred on its events'' mean.']);
    codingDropdown.Layout.Row = 9;
    codingDropdown.Layout.Column = [2 4];

    % Which channels the artefact scan looks at: the toolbox scans them all.
    channelsLabel = uilabel(settings, 'Text', 'Artefact scan on:');
    channelsLabel.Layout.Row = 10;
    channelsLabel.Layout.Column = 1;
    channelsDropdown = uidropdown(settings, 'Tag', 'artifactChannels', ...
        'Items', {'Every channel (the toolbox''s default)', 'The scalp EEG only'}, ...
        'ItemsData', {'all', 'scalp'}, 'Value', choiceSeed(seed.artifactChannels, {'all', 'scalp'}), ...
        'Tooltip', ['An eye channel swings several times as far as the EEG, so scanning it ' ...
         'can leave out much of the recording. The scalp EEG only leaves the eye and other ' ...
         'peripheral channels out of the scan; they are still fitted.']);
    channelsDropdown.Layout.Row = 10;
    channelsDropdown.Layout.Column = [2 4];

    % What happens around a cut in the recording: the toolbox fits across
    % the join as it is; the alternative leaves a response window out on
    % each side.
    cutsLabel = uilabel(settings, 'Text', 'Around cuts:');
    cutsLabel.Layout.Row = 11;
    cutsLabel.Layout.Column = 1;
    cutsDropdown = uidropdown(settings, 'Tag', 'excludeAroundCuts', ...
        'Items', {'Fitted as they are (the toolbox''s default)', ...
                  'A response window left out either side'}, ...
        'ItemsData', {'off', 'on'}, 'Value', onOff(seed.excludeAroundCuts), ...
        'Tooltip', ['A cut (a boundary event) joins two stretches of recording that were ' ...
         'never adjacent. The toolbox fits across it as if they were; leaving out the ' ...
         'longer side of the response window on each side keeps every event''s fitted ' ...
         'window off the join. A trial whose window spans a cut is left out either way.']);
    cutsDropdown.Layout.Row = 11;
    cutsDropdown.Layout.Column = [2 4];

    % An event without a number its formula uses: the toolbox's own
    % uf_imputeMissing, whose default is the median. The model is shown again,
    % since a refusal, or the events 'drop' leaves out, depend on it.
    missingLabel = uilabel(settings, 'Text', 'Missing values:');
    missingLabel.Layout.Row = 12;
    missingLabel.Layout.Column = 1;
    missingDropdown = uidropdown(settings, 'Tag', 'missingValues', ...
        'Items', {'Filled in with the median (the toolbox''s default)', ...
                  'Filled in with the mean', ...
                  'Filled in with a random one of the others (marginal)', ...
                  'Those events left out (drop)', ...
                  'Refused, naming the bin'}, ...
        'ItemsData', {'median', 'mean', 'marginal', 'drop', 'refuse'}, ...
        'Value', choiceSeed(seed.missingValues, {'median', 'mean', 'marginal', 'drop', 'refuse'}), ...
        'ValueChangedFcn', @(~, ~) showModel(), ...
        'Tooltip', ['An event whose formula names a field it has no number for. The toolbox ' ...
         'fills the number in from the bin''s other events, or leaves the event out of the ' ...
         'model, which then does not overlap-correct it either. A factor without a level ' ...
         'is refused whatever is chosen: the toolbox cannot fill in a level.']);
    missingDropdown.Layout.Row = 12;
    missingDropdown.Layout.Column = [2 4];

    % Deconvolution, or the toolbox's other route: a regression on epochs,
    % the same design without the neighbours' overlap taken out. The model
    % is shown again, since only a deconvolution has pairs to warn about.
    overlapLabel = uilabel(fitting, 'Text', 'Overlap correction:');
    overlapLabel.Layout.Row = 1;
    overlapLabel.Layout.Column = 1;
    overlapDropdown = uidropdown(fitting, 'Tag', 'overlapCorrection', ...
        'Items', {'On: fitted against the continuous recording (deconvolution)', ...
                  'Off: a regression on epochs (uf_epoch, uf_glmfit_nodc)'}, ...
        'ItemsData', {'on', 'off'}, 'Value', onOff(seed.overlapCorrection), ...
        'ValueChangedFcn', @(~, ~) onFitChanged(true), ...
        'Tooltip', ['Off fits the same formulas to epochs, each event''s window on its own, ' ...
         'as the toolbox does to compare: with y ~ 1 a bin is the mean of its epochs. An ' ...
         'epoch touching an artefact or a cut, or running off the recording, is left out.']);
    overlapDropdown.Layout.Row = 1;
    overlapDropdown.Layout.Column = [2 4];

    % How each response is expanded over the window: the toolbox's
    % uf_timeexpandDesignmat, one beta per sample unless asked otherwise.
    basisLabel = uilabel(fitting, 'Text', 'Time basis:');
    basisLabel.Layout.Row = 2;
    basisLabel.Layout.Column = 1;
    basisDropdown = uidropdown(fitting, 'Tag', 'timeBasis', ...
        'Items', {'One beta per sample (stick, the toolbox''s default)', 'Cubic splines', ...
                  'A Fourier set'}, ...
        'ItemsData', {'stick', 'splines', 'fourier'}, ...
        'Value', choiceSeed(seed.timeBasis, {'stick', 'splines', 'fourier'}), ...
        'ValueChangedFcn', @(~, ~) onBasisChanged(), ...
        'Tooltip', ['How the toolbox''s uf_timeexpandDesignmat spreads each response over the ' ...
         'window: one beta per sample, or a number of cubic splines or Fourier frequencies, ' ...
         'which share their betas between neighbouring samples and give smoother waveforms ' ...
         'from fewer parameters. With overlap correction only.']);
    basisDropdown.Layout.Row = 2;
    basisDropdown.Layout.Column = [2 3];
    basisCountField = uieditfield(fitting, 'numeric', 'Tag', 'timeBasisCount', ...
        'Value', seed.timeBasisCount, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
        'Tooltip', 'How many splines or frequencies (the toolbox''s default is 30).');
    basisCountField.Layout.Row = 2;
    basisCountField.Layout.Column = 4;
    onBasisChanged();

    solverLabel = uilabel(fitting, 'Text', 'Solver:');
    solverLabel.Layout.Row = 3;
    solverLabel.Layout.Column = 1;
    solverDropdown = uidropdown(fitting, 'Tag', 'solver', ...
        'Items', {'The toolbox''s default (lsmr; pinv on epochs)', ...
                  'Exact (MATLAB''s own solver; needs much memory)', ...
                  'Regularised (glmnet, cross-validated)'}, ...
        'ItemsData', {'toolbox', 'matlab', 'glmnet'}, ...    % not 'default': graphics reads
        'Value', choiceSeed(strrep(char(string(seed.solver)), 'default', 'toolbox'), ...   % that as "reset"
            {'toolbox', 'matlab', 'glmnet'}), ...
        'ValueChangedFcn', @(~, ~) onFitChanged(false), ...
        'Tooltip', ['The exact solver needs no iteration limit, but for a long recording it ' ...
         'can need tens of gigabytes, as the toolbox warns. glmnet shrinks the betas towards ' ...
         'zero by an amount chosen by cross-validation: lasso, ridge or in between.']);
    solverDropdown.Layout.Row = 3;
    solverDropdown.Layout.Column = [2 4];

    % The solver's limit: the toolbox's own 400 unless a fit said it ran out.
    iterationsLabel = uilabel(fitting, 'Text', 'Solver iterations:');
    iterationsLabel.Layout.Row = 4;
    iterationsLabel.Layout.Column = 1;
    iterationsField = uieditfield(fitting, 'numeric', 'Tag', 'solverIterations', ...
        'Value', seed.solverIterations, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
        'Tooltip', ['The most iterations the solver takes per channel (the toolbox''s own ' ...
         'default is 400). A fit that runs out says so; raise this when it does and the ' ...
         'design is otherwise sound.']);
    iterationsField.Layout.Row = 4;
    iterationsField.Layout.Column = 2;
    alphaLabel = uilabel(fitting, 'Text', 'glmnet alpha (1 lasso, 0 ridge):');
    alphaLabel.Layout.Row = 4;
    alphaLabel.Layout.Column = 3;
    alphaField = uieditfield(fitting, 'numeric', 'Tag', 'glmnetAlpha', ...
        'Value', seed.glmnetAlpha, 'Limits', [0 1], ...
        'Tooltip', 'glmnet''s mix of penalties: 1 is lasso (the toolbox''s default), 0 ridge.');
    alphaField.Layout.Row = 4;
    alphaField.Layout.Column = 4;

    % How well the model explains the recording: the toolbox's own
    % uf_checkmodelfit, which it computes only when asked, and the event
    % code its cross-validation cuts the folds at.
    NO_FOLD = '<none>';   % not '': the dropdown's value for "none chosen"
    fitLabel = uilabel(fitting, 'Text', 'Model fit:');
    fitLabel.Layout.Row = 5;
    fitLabel.Layout.Column = 1;
    fitDropdown = uidropdown(fitting, 'Tag', 'modelFit', ...
        'Items', {'Not measured (the toolbox''s default)', 'R2 of the fit', ...
                  'Partial R2 of each term (refits once per term)', ...
                  'Cross-validated R2 (refits once per fold)', ...
                  'Cross-validated partial R2'}, ...
        'ItemsData', {'none', 'R2', 'partialR2', 'crossValR2', 'crossValpartialR2'}, ...
        'Value', choiceSeed(seed.modelFit, {'none', 'R2', 'partialR2', 'crossValR2', 'crossValpartialR2'}), ...
        'ValueChangedFcn', @(~, ~) onFitMeasureChanged(), ...
        'Tooltip', ['How much of the recording the model explains, per channel, by the ' ...
         'toolbox''s uf_checkmodelfit, over the samples it covers. Partial R2 refits the ' ...
         'model without each term in turn; cross-validation fits each fold on the rest. ' ...
         'Needs overlap correction.']);
    fitDropdown.Layout.Row = 5;
    fitDropdown.Layout.Column = [2 4];
    foldLabel = uilabel(fitting, 'Text', 'Folds cut at:');
    foldLabel.Layout.Row = 6;
    foldLabel.Layout.Column = 1;
    codes = unique(arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false));
    foldDropdown = uidropdown(fitting, 'Tag', 'foldEvent', ...
        'Items', [{'(none chosen)'}, codes], 'ItemsData', [{NO_FOLD}, codes], ...
        'Value', choiceSeed(char(strjoin(cellstr(string(seed.foldEvent)), '')), [{NO_FOLD}, codes]), ...
        'Tooltip', ['For cross-validation: the event code the recording is cut into folds at, ' ...
         'such as a break between blocks. The toolbox refuses a code that falls where an ' ...
         'event''s response is modelled.']);
    foldDropdown.Layout.Row = 6;
    foldDropdown.Layout.Column = [2 4];
    onFitMeasureChanged();

    % Collinearity: the toolbox's own uf_vif, which it computes only when
    % asked, of the design or of the time-expanded design.
    vifLabel = uilabel(fitting, 'Text', 'Collinearity (VIF):');
    vifLabel.Layout.Row = 7;
    vifLabel.Layout.Column = 1;
    vifDropdown = uidropdown(fitting, 'Tag', 'collinearity', ...
        'Items', {'Not checked (the toolbox''s default)', ...
                  'Of the time-expanded design (much memory for long windows)'}, ...
        'ItemsData', {'none', 'expanded'}, ...
        'Value', choiceSeed(seed.collinearity, {'none', 'expanded'}), ...
        'Tooltip', ['The variance inflation factor of every column of the time-expanded ' ...
         'design, by the toolbox''s uf_vif: what the overlap does to telling the bins ' ...
         'apart. Splines give high values that are usually no problem.']);
    vifDropdown.Layout.Row = 7;
    vifDropdown.Layout.Column = [2 4];

    % How the overlap-corrected trials are made: Alakazam's arithmetic (the
    % default, first, so choiceSeed falls back to it), which gives the same
    % trials far sooner on a long recording, or the toolbox's uf_erpimage,
    % channel by channel (Unfold.fitBins, TrialsMethod).
    trialsLabel = uilabel(fitting, 'Text', 'Trials made by:');
    trialsLabel.Layout.Row = 8;
    trialsLabel.Layout.Column = 1;
    trialsDropdown = uidropdown(fitting, 'Tag', 'trialsMethod', ...
        'Items', {'Alakazam: less compatible, but much faster', ...
                  'The toolbox''s uf_erpimage, channel by channel'}, ...
        'ItemsData', {'fast', 'toolbox'}, ...
        'Value', choiceSeed(seed.trialsMethod, {'fast', 'toolbox'}), ...
        'Tooltip', ['For overlap-corrected trials. uf_erpimage epochs the whole recording ' ...
         'around every modelled event, three times for each channel, which takes long on a ' ...
         'long recording with many events. Alakazam computes the same trials only at their ' ...
         'own samples and for every channel at once; it agrees with the toolbox to rounding ' ...
         'error, but it is Alakazam''s arithmetic, not the toolbox''s call.']);
    trialsDropdown.Layout.Row = 8;
    trialsDropdown.Layout.Column = [2 4];
    onFitChanged(false);

    % The formulas, the model they make, and what a formula can use share
    % the stretchy row: a formula is only readable next to the columns it
    % turns into and the fields it can name. The right-hand column shares
    % the width rather than taking a fixed strip: its lines (a field with
    % its unit and count, an example formula, an event code with its count)
    % wrapped at every word in a narrow one, and it grows with the window.
    middle = uigridlayout(outer, [1 2], 'ColumnWidth', {'1x', '1x'}, ...
        'Padding', [0 0 0 0], 'ColumnSpacing', 10);
    left = uigridlayout(middle, [4 1], 'RowHeight', {'fit', '1x', 'fit', '1.5x'}, ...
        'Padding', [0 0 0 0], 'RowSpacing', 4);
    formulaHead = uigridlayout(left, [1 2], 'ColumnWidth', {'1x', 130}, 'Padding', [0 0 0 0]);
    uilabel(formulaHead, 'WordWrap', 'on', 'Text', ...
        'What explains each bin''s response, in Unfold''s notation (double-click to edit):');
    uibutton(formulaHead, 'Text', 'Copy to every bin', 'ButtonPushedFcn', @(~, ~) onCopyFormula(), ...
        'Tooltip', 'Give every bin the formula of the row last clicked');
    formulaTable = uitable(left, 'Tag', 'formulas', 'ColumnName', {'Bin', 'Formula'}, ...
        'RowName', {}, 'ColumnEditable', [false true], 'ColumnWidth', {150, 'auto'}, ...
        'Data', cell(0, 2), ...
        'CellEditCallback', @(~, event) onFormulaEdited(event), ...
        'CellSelectionCallback', @(~, event) onFormulaSelected(event));
    uilabel(left, 'Text', 'The model:');
    modelList = uitextarea(left, 'Editable', 'off', 'Value', {''});

    right = uigridlayout(middle, [4 1], 'RowHeight', {'fit', '2x', 'fit', '1x'}, ...
        'Padding', [0 0 0 0], 'RowSpacing', 4);
    uilabel(right, 'WordWrap', 'on', 'Text', 'What a formula can use:');
    uitextarea(right, 'Editable', 'off', 'Tag', 'fields', 'Value', fieldReference(EEG));
    uilabel(right, 'WordWrap', 'on', 'Text', 'Events in no bin, modelled and dropped:');
    % A checkbox tree, not a multi-select list box. A list box deselects only
    % with Ctrl+click, and a plain click on the one selected item leaves it
    % selected, so going back to "none" was out of reach in practice. A tick
    % box toggles on a plain click, which is what an optional choice needs.
    otherTree = uitree(right, 'checkbox', 'CheckedNodesChangedFcn', @(~, ~) onOtherEventsChanged(), ...
        'Tooltip', ['Each ticked code gets its own full response, fitted and then dropped, so ' ...
         'its overlap is taken out of the bins. A code with only a handful of events adds a ' ...
         'whole window of parameters for very little.']);

    uilabel(outer, 'WordWrap', 'on', 'FontColor', [0.35 0.35 0.35], 'Text', [ ...
        'The window should cover the whole response, including anything that precedes the event ' ...
        'itself. Bad stretches are left out by ignoring them in the model rather than by dropping ' ...
        'epochs, so an event beside one keeps the rest of its data; the defaults (150 uV in a ' ...
        '2000 ms window, stepped 100 ms, on every channel) are the toolbox''s own. Events in ' ...
        'no bin are worth modelling: overlap is only removed where it is accounted for, so a ' ...
        'response or a following stimulus left out still overlaps, it just stops being ' ...
        'separated out. ' ...
        'The threshold is peak-to-peak within the moving window, and a window that exceeds it ' ...
        'is left out whole, so in reading or free viewing, where eye movements are the task, ' ...
        'it is usually better off (0) with the saccades and blinks modelled instead. A term in ' ...
        'a formula is fitted and held at the same values for every bin that uses it, so a ' ...
        'difference between bins is not one in that term.']);

    buttons = uigridlayout(outer, [1 3], 'ColumnWidth', {'1x', 90, 90}, 'Padding', [8 6 8 6]);
    cancelButton = uibutton(buttons, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) onCancel());
    cancelButton.Layout.Column = 2;
    okButton = uibutton(buttons, 'Text', 'OK', 'BackgroundColor', accentColor, ...
        'FontColor', [1 1 1], 'ButtonPushedFcn', @(~, ~) onOK());
    okButton.Layout.Column = 3;
    fig.CloseRequestFcn = @(~, ~) onCancel();

    showModel();
    uiwait(fig);

    function onBasisChanged()
    %ONBASISCHANGED  The number of splines or frequencies means nothing for
    %   one beta per sample.
        if ~exist('basisCountField', 'var')
            return;
        end
        basisCountField.Enable = matlab.lang.OnOffSwitchState(~strcmp(basisDropdown.Value, 'stick'));
    end

    function onFitMeasureChanged()
    %ONFITMEASURECHANGED  The fold event only means something for the
    %   cross-validated fits, which, like any fit, need overlap correction.
        if ~exist('foldDropdown', 'var')
            return;
        end
        foldDropdown.Enable = matlab.lang.OnOffSwitchState(startsWith(fitDropdown.Value, 'crossVal'));
    end

    function onOutputChanged()
    %ONOUTPUTCHANGED  The values to evaluate at, and how the other terms are
    %   added, only mean something for the terms, so they are only open for
    %   them.
        forTerms = matlab.lang.OnOffSwitchState(strcmp(outputDropdown.Value, 'terms'));
        evaluateField.Enable = forTerms;
        marginalDropdown.Enable = forTerms;
        if exist('spacingDropdown', 'var')
            spacingDropdown.Enable = forTerms;
        end
        % How the trials are made matters only for overlap-corrected trials:
        % without overlap correction they are the epochs as cut.
        if exist('trialsDropdown', 'var') && exist('overlapDropdown', 'var')
            trialsDropdown.Enable = matlab.lang.OnOffSwitchState( ...
                strcmp(outputDropdown.Value, 'trials') && strcmp(overlapDropdown.Value, 'on'));
        end
    end

    % ---- the formulas ------------------------------------------------- %
    function formulas = tableFormulas()
    %TABLEFORMULAS  The formulas as the table shows them, one per bin, read
    %   from the table itself so OK returns what is on screen however it got
    %   there.
        data = formulaTable.Data;
        formulas = struct('bin', {}, 'formula', {});
        for k = 1:size(data, 1)
            formulas(end + 1) = struct('bin', data{k, 1}, 'formula', normalFormula(data{k, 2})); %#ok<AGROW>
        end
    end

    function seedFormulaTable(tagged)
    %SEEDFORMULATABLE  One row per bin that gets fitted, with the formula
    %   typed for it before, if any. Rebuilt only when the bins themselves
    %   change, so an edit does not reset the table under the pointer.
        labels = ordinaryBinLabels(tagged);
        data = formulaTable.Data;
        if ~isempty(data) && isequal(reshape(data(:, 1), 1, []), labels)
            return;
        end
        remember(tableFormulas());
        legacy = legacyFormulas(tagged);
        data = [reshape(labels, [], 1), repmat({'y ~ 1'}, numel(labels), 1)];
        for k = 1:numel(labels)
            hit = find(strcmp({remembered.bin}, labels{k}), 1);
            if ~isempty(hit)
                data{k, 2} = remembered(hit).formula;
            elseif isKey(legacy, labels{k})
                data{k, 2} = legacy(labels{k});
            end
        end
        formulaTable.Data = data;
        onOff = {'off', 'on'};       % a table takes the words, not an OnOffSwitchState
        formulaTable.Enable = onOff{1 + ~isempty(labels)};
        selectedRow = 1;
    end

    function legacy = legacyFormulas(tagged)
    %LEGACYFORMULAS  The formulas stored covariates always meant, by bin:
    %   what Unfold.binModel makes of them, which adds each one only to the
    %   bins whose events carry it and vary in it. Done once, since from
    %   then on the table holds them.
        legacy = containers.Map('KeyType', 'char', 'ValueType', 'char');
        if isempty(legacyCovariates)
            return;
        end
        try
            plan = Unfold.binModel(tagged, 'OtherEvents', otherSelection, ...
                'Covariates', legacyCovariates);
            for k = 1:numel(plan.typeLabels)
                legacy(plan.typeLabels{k}) = plan.formulas{k};
            end
        catch
            % The bins cannot be fitted as they stand; the model list says why.
        end
        legacyCovariates = {};
    end

    function remember(formulas)
    %REMEMBER  Keep FORMULAS by bin label, replacing any older ones.
        for k = 1:numel(formulas)
            hit = strcmp({remembered.bin}, formulas(k).bin);
            remembered(hit) = [];
            remembered(end + 1) = formulas(k); %#ok<AGROW>
        end
    end

    function onFormulaEdited(event)
    %ONFORMULAEDITED  A formula tidied the way the fit will read it (an
    %   empty one is 'y ~ 1', a missing 'y ~' is added), then the preview.
        row = event.Indices(1);
        formulaTable.Data{row, 2} = normalFormula(event.NewData);
        remember(tableFormulas());
        showModel();
    end

    function onFormulaSelected(event)
        if ~isempty(event.Indices)
            selectedRow = event.Indices(1, 1);
        end
    end

    function onCopyFormula()
    %ONCOPYFORMULA  The formula of the row last clicked, in every row.
        data = formulaTable.Data;
        if isempty(data)
            return;
        end
        row = min(max(selectedRow, 1), size(data, 1));
        data(:, 2) = data(row, 2);
        formulaTable.Data = data;
        remember(tableFormulas());
        showModel();
    end

    % ---- the events in no bin ---------------------------------------- %
    function onOtherEventsChanged()
    %ONOTHEREVENTSCHANGED  The ticks become the choice, then the preview.
        otherSelection = otherTreeSelection();
        showModel();
    end

    function selection = otherTreeSelection()
    %OTHERTREESELECTION  The choice as the ticks show it. Every code ticked is
    %   kept as 'all' rather than as today's list, so a replay on a recording
    %   that has a code this one lacks still models it; anything less is kept
    %   as the codes themselves. Read from the tree, not from a copy, so OK
    %   returns what is on screen whatever route the ticks took to get there.
        selection = otherSelection;
        if strcmpi(otherTree.Enable, 'off')
            return;          % nothing to choose from: the choice stands
        end
        checked = otherTree.CheckedNodes;
        if isempty(checked)
            selection = {};
        elseif numel(checked) == numel(otherTree.Children)
            selection = 'all';
        else
            selection = reshape({checked.NodeData}, 1, []);
        end
    end

    function refreshOtherTree(unbinned)
    %REFRESHOTHERTREE  One tick box per code in no bin, with its count, ticked
    %   when it is modelled. The count is the point: it is what shows that a
    %   code with one event is not worth a whole window of parameters. Rebuilt
    %   only when the codes themselves change, so ticking one does not reset
    %   the list under the pointer.
        if isempty(unbinned)
            delete(otherTree.Children);
            uitreenode(otherTree, 'Text', '(none: every event is in a bin)');
            otherTree.Enable = 'off';
            return;
        end
        labels = arrayfun(@(u) sprintf('%s  (%d event(s))', u.code, u.n), ...
            unbinned, 'UniformOutput', false);
        current = otherTree.Children;
        if strcmpi(otherTree.Enable, 'off') || numel(current) ~= numel(labels) ...
                || ~isequal({current.Text}, reshape(labels, 1, []))
            delete(current);
            for k = 1:numel(unbinned)
                uitreenode(otherTree, 'Text', labels{k}, 'NodeData', unbinned(k).code);
            end
        end
        otherTree.Enable = 'on';
        nodes = otherTree.Children;
        setChecked(otherTree, nodes([unbinned.modelled]));
    end

    % ---- the rest ------------------------------------------------------ %
    function onFitChanged(refresh)
    %ONFITCHANGED  The iteration limit is lsmr's, which only the toolbox's
    %   default solver with overlap correction runs; alpha is glmnet's. With
    %   REFRESH the model is shown again (overlap correction decides which
    %   warnings it has).
        lsmr = strcmp(solverDropdown.Value, 'toolbox') && strcmp(overlapDropdown.Value, 'on');
        iterationsField.Enable = matlab.lang.OnOffSwitchState(lsmr);
        alphaField.Enable = matlab.lang.OnOffSwitchState(strcmp(solverDropdown.Value, 'glmnet'));
        onOutputChanged();
        if refresh
            showModel();
        end
    end

    function onBaselineToggled()
    %ONBASELINETOGGLED  The two fields follow the checkbox, so an unticked box
    %   cannot leave a window behind that looks as if it were being used.
        state = matlab.lang.OnOffSwitchState(baselineBox.Value);
        baseStartField.Enable = state;
        baseStopField.Enable = state;
    end

    function onDefineBins()
    %ONDEFINEBINS  DefineBins' editor, without its epoching half.
        note = ['These bins are fitted against the continuous recording, so there is no epoch ' ...
                'window to set here: the response window in the Deconvolve dialog behind this ' ...
                'one plays that part. Everything else in the language works as it does in ' ...
                'DefineBins, combination bins ("bin 3 = bin 1 - bin 2") included. A combination ' ...
                'bin is not fitted, since it has no events of its own; this step works it out ' ...
                'afterwards by subtracting the bins it refers to, so it arrives with the rest.'];
        result = DefineBinsDialog(binScript, {'', ''}, ...
            struct('showEpoch', false, 'name', 'Deconvolve: define bins', 'note', note));
        if isempty(result); return; end
        binScript = strtrim(result.script);
        showModel();
    end

    function [tagged, why] = binnedDataset()
    %BINNEDDATASET  The recording with its bins on it, or why there are none.
    %   This is the same call Deconvolve makes, so what the preview lists is
    %   what would be fitted rather than an approximation of it.
        why = '';
        if isempty(binScript)
            tagged = EEG;
            if ~hasTags
                tagged = [];
                why = ['No bins yet. Would you press "Define bins..." and say which events ' ...
                       'each bin holds? This recording carries no bin tags of its own, so ' ...
                       'there is nothing to fit until it does.'];
            end
            return;
        end
        try
            tagged = DefineBins(EEG, struct('script', binScript, 'tagsOnly', true));
        catch err
            tagged = [];
            why = err.message;
        end
    end

    function showModel()
    %SHOWMODEL  The model as it stands, or why there is not one.
        [tagged, why] = binnedDataset();
        binsLabel.Text = binSourceText();
        if isempty(tagged)
            formulaTable.Data = cell(0, 2);
            formulaTable.Enable = 'off';
            modelList.Value = [{'No model yet:'}, {''}, {why}];
            refreshOtherTree([]);
            return;
        end
        seedFormulaTable(tagged);
        try
            plan = Unfold.binModel(tagged, 'OtherEvents', otherSelection, ...
                'Formulas', tableFormulas(), 'MissingValues', missingDropdown.Value, ...
                'OverlapCorrection', strcmp(overlapDropdown.Value, 'on'), ...
                'CodingSchema', codingDropdown.Value);
            % A code chosen earlier can stop being "in no bin" when the bins
            % change, and would then only produce a note saying it is absent.
            % The list is the user's view of the choice, so the choice follows
            % it: forget codes that are no longer on offer and ask again.
            if iscell(otherSelection)
                available = {plan.unbinnedCodes.code};
                if ~all(ismember(otherSelection, available))
                    otherSelection = intersect(otherSelection, available, 'stable');
                    plan = Unfold.binModel(tagged, 'OtherEvents', otherSelection, ...
                        'Formulas', tableFormulas(), 'MissingValues', missingDropdown.Value, ...
                        'OverlapCorrection', strcmp(overlapDropdown.Value, 'on'), ...
                        'CodingSchema', codingDropdown.Value);
                end
            end
            refreshOtherTree(plan.unbinnedCodes);
            design = designLines(tagged, plan);
        catch err
            modelList.Value = [{'This design cannot be fitted:'}, {''}, {err.message}];
            return;
        end
        lines = {sprintf('%d bin(s) will be fitted:', numel(plan.binLabels))};
        for k = 1:numel(plan.binLabels)
            lines{end + 1} = sprintf('   %-30s %5d event(s)', plan.binLabels{k}, plan.binCounts(k)); %#ok<AGROW>
        end
        if ~isempty(plan.comboBins)
            lines{end + 1} = '';
            lines{end + 1} = sprintf(['%d difference bin(s) come out of this step too, worked ' ...
                'out here once the fit is done by subtracting the fitted waveforms. Nothing ' ...
                'else has to be run for them.'], numel(plan.comboBins));
        end
        lines = [lines, design];
        if isempty(plan.nuisanceTypes)
            lines{end + 1} = '';
            lines{end + 1} = ['No nuisance events are modelled, so any overlap from events ' ...
                'outside these bins stays in the result.'];
        else
            lines{end + 1} = '';
            lines{end + 1} = sprintf('Nuisance event type(s), fitted and then dropped: %s', ...
                strjoin(strrep(plan.nuisanceTypes, 'evt_', ''), ', '));
        end
        for k = 1:numel(plan.notes)
            lines{end + 1} = ''; %#ok<AGROW>
            lines{end + 1} = ['Note: ' plan.notes{k}]; %#ok<AGROW>
        end
        locked = struct('note', {});
        if strcmp(overlapDropdown.Value, 'on')   % only a deconvolution has pairs to tell apart
            locked = Unfold.timeLockedEvents(plan, tagged.srate, [startField.Value stopField.Value]);
        end
        for k = 1:numel(locked)
            lines{end + 1} = ''; %#ok<AGROW>
            lines{end + 1} = ['Warning: ' locked(k).note]; %#ok<AGROW>
        end
        modelList.Value = lines;
        warnAboutLockedEvents({locked.note});
    end

    function warnAboutLockedEvents(notes)
    %WARNABOUTLOCKEDEVENTS  An alert when the model gains a time-locked pair.
    %   Only when the pairs change, so editing a formula does not raise the
    %   same alert again; unticking the code clears them, and ticking it
    %   again warns again. An alert rather than a question: it leaves the
    %   dialog usable, and the model list keeps the warning in view.
        if ~isempty(notes) && ~isequal(notes, lockedWarned)
            uialert(fig, strjoin(notes, [newline newline]), 'Events locked together', ...
                'Icon', 'warning');
        end
        lockedWarned = notes;
    end

    function text = binSourceText()
    %BINSOURCETEXT  Which of the two possible sources of bins is in use, said
    %   plainly, because the difference decides what the fit means.
        if ~isempty(binScript)
            text = sprintf('Bins: defined here, %d line(s). "Define bins..." to edit them.', ...
                numel(strsplit(binScript, newline)));
        elseif hasTags
            text = sprintf(['Bins: the %d already tagged on this dataset. "Define bins..." ' ...
                'to fit different ones instead.'], numel(EEG.bindesc));
        else
            text = 'Bins: none yet. "Define bins..." to say which events each bin holds.';
        end
    end

    function onOK()
        otherSelection = otherTreeSelection();
        candidate = struct('binScript', binScript, ...
            'formulas', {tableFormulas()}, ...
            'evaluateAt', strtrim(evaluateField.Value), ...
            'windowMs', [startField.Value stopField.Value], ...
            'baselineMs', [], ...
            'otherEvents', {otherSelection}, ...
            'artifactThresholdUv', thresholdField.Value, ...
            'artifactWindowMs', artWindowField.Value, ...
            'artifactStepMs', artStepField.Value, ...
            'artifactChannels', channelsDropdown.Value, ...
            'missingValues', missingDropdown.Value, ...
            'overlapCorrection', strcmp(overlapDropdown.Value, 'on'), ...
            'excludeAroundCuts', strcmp(cutsDropdown.Value, 'on'), ...
            'codingSchema', codingDropdown.Value, ...
            'modelFit', fitDropdown.Value, ...
            'collinearity', vifDropdown.Value, ...
            'timeBasis', basisDropdown.Value, ...
            'timeBasisCount', basisCountField.Value, ...
            'foldEvent', {nonEmpty(foldDropdown.Value, NO_FOLD)}, ...
            'evaluateSpacing', spacingDropdown.Value, ...
            'solver', strrep(solverDropdown.Value, 'toolbox', 'default'), ...
            'glmnetAlpha', alphaField.Value, ...
            'marginal', marginalDropdown.Value, ...
            'solverIterations', iterationsField.Value, ...
            'output', outputDropdown.Value, ...
            'trialsMethod', trialsDropdown.Value);
        if candidate.windowMs(1) >= candidate.windowMs(2)
            uialert(fig, 'Would the window start come before its stop?', 'Check the window');
            return;
        end
        if candidate.overlapCorrection && strcmp(candidate.solver, 'glmnet') ...
                && strcmp(candidate.output, 'trials')
            uialert(fig, Unfold.glmnetTrialsMessage(), 'Trials need another solver');
            return;
        end
        if ~strcmp(candidate.timeBasis, 'stick') && strcmp(candidate.output, 'trials')
            uialert(fig, Unfold.basisTrialsMessage(), 'Trials need one beta per sample');
            return;
        end
        if ~strcmp(candidate.timeBasis, 'stick') && ~candidate.overlapCorrection
            uialert(fig, ['A time basis shapes the deconvolution''s time expansion, which a ' ...
                'regression on epochs does not have. Would you switch Overlap correction on, or ' ...
                'keep one beta per sample?'], 'Time basis');
            return;
        end
        if strcmp(candidate.collinearity, 'expanded') && ~candidate.overlapCorrection
            uialert(fig, ['The collinearity of the time-expanded design needs that design, which ' ...
                'only the deconvolution builds. Would you switch Overlap correction on, or leave ' ...
                'the check out?'], 'Collinearity');
            return;
        end
        if ~strcmp(candidate.modelFit, 'none') && ~candidate.overlapCorrection
            uialert(fig, ['A model fit is measured on the deconvolution''s design, by the ' ...
                'toolbox''s uf_checkmodelfit, so it needs Overlap correction on.'], 'Model fit');
            return;
        end
        if startsWith(candidate.modelFit, 'crossVal') && isempty(candidate.foldEvent)
            uialert(fig, ['Cross-validation cuts the recording into folds at an event code. ' ...
                'Would you choose one under "Folds cut at", such as a break between blocks?'], ...
                'Model fit');
            return;
        end
        if baselineBox.Value
            candidate.baselineMs = [baseStartField.Value baseStopField.Value];
            if candidate.baselineMs(1) >= candidate.baselineMs(2)
                uialert(fig, 'Would the baseline start come before its stop?', ...
                    'Check the baseline');
                return;
            end
            if candidate.baselineMs(1) < candidate.windowMs(1) || ...
                    candidate.baselineMs(2) > candidate.windowMs(2)
                uialert(fig, sprintf(['The baseline window has to sit inside the response ' ...
                    'window (%g to %g ms), since that is the only part that gets fitted.'], ...
                    candidate.windowMs(1), candidate.windowMs(2)), 'Check the baseline');
                return;
            end
        end
        [tagged, why] = binnedDataset();
        if isempty(tagged)
            uialert(fig, why, 'Check the bins');
            return;
        end
        try
            plan = Unfold.binModel(tagged, 'OtherEvents', otherSelection, ...
                'Formulas', candidate.formulas, 'MissingValues', candidate.missingValues, ...
                'OverlapCorrection', candidate.overlapCorrection, ...
                'CodingSchema', candidate.codingSchema);
            designLines(tagged, plan);
        catch err
            uialert(fig, err.message, 'Check the design');
            return;
        end
        % The values to evaluate the terms at, read the way the fit reads
        % them and against the same design, so what is refused here is what
        % the fit would refuse. Without the toolbox the fit checks them.
        if strcmp(candidate.output, 'terms') && ~isempty(candidate.evaluateAt) && Unfold.isAvailable()
            try
                designed = quietDesign(tagged, plan);
                Unfold.predictionValues(candidate.evaluateAt, designed.unfold, plan.typeLabels);
            catch err
                uialert(fig, err.message, 'Check the values to evaluate at');
                return;
            end
        end
        options = candidate;
        uiresume(fig);
        delete(fig);
    end

    function onCancel()
        uiresume(fig);
        delete(fig);
    end
end

% ======================================================================= %
function lines = designLines(tagged, plan)
%DESIGNLINES  The columns the toolbox builds for each bin's formula, from
%   Unfold.designMatrix on the events alone: the same call the fit makes, so
%   a formula it refuses is refused here, naming the bin, and throws. Only
%   bins whose formula is more than 'y ~ 1' are listed, since an intercept
%   is what every bin has anyway. Without the toolbox the formulas have been
%   checked against the events only, which is said.
    lines = {};
    trivial = cellfun(@(f) strcmp(strrep(f, ' ', ''), 'y~1'), plan.formulas);
    if all(trivial)
        return;
    end
    if ~Unfold.isAvailable()
        lines = {'', ['The Unfold toolbox is not installed here, so the formulas are checked ' ...
            'against the events only; the toolbox reads their notation when it fits.']};
        return;
    end
    designed = quietDesign(tagged, plan);
    unfold = designed.unfold;
    lines = {'', 'Each formula, as the toolbox builds it:'};
    nbins = numel(plan.eventTypes) - numel(plan.nuisanceTypes);
    for t = find(~trivial(1:nbins))
        cols = find(unfold.cols2eventtypes == t);
        variableOf = unfold.cols2variablenames(cols);
        parts = {};
        for v = reshape(unique(variableOf, 'stable'), 1, [])
            name = regexprep(unfold.variablenames{v}, '^\d+_', '');
            at = cols(variableOf == v);
            switch unfold.variabletypes{v}
                case 'intercept'
                    parts{end + 1} = 'intercept'; %#ok<AGROW>
                case 'spline'
                    surface = plan.pooled(strcmp({plan.pooled.name}, name) ...
                        & ~cellfun(@isempty, {plan.pooled.pair}));
                    if isempty(surface)
                        parts{end + 1} = sprintf('%s as a spline (%d columns)', name, numel(at)); %#ok<AGROW>
                    else
                        parts{end + 1} = sprintf('%s and %s as a 2D spline (%d columns)', ...
                            surface(1).pair{1}, surface(1).pair{2}, numel(at)); %#ok<AGROW>
                    end
                case 'categorical'
                    parts{end + 1} = sprintf('%s as a factor (%s)', name, ...
                        strjoin(regexprep(unfold.colnames(at), '^\d+_', ''), ', ')); %#ok<AGROW>
                otherwise
                    parts{end + 1} = sprintf('%s (%s)', name, unfold.variabletypes{v}); %#ok<AGROW>
            end
        end
        lines{end + 1} = sprintf('   %s: %s', plan.typeLabels{t}, strjoin(parts, '; ')); %#ok<AGROW>
        lines{end + 1} = sprintf('      %s', plan.formulas{t}); %#ok<AGROW>
    end
end

function designed = quietDesign(tagged, plan) %#ok<INUSD>  plan is read inside evalc
%QUIETDESIGN  Unfold.designMatrix on the events alone, without the progress
%   lines uf_designmat prints on every call: the preview reruns it on every
%   edit, and the command window is for the fit's own narration.
    stub = struct('event', tagged.event, 'srate', tagged.srate, 'pnts', tagged.pnts); %#ok<NASGU>
    designed = [];
    evalc('designed = Unfold.designMatrix(stub, plan);');
end

function formula = normalFormula(formula)
%NORMALFORMULA  A formula as the fit reads it: an empty one is 'y ~ 1', and
%   one without its left-hand side gets 'y ~ ' (Unfold.binModel does the
%   same, so this only makes the table show what will be fitted).
    formula = strtrim(char(string(formula)));
    if isempty(formula)
        formula = 'y ~ 1';
    elseif ~contains(formula, '~')
        formula = ['y ~ ' formula];
    end
end

function formulas = storedFormulas(stored)
%STOREDFORMULAS  Stored formulas as a struct array of .bin and .formula,
%   whatever shape a template gave them.
    formulas = struct('bin', {}, 'formula', {});
    if ~isstruct(stored) || ~all(isfield(stored, {'bin', 'formula'}))
        return;
    end
    for k = 1:numel(stored)
        formulas(end + 1) = struct('bin', char(string(stored(k).bin)), ...
            'formula', normalFormula(stored(k).formula)); %#ok<AGROW>
    end
end

function labels = ordinaryBinLabels(EEG)
%ORDINARYBINLABELS  The bins that are fitted, in order: every bin but the
%   combination ones, which are worked out from the others after the fit.
    labels = {};
    if ~isfield(EEG, 'bindesc') || isempty(EEG.bindesc)
        return;
    end
    combo = false(1, numel(EEG.bindesc));
    if isfield(EEG.bindesc, 'combo')
        combo = arrayfun(@(b) ~isempty(b.combo), EEG.bindesc);
    end
    labels = reshape(cellstr(string({EEG.bindesc(~combo).label})), 1, []);
end

function lines = fieldReference(EEG)
%FIELDREFERENCE  The event fields a formula can name, by what they can be,
%   and the notation, so a formula can be written without leaving the
%   dialog. Numbers come from Unfold.eventCovariates, which knows EYE-EEG's
%   units and which of them are angles; text fields with a handful of
%   values are the factors.
    lines = {};
    found = Unfold.eventCovariates(EEG);
    linear = found(strcmpi({found.kind}, 'linear'));
    circular = found(~strcmpi({found.kind}, 'linear'));
    if ~isempty(linear)
        lines = [lines, {'Numbers (a term, or spl):'}, ...
            arrayfun(@fieldItem, linear, 'UniformOutput', false)];
    end
    if ~isempty(circular)
        lines = [lines, {'', 'Angles (circspl):'}, ...
            arrayfun(@fieldItem, circular, 'UniformOutput', false)];
    end
    factors = textFields(EEG);
    if ~isempty(factors)
        lines = [lines, {'', 'Text (cat):'}, ...
            arrayfun(@(f) sprintf('   %s, %d levels', f.name, f.levels), factors, ...
            'UniformOutput', false)];
    end
    if isempty(lines)
        lines = {'(no event field a formula could use)'};
    end
    lines = [lines, {'', 'The notation:', ...
        '   y ~ 1   the bin''s own waveform', ...
        '   y ~ 1 + rt   plus a straight line in rt', ...
        '   y ~ 1 + spl(rt, 5)   a smooth curve, 5 splines', ...
        '   y ~ 1 + cat(side)   a factor, its first level the reference', ...
        '   y ~ 1 + cat(side) * rt   a factor, rt, and their interaction', ...
        '   circspl(angle, 5, 0, 360)   a curve round a circle', ...
        '   2dspl(x, y, 5)   a smooth surface over two fields', ...
        '', 'An event without a number its formula names is filled in or left out', ...
        'as Missing values says; a factor needs a level on every event.', ...
        'Give a factor''s levels names longer than one letter: the toolbox cannot build', ...
        '"L" and "R" when another event type lacks the field.'}];
end

function text = fieldItem(candidate)
    text = ['   ' candidate.name];
    if ~isempty(candidate.unit)
        text = sprintf('%s (%s)', text, candidate.unit);
    end
    text = sprintf('%s, n = %d', text, candidate.n);
end

function fields = textFields(EEG)
%TEXTFIELDS  Event fields holding text with between 2 and 20 distinct
%   values: the ones cat() can make a factor of. More than that is a label
%   per event (a word, a file name), not a condition. The fields that place
%   and tag the events are left out, as Unfold.eventCovariates leaves them.
    fields = struct('name', {}, 'levels', {});
    if ~isfield(EEG, 'event') || isempty(EEG.event)
        return;
    end
    skip = {'type', 'latency', 'duration', 'urevent', 'bini', 'epoch'};
    for name = reshape(fieldnames(EEG.event), 1, [])
        if any(strcmpi(name{1}, skip))
            continue;
        end
        values = {EEG.event.(name{1})};
        values = values(~cellfun(@isempty, values));
        if isempty(values) || ~all(cellfun(@(v) ischar(v) || (isstring(v) && isscalar(v)), values))
            continue;
        end
        levels = numel(unique(cellfun(@char, values, 'UniformOutput', false)));
        if levels >= 2 && levels <= 20
            fields(end + 1) = struct('name', name{1}, 'levels', levels); %#ok<AGROW>
        end
    end
end

% ======================================================================= %
function setChecked(tree, nodes)
%SETCHECKED  Tick exactly NODES in a checkbox tree. "Nothing ticked" has to
%   be the plain [] rather than an empty array of tree nodes: CheckBoxTree
%   rejects the latter as not being children of the tree, and indexing the
%   children with a mask that matches nothing produces exactly that, so the
%   dialog failed to open whenever nothing was stored.
    if isempty(nodes)
        tree.CheckedNodes = [];
    else
        tree.CheckedNodes = nodes;
    end
end

% ======================================================================= %
function tf = datasetCarriesBins(EEG)
%DATASETCARRIESBINS  Bins already on the recording, tags and all.
    tf = isfield(EEG, 'bindesc') && ~isempty(EEG.bindesc) ...
        && isfield(EEG, 'event') && ~isempty(EEG.event) && isfield(EEG.event, 'bini');
end

function script = lastDefineBinsScript()
%LASTDEFINEBINSSCRIPT  The script last run through DefineBins in this
%   workspace, as a starting point. It is only a seed: it is shown in the
%   editor and can be changed or cleared there before anything is fitted.
    script = '';
    stored = TransformSettings.get('DefineBins');
    if isstruct(stored) && isfield(stored, 'script')
        script = strtrim(char(string(stored.script)));
    end
end

function seed = defaults()
%DEFAULTS  First-run settings: the paper's window and the toolbox's own
%   artefact parameters, channels and marginal effects, with nuisance events
%   modelled because leaving them out quietly weakens the correction the
%   transformation exists for. covariates is only read, from options stored
%   before formulas existed.
    seed = struct('binScript', '', 'covariates', {{}}, ...
        'formulas', {struct('bin', {}, 'formula', {})}, 'evaluateAt', '', ...
        'windowMs', [-200 800], 'artifactThresholdUv', 150, ...
        'artifactWindowMs', 2000, 'artifactStepMs', 100, 'artifactChannels', 'all', ...
        'marginal', 'MEM', 'missingValues', 'median', 'solverIterations', 400, ...
        'solver', 'default', 'glmnetAlpha', 1, 'overlapCorrection', true, ...
        'excludeAroundCuts', false, 'codingSchema', 'reference', ...
        'evaluateSpacing', 'quantile', 'modelFit', 'none', 'foldEvent', {{}}, ...
        'timeBasis', 'stick', 'timeBasisCount', 30, 'collinearity', 'none', ...
        'output', 'average', 'trialsMethod', 'fast');
end

function value = onOff(flag)
%ONOFF  A stored overlapCorrection, true or false, as the dropdown's 'on' or
%   'off'; anything unreadable is the default, on.
    value = 'on';
    if (islogical(flag) || isnumeric(flag)) && isscalar(flag) && ~flag
        value = 'off';
    end
end

function value = choiceSeed(stored, offered)
%CHOICESEED  A stored choice the dropdown offers, compared without regard to
%   case; anything else (a hand-edited template, say) is the first offered,
%   the toolbox's default.
    value = offered{1};
    hit = strcmpi(offered, char(string(stored)));
    if any(hit)
        value = offered{find(hit, 1)};
    end
end

function value = outputSeed(stored)
%OUTPUTSEED  A stored output choice the dropdown can show: anything it does
%   not offer (a hand-edited template, say) falls back to the waveforms,
%   which is what Deconvolve itself does with it.
    value = lower(char(string(stored)));
    if ~any(strcmp(value, {'average', 'trials', 'terms'}))
        value = 'average';
    end
end

function value = nonEmpty(code, none)
%NONEMPTY  A chosen fold event as a list of codes, none when NONE was chosen.
    value = {};
    if ~isempty(code) && ~strcmp(code, none)
        value = {char(code)};
    end
end
