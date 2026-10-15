function text = basisTrialsMessage()
%BASISTRIALSMESSAGE  Why a fit with a time basis other than one beta per
%   sample gives no overlap-corrected trials, for Deconvolve's dialog and for
%   Unfold.fitBins, which refuse the combination.
%
%   The trials are the toolbox's own uf_erpimage (Unfold.fitBins), and
%   uf_erpimage in Unfold 1.3.1 multiplies the time-expanded design, whose
%   columns are then splines or frequencies, by the betas uf_condense has
%   turned back into samples, which do not match ("Incorrect dimensions for
%   matrix multiplication"). A toolbox behaviour kept as it is, not worked
%   around.
%
%   See also UNFOLD.FITBINS, UNFOLD.GLMNETTRIALSMESSAGE, DECONVOLVEDIALOG.
    text = ['The overlap-corrected trials come from the Unfold toolbox''s own ' ...
        'uf_erpimage, which cannot make them with a time basis of splines or ' ...
        'frequencies (Unfold 1.3.1): it multiplies the design, whose columns are then ' ...
        'those functions, by betas already turned back into samples, and the two do not ' ...
        'match. Would you keep one beta per sample for the trials, or ask the spline or ' ...
        'Fourier fit for one waveform per bin or for the terms?'];
end
