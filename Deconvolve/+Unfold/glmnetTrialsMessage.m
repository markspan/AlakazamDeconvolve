function text = glmnetTrialsMessage()
%GLMNETTRIALSMESSAGE  Why a glmnet fit gives no overlap-corrected trials, for
%   Deconvolve's dialog and for Unfold.fitBins, which refuse the combination.
%
%   The trials are the toolbox's own uf_erpimage (Unfold.fitBins), and
%   uf_erpimage in Unfold 1.3.1 cannot make them for a glmnet fit: uf_glmfit
%   adds a column of its own to the design for glmnet, an intercept for the
%   whole recording ('glmnet-DC-Correction'), and uf_erpimage asks for a beta
%   of every column, which uf_condense, leaving that column out, does not
%   return. A toolbox behaviour kept as it is, not worked around.
%
%   See also UNFOLD.FITBINS, DECONVOLVEDIALOG.
    text = ['The overlap-corrected trials come from the Unfold toolbox''s own ' ...
        'uf_erpimage, which cannot make them for a glmnet fit (Unfold 1.3.1): glmnet ' ...
        'adds a column of its own to the design, an intercept for the whole recording, ' ...
        'and uf_erpimage asks for a beta of that column, which uf_condense does not ' ...
        'return. Would you choose the toolbox''s own solver or MATLAB''s exact one for ' ...
        'the trials, or ask glmnet for one waveform per bin or for the terms?'];
end
