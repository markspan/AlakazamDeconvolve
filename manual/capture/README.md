# How the pictures were made

The pictures on the manual page were taken with Alakazam's own capture tool
(`src/help/capture` in Alakazam) while Deconvolve was still built in, on the
authors' face recording (`face_saccades_opendata_fig10.set`, from the
paper's OSF page). The deconvolved nodes come from
`Deconvolve/templates/FaceSaccadesDeconvolution.alztemplate`, with the
**Result** changed for the waveforms and the model's terms as the page
describes; the comparison without deconvolution comes from
`FaceSaccadesPlain.alztemplate` here: the same three bins as plain epochs
from -1500 to 1000 ms, then Baseline, ArtefactDetect and Average.
