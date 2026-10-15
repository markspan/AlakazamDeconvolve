# Deconvolve for Alakazam

A plugin for [Alakazam](https://github.com/markspan/Alakazam): the alternative to epoching and averaging where responses overlap: every bin fitted at once against the continuous recording by the Unfold toolbox, giving overlap-corrected waveforms, single trials or the model's terms.

It was part of Alakazam until it became a plugin of its own, so that it can be
installed, upgraded and documented on its own.

## Manual

The manual is on this repository's pages:
<https://markspan.github.io/AlakazamDeconvolve/>.

## Installing

In Alakazam, **Install plugins** on the Alakazam tab lists this repository:
select it and press **Install**. It installs the latest release. A zip from
the [releases](https://github.com/markspan/AlakazamDeconvolve/releases), or one built with
`build`, installs the same way through **From a link**, and so does the link
to this repository. The tool appears on the Tools tab, as it did before.

It needs an Alakazam that no longer carries `Deconvolve` itself (newer than
V0.4.5.2). Workspaces and templates made with the built-in tool open and
recalculate as before once the plugin is installed; Alakazam offers to
install it when one needs it.

## Building and testing

```matlab
build()                                    % dist/Deconvolve-<version>.zip
runtests('tests')                          % the tool's tests and the install
```

The tests need Alakazam's source: the folder in `ALAKAZAM_ROOT`, or the
Alakazam repository beside this one.

The fit runs the Unfold toolbox (1.3.1, pinned), which Alakazam downloads on first use, with consent. The tests that fit need it, and skip without it; they never download it. How Deconvolve calls the toolbox, step by step, is in [docs/unfold-technical-note.md](docs/unfold-technical-note.md).

Two templates come with the plugin, in `Deconvolve/templates`, and are installed with it: `FaceSaccadesDeconvolution` (Figure 11 of Ehinger and Dimigen, 2019, on the authors' recording) and `ReadingDeconvolution` (fixation-related potentials in reading, on EYE-EEG's test data; it needs the [Eye tracking plugin](https://github.com/markspan/AlakazamEyetracking) too). TemplateReplayTest replays both on their data, under the data folder (ALAKAZAM_DATA, Alakazam's Data folder, or D:\data), and skips without it.

## Licence

GPL-3.0, as Alakazam.
