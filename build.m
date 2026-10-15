function zipFile = build(outFolder)
%BUILD  The plugin as Alakazam's installer takes it: a zip holding the
%   plugin's folder, named after the plugin and its version.
%
%   ZIPFILE = build() writes dist/Deconvolve-<version>.zip beside this file;
%   build(FOLDER) writes it into FOLDER. Install it in Alakazam with Install
%   plugins (it lists this repository), or attach it to a GitHub release.
%
%   The tests, the manual and this file stay out of the zip: only the folder
%   Alakazam runs is in it.
    root = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(outFolder)
        outFolder = fullfile(root, 'dist');
    end
    if ~isfolder(outFolder)
        mkdir(outFolder);
    end
    tools = {'Deconvolve'};
    manifest = jsondecode(fileread(fullfile(root, tools{1}, [tools{1} '.json'])));
    zipFile = fullfile(outFolder, sprintf('%s-%s.zip', tools{1}, manifest.Version));
    if isfile(zipFile)
        delete(zipFile);
    end
    zip(zipFile, tools, root);
end
