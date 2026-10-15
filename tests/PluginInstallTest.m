classdef PluginInstallTest < matlab.unittest.TestCase
%PLUGININSTALLTEST  The built zip installs through Alakazam's own installer,
%   as Install plugins would put it in, and comes out again.
%
%   Run without the plugin's source on the path: the installer refuses a
%   plugin whose name is already a function on the path. Plugins go into a
%   temporary folder (ALAKAZAM_PLUGIN_FOLDER), not the real one.
%
%   Run with: runtests('tests/PluginInstallTest.m').
%
%   See also PLUGINS (Alakazam), BUILD.

    properties (Constant)
        Root = fileparts(fileparts(mfilename('fullpath')))
        Tools = {'Deconvolve'}
    end

    properties
        Alakazam = ''
        Zip = ''
    end

    methods (TestClassSetup)
        function prepare(testCase)
            testCase.Alakazam = getenv('ALAKAZAM_ROOT');
            if isempty(testCase.Alakazam)
                testCase.Alakazam = fullfile(fileparts(PluginInstallTest.Root), 'Alakazam');
            end
            testCase.assumeTrue(isfolder(testCase.Alakazam), ...
                'Alakazam''s source is not found (set ALAKAZAM_ROOT).');
            testCase.assumeEmpty(which('Deconvolve'), ...
                'The plugin is on the path already, which the installer would refuse.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.Alakazam, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(PluginInstallTest.Root));
            out = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            testCase.Zip = build(out);
        end
    end

    methods (TestMethodSetup)
        function ownPluginFolder(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            previous = getenv('ALAKAZAM_PLUGIN_FOLDER');
            testCase.addTeardown(@() setenv('ALAKAZAM_PLUGIN_FOLDER', previous));
            setenv('ALAKAZAM_PLUGIN_FOLDER', folder);
            before = path();
            testCase.addTeardown(@() path(before));
        end
    end

    methods (Test)
        function theInstallerFindsItAndHasNoObjection(testCase)
            candidate = Plugins.prepare(testCase.Zip, ...
                'BuiltInRoot', fullfile(testCase.Alakazam, 'src', 'Transformations'));
            testCase.addTeardown(@() Plugins.discard(candidate));
            testCase.verifyEqual(sort({candidate.Plugins.Name}), sort(PluginInstallTest.Tools));
            for plugin = candidate.Plugins
                testCase.verifyEmpty(plugin.Problems, strjoin(plugin.Problems, ' '));
                testCase.verifyEmpty(plugin.Warnings, strjoin(plugin.Warnings, ' '));
            end
        end

        function itInstallsAndUninstalls(testCase)
            candidate = Plugins.prepare(testCase.Zip, ...
                'BuiltInRoot', fullfile(testCase.Alakazam, 'src', 'Transformations'));
            names = Plugins.install(candidate);
            testCase.verifyEqual(sort(names), sort(PluginInstallTest.Tools));
            for t = PluginInstallTest.Tools
                testCase.verifyEqual(which(t{1}), fullfile(Plugins.folder(), t{1}, [t{1} '.m']));
            end
            for t = PluginInstallTest.Tools
                Plugins.uninstall(t{1});
                testCase.verifyEmpty(which(t{1}));
            end
        end
    end
end
