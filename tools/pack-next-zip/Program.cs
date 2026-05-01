using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Windows.Forms;

internal static class Program
{
    private static readonly string[] DirectoriesToPack =
    {
        "artifacts",
        "cases",
        "docs",
        "src",
        "tests"
        /*ANALYSIS_DIRECTORY*/
    };

    private static readonly string[] RootFilesToPack =
    {
        "AGENTS.md",
        "README.md",
        "CMakeLists.txt",
        "CMakeLists.md"
    };

    [STAThread]
    private static int Main()
    {
        var executableTitle = Path.GetFileNameWithoutExtension(Application.ExecutablePath);

        try
        {
            var rootPath = AppDomain.CurrentDomain.BaseDirectory;
            CreateArchive(rootPath);
            return 0;
        }
        catch (Exception exception)
        {
            MessageBox.Show(
                exception.Message,
                executableTitle,
                MessageBoxButtons.OK,
                MessageBoxIcon.Error);
            return 1;
        }
    }

    private static void CreateArchive(string rootPath)
    {
        var resolvedRoot = Path.GetFullPath(rootPath);
        var availableDirectories = DirectoriesToPack
            .Where(directoryName => Directory.Exists(Path.Combine(resolvedRoot, directoryName)))
            .ToList();
        var availableFiles = RootFilesToPack
            .Where(fileName => File.Exists(Path.Combine(resolvedRoot, fileName)))
            .ToList();

        if (availableDirectories.Count == 0 && availableFiles.Count == 0)
        {
            throw new InvalidOperationException(
                string.Format("No requested directories or files were found in '{0}'.", resolvedRoot));
        }

        var nextNumber = GetNextArchiveNumber(resolvedRoot);
        var archivePath = Path.Combine(resolvedRoot, string.Format("{0}.zip", nextNumber));

        using (var stream = new FileStream(archivePath, FileMode.CreateNew, FileAccess.ReadWrite, FileShare.None))
        using (var archive = new ZipArchive(stream, ZipArchiveMode.Create))
        {
            foreach (var directoryName in availableDirectories)
            {
                AddDirectoryRecursive(archive, resolvedRoot, directoryName);
            }

            foreach (var fileName in availableFiles)
            {
                var fullPath = Path.Combine(resolvedRoot, fileName);
                var normalizedEntryName = NormalizeEntryName(fileName);
                archive.CreateEntryFromFile(fullPath, normalizedEntryName, CompressionLevel.Fastest);
            }
        }
    }

    private static long GetNextArchiveNumber(string rootPath)
    {
        var numericZipNames = Directory.EnumerateFiles(rootPath, "*.zip", SearchOption.TopDirectoryOnly)
            .Select(Path.GetFileNameWithoutExtension)
            .Select(baseName =>
            {
                long value;
                return long.TryParse(baseName, out value) ? (long?)value : null;
            })
            .Where(value => value.HasValue)
            .Select(value => value.Value)
            .ToList();

        return numericZipNames.Count == 0 ? 1L : numericZipNames.Max() + 1L;
    }

    private static void AddDirectoryRecursive(ZipArchive archive, string sourceRoot, string relativeDirectoryPath)
    {
        var fullDirectoryPath = Path.Combine(sourceRoot, relativeDirectoryPath);
        var childDirectories = Directory.GetDirectories(fullDirectoryPath);
        var childFiles = Directory.GetFiles(fullDirectoryPath);

        if (childDirectories.Length == 0 && childFiles.Length == 0)
        {
            archive.CreateEntry(NormalizeEntryName(relativeDirectoryPath) + "/");
            return;
        }

        foreach (var childFile in childFiles)
        {
            var entryName = NormalizeEntryName(GetRelativePath(sourceRoot, childFile));
            archive.CreateEntryFromFile(childFile, entryName, CompressionLevel.Fastest);
        }

        foreach (var childDirectory in childDirectories)
        {
            var childRelativePath = GetRelativePath(sourceRoot, childDirectory);
            AddDirectoryRecursive(archive, sourceRoot, childRelativePath);
        }
    }

    private static string GetRelativePath(string basePath, string fullPath)
    {
        var trimmedBasePath = basePath.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
        if (!trimmedBasePath.EndsWith(Path.DirectorySeparatorChar.ToString(), StringComparison.Ordinal))
        {
            trimmedBasePath += Path.DirectorySeparatorChar;
        }

        var baseUri = new Uri(trimmedBasePath, UriKind.Absolute);
        var targetUri = new Uri(fullPath, UriKind.Absolute);
        return Uri.UnescapeDataString(baseUri.MakeRelativeUri(targetUri).ToString())
            .Replace('/', Path.DirectorySeparatorChar);
    }

    private static string NormalizeEntryName(string path)
    {
        return path.Replace(Path.DirectorySeparatorChar, '/')
            .Replace(Path.AltDirectorySeparatorChar, '/');
    }
}
