using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using Ionic.Zip;
using System.Runtime.CompilerServices;

namespace SplitProteinPrediction {
    class OpenZip {

        public void OpenZipFile(string zipfile, string SaveDir) {
            // Create FileStream for output ZIP archive
            // Files to be added to archive
            string ZipFileDir = zipfile;
            string ExtractTo = SaveDir;

            if (File.Exists(ZipFileDir)) {
                bool WrongFileExtension = false;
                bool FileIsEncrypted = false;
                using (ZipFile zip = ZipFile.Read(ZipFileDir)) {
                    foreach (ZipEntry e in zip) {
                        if (e.UsesEncryption == true) {
                            FileIsEncrypted = true;
                        }
                        if (!e.FileName.EndsWith(".pdb")) {
                            WrongFileExtension = true;
                            break;
                        }
                    }
                    if (FileIsEncrypted == false) {
                        if (WrongFileExtension == false) {
                            if (Directory.Exists(ExtractTo)) {
                                throw new SplitProteinException("Directory already exists...");
                            } else {
                                Directory.CreateDirectory(ExtractTo);
                                zip.ExtractSelectedEntries("name = *.pdb", "", ExtractTo, ExtractExistingFileAction.OverwriteSilently);
                            }
                        } else {
                            throw new SplitProteinException("The archive contains non pdb files");
                        }
                    } else {
                        throw new SplitProteinException("Zip Files are password protected");
                    }
                }
            }
        }
    }
}
