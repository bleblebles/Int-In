using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;
using System.IO;
using System.Configuration;
using System.Globalization;
using System.Threading;

namespace SplitProteinPrediction {
    class RunPDB2PQR {
        public string Runcmd(string input_file) {
            string program_file = ConfigurationManager.AppSettings.Get("Pdb2pqr");
            string filename = Path.GetFileName(input_file);
            string path = input_file.Replace(filename, "");

            string output = "";
            string save_file = path + filename.Replace(".pdb", "_2.pdb");
            if (!File.Exists(save_file)) {
                string command = program_file + " --ff=parse --chain --ph-calc-method=propka --with-ph=7 " + input_file.ToString() + " " + save_file;
                Process bash = new Process();
                string terminal = "/bin/bash";

                //Console.WriteLine(command);
                bash.StartInfo.FileName = terminal;
                bash.StartInfo.RedirectStandardInput = true;
                bash.StartInfo.RedirectStandardOutput = true;
                bash.StartInfo.CreateNoWindow = false;
                bash.StartInfo.UseShellExecute = false;
                bash.Start();
                bash.StandardInput.WriteLine(command);
                bash.StandardInput.Flush();
                bash.StandardInput.Close();
                output = bash.StandardOutput.ReadToEnd();
                bash.WaitForExit();
            }
            if (output.Contains("error:")) {
                throw new SplitProteinException("pdb2pqr encountered an error.");
            }
            if (!File.Exists(save_file)) {
                throw new SplitProteinException("pdb2pqr encountered an error");
            }
            return save_file;
        }
    }
}
