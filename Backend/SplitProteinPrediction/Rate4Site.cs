using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Threading;
using Microsoft.Win32.SafeHandles;
using System.Configuration;

namespace SplitProteinPrediction {

    class Rate4Site {



        public string RunRate4Site(string UniqueID, string savedir) {
            string path_output = savedir + UniqueID + "_Rate.res";

            string strCmdText = "rate4site -s " + savedir + UniqueID + "_MSA_Cluster.clw -o " + savedir + UniqueID + "_Rate.res -a OriginSeq";
            string terminal = "/bin/bash";

            Process bash = new Process();
            bash.StartInfo.FileName = terminal;
            bash.StartInfo.RedirectStandardInput = true;
            bash.StartInfo.RedirectStandardOutput = true;
            bash.StartInfo.CreateNoWindow = false;
            bash.StartInfo.UseShellExecute = false;
            bash.Start();
            bash.StandardInput.WriteLine(strCmdText);
            bash.StandardInput.Flush();
            bash.StandardInput.Close();
            bash.WaitForExit();
            
            return path_output;
        }


        public List<string> ReadRate4Site(string path) {
            // Read file using StreamReader. Reads file line by line    
            List<string> Score = new List<string>();
            using (StreamReader file = new StreamReader(path)) {
                int counter = 0;
                string ln;
                bool ReadScore = false;
                while ((ln = file.ReadLine()) != null) {
                    //ln
                    if (ln.Contains("#Average")) {
                        ReadScore = false;
                    }
                    if (ReadScore == true) {
                        List<string> Values = ln.Split("[")[0].Split(" ").ToList();
                        string splitscore = "";
                        foreach (string single_char in Values) {
                            if (single_char != "") {
                                splitscore = single_char;
                            }
                        }
                        Score.Add(splitscore);
                    }
                    if (ln.Contains("#LL")) {
                        ReadScore = true;
                    }
                }
                file.Close();
            }

            return Score;
        }

        public List<double> NormalizeRate4Site(List<string> Rate4SiteValues) {
            string currentDecSep = Thread.CurrentThread.CurrentCulture.NumberFormat.NumberDecimalSeparator.ToString();
            string ReplaceSeparator = ".";
            if (currentDecSep == ".") {
                ReplaceSeparator = ",";
            }

            List<double> Rate4SiteValuesDoubles = Rate4SiteValues.Select(x => double.Parse(x.Replace(ReplaceSeparator, currentDecSep))).ToList();

            List<double> NormR4SVals = new List<double>();

            List<double> PositiveValues = Rate4SiteValuesDoubles.Where(x => x > 0).ToList();
            List<double> NegativeValues = Rate4SiteValuesDoubles.Where(x => x <= 0).ToList();
            double max_PosValues = PositiveValues.Max();
            double min_PosValues = PositiveValues.Min();
            double max_NegValues = NegativeValues.Max();
            double min_NegValues = NegativeValues.Min();

            NormR4SVals = (from i in Rate4SiteValuesDoubles select i >= 0 ? (1 - (i - min_PosValues) / (max_PosValues - min_PosValues)) * 0.5 : (1 - (i - min_NegValues) / (max_NegValues - min_NegValues)) * 0.5 + 0.5).ToList();

            return NormR4SVals;
        }

    }
}
