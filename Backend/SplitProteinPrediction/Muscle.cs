using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Threading;
using Microsoft.Win32.SafeHandles;
using System.Configuration;

namespace SplitProteinPrediction {

    class Muscle {

       
        public List<string> GetHmmerDataForMSA(string OriginalSeq, List<Hmmer.HmmerAlignment> HmmerAlign) {
            List<string> HitStrings = new List<string>();
            HitStrings.Add(">OriginSeq\n" + OriginalSeq);
            int hitcount = 1;
            foreach (Hmmer.HmmerAlignment SingleAlign in HmmerAlign) {
                hitcount++;
                HitStrings.Add("> " + SingleAlign.Name + "\n" + SingleAlign.MatchSeq);
            }

            return HitStrings;
        }

        public List<string> SampleSeq(List<string> input) {
            List<string> output = new List<string>();
            bool sample = true;
            float select_sequences = 2000f;
            int has_sequences_nbr = input.Count();
            float interval = 6f;

            if (select_sequences > has_sequences_nbr) {
                select_sequences = has_sequences_nbr;
            }

            if (sample == true) {
                interval = (has_sequences_nbr / select_sequences);
                if (interval < 1f) {
                    interval = 1f;
                    return input;
                }
                select_sequences = select_sequences * interval;
            } else {
                interval = 1f;
                return input;
            }
            //select_sequences
            for (float i = 0f; i < select_sequences; i += interval) {
                output.Add(input[(int)MathF.Floor(i)]);
            }
            return output;
        }
        public string RunMuscle(string UniqueID, List<string> HitStrings, string savedir) {
            //generate fasta with unique id:

            string path = savedir + UniqueID + "_muscleinput.fasta";
            string pathresult = savedir + UniqueID + "_MSA.clw";
            string strCmdText = "muscle -in " + path + " -out " + pathresult;
            string terminal = "/bin/bash";

            //write file
            if (!File.Exists(path)) {
                string content = "";
                foreach (string str in HitStrings) {
                    content += str + "\n";
                }
                File.WriteAllText(path, content);
            }
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
            
            return pathresult;
        }

    }
}
