using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Configuration;

namespace SplitProteinPrediction {
    class Run_DSSP {

        public PDBContent DSSP_Cmd(PDBContent cont, string file, string UniqueID, string ResultsDir) {

            string saveFile = ResultsDir + UniqueID + ".dssp";
            while (File.Exists(saveFile)) {
                Console.WriteLine("DSSP file already existed deleting old file");
                File.Delete(saveFile);
            }
            string strCmdText = "dssp -i \"" + file + "\" -o \"" + saveFile + "\"";

            List<string> DSSP_Chars = new List<string>();
            Process bash = new Process();
            string terminal = "/bin/bash";

            bash.StartInfo.FileName = terminal;
            bash.StartInfo.RedirectStandardInput = true;
            bash.StartInfo.RedirectStandardOutput = true;
            bash.StartInfo.CreateNoWindow = true;
            bash.StartInfo.UseShellExecute = false;
            bash.Start();

            bash.StandardInput.WriteLine(strCmdText);
            bash.StandardInput.Flush();
            bash.StandardInput.Close();
            bash.WaitForExit();

            //read output file and then delete it
            if (File.Exists(saveFile)) {
                // Read a text file line by line.
                bool saveStuff = false;
                string[] lines = File.ReadAllLines(saveFile);
                foreach (string line in lines) {
                    char[] charArr = line.ToCharArray();
                    if (saveStuff == true) {
                        string AminoAcid = String.Join("", charArr.Skip(13).Take(1).ToArray()).Trim();
                        string SecStructure = String.Join("", charArr.Skip(16).Take(1).ToArray()).Trim();
                        if (AminoAcid != "!") {// Aminoacid is ! when the chain stops, this is problematic because then it would be interpreted as a loop instead of being nothing
                            if (SecStructure == "S" || SecStructure == "T" || SecStructure == "C" || SecStructure == "") {
                                //Out of simplicity all loops will be denoted "c" and all secondary structures "b"
                                //It's a loop!
                                DSSP_Chars.Add("c");
                            } else {
                                DSSP_Chars.Add("b");
                            }
                        }
                    }
                    if (line.Contains("  #  RESIDUE")) {
                        saveStuff = true;
                    }
                }
                //File.Delete(saveFile);
            } else {
                throw new SplitProteinException("DSSP file doesn't exist!");
            }
            cont.SecondaryStructure = DSSP_Chars;

            return cont;
        }

    }
}
