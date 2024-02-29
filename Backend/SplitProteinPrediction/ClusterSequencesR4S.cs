using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Threading;
using System.Configuration;
using Microsoft.Win32.SafeHandles;

namespace SplitProteinPrediction {
    class ClusterSequencesR4S {

        private List<string> masa_names = new List<string>();
        private Dictionary<string, string> MSADict = new Dictionary<string, string>();//key = name, value=sequence

        public void GetMSAList(string path) {
            masa_names.Clear();
            MSADict.Clear();

            using (StreamReader files = new StreamReader(path)) {

                string name = "NONE";
                string ln;

                while ((ln = files.ReadLine()) != null) {
                    if (ln.Trim() == "") {
                        continue;
                    }
                    if (ln.Contains("CLUSTAL")) {
                        continue;
                    }
                    if (ln.Contains("MUSCLE")) {
                        continue;
                    }
                    if (ln.Contains("*")) {
                        continue;
                    }
                    //Add to the lists:
                    if (string.Join("", ln.Take(1).ToArray()) == ">") {
                        name = ln.Replace(">", "").Trim();
                        if (!masa_names.Contains(name)) {
                            masa_names.Add(name);
                        }
                    } else {
                        if (MSADict.ContainsKey(name)) {
                            MSADict[name] += ln.Trim();
                        } else {
                            MSADict.Add(name, ln.Trim());
                        }

                    }
                }
                files.Close();

            }

        }


        static double StringCompare(string a, string b) {

            //Problem: "-" are falsifying the score but they are needed to compare sequences


            if (a == b) //Same string, no iteration needed.
                return 100;
            if ((a.Length == 0) || (b.Length == 0)) //One is empty, second is not
            {
                return 0;
            }
            double maxLen = a.Length > b.Length ? a.Length : b.Length;
            double minLen = a.Length < b.Length ? a.Length : b.Length;
            int sameCharAtIndex = 0;
            int Comparisions = 0;
            string comparea = "";
            string compareb = "";
            for (int i = 0; i < minLen; i++) //Compare char by char
            {
                if (a[i].ToString() == "-" && b[i].ToString() == "-") {
                    //Do nothing
                } else {

                    Comparisions++;
                    if (a[i] == b[i]) {
                        sameCharAtIndex++;
                    }
                }

            }

            return sameCharAtIndex / (float)Comparisions * 100f;
        }



        List<string> ClusterSequences = new List<string>();
        List<string> ClusterNames = new List<string>();
        public void ClusterSeq() {
            //int SeqMaxNbr = 900;
            ClusterSequences.Clear();
            ClusterNames.Clear();
            //int ClusterCount = 0;
            List<string> CurrentSequences = new List<string>();
            foreach (KeyValuePair<string, string> entry in MSADict)//go through each 
            {
                string SeqName = entry.Key;
                if (SeqName != "OriginSeq") {
                    string Seq = entry.Value;
                    bool ClusterExists = false;
                    int index = 0;
                    foreach (string CurrentSeq in ClusterSequences) {
                        double PercentageSame = StringCompare(Seq, CurrentSeq);
                        if (PercentageSame >= 95f) {
                            ClusterExists = true;
                        }
                        index++;
                    }
                    if (ClusterExists == false) {
                        ClusterNames.Add(SeqName);
                        ClusterSequences.Add(Seq);

                    }
                }
            }
        }
        public void GenerateClusteredAlignment(string UniqueID, string savedir) {
            string path = savedir + UniqueID + "_MSA.clw";
            string path_save = savedir + UniqueID + "_MSA_Cluster.clw";


            GetMSAList(path);//Generate the dictionary with the sequences
            ClusterSeq();//now we cluster all the sequences

            Dictionary<string, string> FinalDictionary = new Dictionary<string, string>();//key = name, value=sequence

            bool sample = true;
            float select_sequences = 150f;
            int has_sequences_nbr = ClusterNames.Count();
            float interval = 6f;

            if (select_sequences > has_sequences_nbr) {
                select_sequences = has_sequences_nbr;
            }

            if (sample == true) {
                interval = has_sequences_nbr / select_sequences;
                if (interval < 1f) {
                    interval = 1f;
                }
                select_sequences = select_sequences * interval;
            } else {
                interval = 1f;
            }
            //select_sequences
            for (float i = 0f; i < select_sequences; i += interval) {
                FinalDictionary.Add(ClusterNames[(int)MathF.Floor(i)], ClusterSequences[(int)MathF.Floor(i)]);
            }
            string FileContent = ">OriginSeq\n" + MSADict["OriginSeq"] + "\n";
            foreach (KeyValuePair<string, string> entry in FinalDictionary)//go through each 
            {
                //make the new file
                string SeqName = entry.Key;
                string Seq = entry.Value;
                FileContent += ">" + SeqName + "\n" + Seq + "\n";
            }

            File.WriteAllText(path_save, FileContent);

        }

    }
}
