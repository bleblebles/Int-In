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
using System.Security;
using System.Configuration;

namespace SplitProteinPrediction {
    static class Async {
        public static Task WaitForExitAsync(this Process process, CancellationToken cancellationToken = default(CancellationToken)) {
            var tcs = new TaskCompletionSource<object>();
            process.EnableRaisingEvents = true;
            process.Exited += (sender, args) => tcs.TrySetResult(null);
            if (cancellationToken != default(CancellationToken))
                cancellationToken.Register(tcs.SetCanceled);
            return tcs.Task;
        }
    }
    class Hmmer {


        public class HmmerAlignment {
            public string Name;
            public string Alignment;
            public string MatchSeq;
        }

        public class DomainAnnotation {
            public string DomainName;
            public float Score;
            public int SeqStart;
            public int SeqEnd;
        }


        public async Task<int> RunHmmer(string UniqueID, string sequence, string dirsave) {
            //generate fasta with unique id:
            string HmmerUniref90 = ConfigurationManager.AppSettings.Get("HmmerUniref90");

            string path = dirsave + UniqueID + "_originseq.fasta";

            string terminal = "/bin/bash";
            string strCmdText = "phmmer -o " + dirsave + UniqueID + "_pairwise.hmmer " + path + " " + HmmerUniref90;
                
            if (!File.Exists(path)) {
                using (StreamWriter sw = File.CreateText(path)) {
                    sw.WriteLine(">MySequence");
                    sw.WriteLine(sequence);
                }
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
            await bash.WaitForExitAsync();
            
            File.Delete(path);
            //now it's finished with the process...
            return 1;
        }

        public List<HmmerAlignment> ReadFile(string UniqueID, List<string> SequenceOriginal, string dirsave) {
            string currentDecSep = Thread.CurrentThread.CurrentCulture.NumberFormat.NumberDecimalSeparator.ToString();
            string ReplaceSeparator = ".";
            if (currentDecSep == ".") {
                ReplaceSeparator = ",";
            }

            string path = dirsave + UniqueID + "_pairwise.hmmer";

            string line;
            bool AlignmentStart = false;
            StreamReader file = new StreamReader(path);

            string HitNameBase = "";
            string HitName = "";
            int line_count = 0;

            int SeqLine = 0;

            bool ScoreFullSeq = false;

            string SeqMatch = "";
            string Alignment = "";

            bool foobar = false;
            int indexfoobar = 0;

            List<string> MatchNameList = new List<string>();
            List<string> MatchList = new List<string>();
            List<string> AlignmentList = new List<string>();
            Dictionary<string, float> EValueHit = new Dictionary<string, float>();

            float eval_cutoff = 0.0001f;
            float Current_eval = 0f;
            int CurrSeq = 0;
            float Score = 0f;
            int GetEvalIndex = 0;
            int GetSequenceNameIndex = 0;
            bool StartGettingEvalues = false;
            string save_text = "";

            while ((line = file.ReadLine()) != null) {
                line_count++;
                if (AlignmentStart == false) {
                    if (String.Join("", line.Take(6).ToArray()) == "Domain") {
                        AlignmentStart = true;
                    } else {
                        //Get the sequence names:
                        if (StartGettingEvalues == true) {
                            if (line.Length >= GetEvalIndex + 7) {
                                float EvalueFloatCurrent;
                                string getEvalue = String.Join("", line.Take(GetEvalIndex + 7).ToArray()).Replace(ReplaceSeparator, currentDecSep);
                                getEvalue = line.Substring(GetEvalIndex, 8); ;
                                //Console.WriteLine(getEvalue);
                                bool res = float.TryParse(getEvalue, out EvalueFloatCurrent);
                                if (res == true) {//is okay..
                                    save_text += "Read E-Value: "+ EvalueFloatCurrent + "\n";
                                    string getSequenceNameplus = String.Join("", line.Skip(GetSequenceNameIndex).ToArray());
                                    string SequenceBaseName = getSequenceNameplus.Split()[0];
                                    EValueHit.Add(SequenceBaseName, EvalueFloatCurrent);
                                    //Console.WriteLine(SequenceBaseName + " -> " + EvalueFloatCurrent);
                                }
                            }
                        }
                        if (line.Contains("E-value")) {
                            GetEvalIndex = line.IndexOf("E-value", line.IndexOf("E-value") + 1) - 1;//-1 because that sometimes happens
                            GetSequenceNameIndex = line.IndexOf("Sequence");
                            StartGettingEvalues = true;
                        }
                    }
                } else {
                    //Get the name of the Hit or if it's multiple domains add it to the name
                    if (String.Join("", line.Take(3).ToArray()) == ">> " || String.Join("", line.Take(4).ToArray()) == "  ==") {//Get 
                        if (Alignment != "") {
                            //Next line in the alignment
                            string NoPlusSpaces = Alignment.Replace(" ", "").Replace("+", "");
                            int IdenticalLength = NoPlusSpaces.Length;
                            int AlignmentLength = SeqMatch.Replace("-", "").Length;
                            Score = (float)IdenticalLength / (float)AlignmentLength;
                            if (ScoreFullSeq == true) {
                                Score = (float)SequenceOriginal.Count() / (float)AlignmentLength;
                            }

                            if (Score >= 0.35f && AlignmentLength >= 70) {//&& Score <= 0.95
                                MatchNameList.Add(HitName);
                                MatchList.Add(SeqMatch);
                                AlignmentList.Add(Alignment);
                                CurrSeq++;
                            }
                        }
                        if (String.Join("", line.Take(3).ToArray()) == ">> ") {
                            //only change hitname when new entry
                            HitNameBase = line.Split()[1];
                            HitName = HitNameBase;
                            Current_eval = EValueHit[HitNameBase];
                        }
                        if (String.Join("", line.Take(4).ToArray()) == "  ==") {
                            //New Domain, so change it's name

                            if (Current_eval <= eval_cutoff) {
                                string domainNbr = String.Join("", line.Skip(12).Take(1).ToArray());
                                HitName = HitNameBase + "_" + domainNbr;
                                foobar = true;
                                indexfoobar = 0;
                            }
                        }

                        SeqMatch = "";
                        Alignment = "";
                    }

                    if (line.Contains("Sequence")) {
                        SeqLine = line_count;
                        //Original seq
                    } else if (line_count == SeqLine + 1) {
                        //Alignment
                        if (foobar == true) {
                            //save the length:
                            indexfoobar = line.TakeWhile(Char.IsWhiteSpace).Count();
                            foobar = false;
                        }
                        if (Current_eval <= eval_cutoff) {
                            Alignment += String.Join("", line.Skip(indexfoobar).ToArray());
                        }
                    } else if (line_count == SeqLine + 2) {
                        if (Current_eval <= eval_cutoff) {
                            //HitSequence
                            string SelectMatch = String.Join("", line.Skip(2).ToArray()).Replace(HitNameBase, "").Replace(" ", "");
                            string NoNumbers = Regex.Replace(SelectMatch, "[0-9]", string.Empty);//delete all numbers
                            SeqMatch += NoNumbers;
                        }
                    }
                }
            }
            int LengthAlignment2 = SeqMatch.Replace("-", "").Length;
            if (Score >= 0.35f && LengthAlignment2 >= 70) {//&& Score <= 0.95
                MatchNameList.Add(HitName);
                MatchList.Add(SeqMatch);
                AlignmentList.Add(Alignment);
            }
            List<HmmerAlignment> HMMERAlign = new List<HmmerAlignment>();

            int select_sequences = MatchNameList.Count();//150
            //select_sequences
            for (int i = 0; i < select_sequences; i += 1) {//do them all
                HMMERAlign.Add(new HmmerAlignment() {
                    Name = MatchNameList[i],
                    MatchSeq = MatchList[i],
                    Alignment = AlignmentList[i]
                });
            }


            //delete file...
            file.Close();
            return HMMERAlign;
        }

    }
}
