using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;
using System.IO;
using System.Globalization;
using System.Threading;
using System.Text.RegularExpressions;

namespace SplitProteinPrediction {

    class PDBParser {

        /*Possible errors: When there are two versions of the same residue (same id) the the program adds them all to the radii instead of only one...*/
        public PDBContent GetPDBData(string file_dir, PDBContent PDBCont, string FilterChainIdent = "none") {
            //Load all the important data inside the PDb Container
            string currentDecSep = Thread.CurrentThread.CurrentCulture.NumberFormat.NumberDecimalSeparator.ToString();
            string ReplaceSeparator = ".";
            if (currentDecSep == ".") {
                ReplaceSeparator = ",";
            }

            PDBCont.FileDir = file_dir;
            bool TakeFirstCa = false;
            AA_Values AAVals = new AA_Values();
            string line;
            int line_count = 0;
            int CurrResNbr = 0;
            int index = 0;
            List<string> Resnbr_Chain = new List<string>();
            Dictionary<string, int> ResSeq_Visited = new Dictionary<string, int>();
            StreamReader file = new StreamReader(file_dir);
            while ((line = file.ReadLine()) != null) {
                char[] charArr = line.ToCharArray();
                if (charArr.Length > 4) {//No ATOM possible
                    string firstname = String.Join("", charArr.Take(4).ToArray());
                    if (firstname == "ATOM" && AAVals.AA_3LetterCodeToSingle.ContainsKey(String.Join("", charArr.Skip(17).Take(3).ToArray()).Trim())) {
                        string ChainIdentNow = String.Join("", charArr.Skip(21).Take(1).ToArray()).Trim().ToUpper();
                        if (FilterChainIdent == "none" || ChainIdentNow == FilterChainIdent.ToUpper()) {
                            line_count++;
                            if (charArr.Length >= 54) {//Need at least 53 chars -> 54
                                string ResSeqNbr = String.Join("", charArr.Skip(22).Take(6).ToArray()).Trim();
                                var matches = Regex.Matches(ResSeqNbr, @"-?\d+");
                                //This because of negative numbers!
                                string GetOnlyNumbers = matches[0].ToString();
                                string GetOnlyChars = string.Join("", ResSeqNbr.Where(c => char.IsLetter(c)).ToArray());
                                CurrResNbr = int.Parse(GetOnlyNumbers);
                                if (index != 0) {
                                    PDBCont.AtomIndexToResidueIndex.Add(index - 1);
                                }
                                //Console.WriteLine(ResSeqNbr + ChainIdentNow);
                                if (!ResSeq_Visited.ContainsKey(ResSeqNbr + ChainIdentNow)) {
                                    ResSeq_Visited.Add(ResSeqNbr + ChainIdentNow, 1);
                                    if (index != 0) {
                                        PDBCont.SplitAtSite.Add(line_count);
                                    }
                                    PDBCont.SingleLetterSequence.Add(AAVals.AA_3LetterCodeToSingle[String.Join("", charArr.Skip(17).Take(3).ToArray()).Trim()]);
                                    index++;
                                    TakeFirstCa = true;
                                }
                                PDBCont.AtomInfo.Add(line);
                                string ResNbrWithChain = CurrResNbr + GetOnlyChars + ":" + ChainIdentNow;
                                Resnbr_Chain.Add(ResNbrWithChain);
                                //Get the right radius:
                                var digits = new[] { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
                                string AtomName = String.Join("", charArr.Skip(12).Take(4).ToArray()).Trim().TrimStart(digits);
                                string ResidueName = String.Join("", charArr.Skip(17).Take(3).ToArray()).Trim();
                                float Radius = 1.87f;
                                string TogetherKey = ResidueName + " " + AtomName;
                                // Console.WriteLine("Key: " + TogetherKey);

                                if (AAVals.AtomVanDerWaalsv2.ContainsKey(TogetherKey)) {
                                    Radius = AAVals.AtomVanDerWaalsv2[TogetherKey];
                                } else {
                                    if (AtomName == "CA" || AtomName == "CB") {
                                        //take 2 first ones
                                        Radius = AAVals.AtomVanDerWaals[String.Join("", AtomName.ToCharArray().Take(2))];
                                    } else {
                                        //take first one
                                        string key = String.Join("", AtomName.ToCharArray().Take(1));
                                        if (AAVals.AtomVanDerWaals.ContainsKey(key)) {
                                            Radius = AAVals.AtomVanDerWaals[key];
                                        } else {
                                            Console.WriteLine(AtomName + " could not be found...");
                                        }
                                    }
                                }

                                //Console.WriteLine(TogetherKey + " -> " + Radius);
                                PDBCont.RadiiList.Add(Radius);
                                string XPos = String.Join("", charArr.Skip(30).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();
                                string YPos = String.Join("", charArr.Skip(38).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();
                                string ZPos = String.Join("", charArr.Skip(46).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();//delete 46 first elements and take next 8
                                                                                                                           //Check if the coordinates can be parsed
                                float XPosFloat;
                                float YPosFloat;
                                float ZPosFloat;
                                bool XParseBool = float.TryParse(XPos, out XPosFloat);
                                bool YParseBool = float.TryParse(YPos, out YPosFloat);
                                bool ZParseBool = float.TryParse(ZPos, out ZPosFloat);
                                if (XParseBool && YParseBool && ZParseBool) {
                                    PDBCont.AtomNames.Add(TogetherKey);
                                    PDBCont.AtomPositions.Add(new Vector3(float.Parse(XPos), float.Parse(YPos), float.Parse(ZPos)));
                                    if (AtomName == "CA" && TakeFirstCa == true) {
                                        int ParseResNbr;
                                        bool ParseResNbrBool = int.TryParse(ResSeqNbr, out ParseResNbr);
                                        if (ParseResNbrBool) {
                                            if (!PDBCont.ResidueIDs.Contains(ParseResNbr)) {
                                                PDBCont.ResidueIDs.Add(ParseResNbr);
                                            }
                                        }
                                        PDBCont.Ca_AtomPositions.Add(new Vector3(float.Parse(XPos), float.Parse(YPos), float.Parse(ZPos)));
                                        TakeFirstCa = false;//Since some pdbs have 2 versions of the Ca for the same Residue
                                    }
                                } else {
                                    throw new SplitProteinException("Coordinates: " + line + " could'nt be converted to floats");
                                }
                            } else {
                                throw new SplitProteinException("The ATOM: " + line + " doesn't contain all the necessary information");
                            }
                        }
                    }
                }
            }
            //Add the last Atom too:
            PDBCont.AtomIndexToResidueIndex.Add(index - 1);
            List<string> ResSeqNbrsUniqueChain = Resnbr_Chain.Distinct().ToList();
            PDBCont.ResNumbersPDBChain = ResSeqNbrsUniqueChain;

            PDBCont.SplitAtSite.Add(line_count + 1); //As a line always gives the starting point of the next residue, so with the last on do one more...
            file.Close();

            return PDBCont;
        }

        public PDBContent GetSplitData(List<string> SplitData, string SplitName, PDBContent PDBCont) {
            //Load all the important data inside the PDb Container
            string currentDecSep = Thread.CurrentThread.CurrentCulture.NumberFormat.NumberDecimalSeparator.ToString();
            string ReplaceSeparator = ".";
            if (currentDecSep == ".") {
                ReplaceSeparator = ",";
            }

            PDBCont.FileDir = SplitName;
            AA_Values AAVals = new AA_Values();
            bool TakeFirstCa = false;
            int line_count = 0;
            int index = 0;
            Dictionary<string, int> ResSeq_Visited = new Dictionary<string, int>();
            foreach (string line in SplitData) {
                char[] charArr = line.ToCharArray();
                if (charArr.Length > 4) {//No ATOM possible
                    string firstname = String.Join("", charArr.Take(4).ToArray());
                    if (firstname == "ATOM") {
                        string ChainIdentNow = String.Join("", charArr.Skip(21).Take(1).ToArray()).Trim().ToUpper();
                        line_count++;
                        if (charArr.Length >= 54) {//the 78th char is the last one...
                            string ResSeqNbr = String.Join("", charArr.Skip(22).Take(4).ToArray()).Trim();
                            if (index != 0) {
                                PDBCont.AtomIndexToResidueIndex.Add(index - 1);
                            }
                            if (!ResSeq_Visited.ContainsKey(ResSeqNbr + ChainIdentNow)) {
                                ResSeq_Visited.Add(ResSeqNbr + ChainIdentNow, 1);
                                if (index != 0) {
                                    PDBCont.SplitAtSite.Add(line_count);
                                }
                                PDBCont.SingleLetterSequence.Add(AAVals.AA_3LetterCodeToSingle[String.Join("", charArr.Skip(17).Take(3).ToArray()).Trim()]);
                                index++;
                                TakeFirstCa = true;
                            }
                            PDBCont.AtomInfo.Add(line);
                            //Get the right radius:
                            var digits = new[] { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
                            string AtomName = String.Join("", charArr.Skip(12).Take(4).ToArray()).Trim().TrimStart(digits); ;
                            string ResidueName = String.Join("", charArr.Skip(17).Take(3).ToArray()).Trim();
                            float Radius = 1.87f;
                            string TogetherKey = ResidueName + " " + AtomName;
                            if (AAVals.AtomVanDerWaalsv2.ContainsKey(TogetherKey)) {
                                Radius = AAVals.AtomVanDerWaalsv2[TogetherKey];
                            } else {
                                if (AtomName == "CA" || AtomName == "CB") {
                                    Radius = AAVals.AtomVanDerWaals[String.Join("", AtomName.ToCharArray().Take(2))];
                                    //Radius = AAVals.AtomVanDerWaals[String.Join("", charArr.Skip(13).Take(2).ToArray()).Trim()];
                                } else {
                                    string key = String.Join("", AtomName.ToCharArray().Take(1));
                                    if (AAVals.AtomVanDerWaals.ContainsKey(key)) {
                                        Radius = AAVals.AtomVanDerWaals[key];
                                    } else {
                                        Console.WriteLine(AtomName + " could not be found...");
                                    }
                                    //Radius = AAVals.AtomVanDerWaals[String.Join("", charArr.Skip(13).Take(1).ToArray()).Trim()];
                                }
                            }
                            PDBCont.RadiiList.Add(Radius);
                            string XPos = String.Join("", charArr.Skip(30).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();
                            string YPos = String.Join("", charArr.Skip(38).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();
                            string ZPos = String.Join("", charArr.Skip(46).Take(8).ToArray()).Replace(ReplaceSeparator, currentDecSep).Trim();
                            //Check if the coordinates can be parsed
                            float XPosFloat;
                            float YPosFloat;
                            float ZPosFloat;
                            bool XParseBool = float.TryParse(XPos, out XPosFloat);
                            bool YParseBool = float.TryParse(YPos, out YPosFloat);
                            bool ZParseBool = float.TryParse(ZPos, out ZPosFloat);
                            if (XParseBool && YParseBool && ZParseBool) {
                                PDBCont.AtomPositions.Add(new Vector3(float.Parse(XPos), float.Parse(YPos), float.Parse(ZPos)));
                                if (AtomName == "CA" && TakeFirstCa == true) {
                                    PDBCont.Ca_AtomPositions.Add(new Vector3(float.Parse(XPos), float.Parse(YPos), float.Parse(ZPos)));
                                    TakeFirstCa = false;
                                }
                            } else {
                                throw new SplitProteinException("Coordinates: " + line + " could'nt be converted to floats");
                            }
                        } else {
                            throw new SplitProteinException("The ATOM: " + line + " doesn't contain all the necessary information");
                        }
                    }
                }
            }
            //Add the last Atom too:
            PDBCont.AtomIndexToResidueIndex.Add(index - 1);
            PDBCont.SplitAtSite.Add(line_count + 1); //As a line always gives the starting point of the next residue, so with the last on do one more...
            return PDBCont;
        }
    }
}
