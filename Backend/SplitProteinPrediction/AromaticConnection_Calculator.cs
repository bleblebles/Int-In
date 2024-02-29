using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;

namespace SplitProteinPrediction
{
    class AromaticConnection_Calculator
    {
        private Dictionary<Vector3, string> GetAromatePositions(int StartIndex, int EndIndex, List<string> AtomNames, List<Vector3> AtomPos)
        {
            float NumberAtoms = 0f;
            Vector3 SumVector = new Vector3();
            Dictionary<string, List<string>> PosChargeProtons = new Dictionary<string, List<string>>() { { "TYR", new List<string>() {"CG","CD1", "CE1", "CZ", "CE2", "CD2" } },
                                                                                                         { "HIS", new List<string>() {"CG","CD2", "NE2", "CE1", "ND1"} },
                                                                                                         { "PHE", new List<string>() {"CG", "CD2","CD1", "CE1", "CZ", "CE2" } },
                                                                                                         { "TRP", new List<string>() {"CG","CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2" } }
                                                                                                         };
            for (int index_curr_res = StartIndex; index_curr_res <= EndIndex; index_curr_res++)
            {
                string currentAtomName = AtomNames[index_curr_res];
                List<string> split_res = currentAtomName.Split(" ").ToList();
                List<string> ImportantAtoms = PosChargeProtons[split_res[0]];
                if (ImportantAtoms.Contains(split_res[1]))
                {
                    SumVector += AtomPos[index_curr_res];
                    NumberAtoms++;
                }
            }
            //now calculate the mean to get thze center of the aromate
            Vector3 ResultVect = SumVector / NumberAtoms;
            Dictionary<Vector3, string> result = new Dictionary<Vector3, string>() { { ResultVect, "CG" } };
            return result;
        }

        private Dictionary<Vector3, string> GetCationPositions(int StartIndex, int EndIndex, List<string> AtomNames, List<Vector3> AtomPos)
        {
            Dictionary<Vector3, string> result = new Dictionary<Vector3, string>();
            List<string> CationAtoms = new List<string>() { "LYS NZ",
                                                            "ARG NH1", "ARG NH2", "ARG NE"};
            for (int index_curr_res = StartIndex; index_curr_res <= EndIndex; index_curr_res++)
            {
                string currentAtomName = AtomNames[index_curr_res];
                List<string> split_res = currentAtomName.Split(" ").ToList();
                if (CationAtoms.Contains(currentAtomName))
                {
                    result.Add(AtomPos[index_curr_res], split_res[1]);
                }
            }
            return result;
        }

        public List<List<string>> GenerateAromaticConnections(PDBContent Content, bool UniqueConnections, float DistancePiPi, float DistanceCatPi)
        {//Also known as ionic bond...
            PDBParser PDBPars = new PDBParser();

            List<List<string>> AromaticConnections = new List<List<string>>();

            List<string> Aromates = new List<string>() { "H", "F", "W", "Y" };// All Aromatic Residues
            List<string> Cations = new List<string>() { "K", "R" };//H is excludet since it's also an aromate

            /*
             * Rules: 1) Aromates make pi-pi bonds when both are neutral or one is positively charged -> Not when both are charged!
             * 2) Aromates bond with cations (if positively charged) when aromates are not charged
            */
            float maxValueAromat = DistancePiPi;
            float CationAromatMaxVal = DistanceCatPi;

            List<string> Sequence = Content.SingleLetterSequence;
            int MaxIndex = Sequence.Count;
            List<List<int>> ResContacts = Content.ResidueContacts;
            List<int> SplitAtSite = Content.SplitAtSite;
            List<Vector3> AtomPositions = Content.AtomPositions;
            List<string> AtomNames = Content.AtomNames;

            for (int curr_index = 0; curr_index < MaxIndex; curr_index++)
            {
                string ResidueCurrent = Sequence[curr_index];
                int CurrentResidueType = 0;
                int StartLineIndexCurrRes = 0;
                if (curr_index != 0)
                {
                    StartLineIndexCurrRes = SplitAtSite[curr_index - 1] - 1;
                }
                int EndLineIndexCurrRes = SplitAtSite[curr_index] - 2;
                Dictionary<Vector3, string> PositionCurrentRes = new Dictionary<Vector3, string>();
                if (Cations.Contains(ResidueCurrent))
                {
                    CurrentResidueType = 1;
                    //Get Cathion Position
                    PositionCurrentRes = GetCationPositions(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames, AtomPositions);

                }
                else if (Aromates.Contains(ResidueCurrent))
                {
                    CurrentResidueType = 2;
                    PositionCurrentRes = GetAromatePositions(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames, AtomPositions);
                    //Get Aromatic Position
                }
                if (CurrentResidueType != 0)
                {
                    foreach (int Contact in ResContacts[curr_index])
                    {//look at each contact of the current residue
                        if (curr_index < Contact || UniqueConnections == false)
                        {
                            string ResidueContact = Sequence[Contact];
                            int StartLineIndexContact = 0;
                            if (Contact != 0)
                            {
                                StartLineIndexContact = SplitAtSite[Contact - 1] - 1;
                            }
                            int EndLineIndexContact = SplitAtSite[Contact] - 2;
                            int ContactResidueType = 0;
                            Dictionary<Vector3, string> PositionContactRes = new Dictionary<Vector3, string>();
                            if (Cations.Contains(ResidueContact))
                            {
                                ContactResidueType = 1;
                                PositionContactRes = GetCationPositions(StartLineIndexContact, EndLineIndexContact, AtomNames, AtomPositions);
                            }
                            else if (Aromates.Contains(ResidueContact))
                            {
                                ContactResidueType = 2;
                                PositionContactRes = GetAromatePositions(StartLineIndexContact, EndLineIndexContact, AtomNames, AtomPositions);
                            }
                            int IsPotentialAromaticConnection = ContactResidueType + CurrentResidueType;
                            if (IsPotentialAromaticConnection >= 3)
                            {//Potential Aromatic connection when 3+ (2 only when one aromate or two cations...)
                             //Get Aromate Charge (If Positive, then no connection with cation but still with aromate if the other one isn't charged)
                             //Tyr may be negatic charge so that one can still have a pi-cation bond but will also have a salt bridge at the same time!!

                                //Get The Closest Position between PositionCurrentRes & PositionContactRes
                                //Now look through the distances between all the possible Negative to Positive Atom Combinations
                                float BestDistance = 1000f;
                                List<string> AtomsPartners = new List<string>() { "CA", "CA" };
                                foreach (KeyValuePair<Vector3, string> PosCurrent in PositionCurrentRes)
                                {
                                    foreach (KeyValuePair<Vector3, string> PosContact in PositionContactRes)
                                    {
                                        float Distance = Vector3.Distance(PosCurrent.Key, PosContact.Key);
                                        if (Distance < BestDistance)
                                        {
                                            BestDistance = Distance;
                                            AtomsPartners[0] = PosCurrent.Value;
                                            AtomsPartners[1] = PosContact.Value;
                                        }
                                    }
                                }
                                if ((BestDistance <= maxValueAromat && IsPotentialAromaticConnection == 4) || (BestDistance <= CationAromatMaxVal && IsPotentialAromaticConnection == 3))
                                {
                                    List<string> AromaticPartners = new List<string>();
                                    AromaticPartners.Add(curr_index + "." + AtomsPartners[0]);
                                    AromaticPartners.Add(Contact + "." + AtomsPartners[1]);
                                    AromaticConnections.Add(AromaticPartners);
                                }
                            }
                        }
                    }
                }
            }

            return AromaticConnections;
        }
    }
}
