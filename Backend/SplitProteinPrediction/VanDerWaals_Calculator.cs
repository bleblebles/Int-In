using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;


namespace SplitProteinPrediction
{
    class VanDerWaals_Calculator
    {
        private Dictionary<Vector3, int> GetAtomPositions(int StartIndex, int EndIndex, List<string> AtomNames, List<Vector3> AtomPos)
        {
            Dictionary<Vector3, int> result = new Dictionary<Vector3, int>();
            List<string> PossibleAtoms = new List<string>() { "GLN NE1", "GLN OE1", "ASN ND2", "ASN OD1" };
            for (int index_curr_res = StartIndex; index_curr_res <= EndIndex; index_curr_res++)
            {
                string currentAtomName = AtomNames[index_curr_res];
                List<string> split_res = currentAtomName.Split(" ").ToList();
                if (PossibleAtoms.Contains(currentAtomName) || split_res[1][0].ToString() == "C" || split_res[1][0].ToString() == "S")
                {
                    result.Add(AtomPos[index_curr_res], index_curr_res); //Only atom name here
                }
            }
            return result;
        }

        public List<List<string>> GenerateVanderwaalsBonds(PDBContent Content, bool UniqueConnections, float Distancevdw)
        {//Also known as ionic bond...
            PDBParser PDBPars = new PDBParser();

            List<List<string>> AromaticConnections = new List<List<string>>();

            List<string> Sequence = Content.SingleLetterSequence;
            int MaxIndex = Sequence.Count;
            List<List<int>> ResContacts = Content.ResidueContacts;
            List<int> SplitAtSite = Content.SplitAtSite;
            List<Vector3> AtomPositions = Content.AtomPositions;
            List<string> AtomNames = Content.AtomNames;

            for (int curr_index = 0; curr_index < MaxIndex; curr_index++)
            {
                string ResidueCurrent = Sequence[curr_index];
                int StartLineIndexCurrRes = 0;
                if (curr_index != 0)
                {
                    StartLineIndexCurrRes = SplitAtSite[curr_index - 1] - 1;
                }
                int EndLineIndexCurrRes = SplitAtSite[curr_index] - 2;
                Dictionary<Vector3, int> PositionCurrentRes = new Dictionary<Vector3, int>();
                PositionCurrentRes = GetAtomPositions(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames, AtomPositions);

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
                        Dictionary<Vector3, int> PositionContactRes = new Dictionary<Vector3, int>();
                        PositionContactRes = GetAtomPositions(StartLineIndexContact, EndLineIndexContact, AtomNames, AtomPositions);

                        List<string> AtomsPartners = new List<string>() { "C", "C" };
                        foreach (KeyValuePair<Vector3, int> PosCurrent in PositionCurrentRes)
                        {
                            float RadiusCurrRes = Content.RadiiList[PosCurrent.Value];
                            foreach (KeyValuePair<Vector3, int> PosContact in PositionContactRes)
                            {
                                //Fix for one C must be involved here!!!
                                string currentAtomName_1 = AtomNames[PosCurrent.Value];
                                List<string> split_res_1 = currentAtomName_1.Split(" ").ToList();
                                string currentAtomName_2 = AtomNames[PosContact.Value];
                                List<string> split_res_2 = currentAtomName_2.Split(" ").ToList();
                                if (split_res_1[1][0].ToString() == "C" || split_res_2[1][0].ToString() == "C") {
                                    float Distance = Vector3.Distance(PosCurrent.Key, PosContact.Key);
                                    //get the radius of both , accessible via the index
                                    float RadiusContactRes = Content.RadiiList[PosContact.Value];
                                    float DistBetweenSurface = Distance - RadiusContactRes - RadiusCurrRes;
                                    if (DistBetweenSurface < Distancevdw) {
                                        List<string> AromaticPartners = new List<string>();
                                        AromaticPartners.Add(curr_index + ".C");
                                        AromaticPartners.Add(Contact + ".C");
                                        //Console.WriteLine(Distance);
                                        AromaticConnections.Add(AromaticPartners);
                                    }
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
