using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;

namespace SplitProteinPrediction {
    class SaltBridge_Calculator {

        //If one 
        public bool CheckNegativeChargedResProtonated(int StartIndex, int EndIndex, List<string> AtomNames) {
            List<string> ProtonationAtoms = new List<string>() { "ASP HD1", "ASP HD2",
                                                                 "GLU HE1", "GLU HE2",
                                                                 "CYS HG",
                                                                 "TYR HH"};//if one of those is there don't form a salt bridge
            for (int index_curr_res = StartIndex; index_curr_res <= EndIndex; index_curr_res++) {
                string currentAtomName = AtomNames[index_curr_res];
                if (ProtonationAtoms.Contains(currentAtomName)) {
                    return true;
                }
            }
            return false;
        }

        public bool CheckPositiveChargedResProtonated(int StartIndex, int EndIndex, List<string> AtomNames) {
            //If they are protonated then everything is okay
            Dictionary<string, List<string>> PosChargeProtons = new Dictionary<string, List<string>>() { { "LYS", new List<string>() {"HZ1","HZ2","HZ3" } },
                                                                                                         { "ARG", new List<string>() { "HH11", "HH12", "HH21", "HH22" } },
                                                                                                         { "HIS", new List<string>() { "HD1", "HE2"} }};
            int HydrogensFound = 0;
            int HydrogensToFind = 0;
            for (int index_curr_res = StartIndex; index_curr_res <= EndIndex; index_curr_res++) {
                string currentAtomName = AtomNames[index_curr_res];
                List<string> split_res = currentAtomName.Split(" ").ToList();
                List<string> HydrogensofResidue = PosChargeProtons[split_res[0]];
                HydrogensToFind = HydrogensofResidue.Count;
                if (HydrogensofResidue.Contains(split_res[1])) {
                    HydrogensFound++;
                }
            }
            if (HydrogensFound >= HydrogensToFind) {
                return true;
            }
            return false;
        }

        public List<List<string>> GenerateSaltBridgesProtein(PDBContent Content, bool UniqueConnections, float DistThreshold) {//Also known as ionic bond...
            //Method to find salt bridges: https://www.ks.uiuc.edu/Research/vmd/plugins/saltbr/ 
            
            PDBParser PDBPars = new PDBParser();
            /*Look at Residue contacts from the other Fragement (but not the residues being sequencially Steps next to the other Fragment)*/

            Dictionary<List<string>, float> PossibleSaltBridges = new Dictionary<List<string>, float>();//It might be possible that multiple salt bridges could be formed with different partners, so one needs to look at all potential ones and choose the one with the lowest distance

            List<string> SaltBridgePartnerOne = new List<string>() { "D", "E", "C", "Y" };// aspartic acid & glutamic acid & Cysteine (neg charge at high pH)
            List<string> SaltBridgePartnerTwo = new List<string>() { "K", "R", "H" };//Lysine & Arginine & histidine (Here it depends on the pKa...)


            List<string> NegativeChargedResidueAtoms = new List<string>() { "ASP OD1", "ASP OD2",
                                                                            "GLU OE1", "GLU OE2",
                                                                            "CYS SG",
                                                                            "TYR OH"};//Check if HD1, HD2 or HE1 and HE2 exist

            List<string> PositiveChargedResidueAtoms = new List<string>() { "LYS NZ",
                                                                            "ARG NH1", "ARG NH2",
                                                                            "HIS ND1", "HIS NE2"};//HIS: if only one: HD1 or  HE2 is present, the charge is lost and it can't be used for a salt bridge

            //Check wether they have a H at their respective Residues, if not then that part is not charged...

            float CutoffValue = DistThreshold;

            List<string> Sequence = Content.SingleLetterSequence;
            int MaxIndex = Sequence.Count;
            List<List<int>> ResContacts = Content.ResidueContacts;
            List<int> SplitAtSite = Content.SplitAtSite;
            List<Vector3> AtomPositions = Content.AtomPositions;
            List<string> AtomNames = Content.AtomNames;

            for (int curr_index = 0; curr_index < MaxIndex; curr_index++) {
                string ResidueCurrent = Sequence[curr_index];
                int CurrentResidueType = 0;
                if (SaltBridgePartnerOne.Contains(ResidueCurrent)) {
                    CurrentResidueType = 1;
                } else if (SaltBridgePartnerTwo.Contains(ResidueCurrent)) {
                    CurrentResidueType = 2;
                }

                if (CurrentResidueType != 0) {
                    foreach (int Contact in ResContacts[curr_index]) {//look at each contact of the current residue
                        if (curr_index < Contact || UniqueConnections == false) {
                            string ResidueContact = Sequence[Contact];
                            int ContactResidueType = 0;
                            if (SaltBridgePartnerOne.Contains(ResidueContact)) {
                                ContactResidueType = 1;
                            } else if (SaltBridgePartnerTwo.Contains(ResidueContact)) {
                                ContactResidueType = 2;
                            }
                            int IsPotentialSaltbridge = ContactResidueType + CurrentResidueType;
                            if (IsPotentialSaltbridge == 3) {//Potential salt bridge
                                int StartLineIndexCurrRes = 0;
                                int StartLineIndexContact = 0;
                                if (curr_index != 0) {
                                    StartLineIndexCurrRes = SplitAtSite[curr_index - 1] - 1;
                                }
                                if (Contact != 0) {
                                    StartLineIndexContact = SplitAtSite[Contact - 1] - 1;
                                }

                                int EndLineIndexCurrRes = SplitAtSite[curr_index] - 2;
                                int EndLineIndexContact = SplitAtSite[Contact] - 2;


                                Dictionary<int, Vector3> NegativelyChargedAtoms = new Dictionary<int, Vector3>();
                                Dictionary<int, Vector3> PositivelyChargedAtoms = new Dictionary<int, Vector3>();

                                for (int index_curr_res = StartLineIndexCurrRes; index_curr_res <= EndLineIndexCurrRes; index_curr_res++) {
                                    string currentAtomName = AtomNames[index_curr_res];
                                    if (NegativeChargedResidueAtoms.Contains(currentAtomName) && !NegativelyChargedAtoms.ContainsKey(index_curr_res)) {
                                        if (CheckNegativeChargedResProtonated(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames) == false) {//check whether HD1, HD2 or HE1 and HE2 exists -> No negative charge!
                                            NegativelyChargedAtoms.Add(index_curr_res, AtomPositions[index_curr_res]);
                                        }
                                    } else if (PositiveChargedResidueAtoms.Contains(currentAtomName) && !PositivelyChargedAtoms.ContainsKey(index_curr_res)) {
                                        if (CheckPositiveChargedResProtonated(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames) == true) {
                                            PositivelyChargedAtoms.Add(index_curr_res, AtomPositions[index_curr_res]);
                                        }
                                    }
                                }

                                for (int index_contact_res = StartLineIndexContact; index_contact_res <= EndLineIndexContact; index_contact_res++) {
                                    string ContactAtomName = AtomNames[index_contact_res];
                                    if (NegativeChargedResidueAtoms.Contains(ContactAtomName) && !NegativelyChargedAtoms.ContainsKey(index_contact_res)) {
                                        if (CheckNegativeChargedResProtonated(StartLineIndexContact, EndLineIndexContact, AtomNames) == false) {
                                            NegativelyChargedAtoms.Add(index_contact_res, AtomPositions[index_contact_res]);
                                        }
                                    } else if (PositiveChargedResidueAtoms.Contains(ContactAtomName) && !PositivelyChargedAtoms.ContainsKey(index_contact_res)) {
                                        if (CheckPositiveChargedResProtonated(StartLineIndexContact, EndLineIndexContact, AtomNames) == true) {
                                            PositivelyChargedAtoms.Add(index_contact_res, AtomPositions[index_contact_res]);
                                        }
                                    }
                                }

                                //Now look through the distances between all the possible Negative to Positive Atom Combinations
                                float BestDistance = 10000f;
                                List<string> AtomsPartners = new List<string>() { "CA", "CA" };
                                foreach (KeyValuePair<int, Vector3> NegativeAtom in NegativelyChargedAtoms) {
                                    foreach (KeyValuePair<int, Vector3> PositiveAtom in PositivelyChargedAtoms) {
                                        float Distance = Vector3.Distance(NegativeAtom.Value, PositiveAtom.Value);
                                        if (Distance < BestDistance) {
                                            string NegAtomName = AtomNames[NegativeAtom.Key];
                                            List<string> split_res_NegativeAtom = NegAtomName.Split(" ").ToList();
                                            AtomsPartners[0] = split_res_NegativeAtom[1];

                                            string PositiveAtomName = AtomNames[PositiveAtom.Key];
                                            List<string> split_res_PositiveAtom = PositiveAtomName.Split(" ").ToList();
                                            AtomsPartners[1] = split_res_PositiveAtom[1];

                                            BestDistance = Distance;
                                        }
                                    }
                                }

                                if (BestDistance <= CutoffValue) {
                                    //It's a salt bridge!!!!
                                    List<string> SaltBridgePartners = new List<string>();
                                    //negative one first
                                    if (CurrentResidueType == 1) {
                                        //curr_res is negatively charged
                                        SaltBridgePartners.Add(curr_index + "." + AtomsPartners[0]);
                                        SaltBridgePartners.Add(Contact + "." + AtomsPartners[1]);
                                    } else {
                                        SaltBridgePartners.Add(Contact + "." + AtomsPartners[0]);
                                        SaltBridgePartners.Add(curr_index + "." + AtomsPartners[1]);
                                    }
                                    PossibleSaltBridges.Add(SaltBridgePartners, BestDistance);
                                }
                            }
                        }
                    }
                }
            }
            //Now we have all the possible salt bridges in PossibleSaltBridges now we order them by BestDistance:
            var orderedPossibleSBridges = PossibleSaltBridges.OrderBy(x => x.Value);
            List<List<string>> SaltBridges = new List<List<string>>();
            foreach (KeyValuePair<List<string>, float> SBridge in orderedPossibleSBridges) {
                string IndexRes1 = SBridge.Key[0];
                string IndexRes2 = SBridge.Key[1];
                //Console.WriteLine("Find sal: " + IndexRes1 + " and " + IndexRes2);
                List<string> SaltBridgePartners = new List<string>();
                SaltBridgePartners.Add(IndexRes1);
                SaltBridgePartners.Add(IndexRes2);
                List<string> SaltBridgePartnersInverse = new List<string>();
                SaltBridgePartnersInverse.Add(IndexRes2);
                SaltBridgePartnersInverse.Add(IndexRes1);
                if (!SaltBridges.Contains(SaltBridgePartners) && !SaltBridges.Contains(SaltBridgePartnersInverse)) {
                    SaltBridges.Add(SaltBridgePartners);
                }
            }
            return SaltBridges;
        }
    }
}
