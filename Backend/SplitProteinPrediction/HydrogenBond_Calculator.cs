using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;

namespace SplitProteinPrediction {
    class HydrogenBond_Calculator {

        public List<Vector3> GetHydrogenPositionsFromDonor(int StartIndex, int EndIndex, List<string> AtomNames, List<Vector3> AtomPositions, string DonorResidueExtension) {
            //Gets the positions of all the Hydrogens linked to a specific donor (filtered via the atom id)
            List<Vector3> output = new List<Vector3>();
            string LookFor = "H" + DonorResidueExtension;
            List<string> LookForVarriations = new List<string>();
            for (int i = 0; i <= 5; i++) {
                if (i == 0) {
                    LookForVarriations.Add(LookFor);
                } else {
                    LookForVarriations.Add(LookFor + i.ToString());
                }
            }
            for (int index = StartIndex; index <= EndIndex; index++) {
                string currentAtomName = AtomNames[index];
                List<string> check_type = currentAtomName.Split(" ").ToList();
                if (LookForVarriations.Contains(check_type[1])) {
                    output.Add(AtomPositions[index]);
                }
            }
            return output;
        }

        public List<int> GetHydrogenIdsFromDonor(int StartIndex, int EndIndex, List<string> AtomNames, List<Vector3> AtomPositions, string DonorResidueExtension) {
            //Gets the positions of all the Hydrogens linked to a specific donor (filtered via the atom id)
            List<int> output = new List<int>();
            string LookFor = "H" + DonorResidueExtension;
            List<string> LookForVarriations = new List<string>();
            for (int i = 0; i <= 5; i++) {
                if (i == 0) {
                    LookForVarriations.Add(LookFor);
                } else {
                    LookForVarriations.Add(LookFor + i.ToString());
                }
            }
            for (int index = StartIndex; index <= EndIndex; index++) {
                string currentAtomName = AtomNames[index];
                List<string> check_type = currentAtomName.Split(" ").ToList();
                if (LookForVarriations.Contains(check_type[1])) {
                    output.Add(index);
                }
            }
            return output;
        }


        public int CountResidueHydrogens(int StartIndex, int EndIndex, List<string> AtomNames, string ResidueExtension) {
            int countHs = 0;
            string LookFor = "H" + ResidueExtension;
            List<string> LookForVarriations = new List<string>();
            for (int i = 0; i <= 4; i++) {
                if (i == 0) {
                    LookForVarriations.Add(LookFor);
                } else {
                    LookForVarriations.Add(LookFor + i.ToString());
                }
            }
            for (int index = StartIndex; index <= EndIndex; index++) {
                string currentAtomName = AtomNames[index];
                List<string> check_type = currentAtomName.Split(" ").ToList();
                if (LookForVarriations.Contains(check_type[1])) {
                    countHs++;
                }
            }
            return countHs;
        }

        class DonorAcceptorConnection {
            public int Acceptor { get; set; }
            public int Donor { get; set; }
            public int DonorHydrogenId { get; set; }
            public float HydrogenDonorDistance { get; set; }
            public float DistanceHydrogenAcceptor { get; set; }
            public float Angle { get; set; }
            public float Distance { get; set; }
        }


        public List<List<string>> GenerateHBonds(PDBContent Content, bool UniqueConnections, float Degthreshold, float Distancethreshold) {
            /*Look at Residue contacts from the other Fragement (but not the residues being sequencially Steps next to the other Fragment)*/
            List<List<string>> ResiduesHyrogenBonds = new List<List<string>>();
            List<float> RadiiList = Content.RadiiList;
            List<int> AtmToResIndex = Content.AtomIndexToResidueIndex;


            List<string> PotentialAcceptorDonorAtoms = new List<string> { "N", "O", "S", "C" };

            float D_A_DistThreshold = Distancethreshold;//3 Angs
            float D_H_A_MinDeg = Degthreshold;

            List<string> Sequence = Content.SingleLetterSequence;
            int MaxIndex = Sequence.Count;
            List<List<int>> ResContacts = Content.ResidueContacts;
            List<int> SplitAtSite = Content.SplitAtSite;
            List<Vector3> AtomPositions = Content.AtomPositions;
            List<string> AtomNames = Content.AtomNames;

            Dictionary<DonorAcceptorConnection, float> AcceptorDonorConnections = new Dictionary<DonorAcceptorConnection, float>();

            for (int curr_index = 0; curr_index < MaxIndex; curr_index++) {
                List<int> CurrResHBonds = new List<int>();

                int StartLineIndexCurrRes = 0;
                if (curr_index != 0) {
                    StartLineIndexCurrRes = SplitAtSite[curr_index - 1] - 1;
                }
                int EndLineIndexCurrRes = SplitAtSite[curr_index] - 2;

                //Get Atom List of all potential donors & acceptors in the current residue
                Dictionary<int, Vector3> DonorsCurrRes = new Dictionary<int, Vector3>();
                Dictionary<int, List<Vector3>> HydrogensFromDonorsCurrRes = new Dictionary<int, List<Vector3>>();
                Dictionary<int, List<int>> HydrogenIdsDonorCurrRes = new Dictionary<int, List<int>>();
                Dictionary<int, Vector3> AcceptorsCurrRes = new Dictionary<int, Vector3>();

                for (int index_curr_res = StartLineIndexCurrRes; index_curr_res <= EndLineIndexCurrRes; index_curr_res++) {
                    string currentAtomName = AtomNames[index_curr_res];
                    List<string> check_type = currentAtomName.Split(" ").ToList();
                    string AtomType = check_type[1].Substring(0, 1);//First Character
                    string ResExtension = check_type[1].Remove(0, 1);//Delete firt character (N/S/O) to get the extenison
                    if (PotentialAcceptorDonorAtoms.Contains(AtomType)) {
                        //Donor: Get The Hydrogen Positions (if none it can't be a donor)
                        List<Vector3> HydrogenPositions = GetHydrogenPositionsFromDonor(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames, AtomPositions, ResExtension);
                        List<int> HydrogenIds = GetHydrogenIdsFromDonor(StartLineIndexCurrRes, EndLineIndexCurrRes, AtomNames, AtomPositions, ResExtension);
                        if (HydrogenPositions.Count() != 0) {//There's a H so it can be a donor, this should actually never be the case because of CanCurrResBeDonor but lets be super save
                            DonorsCurrRes.Add(index_curr_res, AtomPositions[index_curr_res]);
                            HydrogensFromDonorsCurrRes.Add(index_curr_res, HydrogenPositions);
                            HydrogenIdsDonorCurrRes.Add(index_curr_res, HydrogenIds);
                        }

                        //Acceptor (Can always be an acceptor):
                        string LookPosAcceptorString = currentAtomName;
                        if (check_type[1] == "O" || check_type[1] == "OC1" || check_type[1] == "OC2" || check_type[1] == "OXT") {
                            LookPosAcceptorString = "O";
                        }

                        bool exclude_causeNH3 = false;
                        if (HydrogenPositions.Count() == 3 && AtomType == "N") {
                            exclude_causeNH3 = true;
                        }
                        if (AtomType != "C" && !exclude_causeNH3) {
                            AcceptorsCurrRes.Add(index_curr_res, AtomPositions[index_curr_res]);
                        }

                    }
                }

                foreach (int Contact in ResContacts[curr_index]) {
                    if (curr_index < Contact || UniqueConnections == false) {//We only want the contacts once not serveral times, from the persprective of multiple residues...
                        int StartLineIndexContact = 0;
                        if (Contact != 0) {
                            StartLineIndexContact = SplitAtSite[Contact - 1] - 1;
                        }
                        int EndLineIndexContact = SplitAtSite[Contact] - 2;
                        //Get Atom List of all potential donors & acceptors for the neighbouring residue
                        Dictionary<int, Vector3> DonorsContactRes = new Dictionary<int, Vector3>();
                        Dictionary<int, List<Vector3>> HydrogensFromDonorsContactRes = new Dictionary<int, List<Vector3>>();
                        Dictionary<int, Vector3> AcceptorsContactRes = new Dictionary<int, Vector3>();
                        Dictionary<int, List<int>> HydrogenIdsDonorContactRes = new Dictionary<int, List<int>>();

                        for (int index_contact_res = StartLineIndexContact; index_contact_res <= EndLineIndexContact; index_contact_res++) {
                            string ContactAtomName = AtomNames[index_contact_res];
                            List<string> check_type = ContactAtomName.Split(" ").ToList();
                            string AtomType = check_type[1].Substring(0, 1);//First Character
                            string ResExtension = check_type[1].Remove(0, 1);//Delete firt character (N/S/O) to get the extenison
                            if (PotentialAcceptorDonorAtoms.Contains(AtomType)) {
                                //Donor: Get The Hydrogen Positions (if none it can't be a donor)
                                List<Vector3> HydrogenPositions = GetHydrogenPositionsFromDonor(StartLineIndexContact, EndLineIndexContact, AtomNames, AtomPositions, ResExtension);
                                List<int> HydrogenIds = GetHydrogenIdsFromDonor(StartLineIndexContact, EndLineIndexContact, AtomNames, AtomPositions, ResExtension);
                                if (HydrogenPositions.Count() != 0) {//There's a H so it can be a donor, this should actually never be the case because of CanCurrResBeDonor but lets be super save
                                    DonorsContactRes.Add(index_contact_res, AtomPositions[index_contact_res]);
                                    HydrogensFromDonorsContactRes.Add(index_contact_res, HydrogenPositions);
                                    HydrogenIdsDonorContactRes.Add(index_contact_res, HydrogenIds);
                                }
                                //Acceptor (Can always be an acceptor):
                                string LookPosAcceptorString = ContactAtomName;
                                if (check_type[1] == "O" || check_type[1] == "OC1" || check_type[1] == "OC2" || check_type[1] == "OXT") {
                                    LookPosAcceptorString = "O";
                                }
                                bool exclude_causeNH3 = false;
                                if (HydrogenPositions.Count() == 3 && AtomType == "N") {
                                    exclude_causeNH3 = true;
                                }
                                if (AtomType != "C" && !exclude_causeNH3) {//C cannot be an acceptor and also NH3 cannot be one
                                    AcceptorsContactRes.Add(index_contact_res, AtomPositions[index_contact_res]);
                                }
                            }
                        }


                        //Check where if we have HBonds between the two residues (Multiple possible)
                        foreach (KeyValuePair<int, Vector3> CurrentResDonor in DonorsCurrRes) {
                            Vector3 CurrResDonorPos = CurrentResDonor.Value;
                            int CurrResAtomIndex = CurrentResDonor.Key;
                            List<Vector3> HydrogenPositions = HydrogensFromDonorsCurrRes[CurrResAtomIndex];
                            foreach (KeyValuePair<int, Vector3> ContactResAcceptor in AcceptorsContactRes) {
                                Vector3 ContactResAcceptorPos = ContactResAcceptor.Value;
                                int ContactResAtomIndex = ContactResAcceptor.Key;
                                float Dist = Vector3.Distance(CurrResDonorPos, ContactResAcceptorPos);
                                double deg = 0;
                                if (Dist <= D_A_DistThreshold) {
                                    //First condition met:
                                    int hydrogenid_index = 0;
                                    foreach (Vector3 HPos in HydrogenPositions) {
                                        float DistHAcceptor = Vector3.Distance(HPos, ContactResAcceptorPos);
                                        Vector3 Vector_Donor_H = HPos - CurrResDonorPos;
                                        Vector3 Vector_h_Acceptor = ContactResAcceptorPos - HPos;
                                        float DotProd = Vector3.Dot(Vector_Donor_H, Vector_h_Acceptor);
                                        double radians = Math.Acos(DotProd / (Vector_Donor_H.Length() * Vector_h_Acceptor.Length()));
                                        deg = (180f / Math.PI) * radians;
                                        if (deg <= D_H_A_MinDeg) {
                                            int GetDonorHydrogenIndex = HydrogenIdsDonorCurrRes[CurrResAtomIndex][hydrogenid_index];
                                            float DonorDist = Vector3.Distance(CurrResDonorPos, HPos);
                                            DonorAcceptorConnection DonorAcceptorContainer = new DonorAcceptorConnection {
                                                Acceptor = ContactResAtomIndex,
                                                Donor = CurrResAtomIndex,
                                                DonorHydrogenId = GetDonorHydrogenIndex,
                                                HydrogenDonorDistance = DonorDist,
                                                DistanceHydrogenAcceptor = DistHAcceptor,
                                                Distance = Dist,
                                                Angle = (float)deg
                                            };
                                            AcceptorDonorConnections.Add(DonorAcceptorContainer, Dist);
                                        }
                                        hydrogenid_index++;
                                    }
                                }

                            }
                        }

                        //Check where we have HBonds between the two residues (Multiple possible)
                        foreach (KeyValuePair<int, Vector3> CurrentResAcceptor in AcceptorsCurrRes) {
                            Vector3 CurrResAcceptorPos = CurrentResAcceptor.Value;
                            int CurrResAtomIndex = CurrentResAcceptor.Key;
                            foreach (KeyValuePair<int, Vector3> ContactResDonor in DonorsContactRes) {
                                Vector3 ContactResDonorPos = ContactResDonor.Value;
                                int ContactResAtomIndex = ContactResDonor.Key;
                                //get all the hydrogen positions
                                List<Vector3> HydrogenPositions = HydrogensFromDonorsContactRes[ContactResAtomIndex];
                                float Dist = Vector3.Distance(CurrResAcceptorPos, ContactResDonorPos);
                                double deg = 0;
                                if (Dist <= D_A_DistThreshold) {
                                    //First condition met:
                                    int hydrogenid_index = 0;
                                    foreach (Vector3 HPos in HydrogenPositions) {
                                        float DistHAcceptor = Vector3.Distance(HPos, CurrResAcceptorPos);
                                        Vector3 Vector_Donor_H = HPos - ContactResDonorPos;
                                        Vector3 Vector_h_Acceptor = CurrResAcceptorPos - HPos;
                                        float DotProd = Vector3.Dot(Vector_Donor_H, Vector_h_Acceptor);
                                        double radians = Math.Acos(DotProd / (Vector_Donor_H.Length() * Vector_h_Acceptor.Length()));
                                        deg = (180f / Math.PI) * radians;
                                        if (deg <= D_H_A_MinDeg) {
                                            int GetDonorHydrogenIndex = HydrogenIdsDonorContactRes[ContactResAtomIndex][hydrogenid_index];
                                            float DonorDist = Vector3.Distance(ContactResDonorPos, HPos);
                                            DonorAcceptorConnection DonorAcceptorContainer = new DonorAcceptorConnection {
                                                Acceptor = CurrResAtomIndex,
                                                Donor = ContactResAtomIndex,
                                                DonorHydrogenId = GetDonorHydrogenIndex,
                                                HydrogenDonorDistance = DonorDist,
                                                DistanceHydrogenAcceptor = DistHAcceptor,
                                                Distance = Dist,
                                                Angle = (float)deg
                                            };
                                            AcceptorDonorConnections.Add(DonorAcceptorContainer, Dist);

                                        }
                                        hydrogenid_index++;
                                    }
                                }

                            }
                        }
                    }
                }
            }

            //now we have all our potential donor-acceptor pairs with distances
            List<int> AcceptorsFinalResult = new List<int>();
            List<int> DonorsFinalResult = new List<int>();


            List<List<string>> HBonds = new List<List<string>>();
            var OrderedAD = AcceptorDonorConnections.OrderBy(pair => pair.Value);
            foreach (KeyValuePair<DonorAcceptorConnection, float> SingleEntry in OrderedAD) {
                DonorAcceptorConnection InfoAD = SingleEntry.Key;

                int AcceptorAtomIndex = InfoAD.Acceptor;
                float DistanceReal = InfoAD.Distance;
                float Angle = InfoAD.Angle;
                string AcceptorName = AtomNames[AcceptorAtomIndex];
                List<string> AcceptorSplitName = AcceptorName.Split(" ").ToList();
                string AcceptorNameFinal = AcceptorSplitName[1];

                int DonorAtomIndex = InfoAD.Donor;
                int HydrogenDonorAtomIndex = InfoAD.DonorHydrogenId;
                string HydrogenAtomName = AtomNames[HydrogenDonorAtomIndex];
                List<string> HydrogenSplitName = HydrogenAtomName.Split(" ").ToList();
                string HydrogenName = HydrogenSplitName[1];

                int AcceptorResidueIndex = AtmToResIndex[AcceptorAtomIndex];
                int DonorResidueIndex = AtmToResIndex[DonorAtomIndex];
                List<string> Pair = new List<string> { AcceptorResidueIndex + "." + AcceptorNameFinal, DonorResidueIndex + "." + HydrogenName, DistanceReal.ToString(), Angle.ToString() };
                if (!HBonds.Contains(Pair)) {
                    HBonds.Add(Pair);
                    AcceptorsFinalResult.Add(AcceptorResidueIndex);
                    DonorsFinalResult.Add(DonorResidueIndex);
                }
            }
            return HBonds;//residue to residue bridge
        }


    }
}
