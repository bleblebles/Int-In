using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;

using System.Configuration;
namespace SplitProteinPrediction {
    class CFragDocking {

        public List<float> Get_CFragDocking(List<string> HydrogenBonds, List<string> SBridges, List<string> Arom, List<string> vDW_string, int SeqLen) {
            /*------------ HBonds: ------------*/
            int CFragDockingInt = int.Parse(ConfigurationManager.AppSettings.Get("CFragDocking_Residues"));

            List<int> ResultHBonds_x = new List<int>();
            List<int> ResultHBonds_y = new List<int>();

            List<float> ResultsEnergyHBonds = new List<float>();
            foreach (string HBond in HydrogenBonds) {
                List<string> NoDistance = HBond.Split("|").ToList();
                List<string> PartnersinBond_2 = NoDistance[0].Split("-").ToList();
                List<int> PartnersinBond_Indexes = (from i in PartnersinBond_2 select int.Parse(i.Split(".")[0])).ToList();
                List<string> type = (from i in PartnersinBond_2 select i.Split(".")[1][0].ToString()).ToList();
                //if type[0] != "C" and type[1] != "C":
                int PartnerA = PartnersinBond_Indexes[0];
                int PartnerB = PartnersinBond_Indexes[1];

                float Energy = 0f;
                float Distance = float.Parse(NoDistance[1]);
                float Angle = float.Parse(NoDistance[2]);
                if (Angle <= 63f && Distance <= 3.5) {//energy by ring 2.0
                    if (Distance <= 1.5) {
                        Energy = 115f;
                    } else if (Distance <= 2.2) {
                        Energy = 40f;
                    } else {
                        Energy = 17f;
                    }
                    ResultsEnergyHBonds.Add(Energy);
                    if (PartnerA < PartnerB) {
                        ResultHBonds_x.Add(PartnerA);
                        ResultHBonds_y.Add(PartnerB);
                    } else {
                        ResultHBonds_x.Add(PartnerB);
                        ResultHBonds_y.Add(PartnerA);
                    }
                }
            }
            /*------------ SaltBridges: ------------*/
            List<int> ResultSBridges_x = new List<int>();
            List<int> ResultSBridges_y = new List<int>();
            foreach (string SBridge in SBridges) {
                List<string> PartnersinBond_2 = SBridge.Split("-").ToList();
                List<int> PartnersinBond_Indexes = (from i in PartnersinBond_2 select int.Parse(i.Split(".")[0])).ToList();
                List<string> type = (from i in PartnersinBond_2 select i.Split(".")[1][0].ToString()).ToList();
                //if type[0] != "C" and type[1] != "C":
                int PartnerA = PartnersinBond_Indexes[0];
                int PartnerB = PartnersinBond_Indexes[1];
                if (PartnerA < PartnerB) {
                    ResultSBridges_x.Add(PartnerA);
                    ResultSBridges_y.Add(PartnerB);
                } else {
                    ResultSBridges_x.Add(PartnerB);
                    ResultSBridges_y.Add(PartnerA);

                }
            }
            /*------------ Aromatic bonds: ------------*/
            List<int> ResultAromBond_x = new List<int>();
            List<int> ResultAromBond_y = new List<int>();
            List<float> ResultsAromEnergy = new List<float>();
            foreach (string AromInt in Arom) {
                List<string> PartnersinBond_2 = AromInt.Split("-").ToList();
                List<int> PartnersinBond_Indexes = (from i in PartnersinBond_2 select int.Parse(i.Split(".")[0])).ToList();
                List<string> type = (from i in PartnersinBond_2 select i.Split(".")[1][0].ToString()).ToList();
                //if type[0] != "C" and type[1] != "C":
                int PartnerA = PartnersinBond_Indexes[0];
                int PartnerB = PartnersinBond_Indexes[1];

                float Energy = 9.6f;
                if (type[0] == "CA" && type[1] == "CA") {
                    Energy = 9.4f;
                }

                ResultsAromEnergy.Add(Energy);
                if (PartnerA < PartnerB) {
                    ResultAromBond_x.Add(PartnerA);
                    ResultAromBond_y.Add(PartnerB);
                } else {
                    ResultAromBond_x.Add(PartnerB);
                    ResultAromBond_y.Add(PartnerA);

                }
            }

            /*------------ van der waals: ------------*/
            List<int> vDW_x = new List<int>();
            List<int> vDW_y = new List<int>();
            foreach (string vDW in vDW_string) {
                List<string> PartnersinBond_2 = vDW.Split("-").ToList();
                List<int> PartnersinBond_Indexes = (from i in PartnersinBond_2 select int.Parse(i.Split(".")[0])).ToList();
                int PartnerA = PartnersinBond_Indexes[0];
                int PartnerB = PartnersinBond_Indexes[1];
                if (PartnerA < PartnerB) {
                    vDW_x.Add(PartnerA);
                    vDW_y.Add(PartnerB);
                } else {
                    vDW_x.Add(PartnerB);
                    vDW_y.Add(PartnerA);
                }
               
            }

            List<float> EnergiesForSplitSites = new List<float>();
            for (int split_site_index = 0; split_site_index < SeqLen-1; split_site_index++) {
                //check all the residues from before the split:
                float SumEnergySite = 0f;

                int start_before_idx = 0;
                int start_after_idx = split_site_index + CFragDockingInt;
                if (start_after_idx > SeqLen - 1) {
                    start_after_idx = SeqLen - 1;
                }

                foreach (int i in Enumerable.Range(start_before_idx, split_site_index+1)) {
                    if (ResultHBonds_x.Contains(i)) {
                        int index_list = ResultHBonds_x.IndexOf(i);
                        if (ResultHBonds_y[index_list] > split_site_index && ResultHBonds_y[index_list] <= start_after_idx) {
                            SumEnergySite += ResultsEnergyHBonds[index_list];
                            
                        }
                    }
                }
                foreach (var i in Enumerable.Range(start_before_idx, split_site_index+1)) {
                    if (ResultSBridges_x.Contains(i)) {
                        int index_list = ResultSBridges_x.IndexOf(i);
                        if (ResultSBridges_y[index_list] > split_site_index && ResultSBridges_y[index_list] <= start_after_idx) {
                            SumEnergySite += 20f;
                            
                        }
                    }
                }

                foreach (var i in Enumerable.Range(start_before_idx, split_site_index+1)) {
                    if (ResultAromBond_x.Contains(i)) {
                        int index_list = ResultAromBond_x.IndexOf(i);
                        if (ResultAromBond_y[index_list] > split_site_index && ResultAromBond_y[index_list] <= start_after_idx) {
                            SumEnergySite += ResultsAromEnergy[index_list];
                            
                        }
                    }
                }
                
                foreach (var i in Enumerable.Range(start_before_idx, split_site_index+1)) {
                    if (vDW_x.Contains(i)) {
                        int index_list = vDW_x.IndexOf(i);
                        if (vDW_y[index_list] > split_site_index && vDW_y[index_list] <= start_after_idx) {
                            SumEnergySite += 6f;
                        }
                    }
                }
                EnergiesForSplitSites.Add(SumEnergySite);
            }

            return EnergiesForSplitSites;
        }
    }
}
