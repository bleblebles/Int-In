using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.IO;
using System.Linq;

namespace SplitProteinPrediction {

    /*This implementation based on the C implementation found in freesasa, see:https://github.com/mittinatten/freesasa/blob/master/src/sasa_sr.c */

    class Shrake_Rupeley_SASA {
        public PDBContent GetAtomASA(PDBContent PDBCont, List<List<int>> AtomBondsPreCalc , bool RegenerateAdj = true, bool ProtA = true) {
            AA_Values AAVals = new AA_Values();
            ASA_Functions ASA_Functions = new ASA_Functions();
            List<Vector3> ASA_SpherePoints = ASA_Functions.GenerateSpherePoints(100);
            int ASA_SpherePointsCount = ASA_SpherePoints.Count;
            float probe = 1.4f;//The Size of the probe is usually 1.4 A
            float cons = 4f * MathF.PI / ASA_SpherePointsCount;//SASA A^2 per Point on the spere (r^2 = 1 in this case...)

            //Get the important stuff from the pdb
            List<Vector3> AtomPositions = PDBCont.AtomPositions;
            List<float> RadiiList = PDBCont.RadiiList;
            List<int> SplitAtSite = PDBCont.SplitAtSite;
            List<string> Sequence = PDBCont.SingleLetterSequence;
            int AtomCount = AtomPositions.Count();

            float BoxMaxLen = 10 * (probe + RadiiList.Max());//Size of the Box where Atoms are summerized in


            //Here: lead that at the beginning:
            List<List<int>> AtomBondList = new List<List<int>>();
            if (RegenerateAdj == true) {
               AtomBondList = ASA_Functions.AdjacentAtomList(AtomPositions, BoxMaxLen);
            } else {
                AtomBondsPreCalc = ASA_Functions.DeleteSplitAtomIndexes(AtomBondsPreCalc, AtomCount, ProtA);
                AtomBondList = AtomBondsPreCalc;
            }

            List<float> Areas = new List<float>();
            //Go through each atom:
            int a_index = 0;
            //Stuff for the residue ASA:
            int Residue_index = 0;
            float ResidueArea = 0;
            //List<float> AreaResList = new List<float>();
            foreach (Vector3 AtomPos in AtomPositions) {
                List<int> BondedAtoms_of_Atom = AtomBondList[a_index];//Get all the adjacent Atoms of the specific Atom
                List<int> Neighbours = ASA_Functions.RefineNeighbourBonds(AtomPositions, RadiiList, BondedAtoms_of_Atom, probe, a_index);//Not really sure why I can't integrate that in the AdjacentAtomList
                int NumberNeighbours = Neighbours.Count();
                float radius = RadiiList[a_index] + probe;//Atom Radius+ProbeRadius
                int n_accessible_point = 0;
                //Go through the random points on the sphere
                foreach (Vector3 Point in ASA_SpherePoints) {
                    //Check if Point is accessible for the probe
                    bool is_accessible = true;
                    Vector3 TestPoint = Point * radius + AtomPos;//Point is normalized so multiplying it via th radius gives us the Point of the probe when it reaches the atomvdW radius
                    List<int> CycledIndices = new List<int>();

                    foreach (int neighbour_index in Neighbours) {
                        Vector3 Atom_pos_j = AtomPositions[neighbour_index];//Get Position of Neighbour
                        float r = RadiiList[neighbour_index] + probe;//Neighbour Radius + Probe Radius
                        float Distance = Vector3.Distance(Atom_pos_j, TestPoint);//Dist between the Accesspoint and our atom pos
                        if (Distance <= r) {//The Neigbouring Atom touches the probe (as the distance between probe and atom is smaller than the radius of the atom + radius of the probe)
                            is_accessible = false;
                            break;
                        }
                    }

                    if (is_accessible == true) n_accessible_point++;//Point is accessible!
                }
                float area = cons * n_accessible_point * radius * radius;//(4*pi)/points on sphere * nbr of free points on the atoms sphere * r^2
                //Console.WriteLine("Atom area:"+ area);
                PDBCont.AtomAreas.Add(area);
                ResidueArea += area;
                //AreaResList.Add(area);
                a_index++;
                if (Residue_index == 110) {
                    //Console.WriteLine("AtomIndex:" + a_index + "Area: " + area);
                }

                //Add the relatve ASA (For the NIS)
                if (Residue_index < SplitAtSite.Count()) {
                    //We want to include the last atom, so we need to check for that too
                    if ((a_index == SplitAtSite[Residue_index] - 1 && SplitAtSite[Residue_index] != AtomCount)) {
                        //Console.WriteLine("Residue ends at Atom number: {0}, atoms used: {1}", a_index, AreaResList.Count());//a_index = Anzahl der Atome die bis jetzt durchgenommen wurden
                        float MaxASA = AAVals.ASA_MaxResidue[Sequence[Residue_index]];
                        //Console.WriteLine("Res area:" + ResidueArea);
                        float RelASA = ResidueArea / MaxASA;
                        PDBCont.RelativeResidueASA.Add(RelASA);
                        ResidueArea = 0;
                        Residue_index++;
                        //AreaResList.Clear();
                    }
                }
            }
            return PDBCont;
        }

    }



}
