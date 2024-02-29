using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.IO;
using System.Linq;

namespace SplitProteinPrediction {

    /*This implementation based on the C implementation found in freesasa, see: https://github.com/mittinatten/freesasa/blob/master/src/sasa_lr.c */
    class Lee_Richards_SASA {

        public PDBContent GetAtomASA(PDBContent PDBCont) {
            AA_Values AAVals = new AA_Values();
            ASA_Functions ASA_Functions = new ASA_Functions();
            //List<Vector3> ASA_SpherePoints = SphereFreeSASAFunction(20);
            int ns = 20;//slices per atom
            float probe = 1.4f;//The Size of the probe is usually 1.4 A

            //Get the important stuff from the pdb
            List<Vector3> AtomPositions = PDBCont.AtomPositions;
            List<float> RadiiList = PDBCont.RadiiList;
            List<int> SplitAtSite = PDBCont.SplitAtSite;
            List<string> Sequence = PDBCont.SingleLetterSequence;
            int AtomCount = AtomPositions.Count();

            float Distance = 2.5f * (probe + RadiiList.Max());
            List<List<int>> AtomBondList = ASA_Functions.AdjacentAtomList(AtomPositions, Distance);
            List<float> Areas = new List<float>();
            //Go through each atom:
            int a_index = 0;
            //Stuff for the residue ASA:
            int Residue_index = 0;
            float ResidueArea = 0f;


            float sasa = 0f;
            foreach (Vector3 AtomPos in AtomPositions) {
                Vector3 AtomSameZ = new Vector3(AtomPos.X, AtomPos.Y, 0f);
                float zi = AtomPos.Z;
                float Ri = RadiiList[a_index] + probe;
                //ASA_SpherePointsCount
                float delta = 2 * Ri / ns;//Thickness
                float z = zi - Ri - 0.5f * delta;
                float atomASA = 0f;

                //Splices durchgehen:
                for (int islice = 0; islice < ns; islice++) {
                    z += delta;
                    float di = MathF.Abs(zi - z);
                    float Ri_prime2 = Ri * Ri - di * di;
                    if (Ri_prime2 < 0) continue;
                    float Ri_prime = MathF.Sqrt(Ri_prime2);
                    if (Ri_prime <= 0) continue;
                    int is_buried = 0;

                    List<Interval> arc = new List<Interval>();

                    //Alle nachbarAtome durchgehen
                    foreach (int NeighbourAtomIndex in AtomBondList[a_index]) {
                        //Hier die Daten der nachbar Atome:
                        Vector3 NeighbourAtomPos = AtomPositions[NeighbourAtomIndex];
                        float zj = NeighbourAtomPos.Z;
                        float dj = MathF.Abs(zj - z);
                        float Rj = RadiiList[NeighbourAtomIndex] + probe;

                        if (dj < Rj) {
                            float Rj_prime2 = Rj * Rj - dj * dj;
                            float Rj_prime = MathF.Sqrt(Rj_prime2);
                            Vector3 NeigbourSameZ = new Vector3(NeighbourAtomPos.X, NeighbourAtomPos.Y, 0f);
                            float dij = Vector3.Distance(AtomSameZ, NeigbourSameZ);
                            if (dij >= Ri_prime + Rj_prime) {
                                continue;
                            }
                            if (dij + Ri_prime < Rj_prime) {
                                is_buried = 1;
                                break;
                            }
                            if (dij + Rj_prime < Ri_prime) {
                                continue;
                            }
                            if(dij == 0) {//Special case which caused major errors, this is to my knowledge not in freesasa
                                is_buried = 1;
                                break;
                            }

                            float alpha = MathF.Acos((Ri_prime2 + dij * dij - Rj_prime2) / (2f * Ri_prime * dij));

                            float xDist = (NeighbourAtomPos.X - AtomPos.X);
                            float yDist = (NeighbourAtomPos.Y - AtomPos.Y);

                            float beta = MathF.Atan2(yDist, xDist) + MathF.PI;

                            float inf = beta - alpha;
                            float sup = beta + alpha;

                            if (inf < 0) inf += 2f * MathF.PI;
                            if (sup > 2 * MathF.PI) sup -= 2f * MathF.PI;
                            
                            if (sup < inf) {
                                Interval Arc_ = new Interval() { Start = 0f, End = sup };
                                Interval Arc__ = new Interval() { Start = inf, End = 2f * MathF.PI };
                                arc.Add(Arc_);
                                arc.Add(Arc__);
                            } else {
                                Interval Arc_ = new Interval() { Start = inf, End = sup };
                                
                                arc.Add(Arc_);
                            }

                        }
                    }


                    //problem: one of the arcs end is NaN
                    if (is_buried == 0) {
                        float SliceASA = delta * Ri * ExposedArcs(arc);
                        atomASA += SliceASA;
                        sasa += SliceASA;
                    }
                }
                PDBCont.AtomAreas.Add(atomASA);
                ResidueArea += atomASA;
                a_index++;

                //Add the relatve ASA (For the NIS)
                if (Residue_index < SplitAtSite.Count()) {
                    //We want to include the last atom, so we need to check for that too
                    if ((a_index == SplitAtSite[Residue_index] - 1 && SplitAtSite[Residue_index] != AtomCount)) {
                        float MaxASA = AAVals.ASA_MaxResidue[Sequence[Residue_index]];
                        float RelASA = ResidueArea / MaxASA;
                        PDBCont.RelativeResidueASA.Add(RelASA);
                        ResidueArea = 0;
                        Residue_index++;
                    }
                }
            }
            return PDBCont;
        }

        public class Interval {
            public float Start;
            public float End;
        }


        /*This function originates from https://stackoverflow.com/questions/11480031/merging-overlapping-time-intervals, credit goes to Andre Calil and Rexx Magnus*/
        public static List<Interval> Merge(List<Interval> intervals) {

            var mergedIntervals = new List<Interval>();
            var orderedIntervals = intervals.OrderBy<Interval, float>(x => x.Start).ToList<Interval>();

            float start = orderedIntervals.First().Start;
            float end = orderedIntervals.First().End;

            Interval currentInterval;
            for (int i = 1; i < orderedIntervals.Count; i++) {
                currentInterval = orderedIntervals[i];
                if (currentInterval.Start < end) {
                    end = (end > currentInterval.End ? end : currentInterval.End);
                } else {
                    mergedIntervals.Add(new Interval() {
                        Start = start,
                        End = end
                    });

                    start = currentInterval.Start;
                    end = currentInterval.End;
                }
            }

            mergedIntervals.Add(new Interval() {
                Start = start,
                End = end
            });

            return mergedIntervals;
        }

        public float ExposedArcs(List<Interval> Arc) {
            float freeRange = 2f * MathF.PI;
            if (Arc.Count() != 0) {
                
                List<Interval> MergedIntervals = Merge(Arc);
                foreach (Interval SingleInterval in MergedIntervals) {

                    float Start = SingleInterval.Start;//plus
                    float End = SingleInterval.End;//minus

                    float Range = End - Start;
                    freeRange -= Range;
                }
            }
            return freeRange;
        }
        public List<double> sort_arcs(List<double> arc, int n) {
            double[] tmp = new double[2];
            int arcj;
            double end = arc.Count() + (2 * n);
            for (int arci = arc.Count() + 2; arci < end; arci++) {
                tmp[0] = arc[arci];
                tmp[1] = arc[arci + 1];
                arcj = arci;
                while ((arcj > arc.Count()) && (arcj - 2 > tmp[0])) {
                    arcj = (int)arc[arcj - 2];
                    arc[arcj + 1] = (int)arc[arcj - 1];
                    arcj -= 2;
                }
                arc[arcj] = tmp[0];
                arc[arcj + 1] = tmp[1];
            }
            return arc;
        }

        public double exposed_arc_length(List<double> arc, int n) {
            int i2;
            double sum;
            double sup;
            double tmp;
            if (n == 0) return 2 * MathF.PI;
            arc = sort_arcs(arc, n);
            sum = arc[0];
            sup = arc[1];
            /* in the following it is assumed that the arc[i2] <= arc[i2+1] */
            for (i2 = 2; i2 < 2 * n; i2 += 2) {
                if (sup < arc[i2]) sum += arc[i2] - sup;
                tmp = arc[i2 + 1];
                if (tmp > sup) sup = tmp;
            }
            return sum + 2 * MathF.PI - sup;
        }

        public List<Vector3> SphereFreeSASAFunction(int n) {
            List<Vector3> Points = new List<Vector3>();
            float dlong = MathF.PI * (3 - MathF.Sqrt(5));
            float dz = 2f / (float)n;
            float longe = 0f;
            float z = 1 - dz / 2;
            for (int i = 0; i < n; i++) {
                float r = MathF.Sqrt(1 - z * z);
                Vector3 PointVector = new Vector3(MathF.Cos(longe) * r, MathF.Sin(longe) * r, z);
                z -= dz;
                longe += dlong;
                Points.Add(PointVector);
            }
            return Points;
        }
    }
}
