using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.IO;
using System.Linq;

namespace SplitProteinPrediction {
    /*These functions are based upon the excellent work from Simon Mitternacht on freesasa: https://github.com/mittinatten/freesasa*/
    class ASA_Functions {

        public List<Vector3> GenerateSpherePoints(int n) {
            List<Vector3> Points = new List<Vector3>();
            float inc = MathF.PI * (3 - MathF.Sqrt(5));
            float offset = 2f / (float)n;
            for (int i = 0; i < n; i++) {
                float y = i * offset - 1 + (offset / 2);
                float r = MathF.Sqrt(1 - y * y);
                float phi = i * inc;
                Vector3 PointVector = new Vector3(MathF.Cos(phi) * r, y, MathF.Sin(phi) * r);
                Points.Add(PointVector);
            }
            return Points;
        }

        public class DefaultDictionary<TKey, TValue> : Dictionary<TKey, TValue> where TValue : new() {
            public new TValue this[TKey key] {
                get {
                    TValue val;
                    if (!TryGetValue(key, out val)) {
                        val = new TValue();
                        Add(key, val);
                    }
                    return val;
                }
                set { base[key] = value; }
            }
        }

        private Dictionary<Vector3, List<int>> MakeBoxes(List<Vector3> AtomPos, float d_max) {
            //Dict with keys as Positions of boxes with d_max length, values are Lists of the indices of the atoms
            int AtomsCount = AtomPos.Count();
            DefaultDictionary<Vector3, List<int>> BoxesCoords = new DefaultDictionary<Vector3, List<int>>();
            for (int i = 0; i < AtomsCount; i++) {
                Vector3 SingleAtomPos = AtomPos[i];
                Vector3 BoxCoord = new Vector3((int)MathF.Floor(SingleAtomPos.X / d_max), (int)MathF.Floor(SingleAtomPos.Y / d_max), (int)MathF.Floor(SingleAtomPos.Z / d_max));
                //Explanation if we'd only use the X values: All atoms of x < 1*d_max will be in one Box the next box will have x > d_max && x <2*d_max and so on...
                BoxesCoords[BoxCoord].Add(i);
            }
            return BoxesCoords;
        }

        private List<List<int>> AddBond(List<Vector3> AtomPos, int Atom1Index, int Atom2Index, List<List<int>> bonds, float d_max) {
            //List of Lists => First Indice = Atom1, second List gives all the atoms which are connected to him
            Vector3 AtomPos1 = AtomPos[Atom1Index];
            Vector3 AtomPos2 = AtomPos[Atom2Index];
            if (Vector3.Distance(AtomPos1, AtomPos2) <= d_max) {
                //Add the connection via the two sides:
                bonds[Atom1Index].Add(Atom2Index);
                bonds[Atom2Index].Add(Atom1Index);
            }
            return bonds;
        }

        public List<int> NeighbouringAtoms(Dictionary<Vector3, List<int>> Boxes, KeyValuePair<Vector3, List<int>> Box) {
            //Look @ half of the neighbouring boxes, the other half will be considered when we look @ the opposite box...
            List<int> NeighbourAtomList = new List<int>();
            //Current Box:
            float x = Box.Key.X;
            float y = Box.Key.Y;
            float z = Box.Key.Z;
            //Add all the surrounding 9 Boxes which have z+1 
            if (Boxes.ContainsKey(new Vector3(x + 1, y + 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y + 1, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x, y + 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x, y + 1, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x + 1, y, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x, y, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x, y, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x - 1, y + 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x - 1, y + 1, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x + 1, y - 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y - 1, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x, y - 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x, y - 1, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x - 1, y, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x - 1, y, z + 1)]);
            if (Boxes.ContainsKey(new Vector3(x - 1, y - 1, z + 1))) NeighbourAtomList.AddRange(Boxes[new Vector3(x - 1, y - 1, z + 1)]);
            //Add the 4 neighbouring boxes on the right side of our box
            if (Boxes.ContainsKey(new Vector3(x + 1, y + 1, z))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y + 1, z)]);
            if (Boxes.ContainsKey(new Vector3(x, y + 1, z))) NeighbourAtomList.AddRange(Boxes[new Vector3(x, y + 1, z)]);
            if (Boxes.ContainsKey(new Vector3(x + 1, y, z))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y, z)]);
            if (Boxes.ContainsKey(new Vector3(x + 1, y - 1, z))) NeighbourAtomList.AddRange(Boxes[new Vector3(x + 1, y - 1, z)]);
            return NeighbourAtomList;
        }

        public List<List<int>> AdjacentAtomList(List<Vector3> AtomPos, float d_max) {
            //List of Adjacent atoms:
            Dictionary<Vector3, List<int>> Boxes = MakeBoxes(AtomPos, d_max);//Summarize atoms inside Boxes so we don't have to check all of the aloms against each other
            List<List<int>> Bonds = new List<List<int>>();//Empty list of Bonds between close 
                                                          //Fill the Bonds of each atom with an empty list (Initialize:)
            int CountAtoms = AtomPos.Count();
            for (int i = 0; i < CountAtoms; i++) {
                Bonds.Add(new List<int> { });//Add an empty thingy...
            }
            //Go through all the boxes one by one...
            foreach (KeyValuePair<Vector3, List<int>> box in Boxes) {
                Vector3 BoxPos = box.Key;//Coordinates of the current Box
                int LengthBoxPos = Boxes[BoxPos].Count();//Count how many Atoms are encompassed in the Box
                for (int i = 0; i < LengthBoxPos; i++) {//Go through all the Atoms in the Box
                    int AtomIndex = box.Value[i];
                    //Connections inside the box (Since the first loop goes over all the atoms, we don't need to check the atoms we've already checked in this loop:
                    for (int i2 = i + 1; i2 < LengthBoxPos; i2++) {
                        int AtomIndexConn = box.Value[i2];
                        Bonds = AddBond(AtomPos, AtomIndex, AtomIndexConn, Bonds, d_max);//Add a bond between AtomIndex & AtomIndexConn
                    }
                    //Connections in neighbouring Boxes:
                    List<int> Neighbours = NeighbouringAtoms(Boxes, box);
                    foreach (int AtomNeighbourIndex in Neighbours) {
                        Bonds = AddBond(AtomPos, AtomIndex, AtomNeighbourIndex, Bonds, d_max);
                    }
                }
            }
            return Bonds;
        }


        public List<int> RefineNeighbourBonds(List<Vector3> AtomPositions, List<float> AtomRaddii, List<int> indices, float probe, int CurrentAtomIndex) {
            //Give list of indices of atoms wich are within probe distance to the Atom (With current Atom index)
            List<int> NeigbourIndices = new List<int>();
            Vector3 AtomPos = AtomPositions[CurrentAtomIndex];
            float AtomRadius = AtomRaddii[CurrentAtomIndex];
            float Radius = AtomRadius + probe + probe;
            foreach (int i in indices) {
                if (i == CurrentAtomIndex) continue;
                //Console.WriteLine(i);
                float iAtomRadius = AtomRaddii[i];
                Vector3 AtomNeighbourPos = AtomPositions[i];
                float Dist = Vector3.Distance(AtomNeighbourPos, AtomPos);
                if (Dist < (Radius + iAtomRadius)) {
                    NeigbourIndices.Add(i);
                }
            }

            return NeigbourIndices;

        }

        public List<List<int>> DeleteSplitAtomIndexes(List<List<int>> AdjAtomsAll, int AtomsInSplit, bool ProtA) {
            int WholeAtomNbr = AdjAtomsAll.Count();
            //Console.WriteLine("ProtA=" + ProtA.ToString() + "All Atom nbr:"+ WholeAtomNbr+ " - AtomIn splits: " + AtomsInSplit);
            List<List<int>> NewAdjList = new List<List<int>>(AdjAtomsAll);
            if (ProtA == true) {
                NewAdjList.RemoveRange(AtomsInSplit, WholeAtomNbr - AtomsInSplit);
                NewAdjList = NewAdjList.Select(FrstList => FrstList.Where(SecListID => SecListID < AtomsInSplit).ToList()).ToList();
            } else {
                int delete_before = WholeAtomNbr - AtomsInSplit;
                NewAdjList.RemoveRange(0, delete_before);
                NewAdjList = NewAdjList.Select(FrstList => FrstList.Where(SecListID => SecListID >= delete_before).ToList().Select(secListID => secListID - delete_before).ToList()).ToList();
            }
            return NewAdjList;
        }
    }
}
