using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Linq;
using System.IO;


namespace SplitProteinPrediction {
    class PDBContent {
        public List<string> Errors = new List<string>();//Error list
        public string FileDir;//File Title
        public List<string> AtomNames = new List<string>();//The coordinates
        public List<Vector3> AtomPositions = new List<Vector3>();//The coordinates
        public List<Vector3> Ca_AtomPositions = new List<Vector3>();//C alpha Coordinates
        public List<string> AtomInfo = new List<string>();//The Line with atomar info
        public List<float> RadiiList = new List<float>();//Atom Radii
        public List<string> SingleLetterSequence = new List<string>();//Sequence as it is in the PDB (Atom Sequence!)
        public List<int> SplitAtSite = new List<int>();//Marks the last line (index+1 !!!) before the Residue changes
        public List<float> AtomAreas = new List<float>();//ASA of the Atoms
        public List<float> RelativeResidueASA = new List<float>();//ASA of every residue divided by the Maximum possible ASA for that residue
        public List<int> ResidueIDs = new List<int>();
        public List<string> ResNumbersPDBChain = new List<string>();
        public float WholeProteinASA;
        //Residue Contacts:
        public List<List<int>> AtomContacts = new List<List<int>>();//Contacts between atoms (5.5A)
        public List<List<int>> ResidueContacts = new List<List<int>>();//Contacts between residues based on atom contacts
        public List<int> AtomIndexToResidueIndex = new List<int>();//For conversion (Residue contacts)
        public List<float> NIS_ValuesPercent = new List<float>();//Non interacting surface Values
        //Split Values:
        public List<float> BindingAffinities = new List<float>();
        public List<string> SecondaryStructure = new List<string>();
    }
}
