using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.IO;
using System.Linq;

namespace SplitProteinPrediction {

    /*This is a re-implementation in C# of the python code of PRODIGY as found here: https://github.com/haddocking/prodigy/blob/main/prodigy/predict_IC.py */

    class Interface_Contacts {
        public PDBContent GenerateResidueContacts(PDBContent PDBCont, bool UniqueContacts = true, float dist = 5.5f) {
            /*Index=Residue id, value=list with contacts*/
            ASA_Functions ASAFucs = new ASA_Functions();
            List<Vector3> AtomPos = PDBCont.AtomPositions;
            List<string> Sequence = PDBCont.SingleLetterSequence;
            List<int> AtmIndToResIndex = PDBCont.AtomIndexToResidueIndex;
            int SeqLength = Sequence.Count();
            List<int> SplitAtSite = PDBCont.SplitAtSite;
            List<List<int>> AtomContacts = ASAFucs.AdjacentAtomList(AtomPos, dist);//Run it again, as the distace cutoff is different here... (Later maybe integrate it into the ASA)
            List<List<int>> ResidueContactsAdd = new List<List<int>>();
            //Convert AtomContacts to Residue Contacts:
            int charIndex = 0;
            int lastSplitSeqChar = 1;
            foreach (string SeqChar in Sequence) {
                if (charIndex < SeqLength) {
                    List<int> ResidueContacts = new List<int>();
                    //Get the Atom numbers which we consider here...
                    int SplitSeqChar = SplitAtSite[charIndex];//SplitAtSite gives the line number (not index!) after which a new Residue comes
                    int Start = lastSplitSeqChar - 1;
                    int length = SplitSeqChar - lastSplitSeqChar;
                    List<List<int>> AreaAccessPointsResidue = AtomContacts.GetRange(Start, length);
                    //Create a List with all the atoms linked to the Residue (iterate through its own atoms) with char_index
                    List<int> ResidueContactsToAtoms = new List<int>();
                    foreach (List<int> ContactsofAtom in AreaAccessPointsResidue) {
                        ResidueContactsToAtoms.AddRange(ContactsofAtom);
                    }
                    //Convert the Atoms to residues
                    List<int> ResidueContactsToResidues = new List<int>();
                    foreach (int AtomIndex in ResidueContactsToAtoms) {
                        int ResidueIndex = AtmIndToResIndex[AtomIndex];
                        if ((UniqueContacts == true && charIndex < ResidueIndex) || (UniqueContacts == false)) {
                            //Check if the ResidueContactsToResidues has already this contact but in reverse
                            //-> Only accept residues which come after our own residue because we've already checked those...
                            if (!ResidueContactsToResidues.Contains(ResidueIndex) && ResidueIndex != charIndex) {
                                ResidueContactsToResidues.Add(ResidueIndex);
                            }
                        }
                    }
                    lastSplitSeqChar = SplitSeqChar;
                    //Add the Residue contacts to a list
                    ResidueContactsAdd.Add(ResidueContactsToResidues);//If there're no matches, then an empty list will be added...
                }
                charIndex++;
                //Save
            }
            PDBCont.ResidueContacts = ResidueContactsAdd;
            return PDBCont;

        }

        public Dictionary<string, int> CountContactTypes(PDBContent WholeProtein, int SplitSite) {
            //split site = 1 => Cut after 1st residue
            AA_Values AAVals = new AA_Values();
            Dictionary<string, int> Bins = new Dictionary<string, int>() { { "AA", 0 }, { "PP", 0 }, { "CC", 0 }, { "AP", 0 }, { "CP", 0 }, { "AC", 0 } };
            //All the contacts between ProtA and B:
            List<List<int>> IC_Contacts = WholeProtein.ResidueContacts.GetRange(0, SplitSite);
            int ProtA_ResidueIndex = 0;
            foreach (List<int> Contacts in IC_Contacts) {
                string LetterProtA = AAVals.aa_character_ic[WholeProtein.SingleLetterSequence[ProtA_ResidueIndex]];
                foreach (int ProtBIndex in Contacts) {
                    if (ProtBIndex >= SplitSite) {//split site = 1 protbindex = 0
                        string LetterProtB = AAVals.aa_character_ic[WholeProtein.SingleLetterSequence[ProtBIndex]];
                        string type = LetterProtA + LetterProtB;
                        type = String.Concat(type.OrderBy(c => c));//Order the string (A to the front)
                        Bins[type]++;
                    }
                }
                ProtA_ResidueIndex++;
            }
            return Bins;
        }

    }
}
