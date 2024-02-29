using System;
using System.Collections.Generic;
using System.Text;

namespace SplitProteinPrediction {
    class AA_Values {


        public Dictionary<string, float> ASA_MaxResidue = new Dictionary<string, float>() { { "A", 107.95f }, { "R", 238.76f }, { "N", 143.94f }, { "D", 140.39f }, { "C", 134.28f }, { "E", 172.25f }, { "Q", 178.50f }, { "G", 80.10f }, { "H", 182.88f }, { "I", 175.12f }, { "L", 178.63f }, { "K", 200.81f }, { "M", 194.15f }, { "F", 199.48f }, { "P", 136.13f }, { "S", 116.50f }, { "T", 139.27f }, { "W", 249.36f }, { "Y", 212.76f }, { "V", 151.44f } };
        
        public Dictionary<string, string> AA_3LetterCodeToSingle = new Dictionary<string, string>() {{ "CYS", "C"}, { "ASP", "D"}, { "SER", "S"}, { "GLN", "Q"}, { "LYS", "K"},
                                                                                            { "ILE", "I"}, { "PRO", "P"}, { "THR", "T"}, { "PHE", "F"}, { "ASN", "N"},
                                                                                            { "GLY", "G"}, { "HIS", "H"}, { "LEU", "L"}, { "ARG", "R"}, { "TRP", "W"},
                                                                                            { "ALA", "A"}, { "VAL","V"}, { "GLU", "E"}, { "TYR", "Y"}, { "MET", "M"}};

        //The van der Waals radii are from naccess, from the PRODIGY git: https://github.com/haddocking/prodigy/blob/main/prodigy/naccess.config

        public Dictionary<string, float> AtomVanDerWaals = new Dictionary<string, float>() { { "CA", 1.87f }, { "CB", 1.87f }, { "C", 1.76f }, { "S", 1.85f }, { "N", 1.65f }, { "O", 1.4f }, { "H", 0f } };

        public Dictionary<string, float> AtomVanDerWaalsv2 = new Dictionary<string, float>() {
                                                                                            { "ALA CB", 1.87f } ,{ "ARG CG", 1.87f } ,{ "ARG CD", 1.87f } ,{ "ARG NE", 1.65f } ,{ "ARG CZ", 1.76f } ,{ "ARG NH1", 1.65f } ,{ "ARG NH2", 1.65f } ,{ "ASN CG", 1.76f } ,{ "ASN OD1", 1.40f } ,{ "ASN ND2", 1.65f } ,{ "ASP CG", 1.76f } ,{ "ASP OD1", 1.40f } ,{ "ASP OD2", 1.40f } ,{ "CYS SG", 1.85f } ,{ "GLN CG", 1.87f } ,{ "GLN CD", 1.76f } ,{ "GLN OE1", 1.40f } ,{ "GLN NE2", 1.65f } ,{ "GLU CG", 1.87f } ,{ "GLU CD", 1.76f } ,{ "GLU OE1", 1.40f } ,{ "GLU OE2", 1.40f } ,{ "GLY CA", 1.87f } ,{ "HIS CG", 1.76f } ,{ "HIS ND1", 1.65f } ,{ "HIS CD2", 1.76f } ,{ "HIS NE2", 1.65f } ,{ "HIS CE1", 1.76f } ,{ "ILE CG1", 1.87f } ,{ "ILE CG2", 1.87f } ,{ "ILE CD1", 1.87f } ,{ "LEU CG", 1.87f } ,{ "LEU CD1", 1.87f } ,{ "LEU CD2", 1.87f } ,{ "LYS CG", 1.87f } ,{ "LYS CD", 1.87f } ,{ "LYS CE", 1.87f } ,{ "LYS NZ", 1.50f } ,{ "MET CG", 1.87f } ,{ "MET SD", 1.85f } ,{ "MET CE", 1.87f } ,{ "PHE CG", 1.76f } ,{ "PHE CD1", 1.76f } ,{ "PHE CD2", 1.76f } ,{ "PHE CE1", 1.76f } ,{ "PHE CE2", 1.76f } ,{ "PHE CZ", 1.76f } ,{ "PRO CG", 1.87f } ,{ "PRO CD", 1.87f } ,{ "SER OG", 1.40f } ,{ "THR OG1", 1.40f } ,{ "THR CG2", 1.87f } ,{ "TRP CG", 1.76f } ,{ "TRP CD1", 1.76f } ,{ "TRP CD2", 1.76f } ,{ "TRP NE1", 1.65f } ,{ "TRP CE2", 1.76f } ,{ "TRP CE3", 1.76f } ,{ "TRP CZ2", 1.76f } ,{ "TRP CZ3", 1.76f } ,{ "TRP CH2", 1.76f } ,{ "TYR CG", 1.76f } ,{ "TYR CD1", 1.76f } ,{ "TYR CD2", 1.76f } ,{ "TYR CE1", 1.76f } ,{ "TYR CE2", 1.76f } ,{ "TYR CZ", 1.76f } ,{ "TYR OH", 1.40f } ,{ "VAL CG1", 1.87f } ,{ "VAL CG2", 1.87f }
                                                                                            };

        public Dictionary<string, string> aa_character_ic = new Dictionary<string, string>() {
                                                                                                {"A", "A"}, {"C", "A"},{"E", "C"}, {"D", "C"}, {"G", "A"}, {"F", "A"}, {"I", "A"}, {"H", "C"}, {"K", "C"},
                                                                                                {"M", "A"}, {"L", "A"}, {"N", "P"}, {"Q", "P"}, {"P", "A"}, {"S", "P"}, {"R", "C"}, {"T", "P"}, {"W", "A"}, {"V", "A"}, {"Y", "A" }
                                                                                            };

        public Dictionary<string, string> aa_character_protorp = new Dictionary<string, string>() {
                                                                                                {"A", "A"}, {"C", "P"}, {"E", "C"}, {"D", "C"}, {"G", "A"}, {"F", "A"}, {"I", "A"}, {"H", "P"}, {"K", "C"},
                                                                                                {"M", "A"}, {"L", "A"}, {"N", "P"}, {"Q", "P"}, {"P", "A"}, {"S", "P"}, {"R", "C"}, {"T", "P"}, {"W", "P"}, {"V", "A"}, {"Y", "P" }
                                                                                            };
        
    }
}
