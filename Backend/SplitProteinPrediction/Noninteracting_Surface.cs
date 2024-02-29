using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace SplitProteinPrediction {
    class Noninteracting_Surface {

        /*This is a re-implementation in C# of the python code of PRODIGY as found here: https://github.com/haddocking/prodigy/blob/main/prodigy/predict_IC.py */

        public PDBContent GenerateNISVals(PDBContent PDBCont) {
            //Iterate through Residues which are exposed:
            AA_Values AAVals = new AA_Values();
            List<float> RelativeResidueASA = PDBCont.RelativeResidueASA;
            List<string> Sequence = PDBCont.SingleLetterSequence;
            List<float> NISResult = new List<float> { 0f, 0f, 0f };
            Dictionary<string, int> ACP_ToIndex = new Dictionary<string, int>() { { "A", 0 }, { "C", 1 }, { "P", 2 } };
            int ResidueIndex = 0;
            foreach(float RelASA in RelativeResidueASA) {
                if(RelASA >= 0.05f) {//It's an NIS
                    NISResult[ACP_ToIndex[AAVals.aa_character_protorp[Sequence[ResidueIndex]]]]++;
                }
                ResidueIndex++;
            }
            List<float> NISResultPercent = new List<float>();
            float sumNIS = NISResult.Take(3).Sum();
            foreach (int NISRes in NISResult) {
                float percent = (100f * NISRes) / sumNIS;
                NISResultPercent.Add(percent);
            }
            PDBCont.NIS_ValuesPercent = NISResultPercent;
            return PDBCont;
        }
    }
}
