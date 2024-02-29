using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace SplitProteinPrediction {
    class Logit {

        public List<double> GetActivityProbability(List<double> Cons, List<float> BAff, List<float> RelASA, List<string> SecStr, List<float> CDocking, int SeqLen) {
            List<double> Probs = new List<double>();
            /*[[-1.72554112  1.2471668   2.13274471 -0.65159228  0.00367445]] [-0.54106848]*/
            /*1/(1+np.exp(-p_y))*/
            float Coeff_Cons = -1.72554112f;
            float Coeff_BAff = 1.2471668f;
            float Coeff_RelASA = 2.13274471f;
            float Coeff_SecStr = -0.65159228f;
            float Coeff_CDocking = 0.00367445f;
            float Coeff_intercept = -0.54106848f;
            for (int split_site = 1; split_site < SeqLen; split_site++) {
                float issecstr = 0f;
                if (SecStr[split_site - 1] == "b" && SecStr[split_site] == "b"){
                    issecstr = 1f;
                }
                float x_value = Coeff_Cons * (float)Cons[split_site - 1]  + Coeff_BAff * BAff[split_site-1] + Coeff_RelASA * RelASA[split_site - 1] + Coeff_SecStr * issecstr +  Coeff_CDocking * CDocking[split_site - 1] + Coeff_intercept;
                double y_probability = 1 / (1 + Math.Exp(-x_value));
                Probs.Add(y_probability);
            }
            return Probs;
        }

    }
}
