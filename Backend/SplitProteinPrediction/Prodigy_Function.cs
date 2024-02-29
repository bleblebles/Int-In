using System;
using System.Collections.Generic;
using System.Text;

namespace SplitProteinPrediction {

    /*This is a re-implementation in C# of the python code of PRODIGY as found here: https://github.com/haddocking/prodigy/blob/main/ */

    class Prodigy_Function {
        public float PRODIGYv1(int ic_cc, int ic_ca, int ic_pp, int ic_pa, float p_nis_a, float p_nis_c) {
            float Fuct = -0.09459f * ic_cc + -0.10007f * ic_ca + 0.19577f * ic_pp + -0.22671f * ic_pa + 0.18681f * p_nis_a + 0.13810f * p_nis_c + -15.9433f; //+ 11.88802542;
            return Fuct;
        }
    }
}
