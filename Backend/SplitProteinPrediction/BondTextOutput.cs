using System;
using System.Collections.Generic;
using System.Text;

namespace SplitProteinPrediction
{
    class BondTextOutput
    {
        public List<string> GetBondOutput(List<List<string>> Bonds, bool HBond = false)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < Bonds.Count; i++)
            {
                List<string> hBondRes = Bonds[i];
                if (HBond == true)
                {
                    output.Add(hBondRes[0] + "-" + hBondRes[1] + "|" + hBondRes[2] + "|" + hBondRes[3]);
                }
                else
                {
                    output.Add(hBondRes[0] + "-" + hBondRes[1]);
                }
            }

            return output;
        }
    }
}
