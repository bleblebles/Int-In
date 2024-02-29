using System;
using System.Collections.Generic;
using System.Text;

namespace SplitProteinPrediction {
    class SplitProteinException : Exception {

        public SplitProteinException(string exception_string)
            : base(String.Format("Error: {0}", exception_string)) {

        }

    }
}
