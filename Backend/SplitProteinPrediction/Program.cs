using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.IO;
using System.Linq;

using System.Configuration;
using System.Collections.Specialized;

using System.Net.Http;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;

using System.Globalization;
using static System.Net.Mime.MediaTypeNames;


namespace SplitProteinPrediction {
    class Program {

        //The progress bar is from here: https://stackoverflow.com/questions/24918768/progress-bar-in-console-application
        private static void ProgressBar(int progress, int tot) {
            //draw empty progress bar
            //Console.WriteLine("31 / " + tot.ToString() + " " + progress.ToString());
            Console.SetCursorPosition(0, Console.CursorTop);
            Console.Write("["); //start
            Console.CursorLeft = 32;
            Console.Write("]"); //end
            Console.CursorLeft = 1;
            float onechunk = 31.0f / tot;

            //draw filled part
            int position = 1;
            for (int i = 0; i < onechunk * progress; i++) {
                Console.BackgroundColor = ConsoleColor.Green;
                Console.CursorLeft = position++;
                Console.Write(" ");
            }



            //draw totals
            Console.CursorLeft = 35;
            Console.BackgroundColor = ConsoleColor.Black;
            Console.Write(progress.ToString() + " of " + tot.ToString() + "    "); //blanks at the end remove any excess
        
        }

        //This snipped is from here: https://stackoverflow.com/a/33409309
        public static string CalculateMD5Hash(string input) {
            // step 1, calculate MD5 hash from input
            System.Security.Cryptography.MD5 md5 = System.Security.Cryptography.MD5.Create();
            byte[] inputBytes = System.Text.Encoding.ASCII.GetBytes(input);
            byte[] hash = md5.ComputeHash(inputBytes);

            // step 2, convert byte array to hex string
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < hash.Length; i++) {
                sb.Append(hash[i].ToString("X2"));
            }
            return sb.ToString();
        }

        public delegate long StatisticalData();
        static void Main(string[] args) {

            /*Initialization stage*/
            string CurrentInputFile = "";
            string FileType = "";
            string CurrentJobTitle = "";
            string inputfile_beginning = "";

            Thread.CurrentThread.CurrentCulture = new CultureInfo(ConfigurationManager.AppSettings.Get("CurrentCulture"));

            PDBParser PDBPase = new PDBParser();
            ASA_Functions ASAFuncs = new ASA_Functions();
            Shrake_Rupeley_SASA ShrakeRupeleySASA = new Shrake_Rupeley_SASA();
            Lee_Richards_SASA LeeRichSASA = new Lee_Richards_SASA();
            Interface_Contacts InterfaceCont = new Interface_Contacts();
            Noninteracting_Surface NIS = new Noninteracting_Surface();
            Prodigy_Function Prodigy = new Prodigy_Function();

            Run_DSSP DSSP = new Run_DSSP();
            CFragDocking CFragDock = new CFragDocking();
            Logit Logitf = new Logit();



            Hmmer HmmerP = new Hmmer();
            Muscle MucleMSA = new Muscle();
            Rate4Site R4S = new Rate4Site();

            //Bonds:
            HydrogenBond_Calculator HBonds = new HydrogenBond_Calculator();
            SaltBridge_Calculator SBridge = new SaltBridge_Calculator();
            BondTextOutput textOutput = new BondTextOutput();
            AromaticConnection_Calculator Arom = new AromaticConnection_Calculator();
            VanDerWaals_Calculator VdwFunct = new VanDerWaals_Calculator();

            ClusterSequencesR4S ClusterAlign = new ClusterSequencesR4S();

            OpenZip ZipO = new OpenZip();
            RunPDB2PQR pdb2pqr = new RunPDB2PQR();

            string DirOutput = "";
            try {
                //Slant-Relief Font from: Mega-relief by Nick Bryant 12/9, dmc1@st-and.ac.uk for the moment. Modified by Mirko Schmitz
                Console.WriteLine("\r\n          ####                \r\n     **/////((##\\)    #       /\\\\\\\\\\\\\\\\\\\\\\                                       /\\\\\\        \r\n   /*///(((((((####)   ##     \\/////\\\\\\///                                     /\\\\\\//\\\\\\                \r\n  /////((######%%%#)    ##%        \\/\\\\\\                       /\\\\\\            /\\\\\\ /\\\\\\        /\\\\\\            \r\n *///(((##)            ####         \\/\\\\\\      /\\\\/\\\\\\\\\\\\    /\\\\\\\\\\\\\\\\\\\\\\      \\//\\\\\\\\//        \\///   /\\\\/\\\\\\\\\\\\\r\n //*((              #######          \\/\\\\\\     \\/\\\\\\////\\\\\\  \\////\\\\\\////      /\\\\\\///\\\\\\         /\\\\\\ \\/\\\\\\////\\\\\\ \r\n ///(   ///(((((######(###            \\/\\\\\\     \\/\\\\\\  \\//\\\\\\    \\/\\\\\\        /\\\\\\/  \\///\\\\\\/\\\\\\  \\/\\\\\\ \\/\\\\\\  \\//\\\\\\ \r\n  //(   ///(((((##((((#/#              \\/\\\\\\     \\/\\\\\\   \\/\\\\\\    \\/\\\\\\ /\\\\   /\\\\\\      \\//\\\\\\//   \\/\\\\\\ \\/\\\\\\   \\/\\\\\\\r\n   //#   ///((((((((((*              /\\\\\\\\\\\\\\\\\\\\\\ \\/\\\\\\   \\/\\\\\\    \\//\\\\\\\\\\   \\//\\\\\\\\\\\\\\\\\\\\\\//\\\\\\   \\/\\\\\\ \\/\\\\\\   \\/\\\\\\\r\n      ((     /////                   \\///////////  \\///    \\///      \\/////     \\/////////// \\///    \\///  \\///    \\///\r\n        #                    ");
                Console.WriteLine("\r\nIf you use this tool please Cite:");
                Console.WriteLine("");
                if (args.Length >= 1) {
                    string CurrentStatus = "";
                    string FileNameInput = args[0];
                    if (FileNameInput == "--help" || FileNameInput ==  "-help" || FileNameInput == "--h" || FileNameInput == "-h" || FileNameInput == "help" || FileNameInput == "h") {
                        Console.WriteLine("Int&in Standalone tool - Help");
                        Console.WriteLine("Use this tool as follows: dotnet Intin.dll InputFileName DirectoryOutput FileNamePrefix");
                        Console.WriteLine("");
                        Console.WriteLine("Aruments: ");
                        Console.WriteLine("InputFileName: The Name (Extension .pdb or .zip) of the input file, the .zip must contain .pdb files");
                        Console.WriteLine("DirectoryOutput (optional): The directory file of the output");
                        Console.WriteLine("FileNamePrefix (optional): A prefix for given to the output files");
                        Console.WriteLine("");
                        Environment.Exit(0);
                    } else if (!File.Exists(FileNameInput)) {
                        throw new SplitProteinException("The Input file does not exist!");
                    }
                    DirOutput = Path.GetDirectoryName(FileNameInput);
                    string Prefix_FileName = "";
                    string OutputFilesName = "";
                    if (args.Length >= 2) {
                        DirOutput = args[1];
                    }
                    if (args.Length >= 3) {
                        Prefix_FileName = args[2];
                    }
                    DirOutput += "/";
                    Console.WriteLine(DirOutput);
                    OutputFilesName = Prefix_FileName + "_" + Path.GetFileNameWithoutExtension(FileNameInput);
                    if (!Directory.Exists(DirOutput)) {//Create output dir if not already existant
                        Directory.CreateDirectory(DirOutput);
                    }

                    FileType = FileNameInput.Split(".").LastOrDefault();

                    List<string> AllowedFileTypes = new List<string> { "pdb", "zip" };
                    if (AllowedFileTypes.Contains(FileType)) {
                        if (FileType == "pdb") {
                            inputfile_beginning = FileNameInput;
                        } else if (FileType == "zip") {
                            inputfile_beginning = FileNameInput;
                        }
                    }

                    if (File.Exists(inputfile_beginning)) {
                        List<string> Files = new List<string>();
                        //extract zip file & go through each file...
                        if (FileType == "zip") {
                            ZipO.OpenZipFile(inputfile_beginning, DirOutput);
                            //get the files from the zip
                            DirectoryInfo dirinfo = new DirectoryInfo(DirOutput);
                            foreach (var file in dirinfo.GetFiles("*.pdb")) {
                                string filename = file.FullName;
                                if (file.Name.Contains(" ")) {//rename input file if it has spaces:
                                    filename = file.FullName.Replace(" ", "_");
                                    File.Move(file.FullName, filename);
                                }
                                Files.Add(filename);
                            }
                        } else {
                            Files.Add(inputfile_beginning);
                        }

                        string FileOutput = "PDBName\tJobName\tSequence\tResidueIDs\tConservation [Residue]\tNormalized Conservation [Residue]\tSec Structure [Residue]\tBindingAffinity [Split]\tNormalized BindingAffinity [Split]\tRelASA [Residue]\tHydrogenBonds\tSBridges\tArom\tvDW\tCFragDocking\tConservation_Model\tRelASA_Model\tActivityPrediction\n";
                        string Filetitle = "Int&In_Output_data";

                        int files_count = Files.Count();
                        int current_file_nbr = 1;
                        CurrentStatus = "";
                        foreach (string inputfile in Files) {
                            CurrentInputFile = inputfile;
                            string FileName = Path.GetFileNameWithoutExtension(inputfile);
                            if (FileType == "zip") {
                                OutputFilesName = Prefix_FileName + "_" + FileName;
                            }
                            Console.WriteLine(OutputFilesName + " is being analysed. (File " + current_file_nbr + "/" + files_count + ")");

                            string UniqueIdHash = CalculateMD5Hash(System.DateTime.Now.ToString("MM/dd/yyyy HH:mm:ss").ToString() + "_" + FileName);
                            UniqueIdHash = FileName + "_" + UniqueIdHash;

                            string inputfile_hydrogens = pdb2pqr.Runcmd(inputfile);
                            inputfile_hydrogens = inputfile_hydrogens.Replace("\\", "/");

                            PDBContent CurrentPDBContent = new PDBContent();
                            CurrentPDBContent = PDBPase.GetPDBData(inputfile_hydrogens, CurrentPDBContent);//Get The pdb2pqr files as we have hydrogens in here and the structure is fixed
                            CurrentPDBContent = DSSP.DSSP_Cmd(CurrentPDBContent, inputfile, UniqueIdHash, DirOutput);//DSSP has to be run with the original input file, as it doesn't understand pqr files...
                            //File.Move(inputfile_hydrogens, inputfile);//rename the pdb2pqr file to pdb, as for Int&in's parser it doesn't matter

                            CurrentPDBContent = LeeRichSASA.GetAtomASA(CurrentPDBContent);//important that CurrentPDBContent contains the pdb2pqr file!
                            CurrentPDBContent.WholeProteinASA = CurrentPDBContent.AtomAreas.Sum();
                            //Make a Proximity List of all the residues in the Protein (based on the Atoms!)
                            CurrentPDBContent = InterfaceCont.GenerateResidueContacts(CurrentPDBContent);
                            CurrentPDBContent = NIS.GenerateNISVals(CurrentPDBContent);
                            List<float> NISVals = CurrentPDBContent.NIS_ValuesPercent;
                            //Run Some external Programs:

                            List<float> RadiiList = CurrentPDBContent.RadiiList;
                            float BoxMaxLen = 10f * (1.4f + RadiiList.Max());//Size of the Box where Atoms are summerized in
                            List<List<int>> AdjAtoms = ASAFuncs.AdjacentAtomList(CurrentPDBContent.AtomPositions, BoxMaxLen);

                            //Calculate the Bonds
                            List<List<string>> hBonds = HBonds.GenerateHBonds(CurrentPDBContent, true, 80f, 3.5f);// generate unique contacts
                            List<string> HBondOutput_Text = textOutput.GetBondOutput(hBonds, true);

                            List<List<string>> sbridge = SBridge.GenerateSaltBridgesProtein(CurrentPDBContent, true, 4f);// generate unique contacts
                            List<string> SBridgeText = textOutput.GetBondOutput(sbridge);

                            List<List<string>> aromb = Arom.GenerateAromaticConnections(CurrentPDBContent, true, 6.5f, 5f);// generate unique contacts
                            List<string> AromText = textOutput.GetBondOutput(aromb);

                            List<List<string>> vdw = VdwFunct.GenerateVanderwaalsBonds(CurrentPDBContent, true, 0.5f);// generate unique contacts
                            List<string> vdWText = textOutput.GetBondOutput(vdw);

                            Console.WriteLine("Non-covalent bonds generated");

                            Task<int> hmmer_exec = HmmerP.RunHmmer(UniqueIdHash, string.Join("", CurrentPDBContent.SingleLetterSequence), DirOutput);

                            int AtomCount = CurrentPDBContent.AtomInfo.Count();
                            int ResidueCount = CurrentPDBContent.SingleLetterSequence.Count();

                            int ProteinLength = CurrentPDBContent.SingleLetterSequence.Count;
                            Console.WriteLine("Generating split energies");
                            for (int split_site = 1; split_site < ResidueCount; split_site++) {
                                string CurrentStatusSplitSite = "Evaluating split site " + split_site + "/" + ResidueCount;
                                ProgressBar(split_site, ResidueCount);
                                //Make Splits:
                                List<string> ProteinA = new List<string>();
                                List<string> ProteinB = new List<string>();
                                PDBContent ProteinACont = new PDBContent();
                                PDBContent ProteinBCont = new PDBContent();
                                //Grab the important Data
                                int PartALastIndex = CurrentPDBContent.SplitAtSite[split_site - 1] - 1;//Since the index is wanted, do -1 as splitAtsite gives us the line
                                int SizeB = AtomCount - PartALastIndex;
                                ProteinA = CurrentPDBContent.AtomInfo.GetRange(0, PartALastIndex);
                                ProteinB = CurrentPDBContent.AtomInfo.GetRange(PartALastIndex, SizeB);
                                //Lead the two parts in:
                                ProteinACont = PDBPase.GetSplitData(ProteinA, "ProtA", ProteinACont);
                                ProteinBCont = PDBPase.GetSplitData(ProteinB, "ProtB", ProteinBCont);
                                //Get the ASA of the two parts (As well as the RelASA)
                                ProteinACont = ShrakeRupeleySASA.GetAtomASA(ProteinACont, AdjAtoms, false, true);
                                ProteinBCont = ShrakeRupeleySASA.GetAtomASA(ProteinBCont, AdjAtoms, false, false);

                                //Cout the IC's between the two & do the Binding affinity prediction:
                                Dictionary<string, int> Ic_cont = InterfaceCont.CountContactTypes(CurrentPDBContent, split_site);
                                float BindingAffinity = Prodigy.PRODIGYv1(Ic_cont["CC"], Ic_cont["AC"], Ic_cont["PP"], Ic_cont["AP"], NISVals[0], NISVals[1]);
                                CurrentPDBContent.BindingAffinities.Add(BindingAffinity);
                            }

                            Console.WriteLine("All split sites analysed, waiting for HMMER");
                            hmmer_exec.Wait();
                            Console.WriteLine("HMMER done, performing MSA");

                            List<Hmmer.HmmerAlignment> HmmerAlign = HmmerP.ReadFile(UniqueIdHash, CurrentPDBContent.SingleLetterSequence, DirOutput);

                            List<string> SequencesForMSA = MucleMSA.GetHmmerDataForMSA(string.Join("", CurrentPDBContent.SingleLetterSequence), HmmerAlign);

                            List<string> SequencesForMSASample = MucleMSA.SampleSeq(SequencesForMSA);//Sampling of the files

                            string MSA_Result = MucleMSA.RunMuscle(UniqueIdHash, SequencesForMSASample, DirOutput);

                            Console.WriteLine("Muscle done");

                            ClusterAlign.GenerateClusteredAlignment(UniqueIdHash, DirOutput);//generate cluster stuff

                            string RunRate4Site = R4S.RunRate4Site(UniqueIdHash, DirOutput);
                            Console.WriteLine("Rate4Site done");

                            List<string> Rate4SiteScores = R4S.ReadRate4Site(RunRate4Site);

                            List<double> Rate4SiteScoresNormalized = R4S.NormalizeRate4Site(Rate4SiteScores);

                            List<float> NormBAff = new List<float>();
                            float minBAff = CurrentPDBContent.BindingAffinities.Min();
                            if (minBAff == 0) {
                                minBAff = 1;
                            }
                            NormBAff = CurrentPDBContent.BindingAffinities.Select(i => i / minBAff).ToList();

                            List<double> Conservation_Model = new List<double>();
                            List<float> relASAModel = new List<float>();
                            for (int split_site = 1; split_site < ResidueCount; split_site++) {
                                //generate the apporopriate values for each split site
                                int takeIndexes = 4;
                                int divideby = 4;
                                int rest = Rate4SiteScoresNormalized.Count - split_site + 1;
                                if (rest < 4) {
                                    takeIndexes = rest;
                                    divideby = rest;
                                }
                                double minus = 0f;
                                if (takeIndexes >= 3) {
                                    minus = Rate4SiteScoresNormalized[split_site + 1];
                                    divideby -= 1;
                                }
                                double Cons = (Rate4SiteScoresNormalized.GetRange(split_site - 1, takeIndexes).Sum() - minus) / (float)divideby;
                                Conservation_Model.Add(Cons);
                                
                                float relasamdl = CurrentPDBContent.RelativeResidueASA[split_site];
                                if (split_site >= 2) {
                                    relasamdl = (CurrentPDBContent.RelativeResidueASA[split_site - 2] + CurrentPDBContent.RelativeResidueASA[split_site]) / 2f;
                                }
                                relASAModel.Add(relasamdl);

                            }
                            //Use the generated data

                            List<float> Docking = CFragDock.Get_CFragDocking(HBondOutput_Text, SBridgeText, AromText, vdWText, ResidueCount);
                            List<double> Probabilities = Logitf.GetActivityProbability(Conservation_Model, NormBAff, relASAModel, CurrentPDBContent.SecondaryStructure, Docking, ResidueCount);
                           
                            FileOutput += FileName + "\t";
                            FileOutput += CurrentJobTitle + "\t";
                            FileOutput += string.Join("", CurrentPDBContent.SingleLetterSequence) + "\t";
                            FileOutput += string.Join(";", CurrentPDBContent.ResNumbersPDBChain) + "\t";
                            FileOutput += string.Join(";", Rate4SiteScores) + "\t";
                            FileOutput += string.Join(";", Rate4SiteScoresNormalized) + "\t";
                            FileOutput += string.Join("", CurrentPDBContent.SecondaryStructure) + "\t";
                            FileOutput += string.Join(";", CurrentPDBContent.BindingAffinities) + "\t";
                            FileOutput += string.Join(";", NormBAff) + "\t";
                            FileOutput += string.Join(";", CurrentPDBContent.RelativeResidueASA) + "\t";

                            FileOutput += string.Join(";", HBondOutput_Text) + "\t";
                            FileOutput += string.Join(";", SBridgeText) + "\t";
                            FileOutput += string.Join(";", AromText) + "\t";
                            FileOutput += string.Join(";", vdWText) + "\t";
                            FileOutput += string.Join(";", Docking) + "\t";

                            FileOutput += string.Join(";", Conservation_Model) + "\t";
                            FileOutput += string.Join(";", relASAModel) + "\t";

                            FileOutput += string.Join(";", Probabilities) + "\t";

                            FileOutput += "\n";

                            current_file_nbr++;
                        }
                        //save file
                        string ResultsDirSave = DirOutput + Filetitle + ".tsv";
                        File.WriteAllText(ResultsDirSave, FileOutput);
                    }
                } else {
                    throw new SplitProteinException("Input arguments wrong");
                }
            } catch (Exception ex) {
                var get_type = ex.GetType().FullName;
                string ErrorText = "An Error occured: \n Message: \n" + ex.Message + "\n\nStackTrace: \n" + ex.StackTrace + "\n\n\nYou can access the help via entering: dotnet Intin.dll -help";
                Console.WriteLine(ErrorText);
                if (DirOutput != "") {
                    File.WriteAllText(DirOutput + "ErrorText.txt", ErrorText);
                }
            }
        }
    }


}
