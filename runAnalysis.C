// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_ROOT)
//#include <ANALYSIS/macros/AddTaskPIDResponse.C>
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <OADB/macros/AddTaskPhysicsSelection.C>
#endif

#include "AliAnalysisTaskHyperFinder3Body.h"

void runAnalysis()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    // Bool_t local = kTRUE;
    Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    // Bool_t gridTest = kTRUE;
    Bool_t gridTest = kFALSE;
    
    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
     
#if !defined (__CINT__) || defined (__CLING__)
  TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
  AliAnalysisTaskPIDResponse* PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(PIDadd.Exec());
  TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multSelection.Exec());

#else
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask* PIDresponseTask = AddTaskPIDResponse();
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* multSelectionTask = AddTaskMultSelection(kFALSE);
#endif

  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(0,1);

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskHyperFinder3Body.cxx++g");
    AliAnalysisTaskHyperFinder3Body *task = reinterpret_cast<AliAnalysisTaskHyperFinder3Body*>(gInterpreter->ExecuteMacro("AddTask_HyperFinder3Body.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskHyperFinder3Body.cxx++g");
    gROOT->LoadMacro("AddTask_HyperFinder3Body.C");
    AliAnalysisTaskHyperFinder3Body *task = AddTask_HyperFinder3Body();
#endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("~/alice/my_task/data/01/AliESDs.root");
        // chain->Add("~/alice/my_task/data/02/AliESDs.root");
        // chain->Add("~/alice/my_task/data/03/AliESDs.root");
        // chain->Add("~/alice/my_task/data/04/AliESDs.root");
        // chain->Add("~/alice/my_task/data/05/AliESDs.root");
        // chain->Add("~/alice/my_task/data/06/AliESDs.root");
        // chain->Add("~/alice/my_task/data/07/AliESDs.root");
        // chain->Add("~/alice/my_task/data/08/AliESDs.root");
        // chain->Add("~/alice/my_task/data/09/AliESDs.root");
        // chain->Add("~/alice/my_task/data/10/AliESDs.root");
        chain->Add("~/alice/my_task/data/11/AliESDs.root");
        // chain->Add("~/alice/my_task/data/12/AliESDs.root");
        chain->Add("~/alice/my_task/data/13/AliESDs.root");
        chain->Add("~/alice/my_task/data/14/AliESDs.root");
        chain->Add("~/alice/my_task/data/15/AliESDs.root");
        chain->Add("~/alice/my_task/data/16/AliESDs.root");
        chain->Add("~/alice/my_task/data/17/AliESDs.root");
        chain->Add("~/alice/my_task/data/18/AliESDs.root");
        chain->Add("~/alice/my_task/data/19/AliESDs.root");
        chain->Add("~/alice/my_task/data/20/AliESDs.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskHyperFinder3Body.cxx AliAnalysisTaskHyperFinder3Body.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskHyperFinder3Body.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        //alienHandler->SetAliPhysicsVersion("vAN-20210316_ROOT6-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20210810_ROOT6-1");
        alienHandler->SetAliPhysicsVersion("vAN-20221227_O2-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20181028_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
        // define the output folders
        //alienHandler->SetGridWorkingDir("HyperTriton2He3piML");
        alienHandler->SetGridWorkingDir("HyperFinder3Body");
        alienHandler->SetGridOutputDir("myOutputDir");

        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");

        //Choose runnumber Type
    
        alienHandler->SetDataPattern("/pass3/*/AliESDs.root");
        //LHC15o highIR pass1:DPG CentralBarrelTracking 20161130 v6
        int totalRunNumber = 50;
        int workRunNumber[50] = {296935, 296938, 296941, 296966, 297031, 297035, 297085, 297117, 297118, 297119,
            297123, 297124, 297128, 297129, 297132, 297133, 297193, 297195, 297196, 297218, 
            297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 297336, 
            297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 
            297415, 297441, 297442, 297446, 297450, 297451, 297452, 297479, 297483, 297512};
           /* {
            18q pass3:
            295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 
            295668, 295673, 295675, 295676, 295712, 295714, 295717, 295718, 295719, 295721, 
            295723, 295725, 295754, 295755, 295758, 295759, 295762, 295763, 295786, 295788,//zhengqiw 
            
            295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 295853, 
            295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 
            295937, 295941, 295942, 296016, 296060, 296062, 296063, 296065, 296066, 296074, 
            296123, 296132, 296133, 296134, 296135, 296142, 296143, 296191, 296192, 296194, 
            296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 296246, 
            296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312,//sozhang 
            
            296375, 296376, 296377, 296378, 296379, 296380, 296381, 296383, 296414, 296415, 
            296419, 296420, 296423, 296424, 296433, 296472, 296509, 296510, 296511, 296512, 
            296516, 296547, 296548, 296550, 296551, 296552, 296553, 296594, 296615, 296616, 
            296618, 296619, 296621, 296622, 296623                                                                    
            18r pass3:
            296690, 296691, 296693, 296694, 296752, 296781, 296784, 296785, 296786, 296787, 
            296790, 296793, 296794, 296799, 296835, 296836, 296838, 296839, 296848, 296850,
            296851, 296852, 296894, 296899, 296900, 296903, 296930, 296931, 296932, 296934,//qshou 
            
            296935, 296938, 296941, 296966, 297031, 297035, 297085, 297117, 297118, 297119,
            297123, 297124, 297128, 297129, 297132, 297133, 297193, 297195, 297196, 297218, 
            297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 297336, 
            297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 
            297415, 297441, 297442, 297446, 297450, 297451, 297452, 297479, 297483, 297512, 
            297537, 297540, 297541, 297542, 297544, 297558, 297588, 297590, 297595,        //yugang
            
            18r pass3 TimeRangeCut:
            296749, 296750, 296849, 296890, 297029, 297194, 297219, 297481(lagre sub jobs)//zhengqiw
            */
        //Set RunNumber
        int runStart = 297624;
        int runEnd   = 295584;
        // alienHandler->AddRunNumber(246994);

        for (int i=0; i<totalRunNumber; i++){
            if (workRunNumber[i] <= runStart && workRunNumber[i] >=runEnd)
            {
               alienHandler->AddRunNumber(workRunNumber[i]);
            }
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(100);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(36000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(3);
        alienHandler->SetMergeViaJDL(kTRUE);
        // alienHandler->SetMergeViaJDL(kFALSE);

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(3);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            // alienHandler->SetRunMode("full");
            alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}