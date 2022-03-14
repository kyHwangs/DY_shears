// Definition of the Squeezer main class
// Author: Ph. Gras CEA/IRFU Saclay
// Jan. 3, 11. Squeezer code for AnaNaS framework
// Jul. 24, 15. Cloned Squeezer and adapted to the Shears framework

#ifndef PRUNER_H
#define PRUNER_H

#include <algorithm>
#include <iomanip>
#include <memory>
#include <iostream>

#include "TObject.h"
#include "ShearsTChain.h"
#include "TH1F.h"
#define DECLARE_PRUNER(Class, Descr)				\
  extern "C" Pruner::ClassRecord *shearsLoadPlugin() { \
    Pruner::ClassRecord *rcd = new Pruner::ClassRecord; \
    rcd->className = #Class; \
    rcd->description = Descr; \
    rcd->instance = new Class(); \
    return rcd; \
  }
#define NPNLOBINS 3
class TFile;


/**
 * Utility class to copy and skim AnaNaS n-tuple files.  This class is used by
 * the standalone application squeezer which is based. It can also be
 * instantiated directly in custom code.
 *
 * By default the class allows to select a list of branches specified in a text
 * file for a list of events identified by their run number and event number
 * read from a second text file. Advanced selection can be peformed by
 * implementing a custom class as described below.
 *
 * Implementation of a custom Pruner:
 *
 *    //The four following lines are needed to import
 *    //the code generated by TTree::MakeClass()
 *    using namespace std;
 *    #define EventTree_cxx
 *    #include "EventTree.h"
 *    void EventTree::Loop(){} //To make the compiler/linker happy.
 *
 *    class MyPruner: public Pruner, EventTree{
 *
 *    protected:
 *         bool init(TTree* tree);
 *         bool filterEvent();
 *    };
 *
 *    DECLARE_PRUNER(MyPruner, "Description of your filter")
 * 
 *    bool MyPruner::init(TTree* tree){
 *      ...Set here the tree branch addresses. If MyPruner inherits
 *         from the class EventTree produced with TTreeMakeClass()
 *         then following code statement is sufficient...
 *         EventTree::Init(tree);
 *         return true;
 *    }
 *
 *    bool MyPruner::filterEvent(){
 *        bool pass;
 *        ... code to perform the event selection,
 *            which will set the variable pass according
 *            if the event pass the section. The
 *            code can modify the branch contents, for
 *            instance to skim the object collections...
 *        return pass && Pruner::filterEvent();
 *    }
 *
 * In the above example, the EventTree.h code generated with TTree::MakeClass()
 * is used to access to the tree branches.  Alternatively the user can write the
 * code to access the branches (calls to tree->SetBranch() method) in the
 * MyPruner::init(TTree* tree) method. In that case the four first lines of code
 * and the inheritance of MyPruner from EventTree are not needed.
 *
 * The EventTree.h file can be produced from a Boabab file: open the file with
 * ROOT, changing the current directory to tupel and calling
 * EventTree->MakeClass()
 *
 * To support subselections, the method declareSubSelection() should be
 * overridden and the vector subSelections must be filled with the name and
 * description of the offered subselections. The filerEvent() method can then
 * use the iSubselection_ index to determine which subselection should be used. 
 * Example:
 *
 *    void MyPruner:declareSubSelections(){
 *       subSeletions_.push_bach(SubSelection("Signal", "Signal selection"));
 *       subSeletions_.push_bach(SubSelection("ControlSample",
 *                                            "Control sample ""selection"));
 *          ...
 *    }
 *
 *
 * The filterBranch can be overidden for a custom branch selection. The default
 * one can optionally copy all the branches or read the list of branches to copy
 * from a text file.
 *
 * To use the custom class with the pruner application specified the class name
 * with the option --pruner.
 *
 */

class Pruner{
protected:
  /** Data types
   */
  ///@{
  struct SubSelection{
    std::string tag;
    std::string description;
    SubSelection(const std::string& tag_, const std::string& description_):
      tag(tag_), description(description_){}
    SubSelection(){}
  };

  std::string className_;
  
  //Event tree index within the input TChain
  Long64_t treeNum_;
  
  TTree* outEventTree_;
  TTree* outHeaderTree_;
  TTree* outDescriptionTree_;
  TTree* outBitFieldsTree_;
  //};
  ///@}

  /** To be field by the daughter class
   */
  ///@{
  /** List of supported subselection tags. To be
   * filled in the init(TTree*) method if several
   * subselections are supported.
   */
  std::vector<SubSelection> subSelections_;
  ///@}

  
  ShearsTChain chain_;
  std::auto_ptr<TFile> fout_;
  TDirectory* foutDir_;

  //  TreeRcd eventTree_;

  Int_t maxEvents_;
  Int_t skipEvents_;

  Int_t ievent_;

  /** Primary dataset of the processed data.
   * This information is provided as argument
   * of the constructor or create() method
   */
  std::string primaryDataset_;

  /** Type of selection. A same filter (Pruner class daughter)
   * can support several selection flavours, whose list must be
   * declared by the declareSubSelections() hook function. 
   * This index refers to the subSelections_ list.
   */
  size_t iSubSelection_;
  
  /** Pointer the run number of last read event
   */
  UInt_t* runNum_;

  /** Pointer the event number of last read event
   */
  UInt_t* eventNum_;

  /** Pointer to the event weight vector
   */
  std::vector<Double_t>* evtWeights_;

  /** Sums of input event weights
   */
  std::vector<Double_t> evtWeightSums_;

  /** Sums of output event weights
   */
  std::vector<Double_t> passedEvtWeightSums_;

int npNLO_;
  std::vector<float> LHEZPx_;
  std::vector<float> LHEZPy_;
float LHEZPt_; 

  std::vector<float> npNLOBinnedEvtWeightSums_;

  std::vector<float> passednpNLOBinnedEvtWeightSums_;
 
  std::vector<float> LHEZPtBinnedEvtWeightSums_;
  std::vector<float> passedLHEZPtBinnedEvtWeightSums_;
 
  TH1F * hnpNLO_;
  TH1F * hWeightednpNLO_;
  TH1F * hLHEZPt_;
  TH1F * hWeightedLHEZPt_;
  /** Current input file base name
   */
  std::string fileBaseName_;
  
  struct EventRcd{
    EventRcd(unsigned run_, unsigned event_, const std::string& filename_):
      run(run_), event(event_), filename(filename_) {}
    unsigned run;
    unsigned event;
    std::string filename;
    bool operator<(const EventRcd&a){
      if(run == a.run) return event < a.event;
      else return run < a.run;
    }
    bool operator==(const EventRcd& a){
      return run == a.run && event == a.event
	&& (filename.size() == 0 || a.filename.size() == 0 || filename == a.filename);
    }
  };
    
  std::vector<EventRcd> eventList_;

  std::vector<std::string> branchList_;

  static const int RUN_OFFSET = 32;

  bool allEvent_;

  //number of copied events:
  int nCopied_;

  //number of processed events:
  int nRead_;
  
public:
  struct ClassRecord;
  struct Plugin;

  int verbose_;
    
  virtual ~Pruner(){
  }

  /** Constructor
   */
  Pruner(const char* subSelection = 0,
	 const char* primary_dataset =0): treeNum_(-1), outEventTree_(0), outHeaderTree_(0),
					  outDescriptionTree_(0), outBitFieldsTree_(0),
					  chain_("tupel/EventTree"),
					  foutDir_(0), maxEvents_(-1), skipEvents_(0),
					  ievent_(-1), runNum_(0), eventNum_(0),
					  allEvent_(true), nCopied_(0), verbose_(0){
    chain_.SetDirectory(0);
    
    if(subSelection){
      setSubSelection(subSelection);
    }
    
    if(primary_dataset){
      primaryDataset_ = primary_dataset;
    }
  }

  /**
   * Loads a pruner from a shared object.
   * @param path The path of the .so file to load the pruner from.
   */
  static void load(const Plugin &plugin);
  
  /** Finds the ClassRecord that corresponds to the requested className
   * @param className, name of the selection class. It
   * must be a class inherited from Pruner and
   * registered with the macro DECLARE_PRUNER(class, description)
   */
  static ClassRecord *find(const std::string &className){
    //TODO: all selectors are created on registration.
    //Memory footprint can be optimised by creating the
    //instance on demand.
    std::map<std::string, Pruner::ClassRecord>::iterator res = daughtersMap().find(className);
    if (res != daughtersMap().end()) {
      return &res->second;
    } else {
      std::cerr << "Selection " << className << " was not found. Available selections can be listed with the option --list-selections\n";
      std::exit(EXIT_FAILURE);
    }
  }

  /** Creates a Pruner instance
   * @param className, name of the selection class. It
   * must be a class inherited from Pruner and
   * registered with the macro DECLARE_PRUNER(class, description)
   * @param subSelection selection tag to select a paritcular selection flavour
   * of a selection class
   * @param primaryDataset primary dataset name of the data to process, which
   * can also be used by the selection class to use different selections depending
   * on the dataset.
   */
  static Pruner* create(const char* className, const char* subSelection = 0, const char* primary_dataset = 0){
    Pruner *pruner = className != nullptr ? find(className)->instance : new Pruner;

    if(pruner){
      pruner->declareSubSelections();
      if(subSelection){
	if(!pruner->setSubSelection(subSelection)){
	  std::cerr << "Subselection " << subSelection << " was not found. Available subselections can be listed with the option --list-selections\n";
	  pruner = 0;
	}
      }
      if(primary_dataset) pruner->primaryDataset_ = primary_dataset;
    }
    return pruner;
  }

//  /** Limits the number of events to copy.
//   * @param val maximum number of events
//   */
//  void setMaxEvents(int val) { maxEvents_ = val; }
//
//  /** Sets the number of events to skip.
//   * @param val number of event to skip.
//   */
//  void setSkipEvents(int val) { skipEvents_ = val; }

  /** List the events in the format which can be read back
   * the specify the list of events to copy.
   * @param o output stream to write the list.
   * @param nInputFiles number of ROOT files to read.
   * @param inputFiles ROOT files to read the events from.
   */
  void listEvents(std::ostream& o, size_t nInputFiles,
		  const char* const inputFiles[]);


  /** List the events in the format which can be read back
   * the specify the list of events to copy.
   * @param o output stream to write the list.
   * @param catalog Shears catalog of the input files.
   */
  void listEvents(std::ostream& o, const char* catalog);


  /** List the event tree branched in the format which can
   * be read back the specify the list of branch to copy.
   * @param o output stream to write the list.
   * @param inputDataFile ROOT file to read the events from.
   */
  void listBranches(std::ostream& o, const char* inputDataFile);

  /** List the event tree branched in the format which can
   * be read back the specify the list of branch to copy.
   * @param o output stream to write the list.
   * @param catalog Shears catalog. First file of the catalog is
   * used to extract the lis of branches,
   */
  void listBranchesFromCat(std::ostream& o, const char* catalog);

  /** List available selections: Pruner modules and supported
   * subselectoins
   */
  void listSelections(std::ostream& o);
  
  void fillPerInputSummary();

  void fillGlobalSummary();
  
  /** Perform the event copy from multiple files.
   * @param nInputs number of input files.
   * @param inputDataFiles input files the events must be read from.
   * @param outputDataFile ouput file the events must be written to.
   * @maxEvents maximum number of events to process. The value -1 indicates
   * to process all events.
   * @skipEvents number of events to skip. Processing will start from
   * (skipEvents + 1) th event.
   */  
  void run(size_t nInputs, const char* const inputDataFiles[], const char* outputDataFile,
	   int maxEvents = -1, int skipEvents = 0);

  /** Perform the event copy from multiple files.
   * @param catalogFile catalog containing the list of files to process.
   * @param outputDataFile ouput file the events must be written to.
   * @maxEvents maximum number of events to process. The value -1 indicates
   * to process all events.
   * @skipEvents number of events to skip in addition to the events of the
   * first skipFiles files.
   * @maxFiles maximum number of files to process.
   * @skipFiles number of files to skip.
   */
  void run(const char* catalogFile, const char* outputDataFile,
	   int maxEvents = -1, int skipEvents = 0,
	   int maxFiles = -1, int skipFiles = 0);
  
  /** Perform the event copy.
   * @param inputDataFile input file the events must be read from.
   * @param outputDataFile ouput file the events must be written to.
   */
  //  void run(char* inputDataFile, char* outputDataFile){
  //    run(1, &inputDataFile, outputDataFile);
  // }

  /** Read the list of events to copy from a text file.
   * If this method is not called, then all the branches are copied.
   * @fileName path to the text file.
   */
  virtual void readEventList(const char* fileName);

  /** Read the list of branches to not copy from a text file.
   * @fileName path to the text file.
   */
  void readBranchList(const char* fileName);


  /** Sets subselection.
   * @param subselection subselection name
   * @return true iff the subselection was found
   */
  virtual bool setSubSelection(const std::string &subSelection){
    for(size_t i = 0; i < subSelections_.size(); ++i){
      if(subSelection == subSelections_[i].tag){
        iSubSelection_ = i;
        return true;
      }
    }
    return false;
  }

  /** Sets message verbosity level
   * @param val verbosity level, 0 for the quiest mode
   */
  void setVerbosity(int val){ verbose_ = val; chain_.setVerbosity(val);}
  
  struct Plugin
  {
    std::string indexFile;
    std::string path;
  };

  struct ClassRecord{
    ClassRecord(): instance(0){}
    Plugin plugin;
    std::string className;
    std::string description;
    Pruner* instance;
  };

  // Returns the (unique) map of daughters.
  static std::map<std::string, Pruner::ClassRecord> &daughtersMap() {
    static std::map<std::string, Pruner::ClassRecord> instance;
    return instance;
  }
  
protected:
  /** Hoock methods to override in the derived classes.
   */
  ///@{
  
  /** Hook function, where a derived class supporting subselections
   * should fill the subSelections_ field with the list of offered
   * subselections. There is not need to override this method
   * if the class does not provide subselections.
   */
  virtual void declareSubSelections() { /*NOOP*/ }
  
  /** Called before the event loop. When implementing a selector inheriting
   * from Pruner class, the SetBranch() method should be called to set the
   * tree branch address and have access to them in the filterEvent() method.
   * If the class inherits from EventTree class generated with TTree::MakeClass
   * then EventTree::Init(tree) should be called here.
   * @param tree pointer to the EventTree tree.
   * @return true if the initialisation succeeded.
   */
  virtual bool init(TChain* tree){ return true;};

  
  /** Check if curren event must be used. Typically use eventNum_ and runNum_
   * fields which contain the event and the run numbers. Return true if event
   * must be copied
   */
  virtual bool filterEvent(){
    if(allEvent_) return true; //copy every event
    bool selected = (std::find(eventList_.begin(), eventList_.end(),
			       EventRcd(*runNum_, *eventNum_, fileBaseName_))
		     != eventList_.end());
    return selected;
  }

  /** Check if a branch must be copied. Return true if the branch must copied,
   * false otherwise.
   * Can be overridden in a derived class to customize the selection. The
   * default implementation is typically sufficient.
   */
  virtual bool filterBranch( const char* branchName){
    if(branchList_.size() == 0) return true;
    return std::find(branchList_.begin(), branchList_.end(),
		     std::string(branchName))
      != branchList_.end();
  }
  ///@}
  
  /** Utility functions to filter collections.
   */
  ///@{

  /** Filter a collection. This method remove elements from a vector.
   * The mask vector indicates which element to keep.
   * The method makeFilterMask can be used to build the mask vector.
   * @param pointer to the coll collection to filter. If the pointer
   * is null the function has no effect
   * @param mask flags indicating the vector elements to keep
   */
  template<typename T>
  void filter(std::vector<T>* coll, std::vector<bool> mask){
    if(coll==0 || coll->size() == 0) return;
    std::vector<T> tmp;
    tmp.reserve(coll->size());
    for(unsigned i = 0; i < mask.size(); ++i){
      if(mask[i]) tmp.push_back((*coll)[i]);
    }
    //efficient way to move tmp content into coll:
    coll->swap(tmp);
  }

  /** Build a mask vector to be provided to the filter function. A filter method
   * is provided. For each collection element index the filter method is called
   * and should return true iff the element should be kept.
   * @param filter. Pointer to the filter function. The function prototype must
   * be: bool filter(int)
   * @param mask. The vector of flags to set.
   */
  void makeFilterMask(bool (Pruner::*filter)(int), std::vector<bool>& mask){
    for(unsigned i = 0; i < mask.size(); ++i){
      mask[i] = (this->*filter)(i);
    }
  }
  ///@}
  
private:
    
  /** initialized event summary tree. This method is already called by
   * init() method.
   */
  void setEventSummaryTree();

  //  /** Initialize input and output files and trees
  //   */
  //  bool init(size_t nInputDataFiles, const char* const inputDataFiles[],
  //	    const char* outputDataFile = 0);


  /** Sets input files. Called by run(...) methods. Seed also setInput(const char* catalog)
   * @param nInputFiles number of input files
   * @param inputDataFiles list of input files
   * @return false in case of failure
   */
  bool setInput(size_t nInputFiles, const char* const inputDataFiles[]);

  /** Sets input files. Called by run(...) methods. Seed also setInput(const char* catalog)
   * @param catalog input file catalog
   * @param maxFiles maximum number of files to process
   * @param skipFiles number of files to skip. Processing will start from the (skipFiles+1) th file
   * of the catalog
   * @return false in case of failure
   */
  bool setInput(const char* catalog, int maxFiles = -1, int skipFiles = 0);

  /** Sets output file. Called by run(...) methods. Seed also setInput(const char* catalog)
   * @return true in case of success, false otherwise
   */
  bool setOutput(const char* outputDataFile);

  /** Links variables to the branches we need to access
   */
  void setBranchAdd();
  
  /** Copy an event.
   */
  void copyEvent();

  /** Pass to the next event.
   */
  bool nextEvent();


  /** Copies tree entries. Called by run(const char*) and run(size_t, const
   *  char*[], const char*).
   */
  void run();

  /** Lists the event tree branched in the format which can
   * be read back the specify the list of branch to copy.
   * Called by public listBranches(...) methods.
   * @param o output stream to write the list.
   */
  void listBranches(std::ostream& o);

  
  /** List the events in the format which can be read back
   * the specify the list of events to copy. Called by the
   * public listEvents(...) methods.
   * @param o output stream to write the list.
   */
  void listEvents(std::ostream& o);
  
};

#endif //PRUNER_H not defined