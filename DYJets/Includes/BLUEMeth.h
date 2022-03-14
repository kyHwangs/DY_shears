#ifndef _BLUEMETH_H_
#define _BLUEMETH_H_

#include "TMatrixD.h"
#include "TNamed.h"
#include "TString.h"
#include "TVectorD.h"
#include <vector>

class TH1;
class TH2;

using std::vector;

class BLUEMeth : public TNamed
{

  public:
    // Standard methods

    BLUEMeth();                                          // default constructor
    BLUEMeth(const char *name, const char *title);       // named constructor
    BLUEMeth(const TString &name, const TString &title); // named constructor
    BLUEMeth(const BLUEMeth &rhs);                       // copy constructor
    ~BLUEMeth();                                         // destructor
    BLUEMeth &operator=(const BLUEMeth &rhs);            // assignment operator
    BLUEMeth *Clone(const char *newname = 0) const;

    // Special constructors

    BLUEMeth(const vector<TH1 *> &measures,
             const vector<vector<TH2 *>> &covariances,
             const char *name = 0,
             const char *title = 0);
    // BLUEMeth (const vector<TVectorD*> &measures, const vector<vector<TMatrixD*>> &covariances,
    // const char* name= 0, const char* title= 0);

    // Set up an existing object

    BLUEMeth &Setup(const vector<TH1 *> &measurements, const vector<vector<TH2 *>> &covariances);
    void Reset();
    void SetMeasurements(const vector<TH1 *> &measurements);
    void SetCovariances(const vector<vector<TH2 *>> &covariances);

    // Accessors

    const vector<TH1 *> &measurements() const;
    const vector<vector<TH2 *>> &covariances() const;

    TH1 *GetCombination(bool diagCrossChannelCov,
                        bool fullCrossChannelCov,
                        bool fullIndivChannelCov,
                        bool modifiedSWA,
                        vector<TH2 *> &covuxaxb,
                        TH2 *&covxaxb);

  private:
    void Init();
    void Destroy();
    void CopyData(const BLUEMeth &rhs);
    void ComputeFullCovariance();

  protected:
    TString _variable;
    unsigned int _N; // number of observables (number of bins in each channel)
    unsigned int _n; // number of measurements (usually 2*_N, but could be different)
    unsigned int _nMeasurements;
    unsigned int _nCovariances;
    bool _diagCrossChannelCov;
    bool _fullCrossChannelCov;
    bool _fullIndivChannelCov;

    vector<TH1 *> _measurements;        // experimental results for each channel (not owned)
    vector<vector<TH2 *>> _covariances; // covariances on the experimental results for each channel
    TVectorD _yi; // vector of the _n experimental results (usually the cross section measurements)
    // vector<TMatrixD> _covyiyj; // std::vector of the nSyst (_n X _n) TMatrixD covariance of _yi
    vector<TMatrixD> _Muij;     // (_n X _n) covariance matrices of _yi
    TMatrixD _Mij;              // (_n X _n) covariance matrices of _yi
    TVectorD _xa;               // BLUE estimate
    TMatrixD _Uia;              // design matrix
    vector<TMatrixD> _covuxaxb; // covariances for the linear estimate
    TMatrixD _covxaxb;          // covariance for the linear estimate

  protected:
    void Assign(const BLUEMeth &rhs); // implementation of assignement operator

  public:
    // ClassDef (BLUEMeth, 1)
};

//==============================================================================
// Inline method definitions
//==============================================================================

inline BLUEMeth::BLUEMeth() : TNamed()
{
    // Default constructor. Use Setup() to prepare for combination.
    Init();
}

inline BLUEMeth::BLUEMeth(const char *name, const char *title) : TNamed(name, title)
{
    // Basic named construcor. Use Setup() to prepare for combination.
    Init();
}

inline BLUEMeth::BLUEMeth(const TString &name, const TString &title) : TNamed(name, title)
{
    // Basic named construcor. Use Setup() to prepare for combination.
    Init();
}

inline BLUEMeth::~BLUEMeth() { Destroy(); }

inline BLUEMeth &BLUEMeth::operator=(const BLUEMeth &rhs)
{
    // Assignemnt operator or copynig BLUEMeth settings
    Assign(rhs);
    return *this;
}

inline const vector<TH1 *> &BLUEMeth::measurements() const
{
    // Return vector of individual experimental results
    return _measurements;
}

inline const vector<vector<TH2 *>> &BLUEMeth::covariances() const
{
    // Return vector of individual experimental results
    return _covariances;
}

#endif
