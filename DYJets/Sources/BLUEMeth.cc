#include "BLUEMeth.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "TDecompSVD.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TVectorD.h"

using namespace std;

// ClassImp(BLUEMeth);

BLUEMeth::BLUEMeth(const vector<TH1 *> &measurements,
                   const vector<vector<TH2 *>> &covariances,
                   const char *name,
                   const char *title)
    : TNamed(name, title)
{
    _variable = name;
    Init();
    Setup(measurements, covariances);
}

// BLUEMeth::BLUEMeth(const vector<TVectorD*> &meas, const vector<vector<TMatrixD*>> &cov, const
// char* name, const char* title)
//    : TNamed (name, title)
//{
//    Init();
//    Setup(meas, cov);
//}
BLUEMeth::BLUEMeth(const BLUEMeth &rhs) : TNamed(rhs.GetName(), rhs.GetTitle())
{
    // Copy constructor.
    Init();
    CopyData(rhs);
}

BLUEMeth *BLUEMeth::Clone(const char *newname) const
{
    BLUEMeth *blue = new BLUEMeth(*this);
    if (newname) blue->SetName(newname);
    return blue;
}

void BLUEMeth::Destroy() {}

void BLUEMeth::Assign(const BLUEMeth &rhs)
{
    if (this == &rhs) return;
    Reset();
    SetNameTitle(rhs.GetName(), rhs.GetTitle());
    CopyData(rhs);
}

void BLUEMeth::CopyData(const BLUEMeth &rhs) { Setup(rhs.measurements(), rhs.covariances()); }

void BLUEMeth::Reset()
{
    Destroy();
    Init();
}

void BLUEMeth::Init()
{
    _N = _n = _nMeasurements = _nCovariances = 0;
    _diagCrossChannelCov = _fullCrossChannelCov = _fullIndivChannelCov = false;
}

BLUEMeth &BLUEMeth::Setup(const vector<TH1 *> &measurements,
                          const vector<vector<TH2 *>> &covariances)
{
    Reset();
    SetMeasurements(measurements);
    SetCovariances(covariances);
    return *this;
}

void BLUEMeth::SetMeasurements(const vector<TH1 *> &measurements)
{
    _measurements = measurements;
    if (!_nMeasurements) {
        _nMeasurements = _measurements.size();
    } else if (_nMeasurements != _measurements.size()) {
        cerr << "Error: measurements size is different than covariances size." << endl;
    }

    if (!_N) {
        _N = _measurements[0]->GetNbinsX();
        _n = _N * _nMeasurements;
    } else if (_N != (unsigned int)_measurements[0]->GetNbinsX()) {
        cerr << "Error: number of observable size is different in measurement and in covariance."
             << endl;
    }

    //--- now we know the number of experiments (number of channels) : _nMeasurements
    //                the number of observables (number of bins) :     _N
    //                the size of the _yi :                            _n = _N * _nMeasurements
    _xa.ResizeTo(_N);
    _yi.ResizeTo(_n);
    _Uia.ResizeTo(_n, _N);

    //--- build up the full measurement _yi ---
    for (unsigned int i = 0; i < _nMeasurements; ++i) {
        for (unsigned int j = 0; j < _N; ++j) {
            //--- watch out the +1 in the GetBinContent method ---
            _yi(_N * i + j) = measurements[i]->GetBinContent(j + 1);
        }
    }

    //--- build up the design matrix _Uia ---
    for (unsigned int i = 0; i < _N; ++i) {
        _Uia(i, i) = 1.;
        _Uia(i + _N, i) = 1.;
    }
}

void BLUEMeth::SetCovariances(const vector<vector<TH2 *>> &covariances)
{
    _covariances = covariances;
    _nCovariances = _covariances[0].size();
    if (!_nMeasurements) {
        _nMeasurements = _covariances.size();
    } else if (_nMeasurements != _covariances.size()) {
        cerr << "Error: covariances size is different than measurement size." << endl;
    }
    if (!_N) {
        _N = _covariances[0][0]->GetNbinsX();
        _n = _N * _nMeasurements;
    } else if (_N != (unsigned int)_covariances[0][0]->GetNbinsX()) {
        cerr << "Error: number of observable size is different in measurement and in covariance."
             << endl;
    }

    //--- now we know the number of experiments (number of channels) : _nMeasurements
    //                the number of observables (number of bins) :     _N
    //                the size of the _Mij :                           _n = _N * _nMeasurements
    _Muij = vector<TMatrixD>{_nCovariances, TMatrixD(_n, _n)};
    _Mij.ResizeTo(_n, _n);
    _covuxaxb = vector<TMatrixD>{_nCovariances, TMatrixD(_N, _N)};
    _covxaxb.ResizeTo(_N, _N);
}

void BLUEMeth::ComputeFullCovariance()
{
    //--- loop on each source of systematic (each covariance matrix) ---
    for (unsigned int iCov = 0; iCov < _nCovariances; ++iCov) {
        //--- loop on each channel (each measurement of a particular observable) ---
        for (unsigned int iMeas = 0; iMeas < _nMeasurements; ++iMeas) {
            //--- loop on the _N bins (observables) and fill the _n X _n matrix ---
            for (unsigned int i = 0; i < _N; ++i) {
                for (unsigned int j = 0; j < _N; ++j) {
                    //--- only diagonal if not _fullIndivChannelCov
                    if (i == j || _fullIndivChannelCov) {
                        //--- watch out the +1 in the GetBinContent method ---
                        if (_covariances[iMeas][iCov]) {
                            _Muij[iCov](i + iMeas * _N, j + iMeas * _N) =
                                _covariances[iMeas][iCov]->GetBinContent(i + 1, j + 1);
                        } else {
                            _Muij[iCov](i + iMeas * _N, j + iMeas * _N) = 0;
                        }
                    }
                }
            }
        }
        //--- _Muij is now a block diagonal matrix containing the covariance of each channel ---
        // for fully correlated systematics (JES, JER, PU, Lumi, XSec) between channels
        // the off diagonal blocks are filed with sigma(i, j+N) = sign(sigma(i,j)) * sqrt(sigma(i,i)
        // * sigma(j+N,j+N))
        // to exhibit the same behavior than individual channel covariance matrix
        TString syst = _covariances[0][iCov] ? _covariances[0][iCov]->GetName() : "";
        if ((syst.Contains("JES") || syst.Contains("JER") || syst.Contains("PU") ||
             syst.Contains("Lumi") || syst.Contains("XSec")) &&
            (_diagCrossChannelCov || _fullCrossChannelCov)) {
            //--- loop on each pair of channels to add correlation between each pair of measurements
            //---
            for (unsigned int iMeas1 = 0; iMeas1 < _nMeasurements; ++iMeas1) {
                for (unsigned int iMeas2 = iMeas1 + 1; iMeas2 < _nMeasurements; ++iMeas2) {
                    //--- loop on the _N bins (observables) and fill the _n X _n matrix ---
                    for (unsigned int i = 0; i < _N; ++i) {
                        for (unsigned int j = 0; j < _N; ++j) {
                            //--- only diagonal if not _fullCrossChannelCov ---
                            if (i == j || _fullCrossChannelCov) {
                                //--- compute the sign of the correlation
                                int sign = _Muij[iCov](i, j) > 0 ? 1 : -1;
                                //--- compute the covariance ---
                                _Muij[iCov](i + iMeas1 * _N, j + iMeas2 * _N) =
                                    sign * sqrt(_Muij[iCov](i + iMeas1 * _N, i + iMeas1 * _N) *
                                                _Muij[iCov](j + iMeas2 * _N, j + iMeas2 * _N));
                                //--- covariance is symmetric ---
                                _Muij[iCov](j + iMeas2 * _N, i + iMeas1 * _N) =
                                    _Muij[iCov](i + iMeas1 * _N, j + iMeas2 * _N);
                            }
                        }
                    }
                }
            }
        }
        //--- now sum up the different covariance matrices ---
        _Mij += _Muij[iCov];
    }
}

TH1 *BLUEMeth::GetCombination(bool diagCrossChannelCov,
                              bool fullCrossChannelCov,
                              bool fullIndivChannelCov,
                              bool modifiedSWA,
                              vector<TH2 *> &covuxaxb,
                              TH2 *&covxaxb)
{
    //--- it is time to compute the lambda coefficients of the BLUE ---

    TMatrixD *lambdaai = NULL;
    //--- setup the way we want to include the correlation between all measurements ---

    if (modifiedSWA) {
        _diagCrossChannelCov = false;
        _fullCrossChannelCov = false;
        _fullIndivChannelCov = false;
        //--- compute the _Muij and _Mij with given setting ---
        ComputeFullCovariance();
        //--- now calculate the lambda (eq. (8) in Valassi's paper)
        TMatrixD _MInvertij = _Mij;
        _MInvertij.Invert();

        TMatrixD KInvertab =
            TMatrixD(TMatrixD(_Uia, TMatrixD::kTransposeMult, _MInvertij), TMatrixD::kMult, _Uia);
        KInvertab.Invert();

        lambdaai = new TMatrixD(
            KInvertab, TMatrixD::kMult, TMatrixD(_Uia, TMatrixD::kTransposeMult, _MInvertij));
        _xa = (*lambdaai) * _yi;
    } else {
        _diagCrossChannelCov = diagCrossChannelCov;
        _fullCrossChannelCov = fullCrossChannelCov;
        _fullIndivChannelCov = fullIndivChannelCov;
        ComputeFullCovariance();
        //--- now calculate the lambda (eq. (8) in Valassi's paper)
        TMatrixD _MInvertij = _Mij;
        _MInvertij.Invert();

        TMatrixD KInvertab =
            TMatrixD(TMatrixD(_Uia, TMatrixD::kTransposeMult, _MInvertij), TMatrixD::kMult, _Uia);
        KInvertab.Invert();

        lambdaai = new TMatrixD(
            KInvertab, TMatrixD::kMult, TMatrixD(_Uia, TMatrixD::kTransposeMult, _MInvertij));
        _xa = (*lambdaai) * _yi;
    }

    _diagCrossChannelCov = diagCrossChannelCov;
    _fullCrossChannelCov = fullCrossChannelCov;
    _fullIndivChannelCov = fullIndivChannelCov;
    ComputeFullCovariance();

    TMatrixD _MInvertij = _Mij;
    _MInvertij.Invert();

    if (modifiedSWA) {
        _covxaxb = TMatrixD(
            (*lambdaai), TMatrixD::kMult, TMatrixD(_Mij, TMatrixD::kMultTranspose, (*lambdaai)));
    } else {
        _covxaxb =
            TMatrixD(TMatrixD(_Uia, TMatrixD::kTransposeMult, _MInvertij), TMatrixD::kMult, _Uia);
        _covxaxb.Invert();
    }
    covxaxb = new TH2D(_covxaxb);
    covxaxb->SetName("TotCombCov" + _variable);

    //    std::cout << "HERE IS THE LAMBDA" << std::endl;
    //
    //    lambdaai->Print();
    //
    //    std::cout << "END OF THE LAMBDA" << std::endl;

    for (unsigned int iCov = 0; iCov < _nCovariances; ++iCov) {
        if (modifiedSWA) {
            _covuxaxb[iCov] =
                TMatrixD((*lambdaai),
                         TMatrixD::kMult,
                         TMatrixD(_Muij[iCov], TMatrixD::kMultTranspose, (*lambdaai)));
        } else {
            TMatrixD _MuInvertij = _Muij[iCov];
            TDecompSVD svdMuij = _MuInvertij;
            _MuInvertij = svdMuij.Invert();
            _covuxaxb[iCov] = TMatrixD(
                TMatrixD(_Uia, TMatrixD::kTransposeMult, _MuInvertij), TMatrixD::kMult, _Uia);
            TDecompSVD svdcovuxaxb = _covuxaxb[iCov];
            _covuxaxb[iCov] = svdcovuxaxb.Invert();
        }
        if (_covariances[0][iCov]) {
            covuxaxb[iCov] = new TH2D(_covuxaxb[iCov]);
            covuxaxb[iCov]->SetName(_covariances[0][iCov]->GetName() + _variable);
        } else {
            covuxaxb[iCov] = 0;
        }
    }

    TH1 *h = (TH1 *)_measurements[0]->Clone("hCombination");
    for (unsigned int i = 0; i < _N; ++i) {
        h->SetBinContent(i + 1, _xa(i));
        h->SetBinError(i + 1, sqrt(_covuxaxb[0](i, i)));
    }

    delete lambdaai;
    return h;
}
