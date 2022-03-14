#include "TClass.h"
#include "TH1.h"
#include "TLegend.h"
#include "TList.h"
#include "TVirtualPad.h"
#include "math.h"
#include <iostream>
#include <limits>

void fixYscale(double linfact, double logfact, double logMaxRange)
{

    const int verbosity = 0;

    if (logfact <= 0) logfact = linfact;

    if (gPad == 0) return;

    double ymax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double yposmin = std::numeric_limits<double>::max();
    std::vector<TH1 *> hs;

    TLegend *tl = 0;

    if (verbosity > 0)
        std::cout << "Fixing y-scal for the pad with name: " << gPad->GetName() << "\n";
    TListIter it(gPad->GetListOfPrimitives());
    while (it.Next()) {
        TObject *obj = *it;
        if (TClass(obj->ClassName()).InheritsFrom("TLegend")) {
            if (verbosity > 0)
                std::cout << "Found a legend with object name " << obj->GetName() << "\n";
            tl = dynamic_cast<TLegend *>(obj);
        }
        if (!TClass(obj->ClassName()).InheritsFrom("TH1")) continue;
        TH1 *h = dynamic_cast<TH1 *>(obj);
        if (h == 0) continue;
        for (int ibin = 1; ibin <= h->GetXaxis()->GetNbins(); ++ibin) {
            double y = h->GetBinContent(ibin);
            double yerr = h->GetBinError(ibin);
            double yerr_low = 0;
            if (!gPad->GetLogy()) yerr_low = yerr; // ignore error bar for log scales
            if (y + yerr > ymax) ymax = y + yerr;
            if (y - yerr_low < ymin) {
                if (verbosity > 0) {
                    std::cout << ">>> " << h->GetDrawOption() << "\n";
                    std::cout << ">>> " << h->GetName()
                              << ": ibin/nbins, y, yerr, yerr_low = " << ibin << "/"
                              << h->GetNbinsX() << ", " << y << ", " << yerr << ", " << yerr_low
                              << "\n";
                }
                ymin = y - yerr_low;
            }
            if (y > 0 && y < yposmin) yposmin = y;
            hs.push_back(h);
        }
    }

    if (verbosity > 0) std::cout << "ymin = " << ymin << "  ymax = " << ymax;

    if (gPad->GetLogy()) {
        if (ymin <= 0) ymin = yposmin; // yposmin / 10.;
        ymin *= pow(ymax / ymin, -(logfact - 1) / 2.);
        double a = pow(ymax / ymin, logfact);
        ymax = ymin * a;
        if (logMaxRange > 0 && ymax / ymin > logMaxRange) {
            ymin = ymax / logMaxRange;
        }
        if (ymax < ymin / 10) ymax = 15 * ymin;
    } else {
        bool add_ymin_margin = true;
        if (ymin > 0) {
            ymin = 0.001 *
                   (ymax -
                    ymin); // starts at slightly higher than 0 to prevent the display of the label
            add_ymin_margin = false;
        }
        //    ymax *= linfact;
        // ymin *= linfact;
        ymax += (linfact - 1) * (ymax - ymin);
        if (add_ymin_margin) ymin -= (linfact - 1) * (ymax - ymin);
    }

    if (verbosity > 0)
        std::cout << "  ->  axis min = " << ymin << "      axis max = " << ymax << std::endl;

    // Prevent overlapping with the legend by zooming the y-axis.
    // Only legend place on top of the plot is handled.
    bool istopleg = false;
    if (tl &&
        (1 - std::max(tl->GetY1NDC(), tl->GetY2NDC())) < std::min(tl->GetY1NDC(), tl->GetY2NDC())) {
        istopleg = true;
    }
    if (verbosity > 0 && tl) std::cout << "Top legend: " << (istopleg ? "yes" : "no") << "\n";

    unsigned i = 0;
    if (hs.size() > i) {
        TH1 *h = (TH1 *)hs[i];
        h->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    gPad->Paint();

    if (istopleg) {
        double ymax_in_leg_area = -std::numeric_limits<double>::max();
        for (auto h : hs) {
            int lb = h->GetXaxis()->FindBin(tl->GetX1());
            int ub = h->GetXaxis()->FindBin(tl->GetX2());
            if (lb > ub) std::swap(lb, ub);
            if (lb == 0) lb = 1;
            if (ub > h->GetXaxis()->GetNbins()) ub = h->GetXaxis()->GetNbins();
            for (int ibin = lb; ibin <= ub; ++ibin) {
                double y = h->GetBinContent(ibin) + h->GetBinError(ibin);
                if (y > ymax_in_leg_area) ymax_in_leg_area = y;
            }
        }

        double y_leg_min = std::min(tl->GetY1NDC(), tl->GetY2NDC());
        if (gPad->GetLogy()) {
            y_leg_min = ymin * std::pow(ymax / ymin, y_leg_min);
        } else {
            y_leg_min = ymin + (ymax - ymin) * y_leg_min;
        }
        if (verbosity > 0)
            std::cout << "tl->GetY1NDC(), tl->GetY2NDC(): " << tl->GetY1NDC() << ","
                      << tl->GetY2NDC() << "\n";

        //    y_leg_min = gPad->PixeltoY(gPad->VtoPixel(y_leg_min) - gPad->GetWh());

        if (verbosity > 0)
            std::cout << "ymax_in_leg_area = " << ymax_in_leg_area << ", y_leg_min = " << y_leg_min
                      << "\n";

        if (y_leg_min < ymax_in_leg_area) {
            double oldymax = ymax;
            if (gPad->GetLogy()) {
                ymax =
                    ymin * pow(ymax / ymin, log(ymax_in_leg_area / ymin) / log(y_leg_min / ymin));
            } else {
                ymax = ymin + (ymax - ymin) * (ymax_in_leg_area - ymin) / (y_leg_min - ymin);
            }
            double margin = 0.03;
            ymax += margin * (ymax - ymin);
            std::cout << "Legend overlap protection, ymax: " << oldymax << " -> " << ymax << "\n";
        }
    }
    for (unsigned i = 0; i < hs.size(); ++i) {
        // i = 0;
        if (hs.size() > i) {
            TH1 *h = (TH1 *)hs[i];
            h->GetYaxis()->SetRangeUser(ymin, ymax);
        }
    }
    gPad->Paint();
}
