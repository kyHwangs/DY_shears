#include "TLatex.h"
#include "TAxis.h"
#include "TH1.h"
#include "TVirtualPad.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <algorithm>
#include <iostream>

void draw_axis_labels( TAxis* a){
   TLatex* txt_ = new TLatex();
   // draw axis labels
   //  float xlab_ = 65;
   txt_->SetTextFont(a->GetLabelFont());
   txt_->SetTextSize(a->GetLabelSize());
   txt_->SetTextAngle(0);


   
   TH1* h = (TH1*) a->GetParent();
   TAxis* ax = h->GetXaxis();
   TAxis* ay = h->GetYaxis();

   int iV = (a == ax) ? 1 : 0;   

   if( iV==1 ){
     txt_->SetTextAlign(23);
   } else {
     txt_->SetTextAlign(32);
   }

   TVirtualPad* pad = TVirtualPad::Pad();
   

   double ymin = pad->GetUymin();
   double ymax = pad->GetUymax();
   double xmin = pad->GetUxmin();//ax->GetBinLowEdge(ax->GetFirst());
   double xmax = pad->GetUxmax();//ax->GetBinUpEdge(ax->GetLast());
   double yoffset = (ax->GetLabelOffset() + 0.01)*(ymax-ymin);
   double xoffset = (ay->GetLabelOffset() + 0.01)*(xmax-xmin);
   
   //   TAxis* a = (iV == 1)? ax : ay;
   
   int imin = a->GetFirst();
   int imax = a->GetLast();
   
   double tickY0 = gPad->GetUymin();
   double tickY1 = tickY0  + 0.02*(gPad->GetUymax() - tickY0);
   TLine l;

   for(int ii=imin-1; ii<=imax; ++ii){
     double val = (ii==imin-1) ?
       a->GetBinLowEdge(imin)
       : a->GetBinUpEdge(ii);
     float x_, y_;
     if( iV==1 ){ //x-axis
       x_ = val;
       y_ = ymin - yoffset;
     } else {
       y_ = val;
       x_ = xmin - xoffset;
     }
     txt_->SetText(x_,y_,Form("%-2.0f",val));

     l.DrawLine(val, tickY0, val, tickY1);

     txt_->DrawClone();
   }
}

void graph_draw_stairs(TGraphAsymmErrors* g, double ymin, double ymax, bool vertToOut){
  //  TVirtualPad* pad = TVirtualPad::Pad();
  int n = g->GetN();
  //  double ymin = pad->GetUymin();
  //  double ymax = pad->GetUymax();
  std::vector<Double_t> x(3);
  std::vector<Double_t> y(3);
  for(int i = 0 ; i < n; ++i){
    x[0] = g->GetX()[i] - g->GetEXlow()[i];
    x[2] = x[1] = g->GetX()[i] + g->GetEXhigh()[i];
    y[0] = g->GetY()[i];
    bool nextBinOut = false;
    bool out = false;
    if(i + 1 < n 
       && (g->GetY()[i+1] > ymax || g->GetY()[i+1] < ymin)){
      nextBinOut = true;
    }
    if(!vertToOut && (y[0] < ymin || y[0] > ymax)) continue;
    if(y[0] > ymax) { y[0] = ymax; out = true; }
    if(y[0] < ymin) { y[0] = ymin; out = true; }
    y[1] = y[0];
    int np;
    if(i < n - 1 && (!nextBinOut || vertToOut)){
      y[2] = g->GetY()[i+1];
      if(y[2] > ymax) y[2] = ymax;
      if(y[2] < ymin) y[2] = ymin;
      np = 3;
    } else{
      np = 2;
    }

    TPolyLine* pl;
    if(out){
      pl = new TPolyLine(np-1, &x[1], &y[1]);
    } else{
      pl = new TPolyLine(np, &x[0], &y[0]);
    }
    pl->SetLineColor(g->GetLineColor());
    pl->SetLineStyle(g->GetLineStyle());
    pl->SetLineWidth(g->GetLineWidth());
    pl->Draw();    
  }
}

bool alignRanges(const TAxis* axref, TAxis* ax){
  double xmin = axref->GetBinLowEdge(axref->GetFirst());
  int first = 1;
  double eps = 1e-9;
  while(first < ax->GetNbins() + 1
	&& fabs(ax->GetBinLowEdge(first)-xmin)
	/ axref->GetBinWidth(axref->GetFirst()) > eps) ++first;
  if(first > ax->GetNbins() + 1) return false;
  int last = first + axref->GetLast() - axref->GetFirst();
  if(last > ax->GetNbins() + 1) return false;
  ax->SetRange(first, last);
  
  //check bin consistency:
  for(int i = 0; i <= axref->GetLast() - axref->GetFirst() + 1; ++i){
    if(fabs(ax->GetBinLowEdge(ax->GetFirst() + i)
	    - axref->GetBinLowEdge(axref->GetFirst() + i))
       / axref->GetBinWidth(axref->GetFirst() + i) > eps){
      return false;
    }
  }
  
  return true;
}

TCanvas* newTCanvas(const char* name, const char* title, Double_t w, Double_t h){
  TCanvas * c = new TCanvas("tmp", title, w, h);
  c->SetWindowSize(w, 0.5*h);
  c->Draw();
  //recipe of https://root.cern.ch/doc/master/classTCanvas.html to set window size
  //does  not work in batch mode (tested with version 6.06/01). As work around
  //we delete the canvas and instantiate a new one.
  w = w + c->GetWindowWidth() - c->GetWw();
  h = h + c->GetWindowHeight() - c->GetWh();
  delete c;
  
  c = new TCanvas(name, title, w, h);
  return c;
}
