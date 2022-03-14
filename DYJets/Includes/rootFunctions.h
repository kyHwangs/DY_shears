class TAxis;
class TCanvas;

void draw_axis_labels(TAxis* a);

//void graph_draw_stairs(TGraphAsymmErrors* g, int maxPoints = 99999);
void graph_draw_stairs(TGraphAsymmErrors* g, double ymin, double ymax, bool vertToOut = true);

bool alignRanges(const TAxis* axref, TAxis* ax);

TCanvas* newTCanvas(const char* name, const char* title, Double_t w, Double_t h);
