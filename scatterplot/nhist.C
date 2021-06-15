#include <TH1F.h>

void nhist() {
    const int N1 = 24777;
    const int N2 = 85731;
    const int N3 = 9501;
    const float zplus = -0.58740105f;
    const float zminus = 0.58740105f;
    const int bin_count = 3;

    float edges[bin_count + 1] = {-1, zplus, zminus, 1};

    TH1F* hist = new TH1F("a", "N distribution", bin_count, edges);

    for (int i = 0; i < N1; i++) hist->Fill(-0.9f);
    for (int i = 0; i < N2; i++) hist->Fill( 0.0f);
    for (int i = 0; i < N3; i++) hist->Fill( 0.9f);

    hist->GetYaxis()->SetRangeUser(0, 90000);
    hist->Draw();
}
