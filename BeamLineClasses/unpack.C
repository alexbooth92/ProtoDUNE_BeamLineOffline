#include "unpack.h"

void SetHistTitles(TH1I* &h, TString const &title, TString const &xTitle, TString const &yTitle)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);

  return;
}

void SetHistTitles(TH1D* &h, TString const &title, TString const &xTitle, TString const &yTitle)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);

  return;
}

void SetHistTitles(TH2D* &h, TString const &title, TString const &xTitle, TString const &yTitle)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);

  return;
}
