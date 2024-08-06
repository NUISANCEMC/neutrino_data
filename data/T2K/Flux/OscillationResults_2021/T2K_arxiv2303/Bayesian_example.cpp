TH1D *GetCredible(TH1D *h, double interval) {
  if (interval < 0 || interval > 1) {
    std::cerr << "Requested interval (" << interval << ") < 0 or > 1, you're asking for the impossible" << std::endl;
    return (TH1D*)NULL;
  }

  double integral = h->Integral();
  double sum = 0.0;
  TH1D *hInt = (TH1D*)h->Clone(Form("%s_h68", h->GetTitle()));
  hInt->Reset();
  TH1D *hCopy = (TH1D*)h->Clone(Form("%s_hCopy", h->GetTitle()));
  while (sum/integral < interval) {
    int maxbin = hCopy->GetMaximumBin();
    sum += hCopy->GetBinContent(maxbin);
    hInt->SetBinContent(maxbin, hCopy->GetBinContent(maxbin));
    hCopy->SetBinContent(maxbin, 0);
  }

  delete hCopy;
  return hInt;
}

double FindMaximumX(TGraph *g) {
  double maxx = g->GetX()[0];
  for (int i = 0; i < g->GetN(); ++i) {
    if (g->GetX()[i] > maxx) maxx = g->GetX()[i];
  }
  return maxx;
}

double FindMaximumY(TGraph *g) {
  double maxy = g->GetY()[0];
  for (int i = 0; i < g->GetN(); ++i) {
    if (g->GetY()[i] > maxy) maxy = g->GetY()[i];
  }
  return maxy;
}

double FindMinimumX(TGraph *g) {
  double minx = g->GetX()[0];
  for (int i = 0; i < g->GetN(); ++i) {
    if (g->GetX()[i] < minx) minx = g->GetX()[i];
  }
  return minx;
}

double FindMinimumY(TGraph *g) {
  double miny = g->GetY()[0];
  for (int i = 0; i < g->GetN(); ++i) {
    if (g->GetY()[i] < miny) miny = g->GetY()[i];
  }
  return miny;
}

void Bayesian_example() {
  TFile *f1 = new TFile("Bayesian_DataRelease.root");

  TGaxis::SetMaxDigits(3);

  TCanvas *canv = new TCanvas("canv", "canv", 1000, 1000);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.5);
  canv->SetRightMargin(canv->GetRightMargin()*0.9);
  TString output = f1->GetName();
  output.ReplaceAll(".root", "");
  output+="_vis.pdf";
  canv->Print(output+"[");

  TPaveText *p = new TPaveText(0,0,1,1, "NDC");
  p->AddText(f1->GetName());
  p->Draw();
  canv->Print(output);

  TIter next(f1->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)next())) {
    std::string name = std::string(key->GetName());
    std::string classname = key->GetClassName();
    if (classname == "TGraph") {
      TLegend *legend = new TLegend(0.25, 0.89, 0.93, 0.95);
      legend->SetNColumns(5);
      TGraph *g1 = NULL;
      TGraph *g2 = NULL;
      g1 = (TGraph*)f1->Get(name.c_str());
      g1->SetName(Form("%s_Old", name.c_str()));

      g1->SetTitle(Form("%s", name.c_str()));

      TGraph *gbf = NULL;

      TGraph *g68 = NULL;
      TGraph *g90 = NULL;
      TGraph *g955 = NULL;
      TGraph *g997 = NULL;

      // see if we can find the other discontinuous region
      TGraph *g68_2 = NULL;
      TGraph *g90_2 = NULL;
      TGraph *g955_2 = NULL;
      TGraph *g997_2 = NULL;

      if (name.find("cred68_0") != std::string::npos) {
        g68 = g1;
        std::string subname = name.substr(0, name.find("cred68"));
        gbf =  (TGraph*)f1->Get((subname + "bestfit").c_str());
        g90 =  (TGraph*)f1->Get((subname + "cred90_0").c_str());
        g955 = (TGraph*)f1->Get((subname + "cred955_0").c_str());
        g997 = (TGraph*)f1->Get((subname + "cred997_0").c_str());

        g68_2 =  (TGraph*)f1->Get((subname + "cred68_1").c_str());
        g90_2 =  (TGraph*)f1->Get((subname + "cred90_1").c_str());
        g955_2 = (TGraph*)f1->Get((subname + "cred955_1").c_str());
        g997_2 = (TGraph*)f1->Get((subname + "cred997_1").c_str());
      } else {
        continue;
      }

      g997->SetLineColor(kGray+1);
      g997->SetLineStyle(kDashed);
      g955->SetLineColor(kBlack);
      g955->SetLineStyle(kDashed);
      g90->SetLineColor(kGray+1);
      g68->SetLineColor(kBlack);

      double maxx = FindMaximumX(g997);
      double maxy = FindMaximumY(g997);
      double minx = FindMinimumX(g997);
      double miny = FindMinimumY(g997);
      maxx < FindMaximumX(g955) ? maxx = FindMaximumX(g955) : maxx = maxx;
      maxy < FindMaximumY(g955) ? maxy = FindMaximumY(g955) : maxy = maxy;
      minx > FindMinimumX(g955) ? minx = FindMinimumX(g955) : minx = minx;
      miny > FindMinimumY(g955) ? miny = FindMinimumY(g955) : miny = miny;

      maxx < FindMaximumX(g90) ? maxx = FindMaximumX(g90) : maxx = maxx;
      maxy < FindMaximumY(g90) ? maxy = FindMaximumY(g90) : maxy = maxy;
      minx > FindMinimumX(g90) ? minx = FindMinimumX(g90) : minx = minx;
      miny > FindMinimumY(g90) ? miny = FindMinimumY(g90) : miny = miny;

      maxx < FindMaximumX(g68) ? maxx = FindMaximumX(g68) : maxx = maxx;
      maxy < FindMaximumY(g68) ? maxy = FindMaximumY(g68) : maxy = maxy;
      minx > FindMinimumX(g68) ? minx = FindMinimumX(g68) : minx = minx;
      miny > FindMinimumY(g68) ? miny = FindMinimumY(g68) : miny = miny;

      if (g997_2 != NULL) {
        maxx < FindMaximumX(g997_2) ? maxx = FindMaximumX(g997_2) : maxx = maxx;
        maxy < FindMaximumY(g997_2) ? maxy = FindMaximumY(g997_2) : maxy = maxy;
        minx > FindMinimumX(g997_2) ? minx = FindMinimumX(g997_2) : minx = minx;
        miny > FindMinimumY(g997_2) ? miny = FindMinimumY(g997_2) : miny = miny;
      }

      if (g955_2 != NULL) {
        maxx < FindMaximumX(g955_2) ? maxx = FindMaximumX(g955_2) : maxx = maxx;
        maxy < FindMaximumY(g955_2) ? maxy = FindMaximumY(g955_2) : maxy = maxy;
        minx > FindMinimumX(g955_2) ? minx = FindMinimumX(g955_2) : minx = minx;
        miny > FindMinimumY(g955_2) ? miny = FindMinimumY(g955_2) : miny = miny;
      }

      if (g90_2 != NULL) {
        maxx < FindMaximumX(g90_2) ? maxx = FindMaximumX(g90_2) : maxx = maxx;
        maxy < FindMaximumY(g90_2) ? maxy = FindMaximumY(g90_2) : maxy = maxy;
        minx > FindMinimumX(g90_2) ? minx = FindMinimumX(g90_2) : minx = minx;
        miny > FindMinimumY(g90_2) ? miny = FindMinimumY(g90_2) : miny = miny;
      }

      if (g68_2 != NULL) {
        maxx < FindMaximumX(g68_2) ? maxx = FindMaximumX(g68_2) : maxx = maxx;
        maxy < FindMaximumY(g68_2) ? maxy = FindMaximumY(g68_2) : maxy = maxy;
        minx > FindMinimumX(g68_2) ? minx = FindMinimumX(g68_2) : minx = minx;
        miny > FindMinimumY(g68_2) ? miny = FindMinimumY(g68_2) : miny = miny;
      }

      double xarray[2] = {minx, maxx};
      double yarray[2] = {miny, maxy};
      TGraph *g = new TGraph(2, xarray, yarray);

      legend->AddEntry(g68, "68%", "l");
      legend->AddEntry(g90, "90%", "l");
      legend->AddEntry(g955, "95.5%", "l");
      legend->AddEntry(g997, "99.7%", "l");

      if (name.find("dm32") == std::string::npos) {
        legend->SetX1(0.15);
      }

      canv->Clear();
      std::string name2 = name.substr(0, name.find("_cred68"));
      name2 = name2.substr(5, name2.size());
      g->SetTitle(name2.c_str());
      g->SetMarkerSize(0);
      g->Draw("A*");
      // Check x-axis title
      if (name2.find("th23_") == 0) g->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
      if (name2.find("th13_") == 0) g->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
      if (name2.find("th13_") == 0) g->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
      if (name2.find("jarlskog_") == 0) g->GetXaxis()->SetTitle("J=s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta_{CP}");
      // Check y-axis title
      if (name2.find("_dCP_") != std::string::npos || name2.find("_dcp_") != std::string::npos) g->GetYaxis()->SetTitle("#delta_{CP}");
      if (name2.find("_dm32_") != std::string::npos) {
        g->GetYaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2})");
        g->GetYaxis()->SetTitleOffset(g->GetYaxis()->GetTitleOffset()*1.4);
      }
      if (name2.find("_th23_") != std::string::npos) g->GetYaxis()->SetTitle("sin^{2}#theta_{23}");

      g997->Draw("C,same");
      g955->Draw("C,same");
      g90->Draw("C,same");
      g68->Draw("C,same");

      if (g68_2 != NULL) {
        g68_2->SetLineColor(g68->GetLineColor());
        g68_2->SetLineStyle(g68->GetLineStyle());
        g68_2->Draw("C,same");
      }
      if (g90_2 != NULL) {
        g90_2->SetLineColor(g90->GetLineColor());
        g90_2->SetLineStyle(g90->GetLineStyle());
        g90_2->Draw("C,same");
      }
      if (g955_2 != NULL) {
        g955_2->SetLineColor(g955->GetLineColor());
        g955_2->SetLineStyle(g955->GetLineStyle());
        g955_2->Draw("C,same");
      }
      if (g997_2 != NULL) {
        g997_2->SetLineColor(g997->GetLineColor());
        g997_2->SetLineStyle(g997->GetLineStyle());
        g997_2->Draw("C,same");
      }

      if (gbf != NULL) {
        gbf->Draw("*,same");
        gbf->SetMarkerStyle(29);
        gbf->SetMarkerSize(3);
        legend->AddEntry(gbf, "Best fit", "p");
      }
      legend->Draw("same");

      canv->Print(output);

      delete legend;
      // Now draw the TH1D too
    } else if (classname == "TH1D") {

      TH1D *hPost = (TH1D*)f1->Get(name.c_str());
      hPost->SetLineColor(kBlack);
      hPost->Draw("hist");
      hPost->GetYaxis()->SetLabelSize(0.035);
      hPost->GetXaxis()->SetLabelSize(0.035);
      hPost->GetXaxis()->SetTitleSize(0.035);
      hPost->GetYaxis()->SetTitleSize(0.035);
      hPost->GetYaxis()->SetTitleOffset(hPost->GetYaxis()->GetTitleOffset()*1.8);
      hPost->GetXaxis()->SetTitleOffset(hPost->GetXaxis()->GetTitleOffset()*1.5);
      hPost->SetLineWidth(2);
      hPost->GetYaxis()->SetTitle("Posterior probability");
      std::string hpostname = key->GetName();
      if (hpostname.find("h1D_") != std::string::npos) hpostname = hpostname.substr(hpostname.find("h1D_")+4, hpostname.size());
      if (hpostname.find("posterior") != std::string::npos) {
        hpostname = hpostname.substr(0, hpostname.find("posterior"))+hpostname.substr(hpostname.find("posterior")+9, hpostname.size());
      }
      hPost->SetTitle(hpostname.c_str());
      if (name.find("dm32") != std::string::npos) hPost->GetXaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2})");

      // Draw an example 68% interval
      double integral = hPost->Integral();
      double sum = 0.0;
      TH1D *h68 = GetCredible(hPost, 0.68);
      TH1D *h955 = GetCredible(hPost, 0.955);

      h955->SetFillColor(kGray+1);
      h955->SetFillStyle(1001);
      h955->Draw("hist, same");

      h68->SetFillColor(kBlack);
      h68->SetFillStyle(1001);
      h68->Draw("hist, same");

      TLegend *legend2 = new TLegend(0.25, 0.90, 0.93, 0.95);
      legend2->SetNColumns(3);
      legend2->AddEntry(hPost, "Posterior", "l");
      legend2->AddEntry(h68, "68% cred.", "f");
      legend2->AddEntry(h955, "95.5% cred.", "f");
      legend2->Draw("same");

      canv->Print(output);

      delete legend2;
    } else {
      continue;
    }
  }
  canv->Print(output+"]");
}
