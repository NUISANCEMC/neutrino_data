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

void Frequentist_example() {
  TFile *f1 = new TFile("Frequentist_DataRelease.root");

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

      std::string hierarchy;
      if (name.find("conf68_0") != std::string::npos) {
        g68 = g1;
        std::string subname = name.substr(0, name.find("conf68"));
        // The best-fit plot has a slightly different naming convention
        // cut out anything after "RC"
        std::string subname2 = name.substr(0, name.find("RC_")+3);
        // Then suffer some more to find if this is NH or IH
        // for the FC the best fit point is called slightly differently
        if (name.find("_FC_") == std::string::npos) {
          hierarchy = name.substr(name.find("RC_")+3, name.size());
        } else {
          if (name.find("_globalMH_") == std::string::npos) {
            hierarchy = name.substr(name.find("FC_")+3, name.size());
          } else {
            hierarchy = name.substr(name.find("globalMH_")+9, name.size());
          }
        }
        hierarchy = hierarchy.substr(0, hierarchy.find("_conf68"));
        subname2 += "conf"+hierarchy+"_bestfit";
        //gr2D_th23_dCP_wRC_FC_globalMH_IH_conf68_0

        gbf =  (TGraph*)f1->Get((subname2).c_str());
        g90 =  (TGraph*)f1->Get((subname + "conf90_0").c_str());
        g955 = (TGraph*)f1->Get((subname + "conf955_0").c_str());
        g997 = (TGraph*)f1->Get((subname + "conf997_0").c_str());

        g68_2 =  (TGraph*)f1->Get((subname + "conf68_1").c_str());
        g90_2 =  (TGraph*)f1->Get((subname + "conf90_1").c_str());
        g955_2 = (TGraph*)f1->Get((subname + "conf955_1").c_str());
        g997_2 = (TGraph*)f1->Get((subname + "conf997_1").c_str());
      } else {
        continue;
      }

      // The Feldman--Cousins corrected plots only have 68 and 90% computed
      double maxx = -999;
      double maxy = -999;
      double minx = 999;
      double miny = 999;
      if (g997 != NULL) {
        g997->SetLineColor(kGray+1);
        g997->SetLineStyle(kDashed);
        maxx < FindMaximumX(g997) ? maxx = FindMaximumX(g997) : maxx = maxx;
        maxy < FindMaximumY(g997) ? maxy = FindMaximumY(g997) : maxy = maxy;
        minx > FindMinimumX(g997) ? minx = FindMinimumX(g997) : minx = minx;
        miny > FindMinimumY(g997) ? miny = FindMinimumY(g997) : miny = miny;
      }
      if (g955 != NULL) {
        g955->SetLineColor(kBlack);
        g955->SetLineStyle(kDashed);
        maxx < FindMaximumX(g955) ? maxx = FindMaximumX(g955) : maxx = maxx;
        maxy < FindMaximumY(g955) ? maxy = FindMaximumY(g955) : maxy = maxy;
        minx > FindMinimumX(g955) ? minx = FindMinimumX(g955) : minx = minx;
        miny > FindMinimumY(g955) ? miny = FindMinimumY(g955) : miny = miny;
      }
      g90->SetLineColor(kGray+1);
      g68->SetLineColor(kBlack);


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

      if (name.find("_FC_") != std::string::npos) {
        legend->SetNColumns(4);
        legend->AddEntry((TObject*)(NULL), "Feldman-Cousins:", "");
      }
      if (name.find("globalMH") != std::string::npos) {
        legend->SetNColumns(5);
        legend->AddEntry((TObject*)(NULL), "Global #Delta(-2lnL)", "");
      }
      legend->AddEntry(g68, "68%", "l");
      legend->AddEntry(g90, "90%", "l");
      if (g955 != NULL) legend->AddEntry(g955, "95.5%", "l");
      if (g997 != NULL) legend->AddEntry(g997, "99.7%", "l");

      if (name.find("dm2") == std::string::npos) {
        legend->SetX1(0.15);
      }

      canv->Clear();
      std::string name2 = name.substr(0, name.find("_conf68"));
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
      if (name2.find("_dCP_") != std::string::npos) g->GetYaxis()->SetTitle("#delta_{CP}");
      if (name2.find("_dm2_") != std::string::npos) {
        if (hierarchy == "NH") g->GetYaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2})");
        if (hierarchy == "IH") g->GetYaxis()->SetTitle("|#Deltam^{2}_{31}| (eV^{2})");
        g->GetYaxis()->SetTitleOffset(g->GetYaxis()->GetTitleOffset()*1.4);
      }
      if (name2.find("_th23_") != std::string::npos) g->GetYaxis()->SetTitle("sin^{2}#theta_{23}");

      if (g997 != NULL) g997->Draw("C,same");
      if (g955 != NULL) g955->Draw("C,same");
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
      hPost->Draw("hist");
      hPost->GetYaxis()->SetLabelSize(0.035);
      hPost->GetXaxis()->SetLabelSize(0.035);
      hPost->GetXaxis()->SetTitleSize(0.035);
      hPost->GetYaxis()->SetTitleSize(0.035);
      hPost->GetYaxis()->SetTitleOffset(hPost->GetYaxis()->GetTitleOffset()*1.8);
      hPost->GetXaxis()->SetTitleOffset(hPost->GetXaxis()->GetTitleOffset()*1.5);
      hPost->SetLineWidth(2);
      hPost->GetYaxis()->SetTitle("#Delta(-2lnL)");
      hPost->GetYaxis()->SetRangeUser(0, hPost->GetMaximum()*1.05);
      std::string hpostname = key->GetName();
      if (hpostname.find("h1D_") != std::string::npos) hpostname = hpostname.substr(hpostname.find("h1D_")+4, hpostname.size());
      if (hpostname.find("chi2") != std::string::npos) {
        hpostname = hpostname.substr(0, hpostname.find("chi2"))+hpostname.substr(hpostname.find("chi2")+4, hpostname.size());
      }
      hPost->SetTitle(hpostname.c_str());

      TLegend *legend2 = new TLegend(0.75, 0.89, 0.93, 0.95);
      legend2->AddEntry(hPost, "#Delta(-2lnL)", "l");

      // Also draw the critical values if we have them
      TGraph *gCrit68 = NULL;
      TGraph *gCrit90 = NULL;
      TGraph *gCrit955 = NULL;
      TGraph *gCrit997 = NULL;
      if ((name.find("th23") != std::string::npos || name.find("dCP") != std::string::npos) && name.find("_wRC_") != std::string::npos) {
        // Do some name manipulation to get where we want
        std::string graphname = name.substr(name.find("h1D_")+4, name.size());
        graphname = "grCritical_"+graphname;
        gCrit68 = (TGraph*)f1->Get((graphname+"_conf68").c_str());
        gCrit90 = (TGraph*)f1->Get((graphname+"_conf90").c_str());
        gCrit955 = (TGraph*)f1->Get((graphname+"_conf955").c_str());
        gCrit997 = (TGraph*)f1->Get((graphname+"_conf997").c_str());

        gCrit68->SetLineColor(kBlack);
        gCrit90->SetLineColor(kGray+1);
        gCrit955->SetLineColor(kBlack);
        gCrit955->SetLineStyle(kDashed);
        if (gCrit997 != NULL) {
          gCrit997->SetLineColor(kGray+1);
          gCrit997->SetLineStyle(kDashed);
        }

        gCrit68->Draw("L*,same");
        gCrit90->Draw("L*,same");
        gCrit955->Draw("L*,same");
        // 99.7% critical level not calculate for th23
        if (gCrit997 != NULL) gCrit997->Draw("L*,same");

        legend2->SetNColumns(5);
        legend2->SetX1(0.15);
        legend2->AddEntry((TObject*)(NULL), "Feldman-Cousins:", "");
        legend2->AddEntry(gCrit68, "68%", "l");
        legend2->AddEntry(gCrit90, "90%", "l");
        legend2->AddEntry(gCrit955, "95.5%", "l");
        if (gCrit997 != NULL) {
          legend2->SetNColumns(6);
          legend2->AddEntry(gCrit997, "99.7%", "l");
        }
      }

      // theta 23 compute for large range, zoom in
      if (name.find("th23") != std::string::npos) {
        hPost->GetXaxis()->SetRangeUser(0.4, 0.65);
        hPost->GetYaxis()->SetRangeUser(0, 30);
      }
      else if (name.find("dm2") != std::string::npos) {
        hPost->GetYaxis()->SetRangeUser(0, 35);
        if (name.find("_NH") != std::string::npos) hPost->GetXaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2})");
        if (name.find("_IH") != std::string::npos)  hPost->GetXaxis()->SetTitle("|#Deltam^{2}_{31}| (eV^{2})");
      }
      else if (name.find("th13_woRC") != std::string::npos) {
        hPost->GetXaxis()->SetRangeUser(10E-3, 45E-3);
        hPost->GetYaxis()->SetRangeUser(0, 40);
      }
      else if (name.find("th13_woRC") != std::string::npos) {
        hPost->GetXaxis()->SetRangeUser(15E-3, 26E-3);
        hPost->GetYaxis()->SetRangeUser(0, 40);
      }


      legend2->Draw("same");

      canv->Print(output);

      delete legend2;
    } else {
      continue;
    }
  }
  canv->Print(output+"]");
}
