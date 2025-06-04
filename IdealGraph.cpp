#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TColor.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"

std::string NumErrScien(double x, std::string name_par)
{
    double mult = 1;
    std::string um_str; // unità di misura
    switch (name_par[0])
    {
    case 'R':
        um_str = "#Omega"; // ohm
        break;
    case 'L':
        mult = 1000;
        um_str = "mH"; // henry
        break;
    case 'C':
        mult = 1'000'000;
        um_str = "#muF"; // farad
        break;
    default:
        um_str = "V";
        break;
    }

    return std::to_string(mult * x) + um_str;
}

Double_t ampl_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};

    Double_t Rw{par[0]};
    Double_t L{par[1]};
    Double_t Vs{par[2]};
    return (Vs * Rw) / sqrt((Rw) * (Rw) + (w * L) * (w * L));
}

Double_t ampl_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};

    Double_t Rt{par[0]};
    Double_t C1C2{par[1]};
    Double_t Vs{par[2]};
    return (Vs * Rt) / sqrt((Rt) * (Rt) + (1 / (w * C1C2)) * (1 / (w * C1C2)));
}

Double_t phase_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t Rw{par[0]};
    Double_t L{par[1]};
    return -(std::atan(w * L / (Rw)));
}

Double_t phase_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};

    Double_t Rt{par[0]};
    Double_t C1C2{par[1]};
    return -(-std::atan(1 / (w * C1C2 * (Rt))));
}

void IdealGraph()
{
    Double_t const Rw{220};
    Double_t const L{0.047};
    Double_t const V0{2.5};
    Double_t const Rt{220};
    Double_t const C1C2{1.45E-06};

    Double_t Max = 1200;
    TF1 *ampl_func_w{new TF1{"ampl_func_w", ampl_woofer, 0., Max, 3}};
    TF1 *ampl_func_t{new TF1{"ampl_func_t", ampl_tweeter, 0., Max, 3}};
    TF1 *phase_func_w{new TF1{"phase_func_w", phase_woofer, 0., Max, 2}};
    TF1 *phase_func_t{new TF1{"phase_func_t", phase_tweeter, 0., Max, 2}};

    ampl_func_w->SetParName(0, "Rw");
    ampl_func_w->SetParName(1, "L ");
    ampl_func_w->SetParName(2, "V0");

    ampl_func_t->SetParName(0, "Rt");
    ampl_func_t->SetParName(1, "C1C2");
    ampl_func_t->SetParName(2, "V0");

    phase_func_w->SetParName(0, "Rw");
    phase_func_w->SetParName(1, "L ");

    phase_func_t->SetParName(0, "Rt");
    phase_func_t->SetParName(1, "C1C2");

    ampl_func_w->SetParameters(Rw, L, V0);
    ampl_func_t->SetParameters(Rt, C1C2, V0);
    phase_func_w->SetParameters(Rw, L);
    phase_func_t->SetParameters(Rt, C1C2);

    ampl_func_w->SetLineColor(kRed + 2);
    ampl_func_t->SetLineColor(kBlue + 2);
    phase_func_w->SetLineColor(kRed + 2);
    phase_func_t->SetLineColor(kBlue + 2);

    TGraph *w_ampl_gr = new TGraph(100);
    TGraph *t_ampl_gr = new TGraph(100);
    TGraph *w_phase_gr = new TGraph(100);
    TGraph *t_phase_gr = new TGraph(100);
    for (int i = 0; i != 100; ++i)
    {
        w_ampl_gr->SetPoint(i, i * (Max) / 100., ampl_func_w->Eval(i * (Max) / 100.));
        t_ampl_gr->SetPoint(i, i * (Max) / 100., ampl_func_t->Eval(i * (Max) / 100.));
        w_phase_gr->SetPoint(i, i * (Max) / 100., phase_func_w->Eval(i * (Max) / 100.));
        t_phase_gr->SetPoint(i, i * (Max) / 100., phase_func_t->Eval(i * (Max) / 100.));
    }

    w_ampl_gr->SetMarkerColor(kWhite);
    t_ampl_gr->SetMarkerColor(kWhite);
    t_phase_gr->SetMarkerColor(kWhite);
    w_phase_gr->SetMarkerColor(kWhite);
    TMultiGraph *ampl_multi{new TMultiGraph()};
    ampl_multi->Add(w_ampl_gr);
    ampl_multi->Add(t_ampl_gr);
    TMultiGraph *phase_multi{new TMultiGraph()};
    phase_multi->Add(t_phase_gr);
    phase_multi->Add(w_phase_gr);

    TCanvas *no_limits_no_cl_canvas{new TCanvas{"no_limits_no_cl_canvas", "Ideal Amplitude and Phase", 0, 0, 1300, 700}};
    no_limits_no_cl_canvas->Divide(2, 1);
    no_limits_no_cl_canvas->cd(1);

    // for legend
    gPad->SetGridy();
    gPad->SetTopMargin(0.2);

    ampl_multi->SetTitle("Ideal Amplitude - Frequency");
    ampl_multi->GetXaxis()->SetTitle("Frequency (Hz)");
    ampl_multi->GetYaxis()->SetTitle("Amplitude (V)");
    ampl_multi->Draw("AP");
    ampl_func_w->Draw("L, SAME");
    ampl_func_t->Draw("L, SAME");

    TLegend *no_limits_no_cl_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_no_cl_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_no_cl_legend_ampl->SetNColumns(3);
    // RIGA 1
    no_limits_no_cl_legend_ampl->AddEntry(ampl_func_w, "Amplitude V_Woofer", "l");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");

    // RIGA 2
    for (int i = 0; i != ampl_func_w->GetNpar(); ++i)
        no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w->GetParName(i)) + " = " + NumErrScien(ampl_func_w->GetParameter(i), ampl_func_w->GetParName(i))).c_str(), "");
    // RIGA 4
    no_limits_no_cl_legend_ampl->AddEntry(ampl_func_t, "Amplitude V_Tweeter", "l");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 5
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t->GetParName(i)) + " = " + NumErrScien(ampl_func_t->GetParameter(i), ampl_func_t->GetParName(i))).c_str(), "");

    no_limits_no_cl_legend_ampl->Draw();

    no_limits_no_cl_canvas->cd(2);

    gPad->SetGridy();
    gPad->SetTopMargin(0.2);

    phase_multi->SetTitle("Ideal Amplitude - Frequency");
    phase_multi->GetXaxis()->SetTitle("Frequency (Hz)");
    phase_multi->GetYaxis()->SetTitle("Phase shift (rad)");
    TF1 *phase_func_0 = new TF1("aa", "0", 0., 1200);
    phase_func_0->SetLineColor(kBlack);

    phase_multi->Draw("AP");
    phase_func_w->Draw("L, SAME");
    phase_func_t->Draw("L, SAME");
    phase_func_0->Draw("L, SAME");

    TLegend *no_limits_no_cl_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_no_cl_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_no_cl_legend_phase->SetNColumns(3);
    // RIGA 1
    no_limits_no_cl_legend_phase->AddEntry(phase_func_w, "Phase V_Woofer", "l");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");

    // RIGA 2
    for (int i = 0; i != phase_func_w->GetNpar(); ++i)
    {
        no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w->GetParName(i)) + " = " + NumErrScien(phase_func_w->GetParameter(i), phase_func_w->GetParName(i))).c_str(), "");
    }
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");

    // RIGA 4
    no_limits_no_cl_legend_phase->AddEntry(phase_func_t, "Phase V_Tweeter", "l");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 5
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t->GetParName(i)) + " = " + NumErrScien(phase_func_t->GetParameter(i), phase_func_t->GetParName(i))).c_str(), "");

    no_limits_no_cl_legend_phase->Draw();
}