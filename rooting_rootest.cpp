#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TGrid.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TVirtualFFT.h"
#include "TH1D.h"
#include "TFile.h"
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <string>

// returns the voltage read from settings.txt
int GetVoltage()
{
    std::ifstream settings_file{"./input_data/settings.txt"};
    if (!settings_file.is_open())
        std::cout << "Failed to open settings.txt" << std::endl;

    std::string input;
    do
    {
        settings_file >> input;
    } while (input != "vpp:" && !settings_file.eof());

    settings_file >> input;
    std::cout << "Vpp: " << atoi(input.c_str()) << std::endl;
    settings_file.close();
    return atoi(input.c_str());
}
void SamuStyle()
{
    gStyle->SetCanvasPreferGL();
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1);
}

Int_t N_BLOCKS{461};
std::string output_data_block_filename{"./data_root/output_data_block_"};
std::string output_param_filename{"./data_root/output_param_"};

int const N_POINTS_ERR{1000};
int const N_POINTS_ERR_TOT{25000};
int const N_POINTS_ERR_SUB{1000};
int const N_SUB_ERR{N_POINTS_ERR_TOT / N_POINTS_ERR_SUB};
int const N_BLOCKS_ERR{251};
// Double_t const R = 220;
// Double_t const R_L = 120;
// Double_t const R_Elvis = 50;
// Double_t const L = 47E-3;
// Double_t const C = 1.43E-6;

// Double_t const F_iniziale = 47E-3;
// Double_t const step_freq = 47E-3;

// need this as a global variable
TGraphErrors *ampl_graph[3];
TGraphErrors *phase_graph[3];

// par[0] -> Rw
// par[1] -> Rl
// par[2] -> L
Double_t ampl_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // clamp limits
    if (f[0] < ampl_graph[0]->GetPointX(0))
        Vs = ampl_graph[0]->GetPointY(0);
    else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
        Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    return (Vs * Rw) / sqrt((Rl + Rw) * (Rl + Rw) + (w * L) * (w * L));
}

Double_t ampl_woofer_2(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // clamp limits
    if (f[0] < ampl_graph[0]->GetPointX(0))
        Vs = ampl_graph[0]->GetPointY(0);
    else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
        Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Cl{par[3]};
    Double_t Num{(1 - w * w * L * Cl) * (1 - w * w * L * Cl) + (Rl * w * Cl) * (Rl * w * Cl)};
    Double_t Den_A{Rw + Rl - 2 * Rw * w * w * L * Cl + Rw * (w * w * L * Cl) * (w * w * L * Cl) + Rw * (Rl * w * Cl) * (Rl * w * Cl)};
    Double_t Den_B{w * L - w * w * w * L * L * Cl - Rl * Rl * w * Cl};
    return Vs * Rw * Num / sqrt(Den_A * Den_A + Den_B * Den_B);
}

Double_t ampl_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // clamp limits
    if (f[0] < ampl_graph[0]->GetPointX(0))
        Vs = ampl_graph[0]->GetPointY(0);
    else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
        Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);
    Double_t Rt{par[0]};
    Double_t Rl1Rl2{par[1]};
    Double_t C1C2{par[2]};
    return (Vs * Rt) / sqrt((Rt + Rl1Rl2) * (Rt + Rl1Rl2) + (1 / (w * C1C2)) * (1 / (w * C1C2)));
}

Double_t phase_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // clamp limits
    if (f[0] < phase_graph[0]->GetPointX(0))
        phase_s = phase_graph[0]->GetPointY(0);
    else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
        phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    return (+phase_s + std::atan(w * L / (Rl + Rw)));
}

Double_t phase_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // clamp limits
    if (f[0] < phase_graph[0]->GetPointX(0))
        phase_s = phase_graph[0]->GetPointY(0);
    else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
        phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rt{par[0]};
    Double_t Rl1Rl2{par[1]};
    Double_t C1C2{par[2]};
    return (+phase_s - std::atan(1 / (w * C1C2 * (Rt + Rl1Rl2))));
}

Double_t ClampAngle(Double_t angle)
{
    while (angle < -TMath::Pi())
        angle += 2 * TMath::Pi();

    while (angle > TMath::Pi())
        angle -= 2 * TMath::Pi();

    return angle;
}

void DrawBlock(int n_block)
{
    if (n_block < 0 || n_block >= N_BLOCKS)
    {
        return;
    }

    TGraphErrors *V_arr[3]{
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

    TF1 *func_arr[3]{
        new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

    std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
    double frequency;
    file_params_in >> frequency;
    double pulsation = TMath::TwoPi() * frequency;
    // set parameters
    for (int i = 0; i != 3; ++i)
    {
        double amplitude;
        double phase;
        file_params_in >> amplitude;
        file_params_in >> phase;

        phase = phase * TMath::Pi() / 180.; // to rad
        phase -= TMath::Pi() / 2.;          // to account for shifting in converter

        func_arr[i]->SetParameters(amplitude, pulsation, phase);
    }
    file_params_in.close();

    V_arr[0]->SetMarkerColor(kBlack);
    V_arr[1]->SetMarkerColor(kRed);
    V_arr[2]->SetMarkerColor(kBlue);

    V_arr[0]->SetLineColor(kBlack);
    V_arr[1]->SetLineColor(kRed);
    V_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetLineColor(kBlack);
    func_arr[1]->SetLineColor(kRed);
    func_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetNpx(10000);
    func_arr[1]->SetNpx(10000);
    func_arr[2]->SetNpx(10000);

    TMultiGraph *V_multi = new TMultiGraph();
    // V_multi->Add(V_arr[0]);
    V_multi->Add(V_arr[1]);
    V_multi->Add(V_arr[2]);

    TCanvas *test_canva = new TCanvas("test_canva", std::to_string(n_block).c_str(), 0, 0, 800, 600);
    V_multi->Draw("AP");
    // func_arr[0]->Draw("SAME");
    func_arr[1]->Draw("SAME");
    func_arr[2]->Draw("SAME");
}

void DrawSubBlock(int n_block, int n_sub_block)
{
    if (n_block < 0 || n_block >= N_BLOCKS_ERR)
    {
        return;
    }

    TGraphErrors *V_arr[3]{
        new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
        new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
        new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

    TF1 *func_arr[3]{
        new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

    std::ifstream file_params_in{"./mega_sweep_output/output_param_" + std::to_string(n_block) + ".txt"};
    double frequency;
    file_params_in >> frequency;
    double pulsation = TMath::TwoPi() * frequency;
    // set parameters
    for (int i = 0; i != 3; ++i)
    {
        double amplitude;
        double phase;
        file_params_in >> amplitude;
        file_params_in >> phase;

        phase = phase * TMath::Pi() / 180.; // to rad
        phase -= TMath::Pi() / 2.;          // to account for shifting in converter

        func_arr[i]->SetParameters(amplitude, pulsation, phase);
    }
    file_params_in.close();

    V_arr[0]->SetMarkerColor(kBlack);
    V_arr[1]->SetMarkerColor(kRed);
    V_arr[2]->SetMarkerColor(kBlue);

    V_arr[0]->SetLineColor(kBlack);
    V_arr[1]->SetLineColor(kRed);
    V_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetLineColor(kBlack);
    func_arr[1]->SetLineColor(kRed);
    func_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetNpx(100000);
    func_arr[1]->SetNpx(100000);
    func_arr[2]->SetNpx(100000);

    TMultiGraph *V_multi = new TMultiGraph();
    V_multi->Add(V_arr[0]);
    V_multi->Add(V_arr[1]);
    V_multi->Add(V_arr[2]);

    TCanvas *test_canva = new TCanvas("test_canva", std::to_string(n_block).c_str(), 0, 0, 800, 600);
    V_multi->Draw("AP");
    func_arr[0]->Draw("SAME");
    func_arr[1]->Draw("SAME");
    func_arr[2]->Draw("SAME");
}

void FitBlock(int n_block)
{
    if (n_block < 0 || n_block >= N_BLOCKS)
    {
        return;
    }

    TGraphErrors *V_arr[3]{
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

    TF1 *func_arr[3]{
        new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
        new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

    std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
    double frequency;
    file_params_in >> frequency;
    double pulsation = TMath::TwoPi() * frequency;
    // set parameters
    for (int i = 0; i != 3; ++i)
    {
        double amplitude;
        double phase;
        file_params_in >> amplitude;
        file_params_in >> phase;

        phase = phase * TMath::Pi() / 180.; // to rad
        phase -= TMath::Pi() / 2.;          // to account for shifting in converter

        func_arr[i]->SetParameters(amplitude, pulsation, phase);
    }
    file_params_in.close();

    func_arr[0]->SetNpx(10000);
    func_arr[1]->SetNpx(10000);
    func_arr[2]->SetNpx(10000);

    func_arr[0]->SetNumberFitPoints(10000);
    func_arr[1]->SetNumberFitPoints(10000);
    func_arr[2]->SetNumberFitPoints(10000);

    V_arr[0]->SetMarkerColor(kBlack);
    V_arr[1]->SetMarkerColor(kRed);
    V_arr[2]->SetMarkerColor(kBlue);

    V_arr[0]->SetLineColor(kBlack);
    V_arr[1]->SetLineColor(kRed);
    V_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetLineColor(kBlack);
    func_arr[1]->SetLineColor(TColor::GetColor(150, 0, 0));
    func_arr[2]->SetLineColor(kBlue);

    std::cout << "\n\nFitting of Block " << n_block << std::endl;
    for (int i = 0; i != 3; ++i)
    {
        std::cout << "V[" << i << "]:\n";
        V_arr[i]->Fit(func_arr[i]);
        std::cout << "Chi ridotto:" << func_arr[i]->GetChisquare() / func_arr[i]->GetNDF() << std::endl;
    }

    TMultiGraph *V_multi = new TMultiGraph();
    V_multi->Add(V_arr[0]);
    V_multi->Add(V_arr[1]);
    V_multi->Add(V_arr[2]);

    TCanvas *test_canva = new TCanvas("test_canva", std::to_string(n_block).c_str(), 0, 0, 800, 600);
    V_multi->Draw("APE");
}

void FitBlock(int n_block, int only_Vi)
{
    if (n_block < 0 || n_block >= N_BLOCKS || only_Vi < 0 || only_Vi > 2)
    {
        return;
    }

    TGraphErrors *V_arr[3]{
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

    // TF1 *func_arr[3]{
    //     new TF1{"V_s_fit", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5])"},
    //     new TF1{"V_w_fit", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5])"},
    //     new TF1{"V_t_fit", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5])"}};

    TF1 *func_arr[3]{
        new TF1{"V_s_fit", "[0]*cos([1]*x - [2]) + [3]"},
        new TF1{"V_w_fit", "[0]*cos([1]*x - [2]) + [3]"},
        new TF1{"V_t_fit", "[0]*cos([1]*x - [2]) + [3]"}};

    std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
    double frequency;
    file_params_in >> frequency;
    double pulsation = TMath::TwoPi() * frequency;
    // set parameters
    for (int i = 0; i != 3; ++i)
    {
        double amplitude;
        double phase;
        file_params_in >> amplitude;
        file_params_in >> phase;
        phase = phase * TMath::Pi() / 180.; // to rad
        phase -= TMath::Pi() / 2.;          // to account for shifting in converter

        func_arr[i]->SetParName(0, "amplitude");
        func_arr[i]->SetParName(1, "pulsation");
        func_arr[i]->SetParName(2, "phase");

        // amplitude = 1.2735;
        // func_arr[i]->SetParameters(amplitude, pulsation, phase, amplitude / 292, pulsation * 3, phase);
        func_arr[i]->SetParameter(0, amplitude);
        func_arr[i]->SetParameter(1, pulsation);
        func_arr[i]->SetParameter(2, phase);

        // just invented--------------------------------------------------------
        func_arr[i]->SetParName(3, "background");
        double background_shift = -0.0117464;
        func_arr[i]->SetParameter(3, background_shift);
        // just invented--------------------------------------------------------

        // func_arr[i]->SetParLimits(0, func_arr[i]->GetParameter(0) - func_arr[i]->GetParameter(0) / 10, func_arr[i]->GetParameter(0) + func_arr[i]->GetParameter(0) / 10);
        // func_arr[i]->SetParLimits(3, func_arr[i]->GetParameter(3) - func_arr[i]->GetParameter(3) / 2, func_arr[i]->GetParameter(3) + func_arr[i]->GetParameter(3) / 2);
        // func_arr[i]->FixParameter(1, pulsation);
        // func_arr[i]->FixParameter(4, pulsation*2.);
    }
    file_params_in.close();

    func_arr[only_Vi]->SetNpx(10000);
    func_arr[only_Vi]->SetNumberFitPoints(10000);

    V_arr[0]->SetMarkerColor(kBlack);
    V_arr[1]->SetMarkerColor(kRed);
    V_arr[2]->SetMarkerColor(kBlue);

    V_arr[0]->SetLineColor(kBlack);
    V_arr[1]->SetLineColor(kRed);
    V_arr[2]->SetLineColor(kBlue);

    func_arr[0]->SetLineColor(TColor::GetColor(40, 40, 40));
    func_arr[1]->SetLineColor(TColor::GetColor(150, 0, 0));
    func_arr[2]->SetLineColor(TColor::GetColor(0, 0, 150));
    std::cout << "\n\nFitting of Block " << n_block << "V[<<only_Vi" << "]" << std::endl;

    //-----------------------------------------------------------------------------
    //
    TGraphErrors *cleaned_graph{new TGraphErrors(*V_arr[only_Vi])};
    //
    //-----------------------------------------------------------------------------

    // V_arr[only_Vi]->Fit(func_arr[only_Vi], "E, M");
    V_arr[only_Vi]->Fit(func_arr[only_Vi]);
    std::cout << "Chi ridotto:" << func_arr[only_Vi]->GetChisquare() / func_arr[only_Vi]->GetNDF() << std::endl;

    // try removing the cos to see what's left behind
    for (int n_point = 0; n_point != V_arr[only_Vi]->GetN(); ++n_point)
    {
        Double_t x{V_arr[only_Vi]->GetPointX(n_point)};
        Double_t y{V_arr[only_Vi]->GetPointY(n_point) - func_arr[only_Vi]->Eval(x)};
        cleaned_graph->SetPoint(n_point, x, y);
    }

    std::cout << "RMS Cleaned Graph (to be added to error vpp): " << TMath::RMS(cleaned_graph->GetN(), cleaned_graph->GetY()) << std::endl;

    // TF1 *cleaned_func = new TF1("cleaned_func", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5]) + [6]*cos([7]*x - [8]) + [9]*cos([10]*x - [11]) + [12]*cos([13]*x-[14])");
    TF1 *cleaned_func = new TF1("cleaned_func", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5])");
    Double_t base_ampl{func_arr[only_Vi]->GetParameter(0)};
    Double_t base_puls{func_arr[only_Vi]->GetParameter(1)};
    Double_t base_phase{func_arr[only_Vi]->GetParameter(2)};

    // from fft
    cleaned_func->SetParameter(0, base_ampl * 27. / 62'000.);
    cleaned_func->SetParameter(1, base_puls * 22.5 / 11.5);
    cleaned_func->SetParameter(2, base_phase);
    cleaned_func->SetParameter(3, base_ampl * 16. / 62'000);
    cleaned_func->SetParameter(4, base_puls * 33.5 / 11.5);
    cleaned_func->SetParameter(5, base_phase);
    // cleaned_func->SetParameter(6, base_ampl * 9. / 62'000);
    // cleaned_func->SetParameter(7, base_puls * 77.5 / 11.5);
    // cleaned_func->SetParameter(8, base_phase);
    // cleaned_func->SetParameter(9, base_ampl * 10.5 / 62'000);
    // cleaned_func->SetParameter(10, base_puls * 121.5 / 11.5);
    // cleaned_func->SetParameter(11, base_phase);
    // cleaned_func->SetParameter(12, base_ampl * 18. / 62'000);
    // cleaned_func->SetParameter(13, base_puls * 4900. / 11.5);
    // cleaned_func->SetParameter(14, base_phase);

    cleaned_func->SetNpx(10000);
    cleaned_func->SetNumberFitPoints(10000);
    cleaned_func->SetLineColor(kGreen);
    cleaned_graph->Fit(cleaned_func, "M");
    std::cout << "Chi ridotto:" << cleaned_func->GetChisquare() / cleaned_func->GetNDF() << std::endl;

    cleaned_graph->SetLineColor(kTeal);
    cleaned_graph->SetMarkerColor(kTeal);

    // V_arr[only_Vi]->GetXaxis()->SetRangeUser(0.03206, 0.03215);
    // V_arr[only_Vi]->GetYaxis()->SetRangeUser(1.26, 1.28);
    // V_arr[only_Vi]->GetYaxis()
    TMultiGraph *multigraph = new TMultiGraph();
    multigraph->Add(V_arr[only_Vi]);
    multigraph->Add(cleaned_graph);

    TCanvas *test_canva = new TCanvas("test_canva", (std::to_string(n_block) + ", [" + std::to_string(only_Vi) + "]").c_str(), 0, 0, 800, 600);
    test_canva->Divide(2, 1);
    test_canva->cd(1);
    test_canva->SetGridy(1);
    multigraph->Draw("APE");

    // fft-------------------------------------------------------------------------
    TH1D *histo_cleaned = new TH1D("histo", "aa", cleaned_graph->GetN(), cleaned_graph->GetPointX(0), cleaned_graph->GetPointX(cleaned_graph->GetN() - 1));
    for (int i = 1; i <= cleaned_graph->GetN(); ++i)
    {
        histo_cleaned->SetBinContent(i, cleaned_graph->GetPointY(i - 1));
    }

    TH1 *fft_result_histo{nullptr};
    TVirtualFFT::SetTransform(nullptr);
    fft_result_histo = histo_cleaned->FFT(fft_result_histo, "MAG, EX");

    test_canva->cd(2);
    fft_result_histo->Draw();
}

void PhaseShiftError(int n_blocks_input)
{
    N_BLOCKS = n_blocks_input;

    Double_t freq_arr[N_BLOCKS];
    Double_t freq_err_arr[N_BLOCKS];

    Double_t phase_arr[3][N_BLOCKS];
    Double_t phase_err_arr[3][N_BLOCKS];

    for (size_t n_block = 0; n_block != N_BLOCKS; ++n_block)
    {
        TGraphErrors *V_arr[3]{
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

        TF1 *func_arr[3]{
            new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
            new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
            new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

        std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
        double frequency;
        file_params_in >> frequency;

        // set parameters
        for (int i = 0; i != 3; ++i)
        {
            double amplitude;
            double phase;
            file_params_in >> amplitude;
            file_params_in >> phase;
            phase = phase * TMath::Pi() / 180.; // to rad

            phase -= TMath::Pi() / 2.; // to account for shifting in converter

            double pulsation = 2 * TMath::Pi() * frequency;

            func_arr[i]->SetParameter(0, amplitude);
            func_arr[i]->SetParameter(1, pulsation);
            // func_arr[i]->SetParameter(2, phase);

            func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
            func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
            // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);

            func_arr[i]->SetNumberFitPoints(10000);
            func_arr[i]->SetNpx(10000);

            if (((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0)) // i.e.: error
            {
                std::cout << "Invalid Fit! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";

                // only V_s freq
                if (i == 0)
                {
                    freq_arr[n_block] = pulsation / (2. * TMath::Pi());
                    freq_err_arr[n_block] = (pulsation / 10.) / (2. * TMath::Pi());
                }
                phase_arr[i][n_block] = ClampAngle(phase);
                phase_err_arr[i][n_block] = ClampAngle(phase) / 10;
            }
            else
            { // only V_s freq

                // ci è stato riferito (Marco e Fonseca) che root quando calcola gli errori sui parametri
                // restituisce l'autovalore della matrice di covarianza che, in teoria, va moltiplicato
                // per la sqrt(chi quadro ridotto). Root non lo fa perché assume sia circa 1. Per correttezza
                // lo moltiplichiamo noi.
                Double_t root_chi{std::sqrt(func_arr[i]->GetChisquare() / func_arr[i]->GetNDF())};
                if (i == 0)
                {
                    freq_arr[n_block] = func_arr[i]->GetParameter(1) / (2. * TMath::Pi());
                    freq_err_arr[n_block] = (func_arr[i]->GetParError(1) * root_chi) / (2. * TMath::Pi());
                }
                phase_arr[i][n_block] = ClampAngle(func_arr[i]->GetParameter(2));
                phase_err_arr[i][n_block] = func_arr[i]->GetParError(2) * root_chi;
            }
            std::cout << "[" << i << "]" << "[" << n_block << "]" << func_arr[i]->GetChisquare() / func_arr[i]->GetNDF() << '\n';
        }
        file_params_in.close();
        for (int i = 0; i != 3; ++i)
        {
            delete V_arr[i];
            delete func_arr[i];
        }
    }

    TGraphErrors *phase_shift_graph[3];
    TF1 *func_phase_shift[3];
    for (int i = 0; i != 3; ++i)
    {
        phase_shift_graph[i] = new TGraphErrors(N_BLOCKS, freq_arr, phase_arr[i], freq_err_arr, phase_err_arr[i]);
        func_phase_shift[i] = new TF1(("func_phase_shift_" + std::to_string(i)).c_str(), "[0]+[1]*x", 0., 25000.);
        func_phase_shift[i]->SetParName(0, "intercept");
        func_phase_shift[i]->SetParName(1, "slope");
        // func_phase_shift[i]->SetParameter(1, -TMath::PiOver2());
    }
    // func_phase_shift[0]->SetParameter(0, 0);
    // func_phase_shift[1]->SetParameter(0, -9.6E-6);
    // func_phase_shift[2]->SetParameter(0, -2.64E-5);

    // REGRESSIONE LINEARE y = A + Bx
    for (int i = 0; i != 3; ++i)
    {
        Double_t sum_wx2 = 0;
        Double_t sum_wx = 0;
        Double_t sum_w = 0;
        for (int j = 0; j != N_BLOCKS; ++j)
        {
            Double_t w = 1 / (phase_err_arr[i][j] * phase_err_arr[i][j]);
            Double_t x = freq_arr[j];

            sum_wx2 += w * x * x;
            sum_wx += w * x;
            sum_w += w;
        }

        Double_t delta_w = sum_w * sum_wx2 - sum_wx * sum_wx;
        Double_t A{0.};
        Double_t B{0.};
        Double_t A_err{std::sqrt(sum_wx2 / delta_w)};
        Double_t B_err{std::sqrt(sum_w / delta_w)};

        for (int j = 0; j != N_BLOCKS; ++j)
        {
            Double_t w = 1 / (phase_err_arr[i][j] * phase_err_arr[i][j]);
            Double_t x = freq_arr[j];
            Double_t y = phase_arr[i][j];

            Double_t alpha = (w * sum_wx2 - w * x * sum_wx) / delta_w;
            Double_t beta = (w * x * sum_w - w * sum_wx) / delta_w;

            A += alpha * y;
            B += beta * y;
        }

        func_phase_shift[i]->SetParameter(0, A);
        func_phase_shift[i]->SetParameter(1, B);
        func_phase_shift[i]->SetParError(0, A_err);
        func_phase_shift[i]->SetParError(1, B_err);
    }

    func_phase_shift[0]->SetLineColor(kBlack);
    func_phase_shift[1]->SetLineColor(kRed);
    func_phase_shift[2]->SetLineColor(kBlue);

    phase_shift_graph[0]->SetLineColor(kBlack);
    phase_shift_graph[1]->SetLineColor(kRed);
    phase_shift_graph[2]->SetLineColor(kBlue);

    phase_shift_graph[0]->SetMarkerColor(kBlack);
    phase_shift_graph[1]->SetMarkerColor(kRed);
    phase_shift_graph[2]->SetMarkerColor(kBlue);

    // for (int i = 0; i != 3; ++i)
    // {
    //     phase_shift_graph[i]->Fit(func_phase_shift[i]);
    // }

    TMultiGraph *phase_shift_multi = new TMultiGraph();

    for (int i = 0; i != 3; ++i)
    {
        phase_shift_multi->Add(phase_shift_graph[i]);
    }

    TCanvas *phase_shift_canvas = new TCanvas("phase_shift_canvas", "Phase Shift", 0, 0, 800, 600);
    phase_shift_multi->Draw("APE");
    for (int i = 0; i != 3; ++i)
    {
        func_phase_shift[i]->Draw("SAME");
        std::cout << "FUNC_ARRAY[" << i << "]\n";
        std::cout << func_phase_shift[i]->GetParName(0) << ": " << func_phase_shift[i]->GetParameter(0) << " +- " << func_phase_shift[i]->GetParError(0) << '\n';
        std::cout << func_phase_shift[i]->GetParName(1) << ": " << func_phase_shift[i]->GetParameter(1) << " +- " << func_phase_shift[i]->GetParError(1) << '\n';
    }

    // for (int i = 0; i != 3; ++i)
    // {
    //     // Double_t root_chi = std::sqrt(func_phase_shift[i]->GetChisquare() / func_phase_shift[i]->GetNDF());
    //     func_phase_shift[i]->SetParError(0, func_phase_shift[i]->GetParError(0)); //* root_chi);
    //     func_phase_shift[i]->SetParError(1, func_phase_shift[i]->GetParError(1)); //* root_chi);
    // }

    TFile *file_out = new TFile(("./risultati/Shift_error/V_" + std::to_string(GetVoltage()) + ".root").c_str(), "RECREATE");
    for (int i = 0; i != 3; ++i)
    {
        func_phase_shift[i]->Write();
    }
    file_out->Close();
}

void FitError(int N_BLOCKS)
{

    // mean of the 3 RMS in each block (RMS of the residual after the fit)
    TGraph *RMS{new TGraph(N_BLOCKS)};

    for (int n_block = 0; n_block != N_BLOCKS; ++n_block)
    {
        TGraphErrors *V_arr[3]{
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

        TF1 *func_arr[3]{
            new TF1{"V_s_fit", "[0]*cos([1]*x - [2]) + [3]"},
            new TF1{"V_w_fit", "[0]*cos([1]*x - [2]) + [3]"},
            new TF1{"V_t_fit", "[0]*cos([1]*x - [2]) + [3]"}};

        std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
        double frequency;
        file_params_in >> frequency;
        double pulsation = TMath::TwoPi() * frequency;
        // set parameters
        Double_t rms = 0; // we'll sum the 3 rms in this variable then divide by 3 and add it to the RMS Tgraph
        for (int i = 0; i != 3; ++i)
        {
            double amplitude;
            double phase;
            file_params_in >> amplitude;
            file_params_in >> phase;
            phase = phase * TMath::Pi() / 180.; // to rad
            phase -= TMath::Pi() / 2.;          // to account for shifting in converter

            func_arr[i]->SetParName(0, "amplitude");
            func_arr[i]->SetParName(1, "pulsation");
            func_arr[i]->SetParName(2, "phase");

            // amplitude = 1.2735;
            // func_arr[i]->SetParameters(amplitude, pulsation, phase, amplitude / 292, pulsation * 3, phase);
            func_arr[i]->SetParameter(0, amplitude);
            func_arr[i]->SetParameter(1, pulsation);
            func_arr[i]->SetParameter(2, phase);

            // just invented--------------------------------------------------------
            func_arr[i]->SetParName(3, "background");
            double background_shift = -0.0117464;
            func_arr[i]->SetParameter(3, background_shift);
            // just invented-------------------------------------------------------

            func_arr[i]->SetNpx(10000);
            func_arr[i]->SetNumberFitPoints(10000);

            if (((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0)) // i.e.: error
            {
                std::cout << "Invalid Fit! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";
                func_arr[i]->SetParameter(0, amplitude);
                func_arr[i]->SetParameter(1, pulsation);
                func_arr[i]->SetParameter(2, phase);
                func_arr[i]->SetParameter(3, background_shift);
            }

            std::cout << "Chi ridotto: [" << n_block << "][" << i << "]:" << func_arr[i]->GetChisquare() / func_arr[i]->GetNDF() << '\n';

            // try removing the cos to see what's left behind
            for (int n_point = 0; n_point != V_arr[i]->GetN(); ++n_point)
            {
                Double_t x{V_arr[i]->GetPointX(n_point)};
                Double_t y{V_arr[i]->GetPointY(n_point) - func_arr[i]->Eval(x)};
                V_arr[i]->SetPoint(n_point, x, y);
            }

            rms += TMath::RMS(V_arr[i]->GetN(), V_arr[i]->GetY());
        }
        file_params_in.close();

        rms = rms / 3.;
        // mean and pulsation->freq
        Double_t freq = (func_arr[0]->GetParameter(1) + func_arr[1]->GetParameter(1) + func_arr[2]->GetParameter(1)) / (TMath::TwoPi() * 3);
        std::cout << "Block[" << n_block << "]:\n";
        std::cout << "\tFreq:<<" << freq << '\n';
        std::cout << "\tRMS:<<" << rms << '\n';

        for (int i = 0; i != 3; ++i)
        {
            delete V_arr[i];
            delete func_arr[i];
        }

        RMS->SetPoint(n_block, freq, rms);
    }
    TCanvas *rms_freq = new TCanvas("rms_freq", "rms_freq", 0, 0, 800, 600);
    RMS->Draw("APE");

    TFile *file_out = new TFile(("./risultati/Background_error/V_" + std::to_string(GetVoltage()) + ".root").c_str(), "RECREATE");
    RMS->Write();
    file_out->Close();

    std::ofstream file_out_array{"./risultati/Background_error/V_" + std::to_string(GetVoltage()) + "_Freq_RMS.txt", std::ios::trunc};
    file_out_array << std::setprecision(10);
    for (int n_block = 0; n_block != N_BLOCKS; ++n_block)
    {
        if (n_block != N_BLOCKS - 1)
            file_out_array << RMS->GetPointX(n_block) << '\t' << RMS->GetPointY(n_block) << '\n';
        else
            file_out_array << RMS->GetPointX(n_block) << '\t' << RMS->GetPointY(n_block);
    }
    file_out_array.close();
}

// //[which voltage][block number][which subdivision]
// Double_t ampl_arr[3][N_BLOCKS_ERR][N_SUB_ERR];
// Double_t phase_arr[3][N_BLOCKS_ERR][N_SUB_ERR];
// Double_t freq_arr[N_BLOCKS_ERR];

// for (int n_block = 0; n_block != N_BLOCKS_ERR; ++n_block)
// {
//     for (int n_sub_block = 0; n_sub_block != N_SUB_ERR; ++n_sub_block)
//     {
//         TGraphErrors *V_arr[3]{
//             new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
//             new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
//             new TGraphErrors{("./mega_sweep_output/output_data_block_" + std::to_string(n_block) + "_" + std::to_string(n_sub_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

//         TF1 *func_arr[3]{
//             new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
//             new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
//             new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

//         std::ifstream file_params_in{"./mega_sweep_output/output_param_" + std::to_string(n_block) + ".txt"};
//         double frequency;
//         file_params_in >> frequency;

//         ampl_graph[0]->SetMarkerColor(kBlack);
//         ampl_graph[1]->SetMarkerColor(kRed);
//         ampl_graph[2]->SetMarkerColor(kBlue);
//         ampl_graph[0]->SetLineColor(kBlack);
//         ampl_graph[1]->SetLineColor(kRed);
//         ampl_graph[2]->SetLineColor(kBlue);

//         phase_graph[0]->SetMarkerColor(kBlack);
//         phase_graph[1]->SetMarkerColor(kRed);
//         phase_graph[2]->SetMarkerColor(kBlue);
//         phase_graph[0]->SetLineColor(kBlack);
//         phase_graph[1]->SetLineColor(kRed);
//         phase_graph[2]->SetLineColor(kBlue);

//         // set parameters
//         for (int i = 0; i != 3; ++i)
//         {
//             double amplitude;
//             double phase;
//             file_params_in >> amplitude;
//             file_params_in >> phase;
//             phase = phase * TMath::Pi() / 180.; // to rad

//             phase -= TMath::Pi() / 2.; // to account for shifting in converter

//             double pulsation = 2 * TMath::Pi() * frequency;

//             func_arr[i]->SetParameter(0, amplitude);
//             func_arr[i]->SetParameter(1, pulsation);
//             // func_arr[i]->SetParameter(2, phase);

//             func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
//             func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
//             // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);
//             func_arr[i]->SetParLimits(2, -TMath::Pi(), TMath::Pi());

//             // func_arr[i]->FixParameter(0, amplitude);
//             func_arr[i]->FixParameter(1, pulsation);

//             func_arr[i]->SetNumberFitPoints(10000);
//             func_arr[i]->SetNumberFitPoints(10000);

//             if (V_arr[i]->Fit(func_arr[i], "QUIET") != 01) // E: Better parameter error estimation
//             {
//                 std::cout << "Invalid Fit! Block n: " << n_block << "_" << n_sub_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";

//                 // only V_s freq
//                 if (i == 0)
//                 {
//                     // rewrites the same thing for each sub block
//                     freq_arr[n_block] = pulsation / (2. * TMath::Pi());
//                 }

//                 ampl_arr[i][n_block][n_sub_block] = amplitude;
//                 phase_arr[i][n_block][n_sub_block] = ClampAngle(phase);
//             }
//             else
//             { // only V_s freq
//                 if (i == 0)
//                 {
//                     // rewrites the same thing for each sub block
//                     freq_arr[n_block] = func_arr[i]->GetParameter(1) / (2. * TMath::Pi());
//                 }
//                 ampl_arr[i][n_block][n_sub_block] = func_arr[i]->GetParameter(0);
//                 phase_arr[i][n_block][n_sub_block] = ClampAngle(func_arr[i]->GetParameter(2));
//             }
//             std::cout << "[" << i << "]" << "[" << n_block << "]" << func_arr[i]->GetChisquare() / func_arr[i]->GetNDF() << '\n';
//         }

//         file_params_in.close();
//         for (int i = 0; i != 3; ++i)
//         {
//             delete V_arr[i];
//             delete func_arr[i];
//         }
//     }
// }

// TGraphErrors *ampl_graph[3];
// Double_t amplitudes_mean[3][N_BLOCKS_ERR];
// Double_t amplitudes_dev_std[3][N_BLOCKS_ERR];

// TGraphErrors *phase_graph[3];
// Double_t phase_mean[3][N_BLOCKS_ERR];
// Double_t phase_dev_std[3][N_BLOCKS_ERR];

// for (int i = 0; i != 3; ++i)
// {
//     for (int n_block = 0; n_block != N_BLOCKS_ERR; ++n_block)
//     {
//         amplitudes_mean[i][n_block] = TMath::Mean(N_SUB_ERR, ampl_arr[i][n_block]);
//         amplitudes_dev_std[i][n_block] = TMath::RMS(N_SUB_ERR, ampl_arr[i][n_block]); // è la dev standard con N-1

//         phase_mean[i][n_block] = TMath::Mean(N_SUB_ERR, phase_arr[i][n_block]);
//         phase_dev_std[i][n_block] = TMath::RMS(N_SUB_ERR, phase_arr[i][n_block]); // è la dev standard con N-1
//     }
//     ampl_graph[i] = new TGraphErrors(N_BLOCKS_ERR, freq_arr, amplitudes_mean[i], amplitudes_dev_std[i]);
//     phase_graph[i] = new TGraphErrors(N_BLOCKS_ERR, freq_arr, phase_mean[i], phase_dev_std[i]);
// }

// TCanvas *canvas_ampl_fit_err = new TCanvas("canvas_ampl_fit_err", "canvas_ampl_fit_err", 0, 0, 800, 601);
// TMultiGraph *multi_ampl_fit_err = new TMultiGraph();
// for (int i = 0; i != 3; ++i)
// {
//     multi_ampl_fit_err->Add(ampl_graph[i]);
// }
// multi_ampl_fit_err->Draw("APE");

// TCanvas *canvas_phase_fit_err = new TCanvas("canvas_phase_fit_err", "canvas_phase_fit_err", 0, 0, 800, 601);
// TMultiGraph *multi_phase_fit_err = new TMultiGraph();
// for (int i = 0; i != 3; ++i)
// {
//     multi_phase_fit_err->Add(phase_graph[i]);
// }
// multi_phase_fit_err->Draw("APE");
// }

void FastFourierTransform(int n_block, int only_Vi)
{
    TGraphErrors *V_arr[3]{
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
        new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

    TH1D *V_histo{new TH1D("V_histo", (std::to_string(n_block) + "[" + std::to_string(only_Vi) + "]").c_str(),
                           V_arr[only_Vi]->GetN(), V_arr[only_Vi]->GetPointX(0), V_arr[only_Vi]->GetPointX(V_arr[only_Vi]->GetN() - 1))};

    for (int i = 1; i != V_arr[only_Vi]->GetN() + 1; ++i)
    {
        V_histo->SetBinContent(i, V_arr[only_Vi]->GetPointY(i - 1));
    }

    TH1 *fft_result_histo{nullptr};
    TVirtualFFT::SetTransform(nullptr);
    fft_result_histo = V_histo->FFT(fft_result_histo, "MAG, EX");

    TCanvas *fourier_canvas = new TCanvas("fourier_canvas", (std::to_string(n_block) + "[" + std::to_string(only_Vi) + "]").c_str(), 0, 0, 800, 600);
    fft_result_histo->GetXaxis()->SetRangeUser(0, 500);
    fft_result_histo->Draw();
}
void rooting_rootest(Int_t input_n_blocks)
{
    N_BLOCKS = input_n_blocks;
    SamuStyle();
    std::ifstream file_count{"./input_data/ampl_data"};

    Double_t const Rw{2.2032000E+02};
    Double_t const Rt{2.2019000E+02};
    Double_t const Rl{1.2018250E+02};     // Actual L's resistance
    Double_t const Rl1Rl2{1.2183500E+02}; // equivalent L's resistance on tweeter
    Double_t const C1C2{1.4225000E-06};
    Double_t const L{4.7154000E-02};

    Double_t const Rw_err{1.4737707E-01};
    Double_t const Rt_err{1.0639549E-01};
    Double_t const Rl_err{4.8897001E-01};     // Actual L's resistance
    Double_t const Rl1Rl2_err{4.5793013E-01}; // equivalent L's resistance on tweeter
    Double_t const C1C2_err{1.4225000E-08};
    Double_t const L_err{4.7154000E-04};

    Double_t freq_arr[N_BLOCKS];
    Double_t freq_err_arr[N_BLOCKS];

    Double_t ampl_arr[3][N_BLOCKS];
    Double_t ampl_err_arr[3][N_BLOCKS];

    Double_t phase_arr[3][N_BLOCKS];
    Double_t phase_err_arr[3][N_BLOCKS];

    for (size_t n_block = 0; n_block != N_BLOCKS; ++n_block)
    {
        TGraphErrors *V_arr[3]{
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %lg %*lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %lg %*lg %lg %lg"},
            new TGraphErrors{(output_data_block_filename + std::to_string(n_block) + ".txt").c_str(), "%lg %*lg %*lg %lg %lg %lg"}};

        TF1 *func_arr[3]{
            new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
            new TF1{"V_w_fit", "[0]*cos([1]*x - [2])"},
            new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

        // TF1 *func_arr[3]{
        //     new TF1{"V_s_fit", "[0]*cos([1]*x - [2])"},
        //     new TF1{"V_w_fit", "[0]*cos([1]*x - [2]) + [3]*cos([4]*x - [5])"},
        //     new TF1{"V_t_fit", "[0]*cos([1]*x - [2])"}};

        std::ifstream file_params_in{output_param_filename + std::to_string(n_block) + ".txt"};
        double frequency;
        file_params_in >> frequency;

        // set parameters
        for (int i = 0; i != 3; ++i)
        {
            double amplitude;
            double phase;
            file_params_in >> amplitude;
            file_params_in >> phase;
            phase = phase * TMath::Pi() / 180.; // to rad

            phase -= TMath::Pi() / 2.; // to account for shifting in converter

            double pulsation = 2 * TMath::Pi() * frequency;

            // test--------------------------------------------------------------------
            // if (i == 1)
            //     func_arr[i]->SetParameters(amplitude, pulsation, phase, amplitude / 292, pulsation * 3, phase);
            // else
            // {

            func_arr[i]->SetParameter(0, amplitude);
            func_arr[i]->SetParameter(1, pulsation);
            // func_arr[i]->SetParameter(2, phase);

            func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
            func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
            // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);

            // func_arr[i]->FixParameter(0, amplitude);
            // func_arr[i]->FixParameter(1, pulsation);
            // }
            // test--------------------------------------------------------------------

            func_arr[i]->SetNumberFitPoints(10000);
            func_arr[i]->SetNpx(10000);

            if (((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0)) // i.e.: error
            {
                std::cout << "Invalid Fit! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";

                // only V_s freq
                if (i == 0)
                {
                    freq_arr[n_block] = pulsation / (2. * TMath::Pi());
                    freq_err_arr[n_block] = (pulsation / 10.) / (2. * TMath::Pi());
                }
                ampl_arr[i][n_block] = amplitude;
                ampl_err_arr[i][n_block] = amplitude / 10.;
                phase_arr[i][n_block] = ClampAngle(phase);
                phase_err_arr[i][n_block] = ClampAngle(phase) / 10;
            }
            else
            { // only V_s freq

                // ci è stato riferito (Marco e Fonseca) che root quando calcola gli errori sui parametri
                // restituisce l'autovalore della matrice di covarianza che, in teoria, va moltiplicato
                // per la sqrt(chi quadro ridotto). Root non lo fa perché assume sia circa 1. Per correttezza
                // lo moltiplichiamo noi.
                Double_t root_chi{std::sqrt(func_arr[i]->GetChisquare() / func_arr[i]->GetNDF())};
                if (i == 0)
                {
                    freq_arr[n_block] = func_arr[i]->GetParameter(1) / (2. * TMath::Pi());
                    freq_err_arr[n_block] = (func_arr[i]->GetParError(1) * root_chi) / (2. * TMath::Pi());
                }
                ampl_arr[i][n_block] = func_arr[i]->GetParameter(0);
                ampl_err_arr[i][n_block] = func_arr[i]->GetParError(0) * root_chi;
                phase_arr[i][n_block] = ClampAngle(func_arr[i]->GetParameter(2));
                phase_err_arr[i][n_block] = func_arr[i]->GetParError(2) * root_chi;
            }
            std::cout << "[" << i << "]" << "[" << n_block << "]" << func_arr[i]->GetChisquare() / func_arr[i]->GetNDF() << '\n';
        }
        file_params_in.close();
        for (int i = 0; i != 3; ++i)
        {
            delete V_arr[i];
            delete func_arr[i];
        }
    }

    // phase correction
    TFile *phase_correction_file{new TFile(("./risultati/Shift_error/V_" + std::to_string(GetVoltage()) + ".root").c_str(), "READ")};
    TF1 *func_phase_shift[3];
    for (int i = 0; i != 3; ++i)
    {
        func_phase_shift[i] = static_cast<TF1 *>(phase_correction_file->Get(("func_phase_shift_" + std::to_string(i)).c_str()));
    }

    for (int i = 1; i != 3; ++i)
    {
        Double_t A_0{func_phase_shift[0]->GetParameter(0)};
        Double_t B_0{func_phase_shift[0]->GetParameter(1)};
        Double_t A_0_err{func_phase_shift[0]->GetParError(0)};
        Double_t B_0_err{func_phase_shift[0]->GetParError(1)};

        Double_t A{func_phase_shift[i]->GetParameter(0)};
        Double_t B{func_phase_shift[i]->GetParameter(1)};
        Double_t A_err{func_phase_shift[i]->GetParError(0)};
        Double_t B_err{func_phase_shift[i]->GetParError(1)};

        for (int n_block = 0; n_block != N_BLOCKS; ++n_block)
        {
            Double_t f{freq_arr[n_block]};
            Double_t f_err{freq_err_arr[n_block]};
            Double_t phase_0 = A_0 + B_0 * f;
            Double_t phase_i = A + B * f;
            Double_t phase_correction = phase_i - phase_0;
            phase_arr[i][n_block] -= phase_correction;

            // NOTA: SBAGLIATO, MANCA LA COVARIANZA TRA A E B
            Double_t phase_correction_err = std::sqrt(A_0_err * A_0_err + A_err * A_err + (f_err * B_0) * (f_err * B_0) + (f_err * B) * (f_err * B) + ((B - B_0) * f_err) * ((B - B_0) * f_err));
            phase_err_arr[i][n_block] = std::sqrt(phase_err_arr[i][n_block] * phase_err_arr[i][n_block] + phase_correction_err * phase_correction_err);
        }
    }

    // ampl_graph declared as local variable
    ampl_graph[0] = new TGraphErrors{N_BLOCKS, freq_arr, ampl_arr[0], freq_err_arr, ampl_err_arr[0]};
    ampl_graph[1] = new TGraphErrors{N_BLOCKS, freq_arr, ampl_arr[1], freq_err_arr, ampl_err_arr[1]};
    ampl_graph[2] = new TGraphErrors{N_BLOCKS, freq_arr, ampl_arr[2], freq_err_arr, ampl_err_arr[2]};

    phase_graph[0] = new TGraphErrors{N_BLOCKS, freq_arr, phase_arr[0], freq_err_arr, phase_err_arr[0]};
    phase_graph[1] = new TGraphErrors{N_BLOCKS, freq_arr, phase_arr[1], freq_err_arr, phase_err_arr[1]};
    phase_graph[2] = new TGraphErrors{N_BLOCKS, freq_arr, phase_arr[2], freq_err_arr, phase_err_arr[2]};

    // fitting

    TF1 *ampl_func_w{new TF1{"ampl_func_w", ampl_woofer_2, 0., 1000., 4}};
    TF1 *ampl_func_t{new TF1{"ampl_func_t", ampl_tweeter, 0., 1000., 3}};
    TF1 *phase_func_w{new TF1{"phase_func_w", phase_woofer, 0., 1000., 3}};
    TF1 *phase_func_t{new TF1{"phase_func_t", phase_tweeter, 0., 1000., 3}};

    ampl_func_w->SetParName(0, "Rw");
    ampl_func_w->SetParName(1, "Rl");
    ampl_func_w->SetParName(2, "L");
    ampl_func_w->SetParName(3, "Cl");

    ampl_func_t->SetParName(0, "Rt");
    ampl_func_t->SetParName(1, "Rl1Rl2");
    ampl_func_t->SetParName(2, "C1C2");

    phase_func_w->SetParName(0, "Rw");
    phase_func_w->SetParName(1, "Rl");
    phase_func_w->SetParName(2, "L");

    phase_func_t->SetParName(0, "Rt");
    phase_func_t->SetParName(1, "Rl1Rl2");
    phase_func_t->SetParName(2, "C1C2");

    ampl_func_w->SetParameter(0, Rw);
    ampl_func_w->SetParameter(1, Rl);
    ampl_func_w->SetParameter(2, L);
    ampl_func_w->SetParameter(3, 1E-9);

    ampl_func_t->SetParameter(0, Rl);
    ampl_func_t->SetParameter(1, Rl1Rl2);
    ampl_func_t->SetParameter(2, C1C2);

    phase_func_w->SetParameter(0, Rw);
    phase_func_w->SetParameter(1, Rl);
    phase_func_w->SetParameter(2, L);

    phase_func_t->SetParameter(0, Rl);
    phase_func_t->SetParameter(1, Rl1Rl2);
    phase_func_t->SetParameter(2, C1C2);

    // PAR LIMITS-----------------------------------------------------------------
    ampl_func_w->SetParLimits(0, Rw - Rw_err, Rw + Rw_err);
    ampl_func_w->SetParLimits(1, Rl - Rl_err, Rl + Rl_err);
    ampl_func_w->SetParLimits(2, L - L_err, L + L_err);

    ampl_func_t->SetParLimits(0, Rt - Rt_err, Rt + Rt_err);
    ampl_func_t->SetParLimits(1, Rl1Rl2 - Rl1Rl2_err, Rl1Rl2 + Rl1Rl2_err);
    ampl_func_t->SetParLimits(2, C1C2 - C1C2_err, C1C2 + C1C2_err);

    phase_func_w->SetParLimits(0, Rw - Rw_err, Rw + Rw_err);
    phase_func_w->SetParLimits(1, Rl - Rl_err, Rl + Rl_err);
    phase_func_w->SetParLimits(2, L - L_err, L + L_err);

    phase_func_t->SetParLimits(0, Rt - Rt_err, Rt + Rt_err);
    phase_func_t->SetParLimits(1, Rl1Rl2 - Rl1Rl2_err, Rl1Rl2 + Rl1Rl2_err);
    phase_func_t->SetParLimits(2, C1C2 - C1C2_err, C1C2 + C1C2_err);

    // GRAPHICS-----------------------------------------------
    ampl_func_t->SetLineColor(kBlue);
    ampl_func_w->SetLineColor(kRed);

    phase_func_t->SetLineColor(kBlue);
    phase_func_w->SetLineColor(kRed);

    ampl_graph[1]->Fit(ampl_func_w, "M, E");
    std::cout << std::endl;
    ampl_graph[2]->Fit(ampl_func_t, "M, E");
    std::cout << std::endl;

    phase_graph[1]->Fit(phase_func_w, "M, E");
    std::cout << std::endl;
    phase_graph[2]->Fit(phase_func_t, "M, E");
    std::cout << std::endl;

    ampl_graph[0]->SetName("V Elvis");
    ampl_graph[1]->SetName("V Woofer");
    ampl_graph[2]->SetName("V Tweeter");

    ampl_graph[0]->SetTitle("V Elvis");
    ampl_graph[1]->SetTitle("V Woofer");
    ampl_graph[2]->SetTitle("V Tweeter");

    ampl_graph[0]->SetMarkerColor(kBlack);
    ampl_graph[1]->SetMarkerColor(kRed);
    ampl_graph[2]->SetMarkerColor(kBlue);

    phase_graph[0]->SetMarkerColor(kBlack);
    phase_graph[1]->SetMarkerColor(kRed);
    phase_graph[2]->SetMarkerColor(kBlue);

    ampl_graph[0]->SetLineColor(kBlack);
    ampl_graph[1]->SetLineColor(kRed);
    ampl_graph[2]->SetLineColor(kBlue);

    phase_graph[0]->SetLineColor(kBlack);
    phase_graph[1]->SetLineColor(kRed);
    phase_graph[2]->SetLineColor(kBlue);

    auto multi_ampl{new TMultiGraph};
    auto multi_phase{new TMultiGraph};
    for (int i = 0; i != 3; ++i)
    {
        multi_ampl->Add(ampl_graph[i]);
        std::cout << std::endl;

        multi_phase->Add(phase_graph[i]);
    }
    // TCanvas *amplitude_canvas{new TCanvas{"amplitude_canvas", "amplitude", 0, 0, 800, 600}};
    // multi_ampl->Draw("ape");

    // TCanvas *phase_canvas{new TCanvas{"phase_canvas", "phase", 0, 0, 800, 600}};
    // multi_phase->Draw("ape");

    TCanvas *result_canvas{new TCanvas{"result_canvas", "amplitude and phase", 0, 0, 1300, 700}};
    result_canvas->Divide(2, 1);
    result_canvas->cd(1);
    multi_ampl->Draw("ape");
    result_canvas->cd(2);
    multi_phase->Draw("ape");
}
