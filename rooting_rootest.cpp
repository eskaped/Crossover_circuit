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
#include "TLegend.h"
#include "TPaveText.h"
#include "TVirtualFFT.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "TH1D.h"
#include "TFile.h"
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <format>

// da https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
#include <sstream>
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

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

std::string GetSweepRange()
{

    std::ifstream settings_file{"./input_data/settings.txt"};
    if (!settings_file.is_open())
        std::cout << "Failed to open settings.txt" << std::endl;

    std::string input;
    do
    {
        settings_file >> input;
    } while (input != "min:" && !settings_file.eof());

    settings_file >> input;
    int min_freq = atoi(input.c_str());
    std::cout << "min: " << min_freq << "Hz\n";
    settings_file.close();
    std::string output;
    switch (min_freq)
    {
    case 10:
        output = "bassi";
        break;
    case 300:
        output = "medi";
        break;
    case 900:
        output = "alti";
        break;
    case 2000:
        output = "altissimi";
        break;
    case 575:
    case 593:
    case 580:
        output = "crossover";
        break;
    default:
        output = "altri";
        break;
    }

    return output;
}
bool isCrossover()
{
    return GetSweepRange() == "crossover";
}
std::string NumErrScien(double x, double x_err, std::string name_par)
{
    std::string x_err_str = to_string_with_precision(x_err, 15); // s == "1e+05"
    if (x_err_str.find('.') == x_err_str.npos)
        x_err_str += ".0";
    int n_first_non_zero_digit = 0;
    while (n_first_non_zero_digit != x_err_str.size())
    {
        if (x_err_str[n_first_non_zero_digit] != '0' && x_err_str[n_first_non_zero_digit] != '.')
            break;
        ++n_first_non_zero_digit;
    }
    int dot_position = x_err_str.find('.');
    int exp = n_first_non_zero_digit - dot_position;
    if (exp < 0)
        ++exp;

    x_err = x_err * std::pow(10, exp);
    x = x * std::pow(10, exp);

    std::string x_str = to_string_with_precision(x, 15);
    if (x_str.find('.') == x_str.npos)
        x_str += ".0";
    x_err_str = to_string_with_precision(x_err, 15);
    if (x_err_str.find('.') == x_err_str.npos)
        x_err_str += ".0";
    x_str = x_str.substr(0, x_str.find('.') + 2);
    x_err_str = x_err_str.substr(0, 3);

    exp *= -1;
    std::string um_str; // unità di misura
    switch (name_par[0])
    {
    case 'R':
        um_str = "#Omega"; // ohm
        break;
    case 'L':
        um_str = "H"; // henry
        break;
    case 'C':
        um_str = "F"; // farad
        break;
    default:
        um_str = "";
        break;
    }

    std::string um_prefix_str = "";
    if (um_str != "" && exp >= -12 && exp <= 12)
    {
        int prefix = 0;
        while (exp <= -3)
        {
            exp += 3;
            --prefix;
        }
        while (exp >= 3)
        {
            exp -= 3;
            ++prefix;
        }
        switch (prefix)
        {
        case 4:
            um_prefix_str = "T";
            break;
        case 3:
            um_prefix_str = "G";
            break;
        case 2:
            um_prefix_str = "M";
            break;
        case 1:
            um_prefix_str = "k";
            break;
        case 0:
            um_prefix_str = "";
            break;
        case -1:
            um_prefix_str = "m"; // milli
            break;
        case -2:
            um_prefix_str = "#mu"; // micro
            break;
        case -3:
            um_prefix_str = "n"; // nano
            break;
        case -4:
            um_prefix_str = "p";
            break;
        }
    }

    std::string return_str = "(" + x_str + " #pm " + x_err_str + ")";
    if (exp != 0)
        return_str += "E" + std::to_string(exp);
    return return_str + um_prefix_str + um_str;
}

void SamuStyle()
{
    gStyle->SetCanvasPreferGL();
    // gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    // gStyle->SetStatFont(13);
    // gStyle->SetTextFont(13);
    // gStyle->SetTitleFont(13);
    // gStyle->SetLegendFont(13);
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

TF1 *ampl_func_VS;

// par[0] -> Rw
// par[1] -> Rl
// par[2] -> L

Double_t ampl_Vs(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Rt{par[3]};
    Double_t Rl1Rl2{par[4]};
    Double_t C1C2{par[5]};
    Double_t V0{par[6]};
    Double_t Relvis{par[7]};

    Double_t WTA{(Rl + Rw) / ((Rl + Rw) * (Rl + Rw) + (w * L) * (w * L))};
    Double_t WTB{(Rt + Rl1Rl2) / ((Rt + Rl1Rl2) * (Rt + Rl1Rl2) + (1 / (w * C1C2)) * (1 / (w * C1C2)))};
    Double_t WTC1{(-w * L) / ((Rl + Rw) * (Rl + Rw) + (w * L) * (w * L))};
    Double_t WTC2{(1 / (w * C1C2)) / ((Rt + Rl1Rl2) * (Rt + Rl1Rl2) + (1 / (w * C1C2)) * (1 / (w * C1C2)))};

    Double_t ReWT = (WTA + WTB) / ((WTA + WTB) * (WTA + WTB) + (WTC1 + WTC2) * (WTC1 + WTC2));
    Double_t ImWT = -(WTC1 + WTC2) / ((WTA + WTB) * (WTA + WTB) + (WTC1 + WTC2) * (WTC1 + WTC2));

    Double_t ReH = 1 - Relvis * ((Relvis + ReWT) / ((Relvis + ReWT) * (Relvis + ReWT) + ImWT * ImWT));
    Double_t ImH = -Relvis * ((ImWT) / ((Relvis + ReWT) * (Relvis + ReWT) + ImWT * ImWT));

    return V0 * std::sqrt(ReH * ReH + ImH * ImH);
}
Double_t ampl_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < ampl_graph[0]->GetPointX(0))
    //     Vs = ampl_graph[0]->GetPointY(0);
    // else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
    //     Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Vs{ampl_func_VS->Eval(f[0])};
    // // clamp limits
    // if (f[0] < ampl_graph[0]->GetPointX(0))
    //     Vs = ampl_graph[0]->GetPointY(0);
    // else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
    //     Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    return (Vs * Rw) / sqrt((Rl + Rw) * (Rl + Rw) + (w * L) * (w * L));
}

Double_t ampl_woofer_2(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < ampl_graph[0]->GetPointX(0))
    //     Vs = ampl_graph[0]->GetPointY(0);
    // else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
    //     Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Vs{ampl_func_VS->Eval(f[0])};

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Cl{par[3]};
    Double_t Num{(1 - w * w * L * Cl) * (1 - w * w * L * Cl) + (Rl * w * Cl) * (Rl * w * Cl)};
    Double_t Den_A{Rw + Rl - 2 * Rw * w * w * L * Cl + Rw * (w * w * L * Cl) * (w * w * L * Cl) + Rw * (Rl * w * Cl) * (Rl * w * Cl)};
    Double_t Den_B{w * L - w * w * w * L * L * Cl - Rl * Rl * w * Cl};
    return Vs * Rw * Num / sqrt(Den_A * Den_A + Den_B * Den_B);
}

Double_t ampl_woofer_3(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < ampl_graph[0]->GetPointX(0))
    //     Vs = ampl_graph[0]->GetPointY(0);
    // else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
    //     Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Vs{ampl_func_VS->Eval(f[0])};

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Cl{par[3]};
    Double_t Den_A{Rw + Rl};
    Double_t Den_B{(w * L) / (1 - w * w * L * Cl)};
    return Vs * Rw / sqrt(Den_A * Den_A + Den_B * Den_B);
}

Double_t ampl_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t Vs{ampl_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < ampl_graph[0]->GetPointX(0))
    //     Vs = ampl_graph[0]->GetPointY(0);
    // else if (f[0] > ampl_graph[0]->GetPointX(ampl_graph[0]->GetN() - 1))
    //     Vs = ampl_graph[0]->GetPointY(ampl_graph[0]->GetN() - 1);

    Double_t Vs{ampl_func_VS->Eval(f[0])};

    Double_t Rt{par[0]};
    Double_t Rl1Rl2{par[1]};
    Double_t C1C2{par[2]};
    return (Vs * Rt) / sqrt((Rt + Rl1Rl2) * (Rt + Rl1Rl2) + (1 / (w * C1C2)) * (1 / (w * C1C2)));
}

Double_t phase_woofer(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < phase_graph[0]->GetPointX(0))
    //     phase_s = phase_graph[0]->GetPointY(0);
    // else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
    //     phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    return -(std::atan(w * L / (Rl + Rw)));
}

Double_t phase_woofer_2(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < phase_graph[0]->GetPointX(0))
    //     phase_s = phase_graph[0]->GetPointY(0);
    // else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
    //     phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Cl{par[3]};
    Double_t Num{w * L - w * w * w * L * L * Cl - Rl * Rl * w * Cl};
    Double_t Den{(Rw + Rl - Rw * w * w * L * Cl) * (1 - w * w * L * Cl) + (w * L + Rl * Rw * w * Cl) * (Rl * w * Cl)};
    return -(std::atan(Num / Den));
}

Double_t phase_woofer_3(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < phase_graph[0]->GetPointX(0))
    //     phase_s = phase_graph[0]->GetPointY(0);
    // else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
    //     phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rw{par[0]};
    Double_t Rl{par[1]};
    Double_t L{par[2]};
    Double_t Cl{par[3]};
    Double_t Num{-w * L / (1 - w * w * L * Cl)};
    Double_t Den{Rw + Rl};
    return -(std::atan(Num / Den));
}

Double_t phase_tweeter(Double_t *f, Double_t *par)
{
    Double_t w{TMath::TwoPi() * f[0]};
    // Double_t phase_s{phase_graph[0]->Eval(f[0])};
    // // clamp limits
    // if (f[0] < phase_graph[0]->GetPointX(0))
    //     phase_s = phase_graph[0]->GetPointY(0);
    // else if (f[0] > phase_graph[0]->GetPointX(phase_graph[0]->GetN() - 1))
    //     phase_s = phase_graph[0]->GetPointY(phase_graph[0]->GetN() - 1);

    Double_t Rt{par[0]};
    Double_t Rl1Rl2{par[1]};
    Double_t C1C2{par[2]};
    return -(-std::atan(1 / (w * C1C2 * (Rt + Rl1Rl2))));
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
    V_multi->Add(V_arr[0]);
    V_multi->Add(V_arr[1]);
    V_multi->Add(V_arr[2]);

    TCanvas *test_canva = new TCanvas("test_canva", std::to_string(n_block).c_str(), 0, 0, 800, 600);
    V_multi->Draw("AP");
    func_arr[0]->Draw("SAME");
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
            func_arr[i]->SetParameter(2, phase);

            // func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
            // func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
            // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);

            func_arr[i]->SetNumberFitPoints(10000);
            func_arr[i]->SetNpx(10000);

            if ((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0) // failed -> Try setting phase + Pi
            {
                std::cout << "Invalid Fit! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";
                std::cout << "\tTrying again with phase + Pi\n";

                phase += TMath::Pi();

                func_arr[i]->SetParameter(0, amplitude);
                func_arr[i]->SetParameter(1, pulsation);
                func_arr[i]->SetParameter(2, phase);

                // func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
                // func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
                // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);

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
    Double_t MIN_FIT = 500;
    Double_t MAX_FIT = 700;



    N_BLOCKS = input_n_blocks;
    SamuStyle();
    std::ifstream file_count{"./input_data/ampl_data"};

    Double_t const Rw{2.2032000E+02};
    Double_t const Rt{2.2019000E+02};
    Double_t const Rl{1.2018250E+02};     // Actual L's resistance
    Double_t const Rl1Rl2{1.2183500E+02}; // equivalent L's resistance on tweeter
    Double_t const C1C2{1.4225000E-06};
    Double_t const L{4.7154000E-02};
    Double_t const Relvis{5.4406400E+01};

    Double_t V0{(GetVoltage() != 7 ? GetVoltage() : 7.5) / 2.};

    Double_t const Rw_err{1.4737707E-01};
    Double_t const Rt_err{1.0639549E-01};
    Double_t const Rl_err{4.8897001E-01};     // Actual L's resistance
    Double_t const Rl1Rl2_err{4.5793013E-01}; // equivalent L's resistance on tweeter
    Double_t const C1C2_err{1.4225000E-08};
    Double_t const L_err{4.7154000E-04};
    Double_t const Relvis_err{4.8626319E+00};

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
            new TF1{"V_s_fit", "[0]*cos([1]*x - [2]) + [3]"},
            new TF1{"V_w_fit", "[0]*cos([1]*x - [2]) + [3]"},
            new TF1{"V_t_fit", "[0]*cos([1]*x - [2]) + [3]"}};

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
            double background = -0.001;

            double amplitude;
            double phase;
            file_params_in >> amplitude;
            file_params_in >> phase;
            phase = phase * TMath::Pi() / 180.; // to rad

            phase -= TMath::Pi() / 2.; // to account for shifting in converter

            double pulsation = 2 * TMath::Pi() * frequency;

            func_arr[i]->SetParameter(0, amplitude);
            func_arr[i]->SetParameter(1, pulsation);
            func_arr[i]->SetParameter(2, phase);
            func_arr[i]->SetParameter(3, background);

            func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
            func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
            // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);

            func_arr[i]->SetNumberFitPoints(10000);
            func_arr[i]->SetNpx(10000);

            if ((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0) // failed -> Try setting phase + Pi
            {
                std::cout << "Invalid Fit! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";
                std::cout << "\tTrying again with phase + Pi\n";

                phase += TMath::Pi();

                func_arr[i]->SetParameter(0, amplitude);
                func_arr[i]->SetParameter(1, pulsation);
                func_arr[i]->SetParameter(2, phase);
                func_arr[i]->SetParameter(3, background);

                func_arr[i]->SetParLimits(0, amplitude - amplitude / 10., amplitude + amplitude / 10.);
                func_arr[i]->SetParLimits(1, pulsation - pulsation / 10., pulsation + pulsation / 10.);
                // func_arr[i]->SetParLimits(2, phase - phase / 10., phase + phase / 10.);
                if (((V_arr[i]->Fit(func_arr[i], "QUIET")) != 0)) // i.e.: error
                {
                    std::cout << "Invalid Fit AGAIN! Block n: " << n_block << ", graph: " << i << " [0,1,2=V_s,V_w,V_t]\n";

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
                else // second fit went well
                {
                    // only V_s freq

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

    //  subtract phase_Vs from all phase graphs
    for (int n_block = 0; n_block != N_BLOCKS; ++n_block)
    {
        Double_t phase_Vs = phase_graph[0]->GetPointY(n_block);
        Double_t phase_Vs_err = phase_graph[0]->GetErrorY(n_block);
        phase_graph[0]->SetPointY(n_block, -(phase_graph[0]->GetPointY(n_block) - phase_Vs));
        phase_graph[1]->SetPointY(n_block, -(phase_graph[1]->GetPointY(n_block) - phase_Vs));
        phase_graph[2]->SetPointY(n_block, -(phase_graph[2]->GetPointY(n_block) - phase_Vs));
        phase_graph[0]->SetPointError(n_block, phase_graph[0]->GetErrorX(n_block), 0.);
        phase_graph[1]->SetPointError(n_block, phase_graph[1]->GetErrorX(n_block),
                                      std::sqrt(phase_graph[1]->GetErrorY(n_block) * phase_graph[1]->GetErrorY(n_block) + phase_Vs_err * phase_Vs_err));
        phase_graph[2]->SetPointError(n_block, phase_graph[2]->GetErrorX(n_block),
                                      std::sqrt(phase_graph[2]->GetErrorY(n_block) * phase_graph[2]->GetErrorY(n_block) + phase_Vs_err * phase_Vs_err));
    }

    // fitting

    ampl_func_VS = new TF1{"ampl_func_VS", ampl_Vs, 0., 1000., 8};
    // ampl_func_VS = new TF1{"ampl_func_VS", "[0]*x+[1]", 0., 1000.};

    TF1 *ampl_func_w{new TF1{"ampl_func_w", ampl_woofer, 0., 25000., 3}};
    TF1 *ampl_func_t{new TF1{"ampl_func_t", ampl_tweeter, 0., 25000., 3}};
    TF1 *phase_func_w{new TF1{"phase_func_w", phase_woofer, 0., 25000., 3}};
    TF1 *phase_func_t{new TF1{"phase_func_t", phase_tweeter, 0., 25000., 3}};

    // with cl------------------------------------------------
    TF1 *ampl_func_w_2{new TF1{"ampl_func_w", ampl_woofer_2, 0., 25000., 3}};
    TF1 *phase_func_w_2{new TF1{"phase_func_w", phase_woofer_2, 0., 25000., 3}};

    TF1 *ampl_func_w_3{new TF1{"ampl_func_w", ampl_woofer_2, 0., 25000., 3}};
    TF1 *phase_func_w_3{new TF1{"phase_func_w", phase_woofer_2, 0., 25000., 3}};
    // with cl------------------------------------------------

    ampl_func_VS->SetNpx(100000);
    ampl_func_VS->SetNumberFitPoints(100000);

    ampl_func_w->SetNpx(100000);
    ampl_func_w->SetNumberFitPoints(100000);
    ampl_func_t->SetNpx(100000);
    ampl_func_t->SetNumberFitPoints(100000);
    phase_func_w->SetNpx(100000);
    phase_func_w->SetNumberFitPoints(100000);
    phase_func_t->SetNpx(100000);
    phase_func_t->SetNumberFitPoints(100000);

    // with cl------------------------------------------------
    ampl_func_w_2->SetNpx(100000);
    ampl_func_w_3->SetNpx(100000);
    phase_func_w_2->SetNpx(100000);
    phase_func_w_3->SetNpx(100000);
    ampl_func_w_2->SetNumberFitPoints(100000);
    ampl_func_w_3->SetNumberFitPoints(100000);
    phase_func_w_2->SetNumberFitPoints(100000);
    phase_func_w_3->SetNumberFitPoints(100000);
    // with cl------------------------------------------------

    ampl_func_VS->SetParName(0, "Rw");
    ampl_func_VS->SetParName(1, "Rl");
    ampl_func_VS->SetParName(2, "L");
    ampl_func_VS->SetParName(3, "Rt");
    ampl_func_VS->SetParName(4, "Rl1Rl2");
    ampl_func_VS->SetParName(5, "C1C2");
    ampl_func_VS->SetParName(6, "V0");
    ampl_func_VS->SetParName(7, "Relvis");

    ampl_func_w->SetParName(0, "Rw");
    ampl_func_w->SetParName(1, "Rl");
    ampl_func_w->SetParName(2, "L");

    ampl_func_t->SetParName(0, "Rt");
    ampl_func_t->SetParName(1, "Rl1Rl2");
    ampl_func_t->SetParName(2, "C1C2");

    phase_func_w->SetParName(0, "Rw");
    phase_func_w->SetParName(1, "Rl");
    phase_func_w->SetParName(2, "L");

    phase_func_t->SetParName(0, "Rt");
    phase_func_t->SetParName(1, "Rl1Rl2");
    phase_func_t->SetParName(2, "C1C2");

    // with cl------------------------------------------------

    ampl_func_w_2->SetParName(0, "Rw");
    ampl_func_w_2->SetParName(1, "Rl");
    ampl_func_w_2->SetParName(2, "L");
    ampl_func_w_2->SetParName(3, "Cl");

    ampl_func_w_3->SetParName(0, "Rw");
    ampl_func_w_3->SetParName(1, "Rl");
    ampl_func_w_3->SetParName(2, "L");
    ampl_func_w_3->SetParName(3, "Cl");

    phase_func_w_2->SetParName(0, "Rw");
    phase_func_w_2->SetParName(1, "Rl");
    phase_func_w_2->SetParName(2, "L");
    phase_func_w_2->SetParName(3, "Cl");

    phase_func_w_3->SetParName(0, "Rw");
    phase_func_w_3->SetParName(1, "Rl");
    phase_func_w_3->SetParName(2, "L");
    phase_func_w_3->SetParName(3, "Cl");
    // with cl------------------------------------------------

    ampl_func_VS->SetParameter(0, Rw);
    ampl_func_VS->SetParameter(1, Rl);
    ampl_func_VS->SetParameter(2, L);
    ampl_func_VS->SetParameter(3, Rt);
    ampl_func_VS->SetParameter(4, Rl1Rl2);
    ampl_func_VS->SetParameter(5, C1C2);
    ampl_func_VS->SetParameter(6, V0);
    ampl_func_VS->SetParameter(7, Relvis);

    ampl_func_w->SetParameter(0, Rw);
    ampl_func_w->SetParameter(1, Rl);
    ampl_func_w->SetParameter(2, L);

    ampl_func_t->SetParameter(0, Rl);
    ampl_func_t->SetParameter(1, Rl1Rl2);
    ampl_func_t->SetParameter(2, C1C2);

    phase_func_w->SetParameter(0, Rw);
    phase_func_w->SetParameter(1, Rl);
    phase_func_w->SetParameter(2, L);

    phase_func_t->SetParameter(0, Rl);
    phase_func_t->SetParameter(1, Rl1Rl2);
    phase_func_t->SetParameter(2, C1C2);

    // with cl------------------------------------------------
    ampl_func_w_2->SetParameter(0, Rw);
    ampl_func_w_2->SetParameter(1, Rl);
    ampl_func_w_2->SetParameter(2, L);
    ampl_func_w_2->SetParameter(3, 1E-8);

    ampl_func_w_3->SetParameter(0, Rw);
    ampl_func_w_3->SetParameter(1, Rl);
    ampl_func_w_3->SetParameter(2, L);
    ampl_func_w_3->SetParameter(3, 1E-8);

    phase_func_w_2->SetParameter(0, Rw);
    phase_func_w_2->SetParameter(1, Rl);
    phase_func_w_2->SetParameter(2, L);
    phase_func_w_2->SetParameter(3, 1E-8);

    phase_func_w_3->SetParameter(0, Rw);
    phase_func_w_3->SetParameter(1, Rl);
    phase_func_w_3->SetParameter(2, L);
    phase_func_w_3->SetParameter(3, 1E-8);
    // with cl------------------------------------------------

    // GRAPHICS-----------------------------------------------
    ampl_func_VS->SetLineColor(kGray + 3);

    ampl_func_t->SetLineColor(kBlue + 4);
    ampl_func_w->SetLineColor(kRed + 2);

    phase_func_t->SetLineColor(kBlue + 4);
    phase_func_w->SetLineColor(kRed + 2);

    // with cl------------------------------------------------
    ampl_func_w_2->SetLineColor(kRed + 2);
    ampl_func_w_3->SetLineColor(kRed + 2);
    phase_func_w_2->SetLineColor(kRed + 2);
    phase_func_w_3->SetLineColor(kRed + 2);
    // with cl------------------------------------------------

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

    // FITTING VS
    std::cout << "\n\n\n\n";
    // Double_t N_SIGMA_VS = 10;
    // ampl_func_VS->SetParLimits(0, Rw - N_SIGMA_VS * Rw_err, Rw + N_SIGMA_VS * Rw_err);
    // ampl_func_VS->SetParLimits(1, Rl - N_SIGMA_VS * Rl_err, Rl + N_SIGMA_VS * Rl_err);
    // ampl_func_VS->SetParLimits(2, L - N_SIGMA_VS * L_err, L + N_SIGMA_VS * L_err);
    // ampl_func_VS->SetParLimits(3, Rt - N_SIGMA_VS * Rt_err, Rt + N_SIGMA_VS * Rt_err);
    // ampl_func_VS->SetParLimits(4, Rl1Rl2 - N_SIGMA_VS * Rl1Rl2_err, Rl1Rl2_err + N_SIGMA_VS * Rl1Rl2_err);
    // ampl_func_VS->SetParLimits(5, C1C2 - N_SIGMA_VS * C1C2_err, C1C2 + N_SIGMA_VS * C1C2_err);
    // ampl_func_VS->SetParLimits(6, V0 - V0/N_SIGMA_VS, V0 + V0/N_SIGMA_VS);
    // ampl_func_VS->SetParLimits(7, Relvis - N_SIGMA_VS * Relvis_err, Relvis + N_SIGMA_VS * Relvis_err);

    // ampl_func_VS->FixParameter(0, Rw);
    // ampl_func_VS->FixParameter(1, Rl);
    // ampl_func_VS->FixParameter(2, L );
    // ampl_func_VS->FixParameter(3, Rt );
    // ampl_func_VS->FixParameter(4, Rl1Rl2);
    // ampl_func_VS->FixParameter(5, C1C2 );
    // ampl_func_VS->FixParameter(6, V0 );
    // ampl_func_VS->FixParameter(7, Relvis);
    // ampl_func_VS->FixParameter()
    std::cout << "------------------------START FITTING VS--------------------------------\n";
    ampl_graph[0]->Fit(ampl_func_VS, "M, E");
    std::cout << "------------------------END FITTING VS----------------------------------\n";

    // FITTING WITHOUT LIMITS (No Cl)-------------------------
    // 1/(2pi)*sqrt((C*C*(R_w*R_w*(R_t+R_l1l2)*(R_t+R_l1l2) - R_t*R_t*(R_l+R_w)*(R_l+R_w)) + sqrt(C*C*C*C*(R_w*R_w*(R_t+R_l1l2)*(R_t+R_l1l2) - R_t*R_t*(R_l+R_w)*(R_l+R_w))*(R_w*R_w*(R_t+R_l1l2)*(R_t+R_l1l2) - R_t*R_t*(R_l+R_w)*(R_l+R_w)) + 4 *R_w*R_w*R_t*R_t*L*L*C*C ))/(2*L*L*C*C*R_t*R_t))
    //R_w, R_t, R_l, R_l1l2, L, C
    //[220.32000,0.14737707], [220.19000, 0.10639549], [120.18250, 0.48897001], [121.83500, 0.45793013],[0.047154000, 0.00047154000],[0.0000014225000, 0.000000014225000]
    
    //[192.0176023,0.109684324461364], [155.7258872, 0.0430284572390734], [108.6318597, 0.0721233730001343], [87.27564753, 0.0279570867606535],[0.04416795751, .000031079980085238],[0.000001975545041,.000000000671207244439028]
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    int ampl_woofer_fit_res{ampl_graph[1]->Fit(ampl_func_w, "Q, M, E", "", MIN_FIT, MAX_FIT)};
    std::cout << std::endl;
    int ampl_tweeter_fit_res{ampl_graph[2]->Fit(ampl_func_t, "Q, M, E", "", MIN_FIT, MAX_FIT)};
    std::cout << std::endl;

    int phase_woofer_fit_res{phase_graph[1]->Fit(phase_func_w, "Q, M, E","", MIN_FIT, MAX_FIT)};
    std::cout << std::endl;
    int phase_tweeter_fit_res{phase_graph[2]->Fit(phase_func_t, "Q, M, E","", MIN_FIT, MAX_FIT)};
    std::cout << std::endl;

    auto multi_ampl_no_lim_no_Cl{new TMultiGraph};
    auto multi_phase_no_lim_no_Cl{new TMultiGraph};

    // skip Vs if we're in the crossover fitting
    for (int i = (isCrossover() ? 1 : 0); i != 3; ++i)
    {
        multi_ampl_no_lim_no_Cl->Add(ampl_graph[i]);
    }

    for (int i = 0; i != 3; ++i)
    {
        multi_phase_no_lim_no_Cl->Add(phase_graph[i]);
    }

    multi_ampl_no_lim_no_Cl->SetTitle("Amplitude - Frequency (No Limits, No Cl)");
    multi_phase_no_lim_no_Cl->SetTitle("Phase - Frequency (No Limits, No Cl)");

    // multi_ampl_no_lim_no_Cl->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_no_lim_no_Cl->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_no_lim_no_Cl->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_no_lim_no_Cl->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_no_lim_no_Cl->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_no_lim_no_Cl->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_no_lim_no_Cl->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *no_limits_no_cl_canvas{new TCanvas{"no_limits_no_cl_canvas", "Amplitude and Phase No_limits_no_Cl", 0, 0, 1300, 700}};
    no_limits_no_cl_canvas->Divide(2, 1);
    no_limits_no_cl_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_ampl_no_lim_no_Cl->Draw("ape");
    TF1* ampl_func_w_ghost = new TF1(*ampl_func_w);
    TF1* ampl_func_t_ghost = new TF1(*ampl_func_t);
    ampl_func_w_ghost->SetLineColorAlpha(kRed+2,0.1);
    ampl_func_t_ghost->SetLineColorAlpha(kBlue+2,0.1);
    ampl_func_w_ghost->SetRange(0., 25000);
    ampl_func_t_ghost->SetRange(0., 25000);
    ampl_func_w_ghost->Draw("SAME");
    ampl_func_t_ghost->Draw("SAME");
    // ampl_func_w->Draw("Same");
    // ampl_func_t->Draw("Same");

    TPad *no_limits_no_cl_pad_residual_ampl = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_no_cl_pad_residual_ampl->SetTopMargin(0.8);
    no_limits_no_cl_pad_residual_ampl->Draw();
    no_limits_no_cl_pad_residual_ampl->SetFillStyle(0);
    no_limits_no_cl_pad_residual_ampl->cd();
    no_limits_no_cl_pad_residual_ampl->SetGridy();

    TGraphErrors *no_limits_no_cl_pad_ampl_w = new TGraphErrors(ampl_graph[1]->GetN());
    for (Int_t i = 0; i != ampl_graph[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph[1]->GetPointX(i);
        Double_t y = ampl_graph[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w->Eval(x);
        no_limits_no_cl_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_no_cl_pad_ampl_t = new TGraphErrors(ampl_graph[2]->GetN());
    for (Int_t i = 0; i != ampl_graph[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph[2]->GetPointX(i);
        Double_t y = ampl_graph[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t->Eval(x);
        no_limits_no_cl_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_no_cl_pad_ampl_w->SetLineColor(kRed + 2);
    no_limits_no_cl_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_no_cl_pad_multi_ampl = new TMultiGraph();
    no_limits_no_cl_pad_multi_ampl->Add(no_limits_no_cl_pad_ampl_w);
    no_limits_no_cl_pad_multi_ampl->Add(no_limits_no_cl_pad_ampl_t);
    no_limits_no_cl_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    no_limits_no_cl_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    no_limits_no_cl_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_no_cl_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    no_limits_no_cl_pad_multi_ampl->Draw("apl");

    no_limits_no_cl_canvas->cd(1);

    TLegend *no_limits_no_cl_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_no_cl_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_no_cl_legend_ampl->SetNColumns(3);
    // RIGA 1
    no_limits_no_cl_legend_ampl->AddEntry(ampl_graph[0], "Amplitude V_S", "ep");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_no_cl_legend_ampl->AddEntry(ampl_graph[1], "Amplitude V_Woofer", "ep");
    no_limits_no_cl_legend_ampl->AddEntry(ampl_func_w, "Amplitude V_Woofer Fit", "l");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w->GetChisquare() / ampl_func_w->GetNDF())).c_str()), "");
    // RIGA 3
    for (int i = 0; i != ampl_func_w->GetNpar(); ++i)
        no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w->GetParName(i)) + " = " + NumErrScien(ampl_func_w->GetParameter(i), ampl_func_w->GetParError(i), ampl_func_w->GetParName(i))).c_str(), "");
    // RIGA 4
    no_limits_no_cl_legend_ampl->AddEntry(ampl_graph[2], "Amplitude V_Tweeter", "ep");
    no_limits_no_cl_legend_ampl->AddEntry(ampl_func_t, "Amplitude V_Tweeter Fit", "l");
    no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t->GetChisquare() / ampl_func_t->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t->GetParName(i)) + " = " + NumErrScien(ampl_func_t->GetParameter(i), ampl_func_t->GetParError(i), ampl_func_t->GetParName(i))).c_str(), "");

    no_limits_no_cl_legend_ampl->Draw();

    no_limits_no_cl_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_phase_no_lim_no_Cl->Draw("ape");
    TF1* phase_func_w_ghost = new TF1(*phase_func_w);
    TF1* phase_func_t_ghost = new TF1(*phase_func_t);
    phase_func_w_ghost->SetLineColorAlpha(kRed+2,0.1);
    phase_func_t_ghost->SetLineColorAlpha(kBlue+2,0.1);
    phase_func_w_ghost->SetRange(0., 25000);
    phase_func_t_ghost->SetRange(0., 25000);
    phase_func_w_ghost->Draw("SAME");
    phase_func_t_ghost->Draw("SAME");

    TPad *no_limits_no_cl_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_no_cl_pad_residual_phase->SetTopMargin(0.8);
    no_limits_no_cl_pad_residual_phase->Draw();
    no_limits_no_cl_pad_residual_phase->SetFillStyle(0);
    no_limits_no_cl_pad_residual_phase->cd();
    no_limits_no_cl_pad_residual_phase->SetGridy();

    TGraphErrors *no_limits_no_cl_pad_phase_w = new TGraphErrors(phase_graph[1]->GetN());
    for (Int_t i = 0; i != phase_graph[1]->GetN(); ++i)
    {
        Double_t x = phase_graph[1]->GetPointX(i);
        Double_t y = phase_graph[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w->Eval(x);
        no_limits_no_cl_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_no_cl_pad_phase_t = new TGraphErrors(phase_graph[2]->GetN());
    for (Int_t i = 0; i != phase_graph[2]->GetN(); ++i)
    {
        Double_t x = phase_graph[2]->GetPointX(i);
        Double_t y = phase_graph[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t->Eval(x);
        no_limits_no_cl_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_no_cl_pad_phase_w->SetLineColor(kRed + 2);
    no_limits_no_cl_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_no_cl_pad_multi_phase = new TMultiGraph();
    no_limits_no_cl_pad_multi_phase->Add(no_limits_no_cl_pad_phase_w);
    no_limits_no_cl_pad_multi_phase->Add(no_limits_no_cl_pad_phase_t);
    no_limits_no_cl_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    no_limits_no_cl_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    no_limits_no_cl_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_no_cl_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    no_limits_no_cl_pad_multi_phase->Draw("apl");

    no_limits_no_cl_canvas->cd(2);
    // TLegend *no_limits_no_cl_legend_phasel{new TLegend()};
    // no_limits_no_cl_legend_phasel->SetMargin(0.05);
    // no_limits_no_cl_legend_phasel->SetEntrySeparation(0);
    // no_limits_no_cl_legend_phasel->SetHeader("Phase Shift - Frequency", "C"); // option "C" allows to center the header
    // no_limits_no_cl_legend_phasel->AddEntry(phase_graph[0], "Phase Shift V_S", "ep");
    // no_limits_no_cl_legend_phasel->AddEntry(phase_graph[1], "Phase Shift V_Woofer", "ep");
    // no_limits_no_cl_legend_phasel->AddEntry(phase_graph[2], "Phase Shift V_Tweeter", "ep");
    // no_limits_no_cl_legend_phasel->AddEntry(phase_func_w, "Phase Shift V_Woofer Fit", "l");
    // no_limits_no_cl_legend_phasel->AddEntry(phase_func_t, "Phase Shift V_Tweeter Fit", "l");
    // no_limits_no_cl_legend_phasel->Draw();

    TLegend *no_limits_no_cl_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_no_cl_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_no_cl_legend_phase->SetNColumns(3);
    // RIGA 1
    no_limits_no_cl_legend_phase->AddEntry(phase_graph[0], "Phase V_S", "ep");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_no_cl_legend_phase->AddEntry(phase_graph[1], "Phase V_Woofer", "ep");
    no_limits_no_cl_legend_phase->AddEntry(phase_func_w, "Phase V_Woofer Fit", "l");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w->GetChisquare() / phase_func_w->GetNDF())).c_str()), "");
    // RIGA 3
    for (int i = 0; i != phase_func_w->GetNpar(); ++i)
        no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w->GetParName(i)) + " = " + NumErrScien(phase_func_w->GetParameter(i), phase_func_w->GetParError(i), phase_func_w->GetParName(i))).c_str(), "");
    // RIGA 4
    no_limits_no_cl_legend_phase->AddEntry(phase_graph[2], "Phase V_Tweeter", "ep");
    no_limits_no_cl_legend_phase->AddEntry(phase_func_t, "Phase V_Tweeter Fit", "l");
    no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t->GetChisquare() / phase_func_t->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t->GetParName(i)) + " = " + NumErrScien(phase_func_t->GetParameter(i), phase_func_t->GetParError(i), phase_func_t->GetParName(i))).c_str(), "");

    no_limits_no_cl_legend_phase->Draw();

    // print_res
    std::ofstream no_limits_no_Cl_file("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_no_Cl_" + std::to_string(GetVoltage()) + ".txt");
    no_limits_no_Cl_file << std::setprecision(10);
    no_limits_no_Cl_file << "Ampl_woofer:\n";
    no_limits_no_Cl_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    no_limits_no_Cl_file << "\tChi^2: " << ampl_func_w->GetChisquare() << '\n';
    no_limits_no_Cl_file << "\tNdf: " << ampl_func_w->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tCHI RIDOTTO: " << ampl_func_w->GetChisquare() / ampl_func_w->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w->GetNpar(); ++i)
        no_limits_no_Cl_file << "\t\t[" << i << "] - " << ampl_func_w->GetParName(i) << " = " << ampl_func_w->GetParameter(i) << " +- " << ampl_func_w->GetParError(i) << '\n';
    no_limits_no_Cl_file << "Ampl_tweeter:\n";
    no_limits_no_Cl_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    no_limits_no_Cl_file << "\tChi^2: " << ampl_func_t->GetChisquare() << '\n';
    no_limits_no_Cl_file << "\tNdf: " << ampl_func_t->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tCHI RIDOTTO: " << ampl_func_t->GetChisquare() / ampl_func_t->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_no_Cl_file << "\t\t[" << i << "] - " << ampl_func_t->GetParName(i) << " = " << ampl_func_t->GetParameter(i) << " +- " << ampl_func_t->GetParError(i) << '\n';
    no_limits_no_Cl_file << "Phase_woofer:\n";
    no_limits_no_Cl_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    no_limits_no_Cl_file << "\tChi^2: " << phase_func_w->GetChisquare() << '\n';
    no_limits_no_Cl_file << "\tNdf: " << phase_func_w->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tCHI RIDOTTO: " << phase_func_w->GetChisquare() / phase_func_w->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w->GetNpar(); ++i)
        no_limits_no_Cl_file << "\t\t[" << i << "] - " << phase_func_w->GetParName(i) << " = " << phase_func_w->GetParameter(i) << " +- " << phase_func_w->GetParError(i) << '\n';
    no_limits_no_Cl_file << "Phase_tweeter:\n";
    no_limits_no_Cl_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    no_limits_no_Cl_file << "\tChi^2: " << phase_func_t->GetChisquare() << '\n';
    no_limits_no_Cl_file << "\tNdf: " << phase_func_t->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tCHI RIDOTTO: " << phase_func_t->GetChisquare() / phase_func_t->GetNDF() << '\n';
    no_limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_no_Cl_file << "\t\t[" << i << "] - " << phase_func_t->GetParName(i) << " = " << phase_func_t->GetParameter(i) << " +- " << phase_func_t->GetParError(i) << '\n';
    no_limits_no_Cl_file.close();

    // Save Canvas
    no_limits_no_cl_canvas->Update();
    no_limits_no_cl_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_no_cl_" + std::to_string(GetVoltage()) + ".png").c_str());

    // NO LIMITS (EXCEPT CL>=0) - Cl 2 -----------------------------------------------------------------
    TGraphErrors *ampl_graph_no_limits_Cl2{new TGraphErrors(*ampl_graph[1])};
    TGraphErrors *phase_graph_no_limits_Cl2{new TGraphErrors(*phase_graph[1])};

    ampl_func_w_2->SetParLimits(3, 0, 1e-6);
    phase_func_w_2->SetParLimits(3, 0, 1e-6);
    ampl_woofer_fit_res = ampl_graph_no_limits_Cl2->Fit(ampl_func_w_2, "Q, M, E, R");
    std::cout << std::endl;

    phase_woofer_fit_res = phase_graph_no_limits_Cl2->Fit(phase_func_w_2, "Q, M, E, R");
    std::cout << std::endl;

    auto multi_ampl_no_lim_Cl_2{new TMultiGraph};
    auto multi_phase_no_lim_Cl_2{new TMultiGraph};

    if (!isCrossover())
        multi_ampl_no_lim_Cl_2->Add(ampl_graph[0]);
    multi_phase_no_lim_Cl_2->Add(phase_graph[0]);
    multi_ampl_no_lim_Cl_2->Add(ampl_graph_no_limits_Cl2);
    multi_phase_no_lim_Cl_2->Add(phase_graph_no_limits_Cl2);
    multi_ampl_no_lim_Cl_2->Add(ampl_graph[2]);
    multi_phase_no_lim_Cl_2->Add(phase_graph[2]);

    // ampl_woofer_fit_res = ampl_graph[1]->Fit(ampl_func_w_2, "Q, M, E, R");
    // std::cout << std::endl;

    // phase_woofer_fit_res = phase_graph[1]->Fit(phase_func_w_2, "Q, M, E, R");
    // std::cout << std::endl;

    // auto multi_ampl_no_lim_Cl_2{new TMultiGraph};
    // auto multi_phase_no_lim_Cl_2{new TMultiGraph};

    // for (int i = 0; i != 3; ++i)
    // {
    //     multi_ampl_no_lim_Cl_2->Add(ampl_graph[i]);
    //     multi_phase_no_lim_Cl_2->Add(phase_graph[i]);
    // }
    multi_ampl_no_lim_Cl_2->SetTitle("Amplitude - Frequency (No Limits, Cl)");
    multi_phase_no_lim_Cl_2->SetTitle("Phase - Frequency (No Limits, Cl)");

    // multi_ampl_no_lim_Cl_2->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_no_lim_Cl_2->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_no_lim_Cl_2->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_no_lim_Cl_2->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_no_lim_Cl_2->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_no_lim_Cl_2->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_no_lim_Cl_2->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *no_limits_Cl_2_canvas{new TCanvas{"no_limits_Cl_2_canvas", "Amplitude and Phase No_limits_Cl_2", 0, 0, 1300, 700}};
    no_limits_Cl_2_canvas->Divide(2, 1);
    no_limits_Cl_2_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_ampl_no_lim_Cl_2->Draw("ape");

    TPad *no_limits_Cl_2_pad_residual_ampl = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_Cl_2_pad_residual_ampl->SetTopMargin(0.8);
    no_limits_Cl_2_pad_residual_ampl->Draw();
    no_limits_Cl_2_pad_residual_ampl->SetFillStyle(0);
    no_limits_Cl_2_pad_residual_ampl->cd();
    no_limits_Cl_2_pad_residual_ampl->SetGridy();

    TGraphErrors *no_limits_Cl_2_pad_ampl_w = new TGraphErrors(ampl_graph[1]->GetN());
    for (Int_t i = 0; i != ampl_graph[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph[1]->GetPointX(i);
        Double_t y = ampl_graph[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w_2->Eval(x);
        no_limits_Cl_2_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_Cl_2_pad_ampl_t = new TGraphErrors(ampl_graph[2]->GetN());
    for (Int_t i = 0; i != ampl_graph[2]->GetN(); ++i)
    {
        Double_t x = ampl_graph[2]->GetPointX(i);
        Double_t y = ampl_graph[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t->Eval(x);
        no_limits_Cl_2_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_Cl_2_pad_ampl_w->SetLineColor(kRed + 2);
    no_limits_Cl_2_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_Cl_2_pad_multi_ampl = new TMultiGraph();
    no_limits_Cl_2_pad_multi_ampl->Add(no_limits_Cl_2_pad_ampl_w);
    no_limits_Cl_2_pad_multi_ampl->Add(no_limits_Cl_2_pad_ampl_t);
    no_limits_Cl_2_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    no_limits_Cl_2_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    no_limits_Cl_2_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_Cl_2_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    no_limits_Cl_2_pad_multi_ampl->Draw("apl");

    no_limits_Cl_2_canvas->cd(1);

    TLegend *no_limits_cl_2_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_cl_2_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_cl_2_legend_ampl->SetNColumns(3);
    // RIGA 1
    no_limits_cl_2_legend_ampl->AddEntry(ampl_graph[0], "Amplitude V_S", "ep");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_cl_2_legend_ampl->AddEntry(ampl_graph_no_limits_Cl2, "Amplitude V_Woofer", "ep");
    no_limits_cl_2_legend_ampl->AddEntry(ampl_func_w_2, "Amplitude V_Woofer Fit", "l");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w_2->GetChisquare() / ampl_func_w_2->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != ampl_func_w_2->GetNpar(); ++i)
        no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w_2->GetParName(i)) + " = " + NumErrScien(ampl_func_w_2->GetParameter(i), ampl_func_w_2->GetParError(i), ampl_func_w_2->GetParName(i))).c_str(), "");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 5
    no_limits_cl_2_legend_ampl->AddEntry(ampl_graph[2], "Amplitude V_Tweeter", "ep");
    no_limits_cl_2_legend_ampl->AddEntry(ampl_func_t, "Amplitude V_Tweeter Fit", "l");
    no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t->GetChisquare() / ampl_func_t->GetNDF())).c_str()), "");
    // RIGA 6
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_cl_2_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t->GetParName(i)) + " = " + NumErrScien(ampl_func_t->GetParameter(i), ampl_func_t->GetParError(i), ampl_func_t->GetParName(i))).c_str(), "");

    no_limits_cl_2_legend_ampl->Draw();

    no_limits_Cl_2_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_phase_no_lim_Cl_2->Draw("ape");

    TPad *no_limits_Cl_2_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_Cl_2_pad_residual_phase->SetTopMargin(0.8);
    no_limits_Cl_2_pad_residual_phase->Draw();
    no_limits_Cl_2_pad_residual_phase->SetFillStyle(0);
    no_limits_Cl_2_pad_residual_phase->cd();
    no_limits_Cl_2_pad_residual_phase->SetGridy();

    TGraphErrors *no_limits_Cl_2_pad_phase_w = new TGraphErrors(phase_graph[1]->GetN());
    for (Int_t i = 0; i != phase_graph[1]->GetN(); ++i)
    {
        Double_t x = phase_graph[1]->GetPointX(i);
        Double_t y = phase_graph[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w_2->Eval(x);
        no_limits_Cl_2_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_Cl_2_pad_phase_t = new TGraphErrors(phase_graph[2]->GetN());
    for (Int_t i = 0; i != phase_graph[2]->GetN(); ++i)
    {
        Double_t x = phase_graph[2]->GetPointX(i);
        Double_t y = phase_graph[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t->Eval(x);
        no_limits_Cl_2_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_Cl_2_pad_phase_w->SetLineColor(kRed + 2);
    no_limits_Cl_2_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_Cl_2_pad_multi_phase = new TMultiGraph();
    no_limits_Cl_2_pad_multi_phase->Add(no_limits_Cl_2_pad_phase_w);
    no_limits_Cl_2_pad_multi_phase->Add(no_limits_Cl_2_pad_phase_t);
    no_limits_Cl_2_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    no_limits_Cl_2_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    no_limits_Cl_2_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_Cl_2_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    no_limits_Cl_2_pad_multi_phase->Draw("apl");

    no_limits_Cl_2_canvas->cd(2);

    TLegend *no_limits_cl_2_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_cl_2_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_cl_2_legend_phase->SetNColumns(3);
    // RIGA 1
    no_limits_cl_2_legend_phase->AddEntry(phase_graph[0], "Phase V_S", "ep");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_cl_2_legend_phase->AddEntry(phase_graph_no_limits_Cl2, "Phase V_Woofer", "ep");
    no_limits_cl_2_legend_phase->AddEntry(phase_func_w_2, "Phase V_Woofer Fit", "l");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w_2->GetChisquare() / phase_func_w_2->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != phase_func_w_2->GetNpar(); ++i)
        no_limits_cl_2_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w_2->GetParName(i)) + " = " + NumErrScien(phase_func_w_2->GetParameter(i), phase_func_w_2->GetParError(i), phase_func_w_2->GetParName(i))).c_str(), "");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 4
    no_limits_cl_2_legend_phase->AddEntry(phase_graph[2], "Phase V_Tweeter", "ep");
    no_limits_cl_2_legend_phase->AddEntry(phase_func_t, "Phase V_Tweeter Fit", "l");
    no_limits_cl_2_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t->GetChisquare() / phase_func_t->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_cl_2_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t->GetParName(i)) + " = " + NumErrScien(phase_func_t->GetParameter(i), phase_func_t->GetParError(i), phase_func_t->GetParName(i))).c_str(), "");

    no_limits_cl_2_legend_phase->Draw();

    // print_res
    std::ofstream no_limits_Cl_2_file("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_Cl_2_" + std::to_string(GetVoltage()) + ".txt");
    no_limits_Cl_2_file << std::setprecision(10);
    no_limits_Cl_2_file << "Ampl_woofer:\n";
    no_limits_Cl_2_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    no_limits_Cl_2_file << "\tChi^2: " << ampl_func_w_2->GetChisquare() << '\n';
    no_limits_Cl_2_file << "\tNdf: " << ampl_func_w_2->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tCHI RIDOTTO: " << ampl_func_w_2->GetChisquare() / ampl_func_w_2->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w_2->GetNpar(); ++i)
        no_limits_Cl_2_file << "\t\t[" << i << "] - " << ampl_func_w_2->GetParName(i) << " = " << ampl_func_w_2->GetParameter(i) << " +- " << ampl_func_w_2->GetParError(i) << '\n';
    no_limits_Cl_2_file << "Ampl_tweeter:\n";
    no_limits_Cl_2_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    no_limits_Cl_2_file << "\tChi^2: " << ampl_func_t->GetChisquare() << '\n';
    no_limits_Cl_2_file << "\tNdf: " << ampl_func_t->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tCHI RIDOTTO: " << ampl_func_t->GetChisquare() / ampl_func_t->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_Cl_2_file << "\t\t[" << i << "] - " << ampl_func_t->GetParName(i) << " = " << ampl_func_t->GetParameter(i) << " +- " << ampl_func_t->GetParError(i) << '\n';
    no_limits_Cl_2_file << "Phase_woofer:\n";
    no_limits_Cl_2_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    no_limits_Cl_2_file << "\tChi^2: " << phase_func_w_2->GetChisquare() << '\n';
    no_limits_Cl_2_file << "\tNdf: " << phase_func_w_2->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tCHI RIDOTTO: " << phase_func_w_2->GetChisquare() / phase_func_w_2->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w_2->GetNpar(); ++i)
        no_limits_Cl_2_file << "\t\t[" << i << "] - " << phase_func_w_2->GetParName(i) << " = " << phase_func_w_2->GetParameter(i) << " +- " << phase_func_w_2->GetParError(i) << '\n';
    no_limits_Cl_2_file << "Phase_tweeter:\n";
    no_limits_Cl_2_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    no_limits_Cl_2_file << "\tChi^2: " << phase_func_t->GetChisquare() << '\n';
    no_limits_Cl_2_file << "\tNdf: " << phase_func_t->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tCHI RIDOTTO: " << phase_func_t->GetChisquare() / phase_func_t->GetNDF() << '\n';
    no_limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_Cl_2_file << "\t\t[" << i << "] - " << phase_func_t->GetParName(i) << " = " << phase_func_t->GetParameter(i) << " +- " << phase_func_t->GetParError(i) << '\n';
    no_limits_Cl_2_file.close();

    // Save Canvas
    no_limits_Cl_2_canvas->Update();
    no_limits_Cl_2_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_Cl_2_" + std::to_string(GetVoltage()) + ".png").c_str());

    // NO LIMITS - Cl 3 -----------------------------------------------------------------
    TGraphErrors *ampl_graph_no_limits_Cl3{new TGraphErrors(*ampl_graph[1])};
    TGraphErrors *phase_graph_no_limits_Cl3{new TGraphErrors(*phase_graph[1])};

    ampl_func_w_3->SetParLimits(3, 0, 1e-6);
    phase_func_w_3->SetParLimits(3, 0, 1e-6);
    ampl_woofer_fit_res = ampl_graph_no_limits_Cl3->Fit(ampl_func_w_3, "Q, M, E, R");
    std::cout << std::endl;

    phase_woofer_fit_res = phase_graph_no_limits_Cl3->Fit(phase_func_w_3, "Q, M, E, R");
    std::cout << std::endl;

    auto multi_ampl_no_lim_Cl_3{new TMultiGraph};
    auto multi_phase_no_lim_Cl_3{new TMultiGraph};

    if (!isCrossover())
        multi_ampl_no_lim_Cl_3->Add(ampl_graph[0]);
    multi_phase_no_lim_Cl_3->Add(phase_graph[0]);
    multi_ampl_no_lim_Cl_3->Add(ampl_graph_no_limits_Cl3);
    multi_phase_no_lim_Cl_3->Add(phase_graph_no_limits_Cl3);
    multi_ampl_no_lim_Cl_3->Add(ampl_graph[2]);
    multi_phase_no_lim_Cl_3->Add(phase_graph[2]);

    // ampl_woofer_fit_res = ampl_graph[1]->Fit(ampl_func_w_3, "Q, M, E, R");
    // std::cout << std::endl;

    // phase_woofer_fit_res = phase_graph[1]->Fit(phase_func_w_3, "Q, M, E, R");
    // std::cout << std::endl;

    // auto multi_ampl_no_lim_Cl_3{new TMultiGraph};
    // auto multi_phase_no_lim_Cl_3{new TMultiGraph};

    // for (int i = 0; i != 3; ++i)
    // {
    //     multi_ampl_no_lim_Cl_3->Add(ampl_graph[i]);
    //     multi_phase_no_lim_Cl_3->Add(phase_graph[i]);
    // }

    multi_ampl_no_lim_Cl_3->SetTitle("Amplitude - Frequency (No Limits, Cl)");
    multi_phase_no_lim_Cl_3->SetTitle("Phase - Frequency (No Limits, Cl)");

    // multi_ampl_no_lim_Cl_3->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_no_lim_Cl_3->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_no_lim_Cl_3->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_no_lim_Cl_3->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_no_lim_Cl_3->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_no_lim_Cl_3->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_no_lim_Cl_3->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *no_limits_Cl_3_canvas{new TCanvas{"no_limits_Cl_3_canvas", "Amplitude and Phase No_limits_Cl_3", 0, 0, 1300, 700}};
    no_limits_Cl_3_canvas->Divide(2, 1);

    no_limits_Cl_3_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);
    multi_ampl_no_lim_Cl_3->Draw("ape");

    TPad *no_limits_Cl_3_pad_residual_ampl = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_Cl_3_pad_residual_ampl->SetTopMargin(0.8);
    no_limits_Cl_3_pad_residual_ampl->Draw();
    no_limits_Cl_3_pad_residual_ampl->SetFillStyle(0);
    no_limits_Cl_3_pad_residual_ampl->cd();
    no_limits_Cl_3_pad_residual_ampl->SetGridy();

    TGraphErrors *no_limits_Cl_3_pad_ampl_w = new TGraphErrors(ampl_graph[1]->GetN());
    for (Int_t i = 0; i != ampl_graph[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph[1]->GetPointX(i);
        Double_t y = ampl_graph[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w_3->Eval(x);
        no_limits_Cl_3_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_Cl_3_pad_ampl_t = new TGraphErrors(ampl_graph[2]->GetN());
    for (Int_t i = 0; i != ampl_graph[2]->GetN(); ++i)
    {
        Double_t x = ampl_graph[2]->GetPointX(i);
        Double_t y = ampl_graph[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t->Eval(x);
        no_limits_Cl_3_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_Cl_3_pad_ampl_w->SetLineColor(kRed + 2);
    no_limits_Cl_3_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_Cl_3_pad_multi_ampl = new TMultiGraph();
    no_limits_Cl_3_pad_multi_ampl->Add(no_limits_Cl_3_pad_ampl_w);
    no_limits_Cl_3_pad_multi_ampl->Add(no_limits_Cl_3_pad_ampl_t);
    no_limits_Cl_3_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    no_limits_Cl_3_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    no_limits_Cl_3_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_Cl_3_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    no_limits_Cl_3_pad_multi_ampl->Draw("apl");
    no_limits_Cl_3_canvas->cd(1);

    TLegend *no_limits_cl_3_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_cl_3_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_cl_3_legend_ampl->SetNColumns(3);
    // RIGA 1
    no_limits_cl_3_legend_ampl->AddEntry(ampl_graph[0], "Amplitude V_S", "ep");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_cl_3_legend_ampl->AddEntry(ampl_graph_no_limits_Cl3, "Amplitude V_Woofer", "ep");
    no_limits_cl_3_legend_ampl->AddEntry(ampl_func_w_3, "Amplitude V_Woofer Fit", "l");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w_3->GetChisquare() / ampl_func_w_3->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != ampl_func_w_3->GetNpar(); ++i)
        no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w_3->GetParName(i)) + " = " + NumErrScien(ampl_func_w_3->GetParameter(i), ampl_func_w_3->GetParError(i), ampl_func_w_3->GetParName(i))).c_str(), "");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 5
    no_limits_cl_3_legend_ampl->AddEntry(ampl_graph[2], "Amplitude V_Tweeter", "ep");
    no_limits_cl_3_legend_ampl->AddEntry(ampl_func_t, "Amplitude V_Tweeter Fit", "l");
    no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t->GetChisquare() / ampl_func_t->GetNDF())).c_str()), "");
    // RIGA 6
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_cl_3_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t->GetParName(i)) + " = " + NumErrScien(ampl_func_t->GetParameter(i), ampl_func_t->GetParError(i), ampl_func_w_3->GetParName(i))).c_str(), "");

    no_limits_cl_3_legend_ampl->Draw();

    no_limits_Cl_3_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_phase_no_lim_Cl_3->Draw("ape");

    TPad *no_limits_Cl_3_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    no_limits_Cl_3_pad_residual_phase->SetTopMargin(0.8);
    no_limits_Cl_3_pad_residual_phase->Draw();
    no_limits_Cl_3_pad_residual_phase->SetFillStyle(0);
    no_limits_Cl_3_pad_residual_phase->cd();
    no_limits_Cl_3_pad_residual_phase->SetGridy();

    TGraphErrors *no_limits_Cl_3_pad_phase_w = new TGraphErrors(phase_graph[1]->GetN());
    for (Int_t i = 0; i != phase_graph[1]->GetN(); ++i)
    {
        Double_t x = phase_graph[1]->GetPointX(i);
        Double_t y = phase_graph[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w_3->Eval(x);
        no_limits_Cl_3_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *no_limits_Cl_3_pad_phase_t = new TGraphErrors(phase_graph[2]->GetN());
    for (Int_t i = 0; i != phase_graph[2]->GetN(); ++i)
    {
        Double_t x = phase_graph[2]->GetPointX(i);
        Double_t y = phase_graph[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t->Eval(x);
        no_limits_Cl_3_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    no_limits_Cl_3_pad_phase_w->SetLineColor(kRed + 2);
    no_limits_Cl_3_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *no_limits_Cl_3_pad_multi_phase = new TMultiGraph();
    no_limits_Cl_3_pad_multi_phase->Add(no_limits_Cl_3_pad_phase_w);
    no_limits_Cl_3_pad_multi_phase->Add(no_limits_Cl_3_pad_phase_t);
    no_limits_Cl_3_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    no_limits_Cl_3_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    no_limits_Cl_3_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    no_limits_Cl_3_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    no_limits_Cl_3_pad_multi_phase->Draw("apl");

    no_limits_Cl_3_canvas->cd(2);

    TLegend *no_limits_cl_3_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // no_limits_cl_3_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    no_limits_cl_3_legend_phase->SetNColumns(3);
    // RIGA 1
    no_limits_cl_3_legend_phase->AddEntry(phase_graph[0], "Phase V_S", "ep");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    no_limits_cl_3_legend_phase->AddEntry(phase_graph_no_limits_Cl3, "Phase V_Woofer", "ep");
    no_limits_cl_3_legend_phase->AddEntry(phase_func_w_3, "Phase V_Woofer Fit", "l");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w_3->GetChisquare() / phase_func_w_3->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != phase_func_w_3->GetNpar(); ++i)
        no_limits_cl_3_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w_3->GetParName(i)) + " = " + NumErrScien(phase_func_w_3->GetParameter(i), phase_func_w_3->GetParError(i), phase_func_w_3->GetParName(i))).c_str(), "");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 4
    no_limits_cl_3_legend_phase->AddEntry(phase_graph[2], "Phase V_Tweeter", "ep");
    no_limits_cl_3_legend_phase->AddEntry(phase_func_t, "Phase V_Tweeter Fit", "l");
    no_limits_cl_3_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t->GetChisquare() / phase_func_t->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_cl_3_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t->GetParName(i)) + " = " + NumErrScien(phase_func_t->GetParameter(i), phase_func_t->GetParError(i), phase_func_t->GetParName(i))).c_str(), "");

    no_limits_cl_3_legend_phase->Draw();

    // print_res
    std::ofstream no_limits_Cl_3_file("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_Cl_3_" + std::to_string(GetVoltage()) + ".txt");
    no_limits_Cl_3_file << std::setprecision(10);
    no_limits_Cl_3_file << "Ampl_woofer:\n";
    no_limits_Cl_3_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    no_limits_Cl_3_file << "\tChi^2: " << ampl_func_w_3->GetChisquare() << '\n';
    no_limits_Cl_3_file << "\tNdf: " << ampl_func_w_3->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tCHI RIDOTTO: " << ampl_func_w_3->GetChisquare() / ampl_func_w_3->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w_3->GetNpar(); ++i)
        no_limits_Cl_3_file << "\t\t[" << i << "] - " << ampl_func_w_3->GetParName(i) << " = " << ampl_func_w_3->GetParameter(i) << " +- " << ampl_func_w_3->GetParError(i) << '\n';
    no_limits_Cl_3_file << "Ampl_tweeter:\n";
    no_limits_Cl_3_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    no_limits_Cl_3_file << "\tChi^2: " << ampl_func_t->GetChisquare() << '\n';
    no_limits_Cl_3_file << "\tNdf: " << ampl_func_t->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tCHI RIDOTTO: " << ampl_func_t->GetChisquare() / ampl_func_t->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t->GetNpar(); ++i)
        no_limits_Cl_3_file << "\t\t[" << i << "] - " << ampl_func_t->GetParName(i) << " = " << ampl_func_t->GetParameter(i) << " +- " << ampl_func_t->GetParError(i) << '\n';
    no_limits_Cl_3_file << "Phase_woofer:\n";
    no_limits_Cl_3_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    no_limits_Cl_3_file << "\tChi^2: " << phase_func_w_3->GetChisquare() << '\n';
    no_limits_Cl_3_file << "\tNdf: " << phase_func_w_3->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tCHI RIDOTTO: " << phase_func_w_3->GetChisquare() / phase_func_w_3->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w_3->GetNpar(); ++i)
        no_limits_Cl_3_file << "\t\t[" << i << "] - " << phase_func_w_3->GetParName(i) << " = " << phase_func_w_3->GetParameter(i) << " +- " << phase_func_w_3->GetParError(i) << '\n';
    no_limits_Cl_3_file << "Phase_tweeter:\n";
    no_limits_Cl_3_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    no_limits_Cl_3_file << "\tChi^2: " << phase_func_t->GetChisquare() << '\n';
    no_limits_Cl_3_file << "\tNdf: " << phase_func_t->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tCHI RIDOTTO: " << phase_func_t->GetChisquare() / phase_func_t->GetNDF() << '\n';
    no_limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t->GetNpar(); ++i)
        no_limits_Cl_3_file << "\t\t[" << i << "] - " << phase_func_t->GetParName(i) << " = " << phase_func_t->GetParameter(i) << " +- " << phase_func_t->GetParError(i) << '\n';
    no_limits_Cl_3_file.close();

    // Save Canvas
    no_limits_Cl_3_canvas->Update();
    no_limits_Cl_3_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/no_limits_Cl_3_" + std::to_string(GetVoltage()) + ".png").c_str());

    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // FITTING WITH LIMITS-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    TGraphErrors *ampl_graph_limits[3]{new TGraphErrors(*ampl_graph[0]),
                                       new TGraphErrors(*ampl_graph[1]),
                                       new TGraphErrors(*ampl_graph[2])};

    TGraphErrors *phase_graph_limits[3]{new TGraphErrors(*phase_graph[0]),
                                        new TGraphErrors(*phase_graph[1]),
                                        new TGraphErrors(*phase_graph[2])};

    TF1 *ampl_func_w_limits = new TF1(*ampl_func_w);
    TF1 *ampl_func_t_limits = new TF1(*ampl_func_t);

    // TF1 *ampl_func_w_limits = new TF1("amplw", "[0]*x+[1]");
    // TF1 *ampl_func_t_limits = new TF1("amplt", "[0]*x+[1]");
    TF1 *phase_func_w_limits = new TF1(*phase_func_w);
    TF1 *phase_func_t_limits = new TF1(*phase_func_t);

    TF1 *ampl_func_w_2_limits = new TF1(*ampl_func_w_2);
    TF1 *ampl_func_w_3_limits = new TF1(*ampl_func_w_3);
    TF1 *phase_func_w_2_limits = new TF1(*phase_func_w_2);
    TF1 *phase_func_w_3_limits = new TF1(*phase_func_w_3);

    ampl_func_w_limits->SetParameter(0, Rw);
    ampl_func_w_limits->SetParameter(1, Rl);
    ampl_func_w_limits->SetParameter(2, L);

    ampl_func_t_limits->SetParameter(0, Rl);
    ampl_func_t_limits->SetParameter(1, Rl1Rl2);
    ampl_func_t_limits->SetParameter(2, C1C2);

    phase_func_w_limits->SetParameter(0, Rw);
    phase_func_w_limits->SetParameter(1, Rl);
    phase_func_w_limits->SetParameter(2, L);

    phase_func_t_limits->SetParameter(0, Rl);
    phase_func_t_limits->SetParameter(1, Rl1Rl2);
    phase_func_t_limits->SetParameter(2, C1C2);

    // with cl------------------------------------------------
    ampl_func_w_2_limits->SetParameter(0, Rw);
    ampl_func_w_2_limits->SetParameter(1, Rl);
    ampl_func_w_2_limits->SetParameter(2, L);
    ampl_func_w_2_limits->SetParameter(3, 1E-8);

    ampl_func_w_3_limits->SetParameter(0, Rw);
    ampl_func_w_3_limits->SetParameter(1, Rl);
    ampl_func_w_3_limits->SetParameter(2, L);
    ampl_func_w_3_limits->SetParameter(3, 1E-8);

    phase_func_w_2_limits->SetParameter(0, Rw);
    phase_func_w_2_limits->SetParameter(1, Rl);
    phase_func_w_2_limits->SetParameter(2, L);
    phase_func_w_2_limits->SetParameter(3, 1E-8);

    phase_func_w_3_limits->SetParameter(0, Rw);
    phase_func_w_3_limits->SetParameter(1, Rl);
    phase_func_w_3_limits->SetParameter(2, L);
    phase_func_w_3_limits->SetParameter(3, 1E-8);
    // with cl------------------------------------------------

    // // PAR LIMITS-----------------------------------------------------------------
    double N_SIGMA = 5;
    // ampl_func_w_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    // ampl_func_w_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    ampl_func_w_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);

    // ampl_func_t_limits->SetParLimits(0, Rt - N_SIGMA * Rt_err, Rt + N_SIGMA * Rt_err);
    // ampl_func_t_limits->SetParLimits(1, Rl1Rl2 - N_SIGMA * Rl1Rl2_err, Rl1Rl2 + N_SIGMA * Rl1Rl2_err);
    ampl_func_t_limits->SetParLimits(2, C1C2 - N_SIGMA * C1C2_err, C1C2 + N_SIGMA * C1C2_err);

    phase_func_w_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    phase_func_w_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    phase_func_w_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);

    phase_func_t_limits->SetParLimits(0, Rt - N_SIGMA * Rt_err, Rt + N_SIGMA * Rt_err);
    phase_func_t_limits->SetParLimits(1, Rl1Rl2 - N_SIGMA * Rl1Rl2_err, Rl1Rl2 + N_SIGMA * Rl1Rl2_err);
    phase_func_t_limits->SetParLimits(2, C1C2 - N_SIGMA * C1C2_err, C1C2 + N_SIGMA * C1C2_err);

    // cl---------------------
    ampl_func_w_2_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    ampl_func_w_2_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    ampl_func_w_2_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);
    ampl_func_w_2_limits->SetParLimits(3, 0., 1e-6);

    phase_func_w_2_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    phase_func_w_2_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    phase_func_w_2_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);
    phase_func_w_2_limits->SetParLimits(3, 0., 1e-6);

    ampl_func_w_3_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    ampl_func_w_3_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    ampl_func_w_3_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);
    ampl_func_w_3_limits->SetParLimits(3, 0., 1e-6);

    phase_func_w_3_limits->SetParLimits(0, Rw - N_SIGMA * Rw_err, Rw + N_SIGMA * Rw_err);
    phase_func_w_3_limits->SetParLimits(1, Rl - N_SIGMA * Rl_err, Rl + N_SIGMA * Rl_err);
    phase_func_w_3_limits->SetParLimits(2, L - N_SIGMA * L_err, L + N_SIGMA * L_err);
    phase_func_w_3_limits->SetParLimits(3, 0., 1e-6);
    // cl---------------------

    // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    // ampl_func_w_limits->SetRange(ampl_graph_limits[1]->GetPointX(0), ampl_graph_limits[1]->GetPointX(ampl_graph_limits[1]->GetN()) - 1);
    // ampl_func_w_limits->SetMinimum(ampl_graph_limits[1]->GetMinimum());
    // ampl_func_w_limits->SetMaximum(ampl_graph_limits[1]->GetMaximum());
    // ampl_func_t_limits->SetRange(ampl_graph_limits[2]->GetPointX(0), ampl_graph_limits[2]->GetPointX(ampl_graph_limits[2]->GetN()) - 1);
    // FITTING WITH LIMITS (No Cl)-------------------------

    ampl_woofer_fit_res = ampl_graph_limits[1]->Fit(ampl_func_w_limits, "Q, M, E, R");
    std::cout << std::endl;
    ampl_tweeter_fit_res = ampl_graph_limits[2]->Fit(ampl_func_t_limits, "Q, M, E, R");
    std::cout << std::endl;

    phase_woofer_fit_res = phase_graph_limits[1]->Fit(phase_func_w_limits, "Q, M, E, R");
    std::cout << std::endl;
    phase_tweeter_fit_res = phase_graph_limits[2]->Fit(phase_func_t_limits, "Q, M, E, R");
    std::cout << std::endl;

    auto multi_ampl_lim_no_Cl{new TMultiGraph};
    auto multi_phase_lim_no_Cl{new TMultiGraph};
    for (int i = (isCrossover() ? 1 : 0); i != 3; ++i)
    {
        multi_ampl_lim_no_Cl->Add(ampl_graph_limits[i]);
    }
    for (int i = 0; i != 3; ++i)
    {
        multi_phase_lim_no_Cl->Add(phase_graph_limits[i]);
    }
    multi_ampl_lim_no_Cl->SetTitle("Amplitude - Frequency (Limited, No Cl)");
    multi_phase_lim_no_Cl->SetTitle("Phase - Frequency (Limited, No Cl)");

    // multi_ampl_lim_no_Cl->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_lim_no_Cl->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_lim_no_Cl->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_lim_no_Cl->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_lim_no_Cl->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_lim_no_Cl->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_lim_no_Cl->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *limits_no_cl_canvas{new TCanvas{"limits_no_cl_canvas", "Amplitude and Phase Limits_no_Cl", 0, 0, 1300, 700}};
    limits_no_cl_canvas->Divide(2, 1);
    limits_no_cl_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_ampl_lim_no_Cl->Draw("ape");

    TPad *limits_no_cl_pad_residual = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_no_cl_pad_residual->SetTopMargin(0.8);
    limits_no_cl_pad_residual->Draw();
    limits_no_cl_pad_residual->SetFillStyle(0);
    limits_no_cl_pad_residual->cd();
    limits_no_cl_pad_residual->SetGridy();

    TGraphErrors *limits_no_cl_pad_ampl_w = new TGraphErrors(ampl_graph_limits[1]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[1]->GetPointX(i);
        Double_t y = ampl_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w_limits->Eval(x);
        limits_no_cl_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_no_cl_pad_ampl_t = new TGraphErrors(ampl_graph_limits[2]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[2]->GetPointX(i);
        Double_t y = ampl_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t_limits->Eval(x);
        limits_no_cl_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    limits_no_cl_pad_ampl_w->SetLineColor(kRed + 2);
    limits_no_cl_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_no_cl_pad_multi_ampl = new TMultiGraph();
    limits_no_cl_pad_multi_ampl->Add(limits_no_cl_pad_ampl_w);
    limits_no_cl_pad_multi_ampl->Add(limits_no_cl_pad_ampl_t);
    limits_no_cl_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    limits_no_cl_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    limits_no_cl_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_no_cl_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    limits_no_cl_pad_multi_ampl->Draw("apl");

    limits_no_cl_canvas->cd(1);
    TLegend *limits_no_cl_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_no_cl_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_no_cl_legend_ampl->SetNColumns(3);
    // RIGA 1
    limits_no_cl_legend_ampl->AddEntry(ampl_graph_limits[0], "Amplitude V_S", "ep");
    limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    limits_no_cl_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_no_cl_legend_ampl->AddEntry(ampl_graph_limits[1], "Amplitude V_Woofer", "ep");
    limits_no_cl_legend_ampl->AddEntry(ampl_func_w_limits, "Amplitude V_Woofer Fit", "l");
    limits_no_cl_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w_limits->GetChisquare() / ampl_func_w_limits->GetNDF())).c_str()), "");
    // RIGA 3
    for (int i = 0; i != ampl_func_w_limits->GetNpar(); ++i)
        limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_w_limits->GetParameter(i), ampl_func_w_limits->GetParError(i), ampl_func_w_limits->GetParName(i))).c_str(), "");
    // RIGA 4
    limits_no_cl_legend_ampl->AddEntry(ampl_graph_limits[2], "Amplitude V_Tweeter", "ep");
    limits_no_cl_legend_ampl->AddEntry(ampl_func_t_limits, "Amplitude V_Tweeter Fit", "l");
    limits_no_cl_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_no_cl_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_t_limits->GetParameter(i), ampl_func_t_limits->GetParError(i), ampl_func_t_limits->GetParName(i))).c_str(), "");

    limits_no_cl_legend_ampl->Draw();

    limits_no_cl_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_phase_lim_no_Cl->Draw("ape");

    TPad *limits_no_cl_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_no_cl_pad_residual_phase->SetTopMargin(0.8);
    limits_no_cl_pad_residual_phase->Draw();
    limits_no_cl_pad_residual_phase->SetFillStyle(0);
    limits_no_cl_pad_residual_phase->cd();
    limits_no_cl_pad_residual_phase->SetGridy();

    TGraphErrors *limits_no_cl_pad_phase_w = new TGraphErrors(phase_graph_limits[1]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[1]->GetPointX(i);
        Double_t y = phase_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w_limits->Eval(x);
        limits_no_cl_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_no_cl_pad_phase_t = new TGraphErrors(phase_graph_limits[2]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[2]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[2]->GetPointX(i);
        Double_t y = phase_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t_limits->Eval(x);
        limits_no_cl_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    limits_no_cl_pad_phase_w->SetLineColor(kRed + 2);
    limits_no_cl_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_no_cl_pad_multi_phase = new TMultiGraph();
    limits_no_cl_pad_multi_phase->Add(limits_no_cl_pad_phase_w);
    limits_no_cl_pad_multi_phase->Add(limits_no_cl_pad_phase_t);
    limits_no_cl_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    limits_no_cl_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    limits_no_cl_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_no_cl_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    limits_no_cl_pad_multi_phase->Draw("apl");

    limits_no_cl_canvas->cd(2);

    TLegend *limits_no_cl_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_no_cl_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_no_cl_legend_phase->SetNColumns(3);
    // RIGA 1
    limits_no_cl_legend_phase->AddEntry(phase_graph_limits[0], "Phase V_S", "ep");
    limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    limits_no_cl_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_no_cl_legend_phase->AddEntry(phase_graph_limits[1], "Phase V_Woofer", "ep");
    limits_no_cl_legend_phase->AddEntry(phase_func_w_limits, "Phase V_Woofer Fit", "l");
    limits_no_cl_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w_limits->GetChisquare() / phase_func_w_limits->GetNDF())).c_str()), "");
    // RIGA 3
    for (int i = 0; i != phase_func_w_limits->GetNpar(); ++i)
        limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w_limits->GetParName(i)) + " = " + NumErrScien(phase_func_w_limits->GetParameter(i), phase_func_w_limits->GetParError(i), phase_func_w_limits->GetParName(i))).c_str(), "");
    // RIGA 4
    limits_no_cl_legend_phase->AddEntry(phase_graph_limits[2], "Phase V_Tweeter", "ep");
    limits_no_cl_legend_phase->AddEntry(phase_func_t_limits, "Phase V_Tweeter Fit", "l");
    limits_no_cl_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_no_cl_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t_limits->GetParName(i)) + " = " + NumErrScien(phase_func_t_limits->GetParameter(i), phase_func_t_limits->GetParError(i), phase_func_t_limits->GetParName(i))).c_str(), "");

    limits_no_cl_legend_phase->Draw();

    // print_res
    std::ofstream limits_no_Cl_file("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_no_Cl_" + std::to_string(GetVoltage()) + ".txt");
    limits_no_Cl_file << std::setprecision(10);
    limits_no_Cl_file << "Ampl_woofer:\n";
    limits_no_Cl_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    limits_no_Cl_file << "\tChi^2: " << ampl_func_w_limits->GetChisquare() << '\n';
    limits_no_Cl_file << "\tNdf: " << ampl_func_w_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tCHI RIDOTTO: " << ampl_func_w_limits->GetChisquare() / ampl_func_w_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w_limits->GetNpar(); ++i)
        limits_no_Cl_file << "\t\t[" << i << "] - " << ampl_func_w_limits->GetParName(i) << " = " << ampl_func_w_limits->GetParameter(i) << " +- " << ampl_func_w_limits->GetParError(i) << '\n';
    limits_no_Cl_file << "Ampl_tweeter:\n";
    limits_no_Cl_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    limits_no_Cl_file << "\tChi^2: " << ampl_func_t_limits->GetChisquare() << '\n';
    limits_no_Cl_file << "\tNdf: " << ampl_func_t_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tCHI RIDOTTO: " << ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_no_Cl_file << "\t\t[" << i << "] - " << ampl_func_t_limits->GetParName(i) << " = " << ampl_func_t_limits->GetParameter(i) << " +- " << ampl_func_t_limits->GetParError(i) << '\n';
    limits_no_Cl_file << "Phase_woofer:\n";
    limits_no_Cl_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    limits_no_Cl_file << "\tChi^2: " << phase_func_w_limits->GetChisquare() << '\n';
    limits_no_Cl_file << "\tNdf: " << phase_func_w_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tCHI RIDOTTO: " << phase_func_w_limits->GetChisquare() / phase_func_w_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w_limits->GetNpar(); ++i)
        limits_no_Cl_file << "\t\t[" << i << "] - " << phase_func_w_limits->GetParName(i) << " = " << phase_func_w_limits->GetParameter(i) << " +- " << phase_func_w_limits->GetParError(i) << '\n';
    limits_no_Cl_file << "Phase_tweeter:\n";
    limits_no_Cl_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    limits_no_Cl_file << "\tChi^2: " << phase_func_t_limits->GetChisquare() << '\n';
    limits_no_Cl_file << "\tNdf: " << phase_func_t_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tCHI RIDOTTO: " << phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF() << '\n';
    limits_no_Cl_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_no_Cl_file << "\t\t[" << i << "] - " << phase_func_t_limits->GetParName(i) << " = " << phase_func_t_limits->GetParameter(i) << " +- " << phase_func_t_limits->GetParError(i) << '\n';
    limits_no_Cl_file.close();

    // Save Canvas
    limits_no_cl_canvas->Update();
    limits_no_cl_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_no_cl_" + std::to_string(GetVoltage()) + ".png").c_str());

    // LIMITS - Cl 2 -----------------------------------------------------------------
    TGraphErrors *ampl_graph_limits_Cl2{new TGraphErrors(*ampl_graph_limits[1])};
    TGraphErrors *phase_graph_limits_Cl2{new TGraphErrors(*phase_graph_limits[1])};

    ampl_func_w_2_limits->SetParLimits(3, 0, 1e-6);
    phase_func_w_2_limits->SetParLimits(3, 0, 1e-6);
    ampl_woofer_fit_res = ampl_graph_limits_Cl2->Fit(ampl_func_w_2_limits, "Q, M, E, R");
    std::cout << std::endl;

    phase_woofer_fit_res = phase_graph_limits_Cl2->Fit(phase_func_w_2_limits, "Q, M, E, R");
    std::cout << std::endl;

    auto multi_ampl_lim_Cl_2{new TMultiGraph};
    auto multi_phase_lim_Cl_2{new TMultiGraph};

    if (!isCrossover())
        multi_ampl_lim_Cl_2->Add(ampl_graph_limits[0]);
    multi_phase_lim_Cl_2->Add(phase_graph_limits[0]);
    multi_ampl_lim_Cl_2->Add(ampl_graph_limits_Cl2);
    multi_phase_lim_Cl_2->Add(phase_graph_limits_Cl2);
    multi_ampl_lim_Cl_2->Add(ampl_graph_limits[2]);
    multi_phase_lim_Cl_2->Add(phase_graph_limits[2]);

    // ampl_woofer_fit_res = ampl_graph[1]->Fit(ampl_func_w_2, "Q, M, E, R");
    // std::cout << std::endl;

    // phase_woofer_fit_res = phase_graph[1]->Fit(phase_func_w_2_limits, "Q, M, E, R");
    // std::cout << std::endl;

    // auto multi_ampl_lim_Cl_2{new TMultiGraph};
    // auto multi_phase_lim_Cl_2{new TMultiGraph};

    // for (int i = 0; i != 3; ++i)
    // {
    //     multi_ampl_lim_Cl_2->Add(ampl_graph[i]);
    //     multi_phase_lim_Cl_2->Add(phase_graph[i]);
    // }
    multi_ampl_lim_Cl_2->SetTitle("Amplitude - Frequency (Limited, Cl)");
    multi_phase_lim_Cl_2->SetTitle("Phase - Frequency (Limited, Cl)");

    // multi_ampl_lim_Cl_2->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_lim_Cl_2->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_lim_Cl_2->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_lim_Cl_2->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_lim_Cl_2->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_lim_Cl_2->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_lim_Cl_2->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *limits_Cl_2_canvas{new TCanvas{"limits_Cl_2_canvas", "Amplitude and Phase Limits_Cl_2", 0, 0, 1300, 700}};
    limits_Cl_2_canvas->Divide(2, 1);
    limits_Cl_2_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_ampl_lim_Cl_2->Draw("ape");

    TPad *limits_Cl_2_pad_residual_ampl = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_Cl_2_pad_residual_ampl->SetTopMargin(0.8);
    limits_Cl_2_pad_residual_ampl->Draw();
    limits_Cl_2_pad_residual_ampl->SetFillStyle(0);
    limits_Cl_2_pad_residual_ampl->cd();
    limits_Cl_2_pad_residual_ampl->SetGridy();

    TGraphErrors *limits_Cl_2_pad_ampl_w = new TGraphErrors(ampl_graph_limits[1]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[1]->GetPointX(i);
        Double_t y = ampl_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w_2_limits->Eval(x);
        limits_Cl_2_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_Cl_2_pad_ampl_t = new TGraphErrors(ampl_graph_limits[2]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[2]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[2]->GetPointX(i);
        Double_t y = ampl_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t_limits->Eval(x);
        limits_Cl_2_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    limits_Cl_2_pad_ampl_w->SetLineColor(kRed + 2);
    limits_Cl_2_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_Cl_2_pad_multi_ampl = new TMultiGraph();
    limits_Cl_2_pad_multi_ampl->Add(limits_Cl_2_pad_ampl_w);
    limits_Cl_2_pad_multi_ampl->Add(limits_Cl_2_pad_ampl_t);
    limits_Cl_2_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    limits_Cl_2_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    limits_Cl_2_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_Cl_2_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    limits_Cl_2_pad_multi_ampl->Draw("apl");
    limits_Cl_2_canvas->cd(1);

    TLegend *limits_cl_2_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_cl_2_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_cl_2_legend_ampl->SetNColumns(3);
    // RIGA 1
    limits_cl_2_legend_ampl->AddEntry(ampl_graph_limits[0], "Amplitude V_S", "ep");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_cl_2_legend_ampl->AddEntry(ampl_graph_limits_Cl2, "Amplitude V_Woofer", "ep");
    limits_cl_2_legend_ampl->AddEntry(ampl_func_w_2_limits, "Amplitude V_Woofer Fit", "l");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w_2_limits->GetChisquare() / ampl_func_w_2_limits->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != ampl_func_w_2_limits->GetNpar(); ++i)
        limits_cl_2_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w_2_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_w_2_limits->GetParameter(i), ampl_func_w_2_limits->GetParError(i), ampl_func_w_2_limits->GetParName(i))).c_str(), "");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 5
    limits_cl_2_legend_ampl->AddEntry(ampl_graph_limits[2], "Amplitude V_Tweeter", "ep");
    limits_cl_2_legend_ampl->AddEntry(ampl_func_t_limits, "Amplitude V_Tweeter Fit", "l");
    limits_cl_2_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 6
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_cl_2_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_t_limits->GetParameter(i), ampl_func_t_limits->GetParError(i), ampl_func_t_limits->GetParName(i))).c_str(), "");

    limits_cl_2_legend_ampl->Draw();

    limits_Cl_2_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);

    multi_phase_lim_Cl_2->Draw("ape");

    TPad *limits_Cl_2_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_Cl_2_pad_residual_phase->SetTopMargin(0.8);
    limits_Cl_2_pad_residual_phase->Draw();
    limits_Cl_2_pad_residual_phase->SetFillStyle(0);
    limits_Cl_2_pad_residual_phase->cd();
    limits_Cl_2_pad_residual_phase->SetGridy();

    TGraphErrors *limits_Cl_2_pad_phase_w = new TGraphErrors(phase_graph_limits[1]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[1]->GetPointX(i);
        Double_t y = phase_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w_2_limits->Eval(x);
        limits_Cl_2_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_Cl_2_pad_phase_t = new TGraphErrors(phase_graph_limits[2]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[2]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[2]->GetPointX(i);
        Double_t y = phase_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t_limits->Eval(x);
        limits_Cl_2_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    limits_Cl_2_pad_phase_w->SetLineColor(kRed + 2);
    limits_Cl_2_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_Cl_2_pad_multi_phase = new TMultiGraph();
    limits_Cl_2_pad_multi_phase->Add(limits_Cl_2_pad_phase_w);
    limits_Cl_2_pad_multi_phase->Add(limits_Cl_2_pad_phase_t);
    limits_Cl_2_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    limits_Cl_2_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    limits_Cl_2_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_Cl_2_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    limits_Cl_2_pad_multi_phase->Draw("apl");

    limits_Cl_2_canvas->cd(2);

    TLegend *limits_cl_2_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_cl_2_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_cl_2_legend_phase->SetNColumns(3);
    // RIGA 1
    limits_cl_2_legend_phase->AddEntry(phase_graph_limits[0], "Phase V_S", "ep");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_cl_2_legend_phase->AddEntry(phase_graph_limits_Cl2, "Phase V_Woofer", "ep");
    limits_cl_2_legend_phase->AddEntry(phase_func_w_2_limits, "Phase V_Woofer Fit", "l");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w_2_limits->GetChisquare() / phase_func_w_2_limits->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != phase_func_w_2_limits->GetNpar(); ++i)
        limits_cl_2_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w_2_limits->GetParName(i)) + " = " + NumErrScien(phase_func_w_2_limits->GetParameter(i), phase_func_w_2_limits->GetParError(i), phase_func_w_2_limits->GetParName(i))).c_str(), "");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 4
    limits_cl_2_legend_phase->AddEntry(phase_graph_limits[2], "Phase V_Tweeter", "ep");
    limits_cl_2_legend_phase->AddEntry(phase_func_t_limits, "Phase V_Tweeter Fit", "l");
    limits_cl_2_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_cl_2_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t_limits->GetParName(i)) + " = " + NumErrScien(phase_func_t_limits->GetParameter(i), phase_func_t_limits->GetParError(i), phase_func_t_limits->GetParName(i))).c_str(), "");

    limits_cl_2_legend_phase->Draw();

    // print_res
    std::ofstream limits_Cl_2_file("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_Cl_2_" + std::to_string(GetVoltage()) + ".txt");
    limits_Cl_2_file << std::setprecision(10);
    limits_Cl_2_file << "Ampl_woofer:\n";
    limits_Cl_2_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    limits_Cl_2_file << "\tChi^2: " << ampl_func_w_2_limits->GetChisquare() << '\n';
    limits_Cl_2_file << "\tNdf: " << ampl_func_w_2_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tCHI RIDOTTO: " << ampl_func_w_2_limits->GetChisquare() / ampl_func_w_2_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w_2_limits->GetNpar(); ++i)
        limits_Cl_2_file << "\t\t[" << i << "] - " << ampl_func_w_2_limits->GetParName(i) << " = " << ampl_func_w_2_limits->GetParameter(i) << " +- " << ampl_func_w_2_limits->GetParError(i) << '\n';
    limits_Cl_2_file << "Ampl_tweeter:\n";
    limits_Cl_2_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    limits_Cl_2_file << "\tChi^2: " << ampl_func_t_limits->GetChisquare() << '\n';
    limits_Cl_2_file << "\tNdf: " << ampl_func_t_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tCHI RIDOTTO: " << ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_Cl_2_file << "\t\t[" << i << "] - " << ampl_func_t_limits->GetParName(i) << " = " << ampl_func_t_limits->GetParameter(i) << " +- " << ampl_func_t_limits->GetParError(i) << '\n';
    limits_Cl_2_file << "Phase_woofer:\n";
    limits_Cl_2_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    limits_Cl_2_file << "\tChi^2: " << phase_func_w_2_limits->GetChisquare() << '\n';
    limits_Cl_2_file << "\tNdf: " << phase_func_w_2_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tCHI RIDOTTO: " << phase_func_w_2_limits->GetChisquare() / phase_func_w_2_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w_2_limits->GetNpar(); ++i)
        limits_Cl_2_file << "\t\t[" << i << "] - " << phase_func_w_2_limits->GetParName(i) << " = " << phase_func_w_2_limits->GetParameter(i) << " +- " << phase_func_w_2_limits->GetParError(i) << '\n';
    limits_Cl_2_file << "Phase_tweeter:\n";
    limits_Cl_2_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    limits_Cl_2_file << "\tChi^2: " << phase_func_t_limits->GetChisquare() << '\n';
    limits_Cl_2_file << "\tNdf: " << phase_func_t_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tCHI RIDOTTO: " << phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF() << '\n';
    limits_Cl_2_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_Cl_2_file << "\t\t[" << i << "] - " << phase_func_t_limits->GetParName(i) << " = " << phase_func_t_limits->GetParameter(i) << " +- " << phase_func_t_limits->GetParError(i) << '\n';
    limits_Cl_2_file.close();

    // Save Canvas
    limits_Cl_2_canvas->Update();
    limits_Cl_2_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_Cl_2_" + std::to_string(GetVoltage()) + ".png").c_str());

    // NO LIMITS - Cl 3 -----------------------------------------------------------------
    TGraphErrors *ampl_graph_limits_Cl3{new TGraphErrors(*ampl_graph_limits[1])};
    TGraphErrors *phase_graph_limits_Cl3{new TGraphErrors(*phase_graph_limits[1])};

    ampl_func_w_3_limits->SetParLimits(3, 0, 1e-6);
    phase_func_w_3_limits->SetParLimits(3, 0, 1e-6);
    ampl_woofer_fit_res = ampl_graph_limits_Cl3->Fit(ampl_func_w_3_limits, "Q, M, E, R");
    std::cout << std::endl;

    phase_woofer_fit_res = phase_graph_limits_Cl3->Fit(phase_func_w_3_limits, "Q, M, E, R");
    std::cout << std::endl;

    auto multi_ampl_lim_Cl_3{new TMultiGraph};
    auto multi_phase_lim_Cl_3{new TMultiGraph};

    if (!isCrossover())
        multi_ampl_lim_Cl_3->Add(ampl_graph_limits[0]);
    multi_phase_lim_Cl_3->Add(phase_graph_limits[0]);
    multi_ampl_lim_Cl_3->Add(ampl_graph_limits_Cl3);
    multi_phase_lim_Cl_3->Add(phase_graph_limits_Cl3);
    multi_ampl_lim_Cl_3->Add(ampl_graph_limits[2]);
    multi_phase_lim_Cl_3->Add(phase_graph_limits[2]);

    // ampl_woofer_fit_res = ampl_graph[1]->Fit(ampl_func_w_3_limits, "Q, M, E, R");
    // std::cout << std::endl;

    // phase_woofer_fit_res = phase_graph[1]->Fit(phase_func_w_3_limits, "Q, M, E, R");
    // std::cout << std::endl;

    // auto multi_ampl_lim_Cl_3{new TMultiGraph};
    // auto multi_phase_lim_Cl_3{new TMultiGraph};

    // for (int i = 0; i != 3; ++i)
    // {
    //     multi_ampl_lim_Cl_3->Add(ampl_graph[i]);
    //     multi_phase_lim_Cl_3->Add(phase_graph[i]);
    // }

    multi_ampl_lim_Cl_3->SetTitle("Amplitude - Frequency (Limited, Cl)");
    multi_phase_lim_Cl_3->SetTitle("Phase - Frequency (Limited, Cl)");

    // multi_ampl_lim_Cl_3->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_ampl_lim_Cl_3->GetXaxis()->SetLabelColor(1, 0.f);
    multi_ampl_lim_Cl_3->GetYaxis()->SetTitle("Amplitude (V)");
    // multi_phase_lim_Cl_3->GetXaxis()->SetTitle("Frequency (Hz)");
    multi_phase_lim_Cl_3->GetXaxis()->SetLabelColor(1, 0.f);
    multi_phase_lim_Cl_3->GetYaxis()->SetTitle("Phase shift (rad)");
    multi_phase_lim_Cl_3->GetYaxis()->SetTitleOffset(1.5f);

    TCanvas *limits_Cl_3_canvas{new TCanvas{"limits_Cl_3_canvas", "Amplitude and Phase Limits_Cl_3", 0, 0, 1300, 700}};
    limits_Cl_3_canvas->Divide(2, 1);
    limits_Cl_3_canvas->cd(1);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);
    multi_ampl_lim_Cl_3->Draw("ape");

    TPad *limits_Cl_3_pad_residual_ampl = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_Cl_3_pad_residual_ampl->SetTopMargin(0.8);
    limits_Cl_3_pad_residual_ampl->Draw();
    limits_Cl_3_pad_residual_ampl->SetFillStyle(0);
    limits_Cl_3_pad_residual_ampl->cd();
    limits_Cl_3_pad_residual_ampl->SetGridy();

    TGraphErrors *limits_Cl_3_pad_ampl_w = new TGraphErrors(ampl_graph_limits[1]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[1]->GetPointX(i);
        Double_t y = ampl_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = ampl_func_w_3_limits->Eval(x);
        limits_Cl_3_pad_ampl_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_Cl_3_pad_ampl_t = new TGraphErrors(ampl_graph_limits[2]->GetN());
    for (Int_t i = 0; i != ampl_graph_limits[2]->GetN(); ++i)
    {
        Double_t x = ampl_graph_limits[2]->GetPointX(i);
        Double_t y = ampl_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = ampl_func_t_limits->Eval(x);
        limits_Cl_3_pad_ampl_t->SetPoint(i, x, y - y_f0);
    }

    limits_Cl_3_pad_ampl_w->SetLineColor(kRed + 2);
    limits_Cl_3_pad_ampl_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_Cl_3_pad_multi_ampl = new TMultiGraph();
    limits_Cl_3_pad_multi_ampl->Add(limits_Cl_3_pad_ampl_w);
    limits_Cl_3_pad_multi_ampl->Add(limits_Cl_3_pad_ampl_t);
    limits_Cl_3_pad_multi_ampl->GetXaxis()->SetLabelSize(0.04);
    limits_Cl_3_pad_multi_ampl->GetYaxis()->SetLabelSize(0.03);
    limits_Cl_3_pad_multi_ampl->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_Cl_3_pad_multi_ampl->GetYaxis()->SetNdivisions(-4);
    limits_Cl_3_pad_multi_ampl->Draw("apl");
    limits_Cl_3_canvas->cd(1);

    TLegend *limits_cl_3_legend_ampl{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_cl_3_legend_ampl->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_cl_3_legend_ampl->SetNColumns(3);
    // RIGA 1
    limits_cl_3_legend_ampl->AddEntry(ampl_graph_limits[0], "Amplitude V_S", "ep");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_cl_3_legend_ampl->AddEntry(ampl_graph_limits_Cl2, "Amplitude V_Woofer", "ep");
    limits_cl_3_legend_ampl->AddEntry(ampl_func_w_3_limits, "Amplitude V_Woofer Fit", "l");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_w_3_limits->GetChisquare() / ampl_func_w_3_limits->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != ampl_func_w_3_limits->GetNpar(); ++i)
        limits_cl_3_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_w_3_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_w_3_limits->GetParameter(i), ampl_func_w_3_limits->GetParError(i), ampl_func_w_3_limits->GetParName(i))).c_str(), "");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, "", "");
    // RIGA 5
    limits_cl_3_legend_ampl->AddEntry(ampl_graph_limits[2], "Amplitude V_Tweeter", "ep");
    limits_cl_3_legend_ampl->AddEntry(ampl_func_t_limits, "Amplitude V_Tweeter Fit", "l");
    limits_cl_3_legend_ampl->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 6
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_cl_3_legend_ampl->AddEntry((TObject *)0, (std::string(ampl_func_t_limits->GetParName(i)) + " = " + NumErrScien(ampl_func_t_limits->GetParameter(i), ampl_func_t_limits->GetParError(i), ampl_func_t_limits->GetParName(i))).c_str(), "");

    limits_cl_3_legend_ampl->Draw();

    limits_Cl_3_canvas->cd(2);

    // for residual and legend
    gPad->SetGridy();
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.2);
    multi_phase_lim_Cl_3->Draw("ape");

    TPad *limits_Cl_3_pad_residual_phase = new TPad("pad", "pad", 0., 0., 1., 1.);
    limits_Cl_3_pad_residual_phase->SetTopMargin(0.8);
    limits_Cl_3_pad_residual_phase->Draw();
    limits_Cl_3_pad_residual_phase->SetFillStyle(0);
    limits_Cl_3_pad_residual_phase->cd();
    limits_Cl_3_pad_residual_phase->SetGridy();

    TGraphErrors *limits_Cl_3_pad_phase_w = new TGraphErrors(phase_graph_limits[1]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[1]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[1]->GetPointX(i);
        Double_t y = phase_graph_limits[1]->GetPointY(i);
        Double_t y_f0 = phase_func_w_3_limits->Eval(x);
        limits_Cl_3_pad_phase_w->SetPoint(i, x, y - y_f0);
    }

    TGraphErrors *limits_Cl_3_pad_phase_t = new TGraphErrors(phase_graph_limits[2]->GetN());
    for (Int_t i = 0; i != phase_graph_limits[2]->GetN(); ++i)
    {
        Double_t x = phase_graph_limits[2]->GetPointX(i);
        Double_t y = phase_graph_limits[2]->GetPointY(i);
        Double_t y_f0 = phase_func_t_limits->Eval(x);
        limits_Cl_3_pad_phase_t->SetPoint(i, x, y - y_f0);
    }

    limits_Cl_3_pad_phase_w->SetLineColor(kRed + 2);
    limits_Cl_3_pad_phase_t->SetLineColor(kBlue + 2);
    TMultiGraph *limits_Cl_3_pad_multi_phase = new TMultiGraph();
    limits_Cl_3_pad_multi_phase->Add(limits_Cl_3_pad_phase_w);
    limits_Cl_3_pad_multi_phase->Add(limits_Cl_3_pad_phase_t);
    limits_Cl_3_pad_multi_phase->GetXaxis()->SetLabelSize(0.04);
    limits_Cl_3_pad_multi_phase->GetYaxis()->SetLabelSize(0.03);
    limits_Cl_3_pad_multi_phase->GetXaxis()->SetTitle("Frequency (Hz)");
    limits_Cl_3_pad_multi_phase->GetYaxis()->SetNdivisions(-4);
    limits_Cl_3_pad_multi_phase->Draw("apl");
    limits_Cl_3_canvas->cd(2);

    TLegend *limits_cl_3_legend_phase{new TLegend(0.01, 0.8, 0.99, 0.93)};
    // limits_cl_3_legend_phase->SetHeader("Amplitude - Frequency", "C"); // option "C" allows to center the header
    limits_cl_3_legend_phase->SetNColumns(3);
    // RIGA 1
    limits_cl_3_legend_phase->AddEntry(phase_graph_limits[0], "Phase V_S", "ep");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 2
    limits_cl_3_legend_phase->AddEntry(phase_graph_limits_Cl2, "Phase V_Woofer", "ep");
    limits_cl_3_legend_phase->AddEntry(phase_func_w_3_limits, "Phase V_Woofer Fit", "l");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_w_3_limits->GetChisquare() / phase_func_w_3_limits->GetNDF())).c_str()), "");
    // RIGA 3 e 4
    for (int i = 0; i != phase_func_w_3_limits->GetNpar(); ++i)
        limits_cl_3_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_w_3_limits->GetParName(i)) + " = " + NumErrScien(phase_func_w_3_limits->GetParameter(i), phase_func_w_3_limits->GetParError(i), phase_func_w_3_limits->GetParName(i))).c_str(), "");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, "", "");
    // RIGA 4
    limits_cl_3_legend_phase->AddEntry(phase_graph_limits[2], "Phase V_Tweeter", "ep");
    limits_cl_3_legend_phase->AddEntry(phase_func_t_limits, "Phase V_Tweeter Fit", "l");
    limits_cl_3_legend_phase->AddEntry((TObject *)0, (("#tilde{#chi}^{2} = " + std::to_string(phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF())).c_str()), "");
    // RIGA 5
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_cl_3_legend_phase->AddEntry((TObject *)0, (std::string(phase_func_t_limits->GetParName(i)) + " = " + NumErrScien(phase_func_t_limits->GetParameter(i), phase_func_t_limits->GetParError(i), phase_func_t_limits->GetParName(i))).c_str(), "");

    limits_cl_3_legend_phase->Draw();

    // print_res
    std::ofstream limits_Cl_3_file("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_Cl_3_" + std::to_string(GetVoltage()) + ".txt");
    limits_Cl_3_file << std::setprecision(10);
    limits_Cl_3_file << "Ampl_woofer:\n";
    limits_Cl_3_file << "\tFit Status: " + (ampl_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_woofer_fit_res) + ")")) << '\n';
    limits_Cl_3_file << "\tChi^2: " << ampl_func_w_3_limits->GetChisquare() << '\n';
    limits_Cl_3_file << "\tNdf: " << ampl_func_w_3_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tCHI RIDOTTO: " << ampl_func_w_3_limits->GetChisquare() / ampl_func_w_3_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_w_3_limits->GetNpar(); ++i)
        limits_Cl_3_file << "\t\t[" << i << "] - " << ampl_func_w_3_limits->GetParName(i) << " = " << ampl_func_w_3_limits->GetParameter(i) << " +- " << ampl_func_w_3_limits->GetParError(i) << '\n';
    limits_Cl_3_file << "Ampl_tweeter:\n";
    limits_Cl_3_file << "\tFit Status: " + (ampl_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(ampl_tweeter_fit_res) + ")")) << '\n';
    limits_Cl_3_file << "\tChi^2: " << ampl_func_t_limits->GetChisquare() << '\n';
    limits_Cl_3_file << "\tNdf: " << ampl_func_t_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tCHI RIDOTTO: " << ampl_func_t_limits->GetChisquare() / ampl_func_t_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != ampl_func_t_limits->GetNpar(); ++i)
        limits_Cl_3_file << "\t\t[" << i << "] - " << ampl_func_t_limits->GetParName(i) << " = " << ampl_func_t_limits->GetParameter(i) << " +- " << ampl_func_t_limits->GetParError(i) << '\n';
    limits_Cl_3_file << "Phase_woofer:\n";
    limits_Cl_3_file << "\tFit Status: " + (phase_woofer_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_woofer_fit_res) + ")")) << '\n';
    limits_Cl_3_file << "\tChi^2: " << phase_func_w_3_limits->GetChisquare() << '\n';
    limits_Cl_3_file << "\tNdf: " << phase_func_w_3_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tCHI RIDOTTO: " << phase_func_w_3_limits->GetChisquare() / phase_func_w_3_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_w_3_limits->GetNpar(); ++i)
        limits_Cl_3_file << "\t\t[" << i << "] - " << phase_func_w_3_limits->GetParName(i) << " = " << phase_func_w_3_limits->GetParameter(i) << " +- " << phase_func_w_3_limits->GetParError(i) << '\n';
    limits_Cl_3_file << "Phase_tweeter:\n";
    limits_Cl_3_file << "\tFit Status: " + (phase_tweeter_fit_res == 0 ? "OK" : ("ERROR(" + std::to_string(phase_tweeter_fit_res) + ")")) << '\n';
    limits_Cl_3_file << "\tChi^2: " << phase_func_t_limits->GetChisquare() << '\n';
    limits_Cl_3_file << "\tNdf: " << phase_func_t_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tCHI RIDOTTO: " << phase_func_t_limits->GetChisquare() / phase_func_t_limits->GetNDF() << '\n';
    limits_Cl_3_file << "\tParameters: \n";
    for (int i = 0; i != phase_func_t_limits->GetNpar(); ++i)
        limits_Cl_3_file << "\t\t[" << i << "] - " << phase_func_t_limits->GetParName(i) << " = " << phase_func_t_limits->GetParameter(i) << " +- " << phase_func_t_limits->GetParError(i) << '\n';
    limits_Cl_3_file.close();

    // Save Canvas
    limits_Cl_3_canvas->Update();
    limits_Cl_3_canvas->SaveAs(("./risultati_finali/Sweep_" + GetSweepRange() + "/limits_Cl_3_" + std::to_string(GetVoltage()) + ".png").c_str());
    limits_Cl_3_canvas->Update();
}
