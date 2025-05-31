#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>

#define MAX_N_BLOCKS 1000
int GetVpp()
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
    std::cout << "vpp: " << atoi(input.c_str()) << std::endl;
    settings_file.close();
    return atoi(input.c_str());
}

int GetNPoint()
{
    std::ifstream settings_file{"./input_data/settings.txt"};
    if (!settings_file.is_open())
        std::cout << "Failed to open settings.txt" << std::endl;

    std::string input;
    do
    {
        settings_file >> input;
    } while (input != "N_samples:" && !settings_file.eof());

    settings_file >> input;
    std::cout << "N_Points: " << atoi(input.c_str()) << std::endl;
    settings_file.close();
    return atoi(input.c_str());
}

double GetError(double frequency, double vpp_error_elvis, double vpp_err_freq_RMS_arr[2][MAX_N_BLOCKS], int N_VPP_ERR_FREQ_RMS)
{
    int closest_higher_freq_index = 0;
    while (closest_higher_freq_index < N_VPP_ERR_FREQ_RMS)
    {
        if (vpp_err_freq_RMS_arr[0][closest_higher_freq_index] >= frequency)
        {
            break;
        }
        ++closest_higher_freq_index;
    }

    // freq > of every freq we have data of
    if (closest_higher_freq_index == N_VPP_ERR_FREQ_RMS)
        return std::sqrt(vpp_error_elvis * vpp_error_elvis + vpp_err_freq_RMS_arr[1][closest_higher_freq_index - 1] * vpp_err_freq_RMS_arr[1][closest_higher_freq_index - 1]);

    // freq < of every freq we have data of
    if (closest_higher_freq_index == 0)
        return std::sqrt(vpp_error_elvis * vpp_error_elvis + vpp_err_freq_RMS_arr[1][closest_higher_freq_index] * vpp_err_freq_RMS_arr[1][closest_higher_freq_index]);

    // Linear Interpolation
    int closest_lower_freq_index = closest_higher_freq_index - 1;
    double RMS_m = vpp_err_freq_RMS_arr[1][closest_lower_freq_index];
    double RMS_M = vpp_err_freq_RMS_arr[1][closest_higher_freq_index];
    double freq_m = vpp_err_freq_RMS_arr[0][closest_lower_freq_index];
    double freq_M = vpp_err_freq_RMS_arr[0][closest_higher_freq_index];

    double RMS_interp = RMS_m + ((RMS_M - RMS_m) / (freq_M - freq_m)) * (frequency - freq_m);
    return std::sqrt(RMS_interp * RMS_interp + vpp_error_elvis * vpp_error_elvis);
}

int main()
{
    // clear output dir
    if (system("rm -rf ./data_root") != 0)
        std::cout << "Failed to \"rm -rf ./data_root\"\n";
    if (system("mkdir ./data_root") != 0)
        std::cout << "Failed to \"rm -rf ./data_root\"\n";
    int N_POINTS{GetNPoint()};
    double const vpp_error_elvis{1.2287040E-3};

    //[0][i] : frequency at which the rms in [1][i] was calculated
    //[1][i] : RMS calculated at the freq in [0][i]
    double vpp_err_freq_RMS_arr[2][MAX_N_BLOCKS];
    std::ifstream vpp_err_freq_rms_filein{"./risultati/Background_error/V_" + std::to_string(GetVpp()) + "_Freq_RMS.txt"};
    if (!vpp_err_freq_rms_filein.is_open())
    {
        std::cout << "Failed to open vpp_err_freq_rms_filein" << std::endl;
        return 0;
    }
    int N_VPP_ERR_FREQ_RMS = 0;
    while (!vpp_err_freq_rms_filein.eof())
    {
        vpp_err_freq_rms_filein >> vpp_err_freq_RMS_arr[0][N_VPP_ERR_FREQ_RMS];
        vpp_err_freq_rms_filein >> vpp_err_freq_RMS_arr[1][N_VPP_ERR_FREQ_RMS];
        // std::cout << vpp_err_freq_RMS_arr[0][N_VPP_ERR_FREQ_RMS] << '\t' << vpp_err_freq_RMS_arr[1][N_VPP_ERR_FREQ_RMS] << '\n';

        ++N_VPP_ERR_FREQ_RMS;
    }
    // std::cout << N_VPP_ERR_FREQ_RMS << '\n';

    // return 0;

    // // note:
    // // we took the error from block 0 - AI0 for each dataset without the circuit
    // //(a frequenze più alte sembra diventare più piccolo, abbiamo preso il circa max)
    // //  0: 5V           0.00215797
    // //  1: 7.5V (7)     0.00301076
    // //  2: 10V          0.00408209
    // double const vpp_background_error_arr[3]{
    //     0,     //  0: 5V
    //     0,     //  1: 7.5V (7)
    //     0      //  2: 10V
    // };

    // // terribile ma stfu
    // double vpp_background_error;
    // if (VPP == 5)
    //     vpp_background_error = vpp_background_error_arr[0];
    // else if (VPP == 7)
    //     vpp_background_error = vpp_background_error_arr[1];
    // else // (VPP == 10)
    //     vpp_background_error = vpp_background_error_arr[2];

    // std::string vpp_error{std::to_string(std::sqrt(vpp_error_elvis * vpp_error_elvis + vpp_background_error * vpp_background_error))};
    std::string vpp_error{"1.2287040E-3"};
    std::string const time_error{"50E-9"}; // from elvis specs
    std::string output_data_block_filename{"./data_root/output_data_block_"};
    std::string output_param_filename{"./data_root/output_param_"};

    // counting how many blocks of data we have
    // copied from https://stackoverflow.com/questions/3072795/how-to-count-lines-of-a-file-in-c
    std::ifstream file_count{"./input_data/ampl_data"};
    if (!file_count.is_open())
        std::cout << "Failed to open ampl_data" << std::endl;

    std::cout << "started n_blocks count" << std::endl;
    long int N_BLOCKS{std::count(std::istreambuf_iterator<char>(file_count),
                                 std::istreambuf_iterator<char>(), '\n') -
                      1};
    file_count.close();
    std::cout << "N_BLOCKS: " << N_BLOCKS << std::endl;

    std::ifstream file_sweep_in{"./input_data/sweep_data"};
    if (!file_sweep_in.is_open())
        std::cout << "Failed to open sweep_data" << std::endl;
    std::ifstream file_ampl_in{"./input_data/ampl_data"};
    if (!file_ampl_in.is_open())
        std::cout << "Failed to open ampl_data" << std::endl;
    std::ifstream file_phase_in{"./input_data/phase_data"};
    if (!file_phase_in.is_open())
        std::cout << "Failed to open phase_data" << std::endl;

    // skip first row (intestazione)
    char garbage = ' ';
    while (garbage != '\n')
        file_ampl_in >> std::noskipws >> garbage;
    garbage = ' ';
    while (garbage != '\n')
        file_phase_in >> std::noskipws >> garbage;
    file_ampl_in >> std::skipws;
    file_phase_in >> std::skipws;

    for (int block_n = 0; block_n != N_BLOCKS; ++block_n)
    {
        std::ofstream file_data_block_out{output_data_block_filename + std::to_string(block_n) + ".txt", std::ios::out | std::ios::trunc};
        std::ofstream file_param_out{output_param_filename + std::to_string(block_n) + ".txt", std::ios::out | std::ios::trunc};

        // get param
        // freq V_s_ampl V_s_phase V_w_ampl V_w_phase V_t_ampl V_t_phase
        std::string frequency;
        std::string amplitude;
        std::string phase;

        // V_s
        file_ampl_in >> frequency;
        double VS_frequency{std::stod(frequency)};
        file_phase_in >> frequency;
        file_ampl_in >> amplitude;
        file_phase_in >> phase;

        file_param_out << std::setprecision(10) << frequency << '\t' << amplitude << '\t' << phase << '\t';

        // V_w
        file_phase_in >> frequency;
        file_ampl_in >> frequency;
        file_ampl_in >> amplitude;
        file_phase_in >> phase;
        file_param_out << std::setprecision(10) << amplitude << '\t' << phase << '\t';

        // V_t
        file_phase_in >> frequency;
        file_ampl_in >> frequency;
        file_ampl_in >> amplitude;
        file_phase_in >> phase;
        file_param_out << std::setprecision(10) << amplitude << '\t' << phase << '\n';

        // sweep data---------------------------------------

        // first read outside the loop to get initial time
        std::string initial_time;
        std::string initial_V_s;
        std::string initial_V_w;
        std::string initial_V_t;

        file_sweep_in >> initial_time; // data che ignoriamoi
        file_sweep_in >> initial_time; // orario
        file_sweep_in >> initial_V_s;  // V_s
        file_sweep_in >> initial_V_w;  // V_w
        file_sweep_in >> initial_V_t;  // V_t

        // togli ora
        double initial_minutes = std::stod(initial_time.substr(3, 2));
        double initial_seconds = std::stod(initial_time.substr(6)) + initial_minutes * 60.;

        file_data_block_out << std::setprecision(10) << 0. << '\t' << initial_V_s << '\t' << initial_V_w << '\t' << initial_V_t << '\t' << time_error << '\t' << std::to_string(GetError(VS_frequency, vpp_error_elvis, vpp_err_freq_RMS_arr, N_VPP_ERR_FREQ_RMS)) << '\n';
        for (int i = 1; i != N_POINTS; ++i)
        {
            std::string time;
            std::string V_s;
            std::string V_w;
            std::string V_t;

            file_sweep_in >> time; // data che ignoriamoi
            file_sweep_in >> time; // orario
            file_sweep_in >> V_s;  // V_s
            file_sweep_in >> V_w;  // V_w
            file_sweep_in >> V_t;  // V_t

            // togli ora
            double minutes = std::stod(time.substr(3, 2));
            double seconds = std::stod(time.substr(6)) + minutes * 60. - initial_seconds;

            file_data_block_out << std::setprecision(10) << seconds << '\t' << V_s << '\t' << V_w << '\t' << V_t << '\t' << time_error << '\t' << std::to_string(GetError(VS_frequency, vpp_error_elvis, vpp_err_freq_RMS_arr, N_VPP_ERR_FREQ_RMS)) << '\n';
        }

        std::cout << "Done BLOCK " << block_n << '\n';
        file_data_block_out.close();
        file_param_out.close();
    }
    std::cout << "N_BLOCKS: " << N_BLOCKS << '\n';
    return 0;
}
