#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

int main()
{
    // clear output dir
    if (system("rm -rf ./data_root") != 0)
        std::cout << "Failed to \"rm -rf ./data_root\"\n";
    if (system("mkdir ./data_root") != 0)
        std::cout << "Failed to \"rm -rf ./data_root\"\n";
    const int N_POINTS{10'000};
    std::string const vpp_error{"1.2287040E-3"};
    std::string const time_error{"50E-9"}; // from elvis specs
    std::string output_data_block_filename{"./data_root/output_data_block_"};
    std::string output_param_filename{"./data_root/output_param_"};

    // counting how many blocks of data we have
    // copied from https://stackoverflow.com/questions/3072795/how-to-count-lines-of-a-file-in-c
    std::ifstream file_count{"./input_data/ampl_data"};
    if(!file_count.is_open())
        std::cout<<"Failed to open ampl_data"<<std::endl;
    
    std::cout<<"started n_blocks count"<<std::endl;
    long int N_BLOCKS{std::count(std::istreambuf_iterator<char>(file_count),
                                 std::istreambuf_iterator<char>(), '\n') -
                      1};
    file_count.close();
    std::cout<<"N_BLOCKS: "<<N_BLOCKS<<std::endl;

    std::ifstream file_sweep_in{"./input_data/sweep_data"};
    if(!file_sweep_in.is_open())
        std::cout<<"Failed to open sweep_data"<<std::endl;
    std::ifstream file_ampl_in{"./input_data/ampl_data"};
    if(!file_ampl_in.is_open())
        std::cout<<"Failed to open ampl_data"<<std::endl;
    std::ifstream file_phase_in{"./input_data/phase_data"};
    if(!file_phase_in.is_open())
        std::cout<<"Failed to open phase_data"<<std::endl;

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

        file_data_block_out << std::setprecision(10) << 0. << '\t' << initial_V_s << '\t' << initial_V_w << '\t' << initial_V_t << '\t' << time_error << '\t' << vpp_error << '\n';
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

            file_data_block_out << std::setprecision(10) << seconds << '\t' << V_s << '\t' << V_w << '\t' << V_t << '\t' << time_error << '\t' << vpp_error << '\n';
        }

        // get param
        // freq V_s_ampl V_s_phase V_w_ampl V_w_phase V_t_ampl V_t_phase
        std::string frequency;
        std::string amplitude;
        std::string phase;

        // V_s
        file_ampl_in >> frequency;
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

        std::cout<<"Done BLOCK "<<block_n<<'\n';
        file_data_block_out.close();
        file_param_out.close();
    }
    return 0;
}
