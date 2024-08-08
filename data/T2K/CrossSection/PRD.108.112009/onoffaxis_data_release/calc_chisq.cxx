//Example function/script to calculate the chi-sq using the ROOT file in the data release
void calc_chisq(const std::string& filename)
{
    //Open file and get the linear bin histograms for the data and nominal MC,
    //and the 2D histogram containing the covariance matrix inverse
    
    std::cout << "Opening " << filename << std::endl;
    auto input_file = TFile::Open(filename.c_str(), "READ");

    if(!input_file->IsOpen())
    {
        std::cout << "Could not open file." << std::endl;
        return;
    }

    auto h_data = input_file->Get<TH1D>("T2K_OnOffAxis_CC0pi_CH_XSec_2DPcos_joint_data");
    auto h_mc   = input_file->Get<TH1D>("T2K_OnOffAxis_CC0pi_CH_XSec_2DPcos_joint_MC");
    auto invcov = input_file->Get<TH2D>("T2K_OnOffAxis_CC0pi_CH_XSec_2DPcos_joint_INVCOV");

    double chisq = 0.0;
    const unsigned int nbins = h_data->GetNbinsX();

    //Loop over each bin and matrix element to calculate the chi-square
    for(unsigned int i = 0; i < nbins; ++i)
    {
        for(unsigned int j = 0; j < nbins; ++j)
        {
            const double x = h_data->GetBinContent(i+1) - h_mc->GetBinContent(i+1);
            const double y = h_data->GetBinContent(j+1) - h_mc->GetBinContent(j+1);
            const double v = invcov->GetBinContent(i+1, j+1);

            chisq += x * v * y;
        }
    }

    //Print the absolute chi-square and the relative chi-square
    std::cout << "Nbins : " << nbins << std::endl;
    std::cout << "Chi-sq: " << chisq << std::endl;
    std::cout << "Chi-sq / N: " << chisq / nbins << std::endl;
    return;
}
