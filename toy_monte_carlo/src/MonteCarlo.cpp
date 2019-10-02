//  ---------------------------------
//  |                               |
//  |  g-2 fast rotation toy model  |
//  |                               |
//  |       antoine chapelain       |
//  |       atc93@cornell.edu       |
//  |                               |
//  |         tyler barrett         |
//  |       tjb269@cornell.edu      |
//  |                               |
//  ---------------------------------

// Includes. =======================================================================================

// ROOT includes.
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>

// C++ includes.
#include <fstream>
#include <time.h>
#include <vector>
#include <cmath>

// Simulation parameters. ==========================================================================

// output filename
const std::string outputPath = "../data/smoothW-correlation-oneTurn.root";

// number of muons to simulate
const double n_muons = 1E6;

// how many cyclotron revolutions to simulate
const double n_turns = 1;

// (optional) input 2D histogram for momentum offset (y) and time (x)
const bool inputCorrelation = false;
const std::string inputCorrelationPath = "../../tarazona.root";
const std::string inputCorrelationName = "correlation";

// default correlation amplitude between frequency and beam profile
// (if input correlation not provided)
const double correlation = 0.005;

// (optional) input TF1 for beam profile
const bool inputBeam = false;
const std::string inputBeamPath = "beamProfile.root";
const std::string inputBeamName = "profile";

TF1 beamProfile (
  "beam",
    "1929 * exp(-(x + 58.66)^2 / (2 * 11.86^2))"
  "+ 4819 * exp(-(x + 2.598)^2 / (2 * 21.27^2))"
  "+ 1208 * exp(-(x - 47.43)^2 / (2 * 10.63^2))",
  -100,
  100
);

/*
TF1 beamProfile (
  "beam",
  "exp(-x^2 / (2 * 25^2)) / (sqrt(2 * pi) * 25)",
  -100,
  100
);
*/

// default muon beam profile Gaussian spread in nanoseconds
// (if input beam profile not provided)
// const double t_spread = 25;

// (optional) input 1D histogram for frequency distribution
const bool inputFrequency = false;
const std::string inputFrequencyPath = "";
const std::string inputFrequencyName = "";

// default muon frequency distribution Gaussian spread in %
// (if input frequency distribution not provided)
const double f_spread = 0.13;

// Declarations. ===================================================================================

double frequency_mean(double time);
double frequency_to_momentum(double frequency);
double momentum_to_frequency(double momentum);

// Constants. ======================================================================================

// Physical constants.
const double a_mu = 11659208.9E-10;                      // anomalous magnetic moment
const double q_mu = 1.602176565E-19;                     // charge of the muon in coulombs
const double c = 299792458;                              // speed of light in m/s
const double m_mu_GeV = 0.1056583715;                    // muon mass in GeV
const double m_mu_kg = m_mu_GeV * 1E9 * q_mu / (c * c);  // muon mass in kg

// Math constants.
const double pi = TMath::Pi();

// Experimental constants.
const double n = 0.108;                 // electric field gradient index
const double b = 1.4513;                // magnetic field in teslas
const double detector_offset = 0.74;    // t_0 in units of cyclotron periods
const double collimator_low = 6662E3;   // collimator lower limit in Hz
const double collimator_high = 6747E3;  // collimator upper limit in Hz

// Ideal relativstic gamma factor and corresponding speed (natural units).
const double gamma_magic = std::sqrt(1 / a_mu + 1);
const double beta_magic = std::sqrt(1 - 1 / (gamma_magic * gamma_magic));

// Ideal cyclotron frequency (Hz), period (nanoseconds), and radius (meters).
const double f_magic = q_mu * b / (2 * pi * gamma_magic * m_mu_kg);
const double t_magic = (1 / f_magic) * 1E9;
const double r_magic = gamma_magic * beta_magic * m_mu_kg * c / (q_mu * b);
const double p_magic = frequency_to_momentum(f_magic);

// Main function. ==================================================================================

int main() {

  // Initialize random number generator.
  TRandom* random = new TRandom3(0);

  // Padding for printed output.
  std::cout << std::endl;

  // Echo constants.
  std::cout << "gamma_magic: " << gamma_magic << std::endl;
  std::cout << "beta_magic: " << beta_magic << std::endl;
  std::cout << "f_magic: " << f_magic / 1000 << " kHz" << std::endl;
  std::cout << "p_magic: " << p_magic << " GeV/c" << std::endl;
  std::cout << "t_magic: " << t_magic << " ns" << std::endl;
  std::cout << "r_magic: " << r_magic << " m" << std::endl;
  std::cout << std::endl;

  // Echo simulation parameters.
  std::cout << "muons: " << n_muons << std::endl;
  std::cout << "fill length: " << (n_turns * t_magic) * 1E3 << " us" << std::endl;

  if (inputFrequency) {
    std::cout << "frequency: " << inputFrequencyPath << std::endl;
  } else {
    std::cout << "frequency spread (1 sigma): " << f_spread << "%" << std::endl;
  }

  if (inputBeam) {
    std::cout << "beam profile: " << inputBeamPath << std::endl;
  } else {
    std::cout << "beam profile: custom definition" << std::endl;
    // std::cout << "beam length (1 sigma): " << t_spread << " ns" << std::endl;
  }

  if (inputCorrelation) {
    std::cout << "beam-frequency correlation: " << inputCorrelationPath << std::endl;
  } else {
    std::cout << "beam-frequency correlation: " << correlation << std::endl;
  }

  std::cout << std::endl;

  // Open the ROOT output file.
  TFile* outputFile = new TFile(outputPath.c_str(), "recreate");

  // Prepare the histogram for the fast rotation signal.
  TH1D fastRotationSignal (
    "fr",
    "Fast Rotation Signal;Time (us);Intensity",
    1000000,  // number of time bins
    0,        // signal start time (us)
    1000      // signal end time (us)
  );

  // Prepare the frequency distribution.
  TH1D frequencyDistribution (
    "freq",
    "Frequency Distribution;Frequency (kHz)",
    150,   // number of frequency bins
    6630,  // lower frequency limit (kHz)
    6780   // upper frequency limit (kHz)
  );

  // Prepare the momentum distribution.
  TH1D momentumDistribution (
    "p",
    "Momentum Distribution;Momentum (GeV/c)",
    340,     // number of momentum bins
    3.0600,  // lower momentum limit (GeV/c)
    3.1300   // upper momentum limit (GeV/c)
  );

  // Prepare the beam-frequency correlation histogram.
  TH2D correlationDistribution (
    "corr",
    "Frequency-Beam Correlation;Time Offset (ns);Frequency (kHz)",
    150,   // number of time offset bins
    -75,   // lower time offset limit (ns)
    75,    // upper time offset limit (ns)
    150,   // number of frequency bins
    6630,  // lower frequency limit (kHz)
    6780   // upper frequency limit (kHz)
  );

  // Prepare the TTree to hold the muon details.
  TTree* data = new TTree("data", "data");
  double beam_offset = 0;
  double frequency = 0;
  data -> Branch("time", &beam_offset, "time/d");
  data -> Branch("freq", &frequency, "freq/d");

  // Initialize the canvas.
  TCanvas canvas;
  canvas.cd();

  TFile* inputCorrelationFile;
  TH2D h_inputCorrelation;
  if (inputCorrelation) {
    inputCorrelationFile = TFile::Open(inputCorrelationPath.c_str());
    h_inputCorrelation = *((TH2D*) inputCorrelationFile -> Get(inputCorrelationName.c_str()));
  }

  TFile* inputBeamFile;
  TH1D* h_inputBeam;
  if (inputBeam) {
    inputBeamFile = TFile::Open(inputBeamPath.c_str());
    h_inputBeam = (TH1D*) inputBeamFile -> Get(inputBeamName.c_str());
  }

  // E-field correction.
  //double c_e (0);

  // Computation timing.
  clock_t t1, t2;
  t1 = clock();
  int benchmark = n_muons / 100;

  // Main simulation loop. =========================================================================

  // Loop over each individual muon.
  for (int muon = 0; muon < n_muons; ++muon) {

    // If no input correlation is specified, use separate beam and frequency distributions.
    if (!inputCorrelation) {

      // Draw this muon's location in the beam profile.
      if (!inputBeam) {
        beam_offset = beamProfile.GetRandom();
      } else {
        beam_offset = h_inputBeam -> GetRandom() * 1E9;
      }

      double mean = frequency_mean(beam_offset);

      // Draw this muon's cyclotron frequency.
      frequency = random -> Gaus(mean, mean * (f_spread / 100));

    // If an input correlation is specified, use it instead.
    } else {

      // Draw beam offset and momentum offset from histogram.
      double momentum_offset = 0;
      h_inputCorrelation.GetRandom2(beam_offset, momentum_offset);

      // Convert beam offset to nanoseconds.
      beam_offset *= 1E9;

      // Convert momentum offset to raw momentum.
      double momentum = (momentum_offset + 1) * p_magic;

      // Convert momentum to frequency.
      frequency = momentum_to_frequency(momentum);

    }

    // Fill the frequency (kHz) and beam (ns) correlation histogram for this muon.
    correlationDistribution.Fill(beam_offset, frequency / 1000);

    // Check if this muon gets vetoed by the collimators.
    if (frequency < collimator_low || frequency > collimator_high) {
      continue;
    }

    // Compute this muon's momentum and cyclotron period (nanoseconds).
    double momentum = frequency_to_momentum(frequency);
    double cyclotron_period = (1 / frequency) * 1E9;

    // Fill the frequency (kHz) and momentum (GeV/c) histograms for this muon.
    frequencyDistribution.Fill(frequency / 1000);
    momentumDistribution.Fill(momentum);

    data -> Fill();

    // Loop over each turn for the current muon.
    double time = 0;
    for (int turn = 0; turn < n_turns; ++turn) {
      // Calculate the time when this muon hits the detector on the current turn.
      time = (detector_offset + turn) * cyclotron_period + beam_offset;
      fastRotationSignal.Fill(time / 1E3);
    }

    // Monitor progress.
    if(muon % benchmark == 0){
        float percent = muon / n_muons * 100;
        std::cout
          << "--> " << percent << "\% accomplished "
          << "(" << muon << "/" << n_muons << "  muons processed)"
          << std::endl;
    }

  }

  // Notify the user of simulation completion.
  t2 = clock();
  float diff ((float) t2 - (float) t1);
  std::cout << "\n  -- Routine DONE --\n" << std::endl;
  std::cout << "Time elapsed: " << diff / CLOCKS_PER_SEC / 60 << " minutes." << std::endl;

  // Calculate the radial distribution and its mean.
  double sum = 0;
  double mean_x = 0;
  std::vector<double> radii (frequencyDistribution.GetNbinsX());
  std::vector<double> intensities (frequencyDistribution.GetNbinsX());
  for (int iBin = 1; iBin <= frequencyDistribution.GetNbinsX(); ++iBin) {
    radii[iBin - 1] = beta_magic * c / (2 * pi * frequencyDistribution.GetBinCenter(iBin));
    intensities[iBin - 1] = frequencyDistribution.GetBinContent(iBin);
    mean_x += radii[iBin - 1] * intensities[iBin - 1];
    sum += intensities[iBin - 1];
  }
  mean_x /= sum;

  /*
  // Calculate the standard deviation of the radial distribution.
  double std_x = 0;
  for (unsigned int i = 0; i < radii.size(); i++) {
    std_x += intensities[i] * (radii[i] - mean_x) * (radii[i] - mean_x);
  }
  std_x = std::sqrt(std_x / (sum - 1));
  */

  // Calculate the e-field correction.
  //c_e = -2 * n * (1 - n) * (beta_magic * beta_magic) * (mean_x * mean_x + std_x * std_x) / (r_magic * r_magic);

  std::cout << "Mean radius = " << mean_x << " mm" << std::endl;
  std::cout << "Mean frequency = " << frequencyDistribution.GetMean() << " kHz" << std::endl;
  //std::cout << "E-field correction = " << c_e << " ppb" << std::endl;

  TGraph graph (radii.size(), radii.data(), intensities.data());
  graph.SetMarkerStyle(20);

  // Fit a linear function to the correlation distribution.
  TProfile* frequencyMeans = correlationDistribution.ProfileX();
  frequencyMeans -> Fit("pol2");

  fastRotationSignal.Scale(1 / (n_muons / (1000000 / frequencyDistribution.GetMean())));

  outputFile -> cd();

  fastRotationSignal.Write();
  frequencyDistribution.Write();
  momentumDistribution.Write();

  fastRotationSignal.Draw();
  canvas.SaveAs("plot/FastRotation.png");

  frequencyDistribution.Draw();
  canvas.SaveAs("plot/FreqDist.pdf");

  momentumDistribution.Draw();
  canvas.SaveAs("plot/MomentumDist.pdf");

  gStyle -> SetOptFit(111);
  gStyle -> SetStatH(0.1);
  gStyle -> SetStatW(0.1);
  correlationDistribution.SetStats(0);
  correlationDistribution.Draw("COLZ");
  frequencyMeans -> Draw("SAMES");
  canvas.SaveAs("plot/Correlation.pdf");

  correlationDistribution.Write();
  graph.Write("r");
  data -> Write();

  outputFile -> Close();

  return 0;

}

// Helper function. ================================================================================

double frequency_mean(double time) {
  // return f_magic;
  // return f_magic * (1 + correlation * (0 - beam_offset) / t_magic);
  return 1000 * (6703 - 0.04705 * time - 0.001174 * time * time);
}

double frequency_to_momentum(double frequency) {
  double gamma = q_mu * b / (2 * pi * frequency * m_mu_kg);
  double energy = gamma * m_mu_GeV;
  return std::sqrt(energy * energy - m_mu_GeV * m_mu_GeV);
}

double momentum_to_frequency(double momentum) {
  double energy = std::sqrt(momentum * momentum + m_mu_GeV * m_mu_GeV);
  double gamma = energy / m_mu_GeV;
  return q_mu * b / (2 * pi * gamma * m_mu_kg);
}
