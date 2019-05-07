// 	---------------------------------
//	|                               |
//	| g-2 fast rotation toy model   |
//	|				|
//	|       antoine chapelain	|
//	|       atc93@cornell.edu	|
//	|				|
// 	---------------------------------

// root include files
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TTree.h>
#include <time.h>
#include <vector>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>

// c++ include files
#include <fstream>    // for reading files
#include <sstream>

// simulation paramaters
const double N_muon = 1000000;     // number of muons to simulate
const double N_turn = 2700; // 6800;			// how many cyclotron revolution to simulate
const double p_spread = 0.00025; //0.09;      // muon momentum distribution spread in %
const double t_spread = 25;	// muon beam length spread in nano second
const bool   fillTree = false;      // fill root tree with time, momentum...
const double detLocation=0.74;       // half-way 

// physical constants
const double a_mu = 11659208.9e-10;      // anomalous moment
const double M_mu = 0.1056583715 ;      // muon mass in gev
const double M_p = 0.9382720813 ; 		// proton mass in gev
const double M_e = 0.510998928e-3 ;  	// electron mass in gev
const double q_e = 1.602176565e-19 ;	// charge of the positron in coulomb
const double c_light= 299792458 ; 		// speed of light in m/s
const double B = 1.4513; 				// b-field in t
const double p_ptcl_magic = M_mu/sqrt(a_mu); 		// muon magic momentum, in gev/c
const double R0= p_ptcl_magic*1E9/(B*c_light);					// the radius of the ideal orbit in meter
const double mass_factor = 1.782661845e-36*1e9 ;    // kg for 1gev
const double muon_life = 2.1969811e-6; 	// muon life time in seconds
const double inch = 0.0254;

// math constants
const double PI = TMath::Pi();

// debugging mode
const int DEBUG=0;

// Set momentum distribution
const int Gaus		= 1;
const int Landau	= 0;
const int Unif		= 0;
const int Custom 	= 0;
const int VacBounds     = 0;

// Set time distribution
const int t_wShape 	= 0;
const int t_Gaus 	= 1;
const int t_Asymm       = 0;

// namespaces
using namespace std;

// ppre-declaration of functions
Double_t ComputeCycloPeriod(double gamma, double beta, double p_ptcl); 
void DrawATLASLabel(TCanvas* _c);


// ---------------
//  main function
// ---------------

int main() {

    // initializaion of variables
    Double_t time;
    Double_t p_ptcl;
    Double_t E_ptcl;
    Double_t gamma;
    Double_t beta;
    TRandom *r3 = new TRandom3(0);
    gRandom = r3;

    cout << "\n-- start fast rotation routine --\n" << endl;

    cout << " momentum spread (1 sigma): " << p_spread << "%" << endl;
    cout << "     beam length (1 sigma): " << t_spread << " ns" << endl;
    cout << " number of simulated ptcle: " << N_muon << endl;
    cout << "               fill length: " << N_turn * 0.14914 << " us" << endl;

    // draw random numbers
    TFile *f_wShape = new TFile("root/wShape.root", "READ");
    TH1D *h_wShape = (TH1D*)f_wShape->Get("h_wShape");
    vector<double> v_p_spread;
    vector<double> v_t_spread;
    TF1 *f1 = new TF1("f1","((sin(50*(x))/(50*(x))+1)*exp(-0.5*((x)/0.3)**2))",-1,1);
    f1->SetNpx(10000);
    TF1 *f2 = new TF1("f2","(exp(-0.5*((x+0.00003)/(0.0005))**2))/(sqrt(2*3.1415926535)*0.0005)",-0.0025,0.0038);
    TF1 *f4 = new TF1("f4","(exp(-0.5*((x+0.0001)/(0.0005))**2))/(sqrt(2*3.1415926535)*0.0005)",-0.00025,0.001);
    //TF1 *f2 = new TF1("f2","(exp(-0.5*((x)/(0.00112))**2))/(sqrt(2*3.1415926535)*0.00112)",-0.0025,0.0025);
    TF1 *f3 = new TF1("f3","(exp(-0.5*((x-15E-9)/(15E-9))**2))/(sqrt(2*3.1415926535)*15E-9)+(exp(-0.5*((x+40E-9)/(20E-9))**2))/(sqrt(2*3.1415926535)*20E-9)",-95E-9,95E-9);
    f3->SetNpx(100000);
    for (int i=0; i<N_muon; ++i)    {
        if (Gaus)       v_p_spread.push_back(gRandom -> Gaus(0, p_spread/100));
        else if (Landau) v_p_spread.push_back(gRandom -> Landau(-0.005, p_spread/250));
        else if (Unif)       v_p_spread.push_back(gRandom -> Uniform(-5*p_spread/100,5*p_spread/100));
        else if (Custom) v_p_spread.push_back((f1->GetRandom())/150);
        else if (VacBounds) v_p_spread.push_back(4.5*f4->GetRandom()+0.75*f2->GetRandom());

        if (t_Gaus) v_t_spread.push_back(gRandom -> Gaus(0, t_spread*1E-9));
        else if (t_wShape) v_t_spread.push_back(h_wShape->GetRandom()*1E-9);
        else if (t_Asymm) v_t_spread.push_back(f3->GetRandom());
    }

    // root objects
    TFile *rootFile = new TFile("root/FastRotation.root","recreate");
    TTree *tr = new TTree("tr","tr");
    tr->Branch("time",&time,"time/d");
    tr->Branch("p_ptcl",&p_ptcl,"p_ptcl/d");

    TH1D *h_frs = new TH1D("fr", "fast rotation", 1000000, 0, 1000);
    h_frs->GetXaxis()->SetTitle("time [#mus]");
    h_frs->GetXaxis()->CenterTitle();
    h_frs->GetXaxis()->SetTitleOffset(1.1);
    h_frs->GetYaxis()->SetTitle("intensity");
    h_frs->GetYaxis()->CenterTitle();
    h_frs->GetYaxis()->SetTitleOffset(1.5);
    gStyle->SetOptStat(1100);
    TCanvas c; c.cd();
    string _level="internal simulation";
    string _atlas_desc="#font[72]{g-2}";

    TH1D *h_freq;
    h_freq  = new TH1D("freq", "freq dist", 150,6630,6780);
    h_freq->GetXaxis()->SetTitle("Frequency [Hz]");
    TH1D *h_p = new TH1D("p", "momentum dist", 340,3.0600,3.1300);
    h_p->GetXaxis()->SetTitle("Momentum [GeV]");

    TLatex _g;
    _g.SetTextSize(0.035);

    //draw horizontal label  
    string _sh="";
    _sh+=_atlas_desc+" "+_level;
    gPad->Update();

    //drawatlaslabel(c);
    gPad->SetTicks(1);
    gPad->SetTicks(1);

    // processing tim book-kepping
    clock_t t1,t2;
    t1 = clock();
    int progress = N_muon / 10;

    double cyclo_period = 0.;

    for (int i=0; i<N_muon; ++i) {

        // draw a random muon momentum
        p_ptcl = (1+v_p_spread.at(i))*p_ptcl_magic;
        h_p->Fill(p_ptcl);

        // compute particle energy
        E_ptcl = sqrt( p_ptcl * p_ptcl + M_mu * M_mu);

        // compute gamma factor
        gamma = E_ptcl / M_mu; 

        // compute beta factor
        beta = p_ptcl / (gamma * M_mu);

        // compute cyclotron frequency
        cyclo_period = ComputeCycloPeriod(gamma, beta, p_ptcl); 

        // original time offset
        double dt = v_t_spread.at(i)*1E6;

        // loop over tusn
        for (int turn=0; turn<N_turn; ++turn) {

            time = (detLocation + turn) * cyclo_period + dt;
            h_frs->Fill(time);

        } // end turn loop

        h_freq->Fill(1E3/cyclo_period);

        // Monitor progress
        if(i % progress == 0){
            float percent = i / N_muon;
            cout << "--> " << percent*100 << "\% accomplished        (" << i << "/" << N_muon << "  muons processed)" << endl;
        }

    } // end Muon loop

    // Processing time book-keeping
    t2 = clock();
    float diff ((float)t2-(float)t1);
    cout << "\n  -- Routine DONE --\n" << endl;
    cout << "Time elapsed: " << diff/CLOCKS_PER_SEC/60 << " minutes." << endl;

    // Fill radial distribution
    vector<float> radii;
    vector<float> intensity;
    double meanR = 0;
    double sum = 0;
    for ( int iBin = 1; iBin<=h_freq->GetNbinsX(); ++iBin ) {
        radii.push_back( 0.9994174214209439 * c_light / (2*PI*h_freq->GetBinCenter(iBin)) );
        intensity.push_back( h_freq->GetBinContent(iBin) );
        meanR += radii.at(iBin-1)*intensity.at(iBin-1);
        sum += intensity.at(iBin-1);
    }

    cout << "Mean radius = " << meanR/sum << std::endl;
    cout << "Mean freq = " << h_freq->GetMean() << std::endl;

    TGraph g(radii.size(),radii.data(),intensity.data());
    g.SetMarkerStyle(20);
    h_frs->Scale(1/ (N_muon/(1000000/h_freq->GetMean())));
    h_frs->Write();
    h_freq->Write();
    h_p->Write();
    g.Write("r");

    h_frs->Draw();
    c.SaveAs("plot/FastRotation.eps");
    c.SaveAs("plot/FastRotation.png");

    h_freq->Draw();
    c.SaveAs("plot/FreqDist.eps");
    c.SaveAs("plot/FreqDist.png");

    h_p->Draw();
    c.SaveAs("plot/MomentumDist.eps");
    c.SaveAs("plot/MomentumDist.png");

    //    tr->Write();
    rootFile->Close();

    return 0;

}

// function to compute travel time for the case:
// - longitudinal initial beam width
// - momentum spread

Double_t ComputeCycloPeriod(double gamma, double beta, double p_ptcl) {

    double rho = p_ptcl / (c_light * 1E-9 * B) *1000;
    double t = 2*PI*rho / ( beta*c_light);
    return t*1E3;

}

