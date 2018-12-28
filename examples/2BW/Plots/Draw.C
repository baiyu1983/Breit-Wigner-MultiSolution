#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TComplex.h"
#include "Draw.h"
#include "TMath.h"
#include "TAxis.h"
#define MJPsi 3.0969
#define MPi 0.1395706
#define MELE 0.000511
#define UNIT 0.38936951 // cross-section in pb
using namespace std;

void Draw(){
    TCanvas *c1 = new TCanvas("c1","Combination of two Breit-Weigner functions",150,150,600,600);
    TCanvas *c2 = new TCanvas("c2","Combination of two Breit-Weigner functions",150,150,600,600);

    double Mass[2] = {4.0,4.2};
    double Width[2] = {0.1,0.15};

    double PhaseAngle[2][2]={{0,0.78539816},{1.462686,-0.600191}};
    double BranchWidth[2][2]={{1.0,0.7},{2.434129,3.026686}};

//    double PhaseSpace[2] = {Phi23JPiPiPsi(Mass[0]),Phi23JPiPiPsi(Mass[1])};

    vector<resonance> states_single;       states_single.clear();
    vector<resonance> states_combine_1;    states_combine_1.clear();
    vector<resonance> states_combine_2;    states_combine_2.clear();

    TF1 * resonance_1_combine_1 = new TF1("resonance_1_combine_1",amp2TF,3.7,4.6,6);
    TF1 * resonance_1_combine_2 = new TF1("resonance_1_combine_2",amp2TF,3.7,4.6,6);
    TF1 * resonance_2_combine_1 = new TF1("resonance_2_combine_1",amp2TF,3.7,4.6,6);
    TF1 * resonance_2_combine_2 = new TF1("resonance_2_combine_2",amp2TF,3.7,4.6,6);
    TF1 * resonance_combine_1 = new TF1("resonance_combine_1",amp2TF,3.7,4.6,10);
    TF1 * resonance_combine_2 = new TF1("resonance_combine_2",amp2TF,3.7,4.6,10);
 
//    cout<<"Mark A"<<endl;
    resonance_1_combine_1->SetParameter(0,1);     resonance_2_combine_1->SetParameter(0,1);
    resonance_1_combine_2->SetParameter(0,1);     resonance_2_combine_2->SetParameter(0,1);
    resonance_1_combine_1->SetParameter(1,Mass[0]);     resonance_2_combine_1->SetParameter(1,Mass[1]);
    resonance_1_combine_2->SetParameter(1,Mass[0]);     resonance_2_combine_2->SetParameter(1,Mass[1]);
    resonance_1_combine_1->SetParameter(2,Width[0]);    resonance_2_combine_1->SetParameter(2,Width[1]);
    resonance_1_combine_2->SetParameter(2,Width[0]);    resonance_2_combine_2->SetParameter(2,Width[1]);
    resonance_1_combine_1->SetParameter(4,BranchWidth[0][0]); resonance_2_combine_1->SetParameter(4,BranchWidth[0][1]);
    resonance_1_combine_2->SetParameter(4,BranchWidth[1][0]); resonance_2_combine_2->SetParameter(4,BranchWidth[1][1]);
    resonance_1_combine_1->SetParameter(3,PhaseAngle[0][0]); resonance_2_combine_1->SetParameter(3,PhaseAngle[0][1]);
    resonance_1_combine_2->SetParameter(3,PhaseAngle[1][0]); resonance_2_combine_2->SetParameter(3,PhaseAngle[1][1]);
    resonance_1_combine_1->SetParameter(6,0); resonance_2_combine_1->SetParameter(6,0);
    resonance_1_combine_2->SetParameter(6,0); resonance_2_combine_2->SetParameter(6,0);
    resonance_1_combine_1->SetParameter(5,0); resonance_2_combine_1->SetParameter(5,0);
    resonance_1_combine_2->SetParameter(5,0); resonance_2_combine_2->SetParameter(5,0);

    resonance_combine_1->SetParameter(0,2);                   resonance_combine_2->SetParameter(0,2);
    resonance_combine_1->SetParameter(1,Mass[0]);             resonance_combine_2->SetParameter(1,Mass[0]);
    resonance_combine_1->SetParameter(2,Width[0]);             resonance_combine_2->SetParameter(2,Width[0]);
    resonance_combine_1->SetParameter(3,PhaseAngle[0][0]);    resonance_combine_2->SetParameter(3,PhaseAngle[1][0]);
    resonance_combine_1->SetParameter(4,BranchWidth[0][0]);   resonance_combine_2->SetParameter(4,BranchWidth[1][0]);
    resonance_combine_1->SetParameter(5,Mass[1]);             resonance_combine_2->SetParameter(5,Mass[1]);
    resonance_combine_1->SetParameter(6,Width[1]);             resonance_combine_2->SetParameter(6,Width[1]);
    resonance_combine_1->SetParameter(7,PhaseAngle[0][1]);    resonance_combine_2->SetParameter(7,PhaseAngle[1][1]);
    resonance_combine_1->SetParameter(8,BranchWidth[0][1]);   resonance_combine_2->SetParameter(8,BranchWidth[1][1]);
    resonance_combine_1->SetParameter(9,0);                   resonance_combine_2->SetParameter(9,0);
    resonance_combine_1->SetParameter(10,0);                  resonance_combine_2->SetParameter(10,0);
//    cout<<"Mark B"<<endl;

    c1->cd();

    cout<<"Para1 = "<<resonance_1_combine_1->GetParameter(0)<<endl;
    cout<<"Para2 = "<<resonance_1_combine_1->GetParameter(1)<<endl;

    cout<<"Image unit = "<<TComplex::I()<<",Amp at 3.9 = "<<resonance_combine_1->Eval(3.9)<<endl;
    cout<<"Image unit = "<<TComplex::I()<<",Amp at 3.9 = "<<resonance_combine_2->Eval(3.9)<<endl;
    cout<<"Image unit = "<<TComplex::I()<<",Amp at 4.05 = "<<resonance_combine_1->Eval(4.05)<<endl;
    cout<<"Image unit = "<<TComplex::I()<<",Amp at 4.05 = "<<resonance_combine_2->Eval(4.05)<<endl;
    resonance_1_combine_1->SetLineColor(2);
    resonance_2_combine_1->SetLineColor(4);
    resonance_1_combine_2->SetLineColor(2);
    resonance_2_combine_2->SetLineColor(4);
    resonance_1_combine_1->SetLineStyle(2);
    resonance_2_combine_1->SetLineStyle(2);
    resonance_1_combine_2->SetLineStyle(2);
    resonance_2_combine_2->SetLineStyle(2);
    resonance_combine_1->SetLineColor(kBlack);
    resonance_combine_2->SetLineColor(kBlack);

    resonance_1_combine_1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_1_combine_1->GetYaxis()->SetTitle("#sigma [pb]");
    resonance_1_combine_2->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_1_combine_2->GetYaxis()->SetTitle("#sigma [pb]");
    resonance_1_combine_1->GetYaxis()->SetRangeUser(0,25);
    resonance_2_combine_1->GetYaxis()->SetRangeUser(0,25);
    resonance_1_combine_2->GetYaxis()->SetRangeUser(0,25);
    resonance_2_combine_2->GetYaxis()->SetRangeUser(0,25);
    resonance_1_combine_1->GetXaxis()->SetLabelSize(0.05);
    resonance_1_combine_1->GetYaxis()->SetLabelSize(0.05);
    resonance_1_combine_2->GetXaxis()->SetLabelSize(0.05);
    resonance_1_combine_2->GetYaxis()->SetLabelSize(0.05);
    resonance_1_combine_1->GetXaxis()->SetTitleSize(0.05);
    resonance_1_combine_1->GetYaxis()->SetTitleSize(0.05);
    resonance_1_combine_2->GetXaxis()->SetTitleSize(0.05);
    resonance_1_combine_2->GetYaxis()->SetTitleSize(0.05);

    resonance_1_combine_1->GetXaxis()->SetTitleOffset(0);
    resonance_1_combine_1->GetYaxis()->SetTitleOffset(0.97);
    resonance_1_combine_2->GetXaxis()->SetTitleOffset(0);
    resonance_1_combine_2->GetYaxis()->SetTitleOffset(0.97);
    
    resonance_1_combine_1->GetXaxis()->SetDecimals(2);
    resonance_2_combine_1->GetXaxis()->SetDecimals(2);
    resonance_1_combine_2->GetXaxis()->SetDecimals(2);
    resonance_2_combine_2->GetXaxis()->SetDecimals(2);
    resonance_combine_1->GetXaxis()->SetDecimals(2); 
    resonance_combine_2->GetXaxis()->SetDecimals(2);
    resonance_1_combine_1->GetXaxis()->CenterTitle(kTRUE);
    resonance_2_combine_1->GetXaxis()->CenterTitle(kTRUE);
    resonance_1_combine_2->GetXaxis()->CenterTitle(kTRUE);
    resonance_2_combine_2->GetXaxis()->CenterTitle(kTRUE);
    resonance_1_combine_1->GetYaxis()->CenterTitle(kTRUE);
    resonance_2_combine_1->GetYaxis()->CenterTitle(kTRUE);
    resonance_1_combine_2->GetYaxis()->CenterTitle(kTRUE);
    resonance_2_combine_2->GetYaxis()->CenterTitle(kTRUE);

    c1->cd();
    resonance_1_combine_1->Draw();
    resonance_2_combine_1->Draw("same");
    resonance_combine_1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_combine_1->GetYaxis()->SetTitle("#sigma [pb]");
    resonance_combine_1->GetXaxis()->SetLabelSize(0.05);
    resonance_combine_1->GetYaxis()->SetLabelSize(0.05);
    resonance_combine_1->GetXaxis()->SetTitleSize(0.05);
    resonance_combine_1->GetYaxis()->SetTitleSize(0.05);
    resonance_combine_1->GetXaxis()->SetTitleOffset(0);
    resonance_combine_1->GetYaxis()->SetTitleOffset(0.97);
    resonance_combine_1->Draw("same");
    c1->Draw();
    c1->Print("Spectrum_1.pdf");    

    c2->cd();
    resonance_1_combine_2->Draw();
    resonance_2_combine_2->Draw("same");
    resonance_combine_2->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_combine_2->GetYaxis()->SetTitle("#sigma [pb]");
    resonance_combine_2->GetXaxis()->SetLabelSize(0.05);
    resonance_combine_2->GetYaxis()->SetLabelSize(0.05);
    resonance_combine_2->GetXaxis()->SetTitleSize(0.05);
    resonance_combine_2->GetYaxis()->SetTitleSize(0.05);
    resonance_combine_2->GetXaxis()->SetTitleOffset(0);
    resonance_combine_2->GetYaxis()->SetTitleOffset(0.97);
    resonance_combine_2->Draw("same");

    c2->Draw();
    c2->Print("Spectrum_2.pdf");

}

double Phi23(double sqrts,double M1,double M2,double M3,double M4,double M5){
      double nn = 100000;
      double E5min = M5;
      double E5max = sqrts/2 - ((M3+M4)*(M4+M3)-M5*M5)/(2*sqrts);
      double step = (E5max-E5min)/(nn);
      double P23 = 0;
//      gRandom->SetSeed(0);

      for(double i=1;i<=nn;i++){
           double E5 = M5+i*step;
           double sigma = sqrts - E5;
           double tau  = sigma*sigma-(E5*E5-M5*M5);
           double tau0 = fabs((tau-(M3+M4)*(M3+M4))*(tau-(M3-M4)*(M3-M4)));
           double E3max = 0.5/tau*(sigma*(tau+M3*M3-M4*M4)+sqrt(E5*E5-M5*M5)*sqrt(tau0));
           double E3min = 0.5/tau*(sigma*(tau+M3*M3-M4*M4)-sqrt(E5*E5-M5*M5)*sqrt(tau0));
           P23 = P23 + (E3max-E3min)*step;
      }

      double K1 = sqrt((sqrts*sqrts-(M1+M2)*(M1+M2))*(sqrts*sqrts-(M1-M2)*(M1-M2)))/(2*sqrts);
      P23 /= (128*TMath::Pi()*TMath::Pi()*TMath::Pi()*K1*sqrts);
//      cout<<"P23 = "<<P23<<endl;
      return P23;
}

double Phi23PiPiJPsi(double sqrts){
      return Phi23(sqrts,MELE,MELE,MPi,MPi,MJPsi);
}

TComplex amp(vector<resonance> states, double sqrts){
       TComplex result;
//       cout<<"sqrts = "<<sqrts<<endl;
       for(size_t i =0; i<states.size();i++ ){
            TComplex this_amp = TComplex::One();
            double M = states.at(i).mass;
            double WT = states.at(i).width;
            double Phase = states.at(i).phase;
            double BR = states.at(i).branchratio;
/*
            cout<<"M = "<<M<<endl;
            cout<<"WT = "<<WT<<endl;
            cout<<"Phase = "<<Phase<<endl;
            cout<<"BR = "<<BR<<endl;
*/
            this_amp *= (M/sqrts);
//            cout<<"this_amp = "<<this_amp<<endl;
            this_amp *= sqrt(12*TMath::Pi()*WT*BR);
//            cout<<"this_amp = "<<this_amp<<endl;
            this_amp *= TComplex::Exp(TComplex::I()*Phase);
//            cout<<"this_amp = "<<this_amp<<endl;
            this_amp *= sqrt(Phi23PiPiJPsi(sqrts)/Phi23PiPiJPsi(M));
//            cout<<"this_amp = "<<this_amp<<endl;
            TComplex BW = TComplex::One()*(sqrts*sqrts-M*M)+TComplex::I()*(M*WT);
//            cout<<"BW = "<<BW<<endl;
//            cout<<"this_amp = "<<this_amp<<endl;
            result += this_amp/BW;
       }
//       cout<<"result = "<<result<<endl;
       return result;
}

double amp2TF(double* x, double * mypar){
    int nres=int(mypar[0]);
    vector<resonance> v_state; v_state.clear();
    for(size_t i=0; i<nres;i++){
        resonance tmp_state;
        tmp_state.mass = mypar[4*i+1];
        tmp_state.width = mypar[4*i+2];
        tmp_state.phase = mypar[4*i+3];
        tmp_state.branchratio = mypar[4*i+4];
        v_state.push_back(tmp_state);
    }
//    cout<<"x[0] = "<<x[0]<<endl;
    return UNIT*amp(v_state,x[0]).Rho2()+x[0]*mypar[4*nres+2]-mypar[4*nres+1];
//    return amp(v_state,x[0]);
}
