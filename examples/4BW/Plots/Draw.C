#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TLatex.h"
#include "TPad.h"
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
    TCanvas *c1 = new TCanvas("c1","Combination of two Breit-Weigner functions",150,10,1500,1000);
    c1->Divide(4,2,0.0,-0.02);

    double Mass[4] = {3.7,3.9,4.2,4.4};
    double Width[4] = {0.1,0.2,0.15,0.25};

    double PhaseAngle[8][4]={{0,0.7,1.4,2.1},{0,3.862422,5.304776,6.187700},{0,1.078305,4.432174,0.309779},{0,0.948158,2.840070,5.739634},{0,4.240727,2.053764,4.397480},{0, 
    4.110580,0.461660,3.544149},{0,1.326463,5.872243,3.949414},{0,4.488885,3.493834,1.753929}};
    double BranchWidth[8][4]={{1.0,1.0,1.0,1.0},{2.888769,4.679042,1.214551,1.163425},{1.154217,1.981055,5.500061,2.371476},{1.103678,1.464017,3.176030,7.674509},{3.334267,9.269441,6.680104,2.759035},{3.188271,6.850195,
    3.857449,8.928718},{1.273885,2.900298,17.468357,18.199914},{3.679958,13.570615,21.216208,21.174242}};
    string s_index[8]={"(I)","(II)","(III)","(IV)","(V)","(VI)","(VII)","(VIII)"};
//    double PhaseSpace[2] = {Phi23JPiPiPsi(Mass[0]),Phi23JPiPiPsi(Mass[1])};

/*
    vector<resonance> states_combine[8];    
    for(size_t i=0;i<8;i++){
        states_combine[i].clear();
    }

    TF1 * resonance_1_combine_1 = new TF1("resonance_1_combine_1",amp2TF,3.7,4.6,6);
    TF1 * resonance_1_combine_2 = new TF1("resonance_1_combine_2",amp2TF,3.7,4.6,6);
    TF1 * resonance_2_combine_1 = new TF1("resonance_2_combine_1",amp2TF,3.7,4.6,6);
    TF1 * resonance_2_combine_2 = new TF1("resonance_2_combine_2",amp2TF,3.7,4.6,6);
*/
    TF1 * resonance[8][4];
    TF1 * resonance_combine[8];
    for(size_t i=0; i<8;i++){
       string tmp_name("resonance_combine");
       tmp_name+="_";
       tmp_name+=to_string(i+1);
       resonance_combine[i] = new TF1(tmp_name.c_str(),amp2TF,3.51,5.1,19);
       resonance_combine[i]->SetParameter(17,0);  resonance_combine[i]->SetParameter(18,0);
       resonance_combine[i]->SetParameter(0,4);
       for(size_t k=1; k<=16;k++){
             int ir = (k-1)/4;
             int ip = k%4;
             if(ip==1){
                  resonance_combine[i]->SetParameter(k,Mass[ir]);
             }else if(ip==2){
                  resonance_combine[i]->SetParameter(k,Width[ir]);
             }else if(ip==3){
                  resonance_combine[i]->SetParameter(k,PhaseAngle[i][ir]);
             }else{
                  resonance_combine[i]->SetParameter(k,BranchWidth[i][ir]);
             }
       }
       resonance_combine[i]->SetLineColor(kBlack);

       if(i<4)resonance_combine[i]->GetYaxis()->SetRangeUser(0,33);
       else resonance_combine[i]->GetYaxis()->SetRangeUser(0,83);

//       resonance_combine[i]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
//       resonance_combine[i]->GetYaxis()->SetTitle("#sigma [pb]");
       resonance_combine[i] -> GetXaxis()->SetLabelSize(0.07);
       resonance_combine[i] -> GetYaxis()->SetLabelSize(0.07); 
       resonance_combine[i]->GetXaxis()->SetTitleOffset(0);
       resonance_combine[i]->GetYaxis()->SetTitleOffset(0.90);
       resonance_combine[i]->GetXaxis()->SetTitleSize(0.05);
       resonance_combine[i]->GetYaxis()->SetTitleSize(0.05);

       if(i>=4) resonance_combine[i]->GetYaxis()->SetTitleOffset(1.3);

       for(size_t j=0;j<4;j++){
             string name_resonance("resonance");
             name_resonance += "_";
             name_resonance += to_string(i+1);
             name_resonance += "_";
             name_resonance += to_string(j+1);
             resonance[i][j] = new TF1(name_resonance.c_str(),amp2TF,3.51,5.1,7);
             resonance[i][j] ->SetParameter(0,1);
             resonance[i][j] ->SetParameter(1,Mass[j]);
             resonance[i][j] ->SetParameter(2,Width[j]);
             resonance[i][j] ->SetParameter(3,PhaseAngle[i][j]);
             resonance[i][j] ->SetParameter(4,BranchWidth[i][j]);
             resonance[i][j] ->SetParameter(5,0);
             resonance[i][j] -> SetParameter(6,0);
             resonance[i][j] -> SetLineStyle(2);            
             resonance[i][j] -> GetXaxis()->SetLabelSize(0.07);
             resonance[i][j] -> GetYaxis()->SetLabelSize(0.07);   
             resonance[i][j]->GetXaxis()->SetTitleSize(0.055);
             resonance[i][j]->GetYaxis()->SetTitleSize(0.055);
             resonance[i][j]->GetXaxis()->SetTitleOffset(0);
             resonance[i][j]->GetYaxis()->SetTitleOffset(0.90);
             if(i<4)resonance[i][j]->GetYaxis()->SetRangeUser(0,33);
             else resonance[i][j]->GetYaxis()->SetRangeUser(0,83);
//             resonance[i][j]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
//             resonance[i][j]->GetYaxis()->SetTitle("#sigma [pb]");
             if(i>=4) resonance[i][j]->GetYaxis()->SetTitleOffset(0.90);
       }
       resonance[i][0]->SetLineColor(kBlue);
       resonance[i][1]->SetLineColor(kRed);
       resonance[i][2]->SetLineColor(kGreen);
       resonance[i][3]->SetLineColor(kViolet);
    }
 
//    cout<<"Mark A"<<endl;
//    cout<<"Mark B"<<endl;

//    c1->cd();
    TPad *mypad[8];
    for(size_t i=0; i<8;i++){
       string tmp_name("c1_");
       tmp_name+=to_string(i+1);
       mypad[i] = (TPad*) c1->GetListOfPrimitives()->FindObject(tmp_name.c_str());
       mypad[i]->cd();
//    resonance_combine[1]->Draw();
       resonance[i][0] -> GetXaxis()->SetLabelSize(0.07);
       resonance[i][0] -> GetYaxis()->SetLabelSize(0.07);
       resonance[i][0] ->GetXaxis()->SetTitleSize(0.055);
       resonance[i][0] ->GetYaxis()->SetTitleSize(0.055);
       resonance[i][1] -> GetXaxis()->SetLabelSize(0.07);
       resonance[i][1] -> GetYaxis()->SetLabelSize(0.07);
       resonance[i][1] ->GetXaxis()->SetTitleSize(0.055);
       resonance[i][1] ->GetYaxis()->SetTitleSize(0.055);
       resonance[i][2] -> GetXaxis()->SetLabelSize(0.07);
       resonance[i][2] -> GetYaxis()->SetLabelSize(0.07);
       resonance[i][2] ->GetXaxis()->SetTitleSize(0.055);
       resonance[i][2] ->GetYaxis()->SetTitleSize(0.055);
       resonance[i][3] -> GetXaxis()->SetLabelSize(0.07);
       resonance[i][3] -> GetYaxis()->SetLabelSize(0.07);
       resonance[i][3] ->GetXaxis()->SetTitleSize(0.055);
       resonance[i][3] ->GetYaxis()->SetTitleSize(0.055);
       resonance[i][0] ->GetXaxis()->SetDecimals(2);
       resonance[i][1] ->GetXaxis()->SetDecimals(2);
       resonance[i][2] ->GetXaxis()->SetDecimals(2);
       resonance[i][3] ->GetXaxis()->SetDecimals(2);
       resonance[i][0]->Draw();
       resonance[i][1]->Draw("same");
       resonance[i][2]->Draw("same");
       resonance[i][3]->Draw("same");
       resonance_combine[i]->GetXaxis()->SetTitleSize(0.055);
       resonance_combine[i]->GetYaxis()->SetTitleSize(0.055);
       resonance_combine[i]->GetXaxis()->SetLabelSize(0.07);
       resonance_combine[i]->GetYaxis()->SetLabelSize(0.07);
       resonance_combine[i]->GetXaxis()->SetDecimals(2);
       resonance_combine[i]->Draw("same");
       TLatex * t_index = new TLatex(0.92,0.94,s_index[i].c_str());
       t_index->SetNDC();
       t_index->SetTextFont(22);
       t_index->SetTextSize(0.075);
       t_index->SetTextAlign(33);
       t_index->Draw();
       c1->Update();
    }
    c1->Print("test.pdf");
/*
    cout<<"Para1 = "<<resonance_1_combine_1->GetParameter(0)<<endl;
    cout<<"Para2 = "<<resonance_1_combine_1->GetParameter(1)<<endl;
*/

    for(size_t r =0; r<8;r++){
       cout<<"Solution["<<r<<"] :Image unit = "<<TComplex::I()<<",Amp at 4.95 = "<<resonance_combine[r]->Eval(4.95)<<endl;
    }
     
    
/*
    resonance_1_combine_1->SetLineColor(2);
    resonance_2_combine_1->SetLineColor(4);
    resonance_1_combine_2->SetLineColor(2);
    resonance_2_combine_2->SetLineColor(4);
    resonance_1_combine_1->SetLineStyle(2);
    resonance_2_combine_1->SetLineStyle(2);
    resonance_1_combine_2->SetLineStyle(9);
    resonance_2_combine_2->SetLineStyle(9);
    resonance_combine_1->SetLineColor(kBlack);
    resonance_combine_2->SetLineColor(kViolet);

    resonance_1_combine_1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_1_combine_1->GetYaxis()->SetTitle("#sigma");
    resonance_1_combine_1->GetYaxis()->SetRangeUser(0,70);
    resonance_2_combine_1->GetYaxis()->SetRangeUser(0,70);
    resonance_1_combine_2->GetYaxis()->SetRangeUser(0,70);
    resonance_2_combine_2->GetYaxis()->SetRangeUser(0,70);
     
    resonance_1_combine_1->Draw();
    resonance_2_combine_1->Draw("same");
    resonance_1_combine_2->Draw("same");
    resonance_2_combine_2->Draw("same");

    resonance_combine_1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    resonance_combine_1->GetYaxis()->SetTitle("#sigma");
    resonance_combine_1->Draw("same");

//    resonance_combine_2->Draw("same");
    c1->Draw();
    c1->Print("Spectrum.pdf");
*/
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
