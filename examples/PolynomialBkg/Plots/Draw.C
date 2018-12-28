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
#define NRES 1
#define NCOMB 8
#define KR 389000000.0
#define MJPsi 3.0969
#define MPi 0.1395706
#define MELE 0.000511

double Mass = 4.0;
double Width = 0.1;
TComplex w2 = TComplex(0.0004,0);
TComplex Pole= Mass*Mass-TComplex::I()*Mass*Width;


TComplex z1[NCOMB] = {0.1*TComplex(0.125,0.25),0.1*TComplex(0.142075363178,-0.227036458143),0.1*TComplex(-0.291249292774,-0.016199639553),0.1*TComplex(0.132140829691,0.246299446179),0.1*TComplex(0.176812567818,-1.08645327767),0.1*TComplex(-0.976873728988,0.397738312753),0.1*TComplex(0.89285153752,0.722812516039),0.1*TComplex(0.14548,-1.091088)};
TComplex w0[NCOMB] = {0.1*TComplex(3.1,0),0.1*TComplex(1.47453349,-2.73741274),0.1*TComplex(1.52684305,2.67757496),0.1*TComplex(3.09957224,-0.05983778),0.1*TComplex(3.09681234,0.16223778),0.1*TComplex(1.61531027,-2.57517496),0.1*TComplex(1.38202636,2.83981274),0.1*TComplex(3.09872,0.1024)};
TComplex w1[NCOMB] = {0.1*TComplex(-0.192,0),0.1*TComplex(-0.192,0.11403391),0.1*TComplex(-0.192,-0.11218403),0.1*TComplex(-0.192,0.00184988),0.1*TComplex(-0.192,-0.00504988),0.1*TComplex(-0.192,0.10898403),0.1*TComplex(-0.192,-0.11723391),0.1*TComplex(-0.192,-0.0032)};

void Draw(){     
  
      TCanvas *c1 = new TCanvas("c1","Combination of two Breit-Weigner functions",150,10,1500,1000);
      c1->Divide(4,2,0.0,-0.02);
      
      TF1 * bw[8];
      TF1 * combine[8];
      TF1 * bkg[8];

      string s_index[8]={"(I)","(II)","(III)","(IV)","(V)","(VI)","(VII)","(VIII)"};

      for(size_t i=0;i<8;i++){
          cout<<"fR["<<i<<"] = "<<1.0E9*z1[i].Rho2()*Phi23PiPiJPsi(4.0)/(4.0*4.0*0.1)<<endl;
      }

      TPad *mypad[8];
      for(size_t i=0; i<8;i++){
         string tmp_name("c1_");
         tmp_name+=to_string(i+1);
         mypad[i] = (TPad*) c1->GetListOfPrimitives()->FindObject(tmp_name.c_str());
         mypad[i]->cd();
         string bw_name("bw");
         string bkg_name("bkg");
         string combine_name("combine");
//         cout<<"Mark A"<<endl;
         bw[i] = new TF1((bw_name+to_string(i)).c_str(),amp2_BW,3.51,5.1,1);
         bkg[i] = new TF1((bkg_name+to_string(i)).c_str(),amp2_BKG,3.51,5.1,1);
         combine[i] = new TF1((combine_name+to_string(i)).c_str(),amp2_COMB,3.51,5.1,1);
         bw[i]->SetNpx(1000);
         bkg[i]->SetNpx(1000);
         combine[i]->SetNpx(1000);
//         cout<<"Mark B"<<endl;
         bw[i]->SetParameter(0,double(i));
//         cout<<"Mark B1"<<endl;
         bkg[i]->SetParameter(0,double(i));
//         cout<<"Mark B2"<<endl;
         combine[i]->SetParameter(0,double(i));
//         cout<<"Mark B3"<<endl;
//         cout<<"Mark C"<<endl;   
         bw[i]->SetLineColor(kRed);
         bkg[i]->SetLineColor(kBlue);
         combine[i]->SetLineColor(kBlack);
         if(i<=3){
           bw[i]->GetYaxis()->SetRangeUser(0,44);
           bkg[i]->GetYaxis()->SetRangeUser(0,44);
           combine[i]->GetYaxis()->SetRangeUser(0,44);
         }else{
           bw[i]->GetYaxis()->SetRangeUser(0,95);
           bkg[i]->GetYaxis()->SetRangeUser(0,95);
           combine[i]->GetYaxis()->SetRangeUser(0,95);
         }
/*
         bw[i]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
         bw[i]->GetYaxis()->SetTitle("#sigma [pb]");
         bkg[i]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
         bkg[i]->GetYaxis()->SetTitle("#sigma [pb]");
         combine[i]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
         combine[i]->GetYaxis()->SetTitle("#sigma [pb]");
*/
         bw[i]->SetLineStyle(2);
         bkg[i]->SetLineStyle(2);
         bw[i]->GetXaxis()->SetLabelSize(0.07);
         bw[i]->GetYaxis()->SetLabelSize(0.07);
         bw[i]->GetXaxis()->SetTitleSize(0.055);
         bw[i]->GetYaxis()->SetTitleSize(0.055);
         bw[i]->GetXaxis()->SetTitleOffset(0.90);
         bw[i]->GetYaxis()->SetTitleOffset(-0.6);
         bkg[i]->GetXaxis()->SetLabelSize(0.07);
         bkg[i]->GetYaxis()->SetLabelSize(0.07);
         bkg[i]->GetXaxis()->SetTitleOffset(0.90);
         bkg[i]->GetYaxis()->SetTitleOffset(-0.6);
         bkg[i]->GetXaxis()->SetTitleSize(0.055);
         bkg[i]->GetYaxis()->SetTitleSize(0.055);
         combine[i]->GetXaxis()->SetLabelSize(0.07);
         combine[i]->GetYaxis()->SetLabelSize(0.07);
         combine[i]->GetXaxis()->SetTitleSize(0.055);
         combine[i]->GetYaxis()->SetTitleSize(0.055);
         combine[i]->GetXaxis()->SetTitleOffset(0.90);
         combine[i]->GetYaxis()->SetTitleOffset(-0.6);
         bw[i]->GetXaxis()->SetDecimals(2);
         bkg[i]->GetXaxis()->SetDecimals(2);
         combine[i]->GetXaxis()->SetDecimals(2);
         TLatex * t_index = new TLatex(0.92,0.94,s_index[i].c_str());
         t_index->SetNDC();
         t_index->SetTextFont(22);
         t_index->SetTextSize(0.075);
         t_index->SetTextAlign(33);
         mypad[i]->Update();
         cout<<"Solution["<<i<<"] "<<"at 4 GeV, bw = "<<bw[i]->Eval(4.0)<<", bkg = "<<bkg[i]->Eval(4.0)<<", combined = "<<combine[i]->Eval(4.0)<<endl;
         bw[i]->Draw();
//         cout<<"Mark D"<<endl;
         bkg[i]->Draw("same");
         combine[i]->Draw("same");
         t_index->Draw();
         c1->Update();
      }
      c1->Print("test.pdf");
}

double amp2_COMB(double *x,double *par){
    int index = int(par[0]);
    return KR*(w2*x[0]*x[0]*x[0]*x[0]+w1[index]*x[0]*x[0]+w0[index]+z1[index]/(x[0]*x[0]-Pole)).Rho2()*12*TMath::Pi()*Phi23PiPiJPsi(x[0])/(x[0]*x[0]);
}

double amp2_BW(double *x,double *par){
    int index = int(par[0]);
    return KR*(z1[index]/(x[0]*x[0]-Pole)).Rho2()*12*TMath::Pi()*Phi23PiPiJPsi(x[0])/(x[0]*x[0]);
}

double amp2_BKG(double *x,double *par){
    int index = int(par[0]);
    return KR*(w2*x[0]*x[0]*x[0]*x[0]+w1[index]*x[0]*x[0]+w0[index]).Rho2()*12*TMath::Pi()*Phi23PiPiJPsi(x[0])/(x[0]*x[0]);
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
