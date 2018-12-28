#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TComplex.h"
#include "Draw.h"
#include "TMath.h"
#include "TAxis.h"
#define MJPsi 3.097
#define MPi 0.139
#define MELE 0.000511
#define M0 3.0
#define M 4.0
#define Gamma 0.1
#define SIGMA 389000000.0

double Phi23(double sqrts,double M1,double M2,double M3,double M4,double M5){
      double nn = 1000;
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

void Draw(){
     TComplex z0 = 0.15;
     TComplex z1 = (0.5+TComplex::I())/40.0;
     TComplex zero = 15.414427718834121-1.2510703163162662*TComplex::I();
     TComplex Pole = M*M-M*Gamma*TComplex::I();
     
     TF1 * f_BW1 = new TF1("BW1",amp2TF_BW,3.4,5.0,2);
     TF1 * f_BW2 = new TF1("BW2",amp2TF_BW,3.4,5.0,2);
     TF1 * f_BKG1 = new TF1("BKG1",amp2TF_BKG1,3.4,5,2);
     TF1 * f_BKG2 = new TF1("BKG2",amp2TF_BKG2,3.4,5,4);
     TF1 * f_Comb1 = new TF1("Combiniation_I",amp2TF_Comb1,3.4,5,4);
     TF1 * f_Comb2 = new TF1("Combiniation_II",amp2TF_Comb2,3.4,5,6);

     f_BW1->SetParameter(0,z1.Re()); f_BW1->SetParameter(1,z1.Im());
     TComplex z1_prime = z1*(TComplex::Conjugate(zero)-Pole)/(zero-Pole);
     f_BW2->SetParameter(0,z1_prime.Re()); f_BW2->SetParameter(1,z1_prime.Im());
     f_BKG1->SetParameter(0,z0.Re());  f_BKG1->SetParameter(1,z0.Im());
     f_BKG2->SetParameter(0,z0.Re());  f_BKG2->SetParameter(1,z0.Im());  f_BKG2->SetParameter(2,zero.Re()); f_BKG2->SetParameter(3,zero.Im());
     f_Comb1->SetParameter(0,z0.Re()); f_Comb1->SetParameter(1,z0.Im()); f_Comb1->SetParameter(2,z1.Re()); f_Comb1->SetParameter(3,z1.Im());
     f_Comb2->SetParameter(0,z0.Re()); f_Comb2->SetParameter(1,z0.Im()); f_Comb2->SetParameter(2,zero.Re()); f_Comb2->SetParameter(3,zero.Im()); f_Comb2->SetParameter(4,z1_prime.Re()); f_Comb2->SetParameter(5,z1_prime.Im());

     TCanvas *c1 = new TCanvas("c1","Combination of two Breit-Weigner functions",150,150,600,600);
     TCanvas *c2 = new TCanvas("c2","Combination of two Breit-Weigner functions",150,150,600,600);


     f_BW1->GetXaxis()->SetLabelSize(0.07);
     f_BW1->GetYaxis()->SetLabelSize(0.07);
     f_BW1->GetXaxis()->SetTitleSize(0.07);
     f_BW1->GetYaxis()->SetTitleSize(0.07);
     f_BW1->GetXaxis()->SetTitleOffset(0);
     f_BW1->GetYaxis()->SetTitleOffset(0.85);
     f_BW1->GetXaxis()->SetDecimals(2);
     f_BW1->GetXaxis()->CenterTitle(kTRUE);
     f_BW2->GetXaxis()->SetLabelSize(0.07);
     f_BW2->GetYaxis()->SetLabelSize(0.07);
     f_BW2->GetXaxis()->SetTitleSize(0.07);
     f_BW2->GetYaxis()->SetTitleSize(0.07);
     f_BW2->GetXaxis()->SetTitleOffset(0);
     f_BW2->GetYaxis()->SetTitleOffset(0.85);
     f_BW2->GetXaxis()->SetDecimals(2);
     f_BW2->GetXaxis()->CenterTitle(kTRUE);
     f_BKG1->GetXaxis()->SetLabelSize(0.07);
     f_BKG1->GetYaxis()->SetLabelSize(0.07);
     f_BKG1->GetXaxis()->SetTitleSize(0.07);
     f_BKG1->GetYaxis()->SetTitleSize(0.07);
     f_BKG1->GetXaxis()->SetTitleOffset(0);
     f_BKG1->GetYaxis()->SetTitleOffset(0.85);
     f_BKG1->GetXaxis()->SetDecimals(2);
     f_BKG1->GetXaxis()->CenterTitle(kTRUE);
     f_BKG2->GetXaxis()->SetLabelSize(0.07);
     f_BKG2->GetYaxis()->SetLabelSize(0.07);
     f_BKG2->GetXaxis()->SetTitleSize(0.07);
     f_BKG2->GetYaxis()->SetTitleSize(0.07);
     f_BKG2->GetXaxis()->SetTitleOffset(0);
     f_BKG2->GetYaxis()->SetTitleOffset(0.85);
     f_BKG2->GetXaxis()->SetDecimals(2);
     f_BKG2->GetXaxis()->CenterTitle(kTRUE);
     f_Comb1->GetXaxis()->SetLabelSize(0.07);
     f_Comb1->GetYaxis()->SetLabelSize(0.07);
     f_Comb1->GetXaxis()->SetTitleSize(0.07);
     f_Comb1->GetYaxis()->SetTitleSize(0.07);
     f_Comb1->GetXaxis()->SetTitleOffset(0);
     f_Comb1->GetYaxis()->SetTitleOffset(0.85);
     f_Comb1->GetXaxis()->SetDecimals(2);
     f_Comb1->GetXaxis()->CenterTitle(kTRUE);
     f_Comb2->GetXaxis()->SetLabelSize(0.07);
     f_Comb2->GetYaxis()->SetLabelSize(0.07);
     f_Comb2->GetXaxis()->SetTitleSize(0.07);
     f_Comb2->GetYaxis()->SetTitleSize(0.07);
     f_Comb2->GetXaxis()->SetTitleOffset(0);
     f_Comb2->GetYaxis()->SetTitleOffset(0.85);
     f_Comb2->GetXaxis()->SetDecimals(2);
     f_Comb2->GetXaxis()->CenterTitle(kTRUE);


     f_BW1->SetLineColor(kRed);
     f_BKG1->SetLineColor(kBlue);
     f_Comb1->SetLineColor(kBlack);
     f_BW1->SetLineStyle(2);
     f_BKG1->SetLineStyle(2);

     f_BW1->GetYaxis()->SetRangeUser(0,100);
     f_BKG1->GetYaxis()->SetRangeUser(0,100);
     f_Comb1->GetYaxis()->SetRangeUser(0,100);

     c1->cd();
     f_BW1->SetNpx(1000);
     f_BKG1->SetNpx(1000);
     f_Comb1->SetNpx(1000); 
     f_BW1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_BW1->GetYaxis()->SetTitle("#sigma [pb]");   
     f_BW1->GetXaxis()->CenterTitle(kTRUE);
     f_BW1->GetYaxis()->CenterTitle(kTRUE);
     f_BW1->GetYaxis()->SetRangeUser(0,20);
     f_BW1->Draw();
     f_BKG1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_BKG1->GetYaxis()->SetTitle("#sigma [pb]"); 
     f_BKG1->GetXaxis()->CenterTitle(kTRUE);
     f_BKG1->GetYaxis()->CenterTitle(kTRUE);
     f_BKG1->GetYaxis()->SetRangeUser(0,20);
     f_BKG1->Draw("same");
     f_Comb1->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_Comb1->GetYaxis()->SetTitle("#sigma [pb]");
     f_Comb1->GetXaxis()->CenterTitle(kTRUE);
     f_Comb1->GetYaxis()->CenterTitle(kTRUE);
     f_Comb1->GetYaxis()->SetRangeUser(0,20);
     f_Comb1->Draw("same");
     TLatex * t_index_1 = new TLatex(0.82,0.84,"(I)");
     t_index_1->SetNDC();
     t_index_1->SetTextFont(22);
     t_index_1->SetTextSize(0.075);
     t_index_1->SetTextAlign(33);
//     t_index->Draw();
     c1->Print("ExpBkg_1.pdf");

     c2->cd();
     f_BW2->SetLineColor(kRed);
     f_BKG2->SetLineColor(kBlue);
     f_Comb2->SetLineColor(kBlack);
     f_BW2->SetLineStyle(2);
     f_BKG2->SetLineStyle(2);

     f_BW2->GetYaxis()->SetRangeUser(0,100);
     f_BKG2->GetYaxis()->SetRangeUser(0,100);
     f_Comb2->GetYaxis()->SetRangeUser(0,100);
  
     f_BW2->SetNpx(1000);
     f_BKG2->SetNpx(1000);
     f_Comb2->SetNpx(1000); 
     f_BW2->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_BW2->GetYaxis()->SetTitle("#sigma [pb]");   
     f_BW2->GetXaxis()->CenterTitle(kTRUE);
     f_BW2->GetYaxis()->CenterTitle(kTRUE);
     f_BW2->GetYaxis()->SetRangeUser(0,20);
     f_BW2->Draw();
     f_BKG2->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_BKG2->GetYaxis()->SetTitle("#sigma [pb]"); 
     f_BKG2->GetXaxis()->CenterTitle(kTRUE);
     f_BKG2->GetYaxis()->CenterTitle(kTRUE);
     f_BKG2->GetYaxis()->SetRangeUser(0,20);
     f_BKG2->Draw("same");
     f_Comb2->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
     f_Comb2->GetYaxis()->SetTitle("#sigma [pb]");
     f_Comb2->GetXaxis()->CenterTitle(kTRUE);
     f_Comb2->GetYaxis()->CenterTitle(kTRUE);
     f_Comb2->GetYaxis()->SetRangeUser(0,20);
     f_Comb2->Draw("same");

     TLatex * t_index_2 = new TLatex(0.82,0.84,"(I)");
     t_index_2->SetNDC();
     t_index_2->SetTextFont(22);
     t_index_2->SetTextSize(0.075);
     t_index_2->SetTextAlign(33);
//     t_index->Draw();
     c2->Print("ExpBkg_2.pdf");    
}

TComplex BW(double *x){
    TComplex BW = x[0]-M*M+TComplex::I()*M*Gamma;
    TComplex Coeff= x[1]+TComplex::I()*x[2];
    return Coeff/BW;
}

TComplex BKG1(double *x){
    double bkg = TMath::Exp(-x[0]/(M0*M0));
    TComplex Coeff = x[1]+TComplex::I()*x[2];
    return bkg*Coeff;
}

TComplex BKG2(double *x){
    double bkg = TMath::Exp(-x[0]/(M0*M0));
    TComplex Coeff = x[1]+TComplex::I()*x[2];
    return bkg*Coeff+Coeff*(TMath::Exp(-x[0]/(M0*M0))-TComplex::Exp(-(x[3]+TComplex::I()*x[4])/(M0*M0)))*(2.0)*TComplex::I()*x[4]/(x[0]-x[3]-TComplex::I()*x[4]);
}

double amp2TF_BW(double *x, double *par){ //x[0] is sqrt{s}, par[0] par[1] are the coefficients' real and image part
    double y[3]={x[0]*x[0],par[0],par[1]};
//    cout<<"s = "<<x[0]<<", bw amplitude = "<<BW(y).Rho2()<<endl;
    return 12*SIGMA*TMath::Pi()*Phi23PiPiJPsi(x[0])*BW(y).Rho2()/(x[0]*x[0]);
}

double amp2TF_BKG1(double *x, double *par){ //x[0] is sqrt{s}, par[0] par[1] are the coefficients' real and image part
    double y[3]={x[0]*x[0],par[0],par[1]};
    return 12*SIGMA*TMath::Pi()*Phi23PiPiJPsi(x[0])*BKG1(y).Rho2()/(x[0]*x[0]);
}
double amp2TF_BKG2(double *x, double *par){ //x[0] is sqrt{s}, par[0] par[1] are the coefficients' real and image part, par[2] and par[3] are zero's real and image part
    double y[5]={x[0]*x[0],par[0],par[1],par[2],par[3]};
    return 12*SIGMA*TMath::Pi()*Phi23PiPiJPsi(x[0])*BKG2(y).Rho2()/(x[0]*x[0]);
}


double amp2TF_Comb1(double *x, double *par){
    double v_bkg[3]  = {x[0]*x[0],par[0],par[1]}; //par[0][1] are z0's real and image, par[2][3] are z1's real and image
    double v_bw[3] = {x[0]*x[0],par[2],par[3]}; //par[0][1] are z0's real and image, par[2][3] are z1's real and image 
//    cout<<"par[0] = "<<par[0]<<", par[1]"<<par[1]<<", par[2] = "<<par[2]<<", par[3] = "<<par[3]<<endl;
//    cout<<"sqrt{s} = "<<x[0]<<", bw amplitude = "<<BW(v_bw).Rho2()<<", bkg amplitude = "<<BKG1(v_bkg).Rho2()<<endl;
    return 12*SIGMA*TMath::Pi()*Phi23PiPiJPsi(x[0])*(BW(v_bw)+BKG1(v_bkg)).Rho2()/(x[0]*x[0]);
}

double amp2TF_Comb2(double *x, double *par){
   double v_bkg[5] = {x[0]*x[0],par[0],par[1],par[2],par[3]};//par[0][1] are z0's real and image, par[2][3] are zero's real and image
   double v_bw[3] = {x[0]*x[0],par[4],par[5]}; //par[4][5] are z1's real and image part

   return 12*SIGMA*TMath::Pi()*Phi23PiPiJPsi(x[0])*(BW(v_bw)+BKG2(v_bkg)).Rho2()/(x[0]*x[0]);
}
