import numpy as np
import math
import os
PI=3.1415926535

def Phi23(sqrts,M1=0.000511,M2=0.000511,M3=0.1395706,M4=0.1395706,M5=3.0969): #2->3 phase space, M1+M2->M3+M4+M5
    nn=100000
    E5min=M5
    E5max=sqrts/2.0-((M3+M4)*(M4+M3)-M5*M5)/(2*sqrts)
    step = (E5max-E5min)/(float(nn))
    P23 = 0.0

    for i in range(1,nn+1):
        E5 = M5+float(i)*step
        sigma = sqrts-E5
        tau = sigma*sigma-(E5*E5-M5*M5)
        tau0 = math.fabs((tau-(M3+M4)*(M3+M4))*(tau-(M3-M4)*(M3-M4)))
        E3max = 0.5/tau*(sigma*(tau+M3*M3-M4*M4)+np.sqrt(E5*E5-M5*M5)*np.sqrt(tau0))
        E3min = 0.5/tau*(sigma*(tau+M3*M3-M4*M4)-np.sqrt(E5*E5-M5*M5)*np.sqrt(tau0))
        P23 = P23 + (E3max-E3min)*step

    K1 = np.sqrt((sqrts*sqrts-(M1+M2)*(M1+M2))*(sqrts*sqrts-(M1-M2)*(M1-M2)))/(2*sqrts)
    P23 /= (128*PI*PI*PI*K1*sqrts)
    return P23


def MultiSolution(bkginputfile='bkg.txt',bwinputfile='bw.txt',outputfile='output.txt'):
    #reading Breit-Weigner
    if not os.path.isfile(bwinputfile):
       print('Resonance input file not exists!')
       return

    if not os.path.isfile(bkginputfile):
       print('Background input file not exists, set background to 0')
      
    Mass  = np.genfromtxt(bwinputfile,delimiter=',',usecols = (0))
    Width = np.genfromtxt(bwinputfile,delimiter=',',usecols = (1))
    Amp   = np.genfromtxt(bwinputfile,delimiter=',',usecols = (2))
    Phase = np.genfromtxt(bwinputfile,delimiter=',',usecols = (3))
    PhaseSpace = np.zeros(Mass.size,dtype=float)
    if Mass.size > 1:
       for i in range(0,Mass.size):
           PhaseSpace[i] = Phi23(Mass[i])
    else:
       PhaseSpace = Phi23(Mass)
    Pole = Mass*Mass - 1j*Mass*Width
    BWCoeff = 0.001*np.sqrt(0.001*Width*Amp/PhaseSpace)*Mass*np.exp(1j*Phase)
    #reading background

    BKGCoeff = np.zeros(10,dtype=complex)
    if os.path.isfile(bkginputfile):
        BKGCoeff_R = np.genfromtxt(bkginputfile,delimiter=',',usecols = (0))
        BKGCoeff_I = np.genfromtxt(bkginputfile,delimiter=',',usecols = (1))
        BKGCoeff = np.zeros(BKGCoeff_R.size,dtype=complex)
        if BKGCoeff_R.size>1:
            BKGCoeff = BKGCoeff_R+1j*BKGCoeff_I
        else:
            BKGCoeff[0] = BKGCoeff_R+1j*BKGCoeff_I
    if not os.path.isfile(bkginputfile):
        BKGCoeff = np.zeros(1,dtype=complex)

    #PolyCoeff=[-49.475 + 1.49*1j, 6.172 - 0.0768*1j,-0.256 + 0.0016*1j,0.004]
    PolyBKGExpand=np.polynomial.polynomial.Polynomial(BKGCoeff)
    PolyBKG=np.polynomial.polynomial.Polynomial(BKGCoeff)
    
    PolyPole = np.polynomial.polynomial.Polynomial([0])
    #print Pole
    #print BWCoeff_R
    #print BWCoeff_I
    
    #BWCoeff = BWCoeff_R+1j*BWCoeff_I
    if Pole.size>1:
        for i in range(0,Pole.size):
            PolyPoleProduction = np.polynomial.polynomial.Polynomial([1])
            for j in range(0,Pole.size):
                PolyTemp=np.polynomial.polynomial.Polynomial([-Pole[j],1])
                if i!=j:
                    PolyPoleProduction = PolyPoleProduction*PolyTemp
            PolyPole = PolyPole+(BWCoeff[i])*PolyPoleProduction
            PolyBKGExpand = PolyBKGExpand*np.polynomial.polynomial.Polynomial([-Pole[i],1])
    elif Pole.size == 1:
         PolyTemp = np.polynomial.polynomial.Polynomial([-Pole,1])
         PolyBKGExpand = PolyBKGExpand*PolyTemp
         PolyPole = np.polynomial.polynomial.Polynomial([BWCoeff])
    else:
         PolyPole= np.polynomial.polynomial.Polynomial([0])
    
    PolyCombine=PolyPole+PolyBKGExpand
    print PolyPole
    print PolyBKGExpand
    print PolyCombine
    Solution=PolyCombine.roots()
    

    print ('========================Welcome to Use the Multiple Solution Finder========================')
    print ('===============================Input Resonance Parameter===================================')
    for i in range(0,BWCoeff.size):
        if BWCoeff.size == 1:
             print('Resonance[%i] : Mass = %f, Width = %f, fR = %f, Phi = %f' % (i+1,Mass,Width, Amp,Phase))
        else:
             print('Resonance[%i] : Mass = %f, Width = %f, fR = %f, Phi = %f' % (i+1,Mass[i],Width[i], Amp[i],Phase[i]))
    print ('===============================Input Background Parameter==================================')
    for i in range(0,BKGCoeff.size):
        if i == BKGCoeff.size-1:
             if BKGCoeff[0].imag>0 and BKGCoeff[0].real>0:
                 print('+%f+%f i'%(BKGCoeff[0].real, BKGCoeff[0].imag))
             elif BKGCoeff[0].imag<0 and BKGCoeff[0].real>0:
                 print('+%f-%f i'%(BKGCoeff[0].real, math.fabs(BKGCoeff[0].imag)))
             elif BKGCoeff[0].imag>0 and BKGCoeff[0].real<0:
                 print('-%f+%f i'%(math.fabs(BKGCoeff[0].real), BKGCoeff[0].imag))
             else:
                 print('-%f-%f i'%(math.fabs(BKGCoeff[0].real), math.fabs(BKGCoeff[0].imag)))
        else:
             if BKGCoeff[BKGCoeff.size-1-i].imag>0 and BKGCoeff[BKGCoeff.size-1-i].real>0:
                 print(('(+%f+%f i)s^(%i)+')%(BKGCoeff[BKGCoeff.size-1-i].real, BKGCoeff[BKGCoeff.size-1-i].imag,BKGCoeff.size-1-i))
             elif BKGCoeff[BKGCoeff.size-1-i].imag<0 and BKGCoeff[BKGCoeff.size-1-i].real>0:
                 print('(+%f-%f i)s^(%i)+'%(BKGCoeff[BKGCoeff.size-1-i].real, math.fabs(BKGCoeff[BKGCoeff.size-1-i].imag),BKGCoeff.size-1-i))
             elif BKGCoeff[BKGCoeff.size-1-i].imag>0 and BKGCoeff[BKGCoeff.size-1-i].real<0:
                 print('(-%f+%f i)s^(%i)+'%(math.fabs(BKGCoeff[BKGCoeff.size-1-i].real), BKGCoeff[BKGCoeff.size-1-i].imag,BKGCoeff.size-1-i))
             else:
                 print('(-%f-%f i)s^(%i)+'%(math.fabs(BKGCoeff[BKGCoeff.size-1-i].real), math.fabs(BKGCoeff[BKGCoeff.size-1-i].imag),BKGCoeff.size-1-i))
    print ('====================================== Printing Zero ======================================')
    print Solution
    print ('==================================== Zero Printing Done ===================================')
    
    BWCAeff_all = np.ndarray(shape=(np.power(2,Solution.size),BWCoeff.size),dtype=complex)
    BKGCAeff_all = np.ndarray(shape=(np.power(2,Solution.size),BKGCoeff.size),dtype=complex)
    
    for k in range(0,np.power(2,Solution.size)):
#         print('input BW Coeff:')
#         print BWCoeff
         tmpBWCoeff = BWCoeff
         tmpPolyDeltaBKG = np.polynomial.polynomial.Polynomial([0])
         tmpBKGCoeff = BKGCoeff
         for l in range(Solution.size):
             if((k&np.power(2,l))>>l == 1):
                   Ratio = (Solution[l].conjugate()-Pole)/(Solution[l]-Pole)
                   tmpBWCoeff = tmpBWCoeff*Ratio
                   if BKGCoeff.size>1: 
                       for n in range(0,BKGCoeff.size-1):
    #                       BKGCoeffTemp = np.zeros(n+1,dtype=complex)
                           BKGCoeffTemp = np.zeros(BKGCoeff.size,dtype=complex)
                           if n ==0:
                               BKGCoeffTemp[0] = 1
                           else:
                               for m in range(0,n+1):
                                   BKGCoeffTemp[m]=np.power(Solution[l],n-m)
                           for m in range(n+1,BKGCoeff.size):
                               BKGCoeffTemp[m]=0
    
                           BKGCoeffTemp = BKGCoeffTemp*tmpBKGCoeff[n+1]
                           BKGCoeffTemp = BKGCoeffTemp*(Solution[l]-Solution[l].conjugate())
                           tmpBKGCoeff = tmpBKGCoeff+BKGCoeffTemp
    #         print tmpBWCoeff
    #         print('k = ')
    #         print k
    #         print('l = ')
    #         print l
    #         print('BKG in this stage = ')
    #         print tmpBKGCoeff
    #     print PhaseSpace/(Mass*Mass*Width)
   #      print np.power(tmpBWCoeff,2)
   #      outputBWCoeff = 1000000000.0*(np.power(tmpBWCoeff,2))*PhaseSpace/(Mass*Mass*Width) //CAUSE BUG!
         outputBWCoeff = np.sqrt(1000000000.0*PhaseSpace/Width)*tmpBWCoeff/Mass
         BKGCAeff_all[k]=tmpBKGCoeff
         BWCAeff_all[k]=outputBWCoeff


    print('Printing Output'.center(140,' '))
    for k in range(0,np.power(2,Solution.size)):
        print(('Solution %i'.center(141,'='))%(k+1))
        print('Resoance'.center(20,' ')+'\t|'+'Mass'.center(30,' ')+'\t|'+'Width'.center(30,' ')+'\t|'+'fR'.center(30,' ')+'\t|'+'Phi'.center(30,' '))
        for i in range(0,BWCAeff_all[k].size):
            if BWCoeff.size == 1:
                print(('%i'.center(20,' ')+'\t|'+'%-.4f'.center(26)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24))%(i+1,Mass,Width,BWCAeff_all[k].real*BWCAeff_all[k].real+BWCAeff_all[k].imag*BWCAeff_all[k].imag,np.angle(BWCAeff_all[k])))
            else:
                print(('%i'.center(20,' ')+'\t|'+'%-.4f'.center(26)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24))%(i+1,Mass[i],Width[i],BWCAeff_all[k][i].real*BWCAeff_all[k][i].real+BWCAeff_all[k][i].imag*BWCAeff_all[k][i].imag,np.angle(BWCAeff_all[k][i])))
        print('-'.center(140,'-'))             
        
        print('s power'.center(30,' ')+'\t|'+'Coefficients'.center(60))
        for i in range(0,BKGCAeff_all[k].size):
            if BKGCAeff_all[k][i].imag <0:
               print(('s^%i'.center(30,' ')+'\t|'+'%f%f i'.center(40))%(i,BKGCAeff_all[k][i].real,BKGCAeff_all[k][i].imag))
            else:
               print(('s^%i'.center(30,' ')+'\t|'+'%f+%f i'.center(40))%(i,BKGCAeff_all[k][i].real,BKGCAeff_all[k][i].imag))
        print('='.center(140,'='))
    #     print('BW Coefficients:')
    #     print tmpBWCoeff
    #     print(np.sqrt(outputBWCoeff.real*outputBWCoeff.real+outputBWCoeff.imag*outputBWCoeff.imag),np.angle(outputBWCoeff))
    #     print('Background:')
    #     print tmpBKGCoeff
#     if BWCoeff.size>1:
#         for i in range(0,BWCoeff.size):
#              BWCAeff_all[k][i] = tmpBWCoeff[i]  
#     else:
#         BWCAeff_all[k][0] = tmpBWCoeff 
#
#     for j in range(0,BKGCoeff.size):
#         BKGCAeff_all[k] = (np.poly1d(tmpPolyDeltaBKG).c)
        
#print BKGCAeff_all
#print BWCAeff_all         
#Solution=Poly.roots()
#PolyCoeff_BKG=[775,-48,1]
#w1=0.004
#z1=(1j+0.5)/4.0
#CAeff_all = np.ndarray(shape=(int(np.power(2,Solution.size)),2),dtype=complex)
    fout = open(outputfile,"w")
    print >> fout,'Printing Output'.center(140,' ')
    for k in range(0,np.power(2,Solution.size)):
        print >> fout,(('Solution %i'.center(141,'='))%(k+1))
        print >> fout,('Resoance'.center(20,' ')+'\t|'+'Mass'.center(30,' ')+'\t|'+'Width'.center(30,' ')+'\t|'+'fR'.center(30,' ')+'\t|'+'Phi'.center(30,' '))
        for i in range(0,BWCAeff_all[k].size):
            if BWCoeff.size == 1:
                print >> fout,(('%i'.center(20,' ')+'\t|'+'%-.4f'.center(26)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24))%(i+1,Mass,Width,BWCAeff_all[k].real*BWCAeff_all[k].real+BWCAeff_all[k].imag*BWCAeff_all[k].imag,np.angle(BWCAeff_all[k])))
            else:
                print >> fout,(('%i'.center(20,' ')+'\t|'+'%-.4f'.center(26)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24)+'\t|'+'%-.6f'.center(24))%(i+1,Mass[i],Width[i],BWCAeff_all[k][i].real*BWCAeff_all[k][i].real+BWCAeff_all[k][i].imag*BWCAeff_all[k][i].imag,np.angle(BWCAeff_all[k][i])))
        print >> fout,('-'.center(140,'-'))

        print >> fout,('s power'.center(30,' ')+'\t|'+'Coefficients'.center(60))
        for i in range(0,BKGCAeff_all[k].size):
            if BKGCAeff_all[k][i].imag <0:
               print >> fout,(('s^%i'.center(30,' ')+'\t|'+'%f%f i'.center(40))%(i,BKGCAeff_all[k][i].real,BKGCAeff_all[k][i].imag))
            else:
               print >> fout,(('s^%i'.center(30,' ')+'\t|'+'%f+%f i'.center(40))%(i,BKGCAeff_all[k][i].real,BKGCAeff_all[k][i].imag))
        print >> fout,('='.center(140,'='))
    fout.close

if __name__ == "__main__":
    print ('Finding all solutions....')
    MultiSolution()
