# Breit-Weigner-MultiSolution
# author: Yu Bai(baiy@seu.edu.cn)
# This program can be implemented to find all possible combination between Breit-Weigners and Polynomail background(Optionnal). You must indicate the input paramter of the resonance. 
0. Make sure python(https://www.python.org/downloads) and numpy(https://scipy.org/install.html) are installed.
1. Write the resonance input file in the form that each line for each resonance. You need to specify the mass(GeV), width(GeV), fR(eV) and phase angle(rad). These variables are separated by command',' like:
     4.0,0.1,1.0,0.785
2. Write the input file for the polynomail background if you need. If you don't create a background input file, the background will be ignored. In the file, each line includes the coefficient for each power term of the polynomial. The coefficients should be put as the increase power. NOTE: EACH COEFFICENTS IS A COMPLEX NUMBER, THE REAL AND IMAGE PART ARE SEPARATED BY COMMAND ','. For example, a background file with the following lines:
     0.31,0.15
    -0.01,0.005
     0.004
    represents a background 0.004*s^2+(-0.01+0.005i)s^2+0.31+0.15i 
3. Specify the input and output file into the python code 'multibw.py' and run as 'python multibw.py'
