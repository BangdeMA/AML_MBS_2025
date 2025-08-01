//
//  main.cpp
//  Solve SDE by Euler method
//  Author: Shizhao Ma
//  Mail: shizhaoma@sjtu.edu.cn
//  Date: 2025/08/01.

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

using namespace std;
/* function declaration */
void Run(int Num, double beta_L0,double beta_L2, double d_L2);                    // case 1 for more patients
void RunSA(int Num,
           double beta_H0, double beta_H2, double beta_B1, double beta_B2, double beta_L0, double beta_L2,double kappa_H0, double kappa_H2, double kappa_B1, double kappa_B2, double kappa_L1, double d_H1, double d_H2, double d_H3, double d_B1, double d_B2, double d_B3, double d_L1, double d_L2, double K_H1B1, double K_H2B3, double K_B1L2, double K_L2B2); // case 2 for sensitive analysis
void RunSingle();                     // case 0 for one patients
void RunMore(int Num, double beta_L0,double beta_L2, double d_L2);
void RunLC(int Num, double A);
void RunOS(int Num, double beta_L0,double beta_L2, double d_L2, double UST);
void RunTherapy(int Num, double beta_L0,double beta_L2, double d_L2,double UST, double UST2);
void RunStoTherapy(int Num, double beta_L0,double beta_L2, double d_L2);
void RunTherapyNoDeath(int Num, double beta_L0,double beta_L2, double d_L2);
void RunTherapyNoStoStop(int Num, double beta_L0,double beta_L2, double d_L2);
void RunTherapyNoStoStopOS(int Num, double beta_L0,double beta_L2, double d_L2);
void RunCombTherapy(int Num, double beta_L0,double beta_L2, double d_L2);
void RunDeathFunction(int Num, double beta_L0,double beta_L2, double d_L2);

/* global variable */
const double  T = 3*365;
//const double  T = 1*365; //more patients without treatment
//const double  T = 0.5*365; // SA time
const double dt = 0.001;
const int     N = T/dt;
double H1[N+1];
double H2[N+1];
double H3[N+1];
double B1[N+1];
double B2[N+1];
double B3[N+1];
double L1[N+1];
double L2[N+1];
double  R[N+1];
double  P[N+1];
double  sigma[N+1];
double lambda[N+1];
double D1[N+1];
double D2[N+1];


/* New global parameters of Michaelis constants*/
double epsilon_1 = 1.0;
double epsilon_2 = 1.0;
double epsilon_3 = 1.0;
double epsilon_4 = 1.0;

double    K_H1B1 = 5e4;//1e7;//4e5;
double    K_H2B3 = 2.6444e+08;
double    K_B1L2 = 7.5e6;
double    K_L2B2 = 2.9042e10;//1e8;//3e10;//

double theta_H1H1 = 6.1566e-07;
double theta_H1L1 = 5e-3;
double theta_H1H3 = 3.1498e-12;
double theta_B1B1 = 3.8095e-05;
double theta_B2B2 = 2.8154e-09;
double theta_B2L2 = 1e-5;
double theta_L1H1 = 1e-10;
double theta_L1L1 = 2e-7;


/* random number generator */
default_random_engine generator;

/****************************************************************************************************************************************************************************************************/
int main()
{
    //const int   Num = 100;
          int index = 0;
    cout << "Enter a number: " << endl;
    cout << "0 means single patient" << endl;
    cout << "1 means more patients without death and therapy" << endl;
    cout << "2 means sensitive analysis" << endl;
    cout << "3 means calculating LCalpha" << endl;
    cout << "4 means more patients for overall survival" << endl;
    cout << "5 means more patients with stochastic therapy, output r_{start}" << endl;
    cout << "6 means more patients with stochastic therapy but without death" << endl;
    cout << "7 means more patients with stochastic therapy r_{begin} and combination therapy" << endl;
    cout << "8 means more patients with stochastic therapy r_{begin} and r{stop}" << endl;
    cin  >> index;
    cout << "Wait a minute" << endl;
    
    switch (index)
    {
        case 0:{// simulate single patient without death and therapy, print out "patient.txt"===============================================
            RunSingle();
            break;
        }
        case 1:{// simulate more patients without death, Num indicates the patients number, print out "Vpatient-Num.txt"====================
            const int Num = 1000;
            double D[Num][3]={0.0};
            ifstream file;
            file.open("vpatients.txt");
            if(!file)
            {
                cout << "can not open vpatients.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            {
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                RunMore(i,D[i][0],D[i][1],D[i][2]);
            }
            break;
        }
        case 2:{// sensitive analysis, Num indicates the patients number, print out "patient-Num.txt"
            const int Num = 1000;
            
            /* read parameter file of sensitive analysis */
            double D[Num][19+4]={0.0};
            ifstream file;
            file.open("parameter.txt");
            if(!file)
            {
                cout << "can not open parameter.txt\n" << endl;
                return 0;
            }
            for (int i=0; i < Num; i++)
            {
                for (int j=0; j < 19+4; j++)
                {
                    file >> D[i][j];
                }
                RunSA(i,D[i][0],D[i][1],D[i][2],D[i][3],D[i][4],D[i][5],D[i][6],D[i][7],D[i][8],D[i][9],D[i][10],D[i][11],D[i][12],D[i][13],D[i][14],D[i][15],D[i][16],D[i][17],D[i][18],D[i][19],D[i][20],D[i][21],D[i][22]);
            }
            break;
        }
        case 3:{// simulate the relationship between drug dosage D with L2star, print out "LCalpha.txt"
            int Num = 101;
            for (int i=0; i<Num; i++)
            {
                RunLC(i,0.001*i);
            }
            break;
        }
        case 4:{// simulate more patients with therapy for OS, Num indicates the patients number, print out "VpatientTherapy-Num.txt"
            //const int Num = 1000;
            const int Num = 950;
            double D[Num][3]={0.0};
            double UST[1000]={0.0};
            ifstream file,file2;
            file.open("vpatients2.txt");
            file2.open("vtherapy.txt");
            if(!file && !file2 )
            {
                cout << "can not open vpatients.txt or vtherapy.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            //for(int i=998; i<999; i++)
            {
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                file2 >> UST[i];
                RunOS(i,D[i][0],D[i][1],D[i][2],UST[i]); //UST = r{begin}
                //RunTherapy(i,D[i][0],D[i][1],D[i][2]);
            }
            break;
        }
        case 5:{// simulate more patients with stochastic therapy, output r_{begin}
            const int Num = 1000;
            double D[Num][3]={0.0};
            ifstream file;
            file.open("vpatients.txt");
            if(!file)
            {
                cout << "can not open vpatients.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            {
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                RunStoTherapy(i,D[i][0],D[i][1],D[i][2]);
            }
            break;
        }
        case 6:{// simulate more patients with stochastic therapy
            const int Num = 1000;
            double D[Num][3]={0.0};
            ifstream file;
            file.open("vpatients.txt");
            if(!file)
            {
                cout << "can not open vpatients.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            {
                cout<< i << endl;
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                RunTherapyNoDeath(i,D[i][0],D[i][1],D[i][2]);
            }
            break;
        }
        case 7:{// simulate more patients with stochastic therapy
            const int Num = 1000;
            double D[Num][3]={0.0};
            ifstream file;
            file.open("vpatients.txt");
            if(!file)
            {
                cout << "can not open vpatients.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            {
                cout<< i << endl;
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                //RunTherapy(i,D[i][0],D[i][1],D[i][2]);
                //RunTherapyNoStoStop(i,D[i][0],D[i][1],D[i][2]);
                //RunTherapyNoStoStopOS(i,D[i][0],D[i][1]*0.546,D[i][2]);
                //RunTherapyNoStoStopOS(i,D[i][0],D[i][1],D[i][2]);
                RunCombTherapy(i,D[i][0],D[i][1],D[i][2]);
            }
            break;
        }
        case 8:{
            const int Num = 1000;
            double D[Num][3]={0.0};
            double UST[Num][2]={0.0};
            ifstream file,file2;
            file.open("vpatients.txt");
            file2.open("vtherapy2.txt");
            if(!file && !file2 )
            {
                cout << "can not open vpatients.txt or vtherapy.txt\n" << endl;
                return 0;
            }
            for(int i=0; i<Num; i++)
            //for(int i=998; i<999; i++)
            {
                for (int j=0; j < 3; j++)
                {
                    file >> D[i][j];
                }
                for (int j=0; j < 2; j++)
                {
                    file2 >> UST[i][j];
                }
                cout << UST[i][0]<<" " << UST[i][1] << endl;
                RunTherapy(i,D[i][0],D[i][1],D[i][2],UST[i][0],UST[i][1]);
            }
            break;
        }
    }
    return 0;
}

/****************************************************************************************************************************************************************************************************/
/* CASE 2 sub function solving SDE for sensitive analysis*/
void RunSA(int Num,
           double beta_H0, double beta_H2, double beta_B1, double beta_B2, double beta_L0, double beta_L2,double kappa_H0, double kappa_H2, double kappa_B1, double kappa_B2, double kappa_L1, double d_H1, double d_H2, double d_H3, double d_B1, double d_B2, double d_B3, double d_L1, double d_L2, double K_H1B1, double K_H2B3, double K_B1L2, double K_L2B2)
{
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    /*double    mu_p = 0.50;
    double sigma_p = 0.03;
    double death;*/
    
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "patient-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << endl;

    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*d_L2*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
         // death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}

/****************************************************************************************************************************************************************************************************/
/* CASE 1 sub function solving SDE for more patients without death and therapy*/
void RunMore(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
        
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "Vpatient-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << endl;

    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2)*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}
/*****************************************************************************************************************************************************************************************************/
/* CASE 0 sub function solving SDE for single patient*/
void RunSingle()
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    double   beta_L0 = 0.5000;
    double   beta_L2 = 0.1240;//*0.546;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5000;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    double      d_L2 = 0.1330;//0.0872;

    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    
    
    /* parameters of chemotherapy
    double t[2048] = {0.0}; // save the beginning time of chemo
    int    k = 0;*/
    
    /* open files and write into initial value */
    ofstream fout;
    fout.open("patient.txt",ios::out);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << endl;

    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function
         if (abs(R[i] - 0.2)<1e-5)
        {
            t[k] = i;
            k++;
        }
        if (i > t[0] && t[0] != 0)
        {
            if(R[i]>0.001)
            {
                D1[i] = 0.0;
                D2[i] = 1.0;
                //cout <<"D1[i]=" << D1[i] << endl;
            }
            else
            {
                D1[i] = 0.0;
                D2[i] = 0.0;
            }
            //cout << "D=" << D1[i] << " " << "t0= " << t[0] << endl;
        }*/
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*d_L2*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
        
        /* write data into patient.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << endl;
            //cout <<"D2[i]=" << D2[i] << endl;
        }
    }
    fout.close();
}

/* CASE 3 simulate the relationship between drug dosage and L2(T), print out "LCalpha.txt"*/
void RunLC(int Num, double A)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    double   beta_L0 = 0.5000;
    double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    double      d_L2 = 0.1330;

    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    
    
    /* parameters of chemotherapy */
    double t[2048] = {0.0}; // save the beginning time of chemo
    int    k = 0;
    
    /* open files and write into initial value */
    ofstream fout;
    fout.open("LCalpha.txt",ios::app);

    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
         if (abs(R[i] - 0.2)<1e-5)
        {
            t[k] = i;
            k++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]<<endl;
        }
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = A;
            D2[i] = D1[i];
            //cout<<Num<<" "<<R[i]<<i<<" "<<D2[i]<<endl;
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
        }
        //cout<< A << endl;
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D2[i])  + mutation;
        //cout<< i << endl;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
    }
    /* write data into LCalpha.txt */
    fout << Num << " " << A << " " << L2[N] << " " << endl;
    fout.close();
}

/* CASE 4 simulate the virtual patient from CASE 1 for overall survival, with therapy */
void RunOS(int Num, double beta_L0,double beta_L2, double d_L2, double UST)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.50;//0.35;//0.26
    double sigma_p = 0.0675;//0.050;//0.03
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    //double  t[2048] = {0.0}; // save the beginning time of chemo
    double  t[16384] = {0.0};
    //double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    //int    k2 = 0;
    
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "VTherapy-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);

    uniform_real_distribution<double> unistoTherapy(0.05,0.20);
    double UST2 = unistoTherapy(generator);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << D2[0] << " "  << UST << " " << t[0] <<endl;
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        //if (abs(R[i] - UST2)<1e-5)
        if (abs(R[i] - UST2)<1e-3)
        //if (abs(R[i] - 0.30)<1e-5)
        {
            t[k] = i;
            k++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = 0.30;//0.2620911527377836;
            D2[i] = D1[i];
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
        }
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            break;
        }
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << D2[i] << " " << UST << " " << t[0] <<endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}

/* CASE X simulate the virtual patient from CASE 1, with therapy */
void RunTherapy(int Num, double beta_L0,double beta_L2, double d_L2, double UST, double UST2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.50;
    double sigma_p = 0.07;//0.065;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    double  t[2048] = {0.0}; // save the beginning time of chemo
    double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    int    k2 = 0;
    
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "VTherapy-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);

    //uniform_real_distribution<double> unistoTherapy(0.01,0.30);
    //uniform_real_distribution<double> unistoTherapy2(0.0001,0.0010);
    //double UST = unistoTherapy(generator);
    //double UST2 = unistoTherapy2(generator);
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << D2[0] << " "  << UST << " " << UST2 << endl;
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        if (abs(R[i] - UST)<1e-5)
        {
                t[k] = i;
                k++;
                //cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
            //if (abs(R[i] - 0.001)<1e-5 && i > t[0] && t[0] != 0)
            if (abs(R[i] - UST2)<1e-5 && i > t[0] && t[0] != 0)
            {
                t2[k2] = i;
                k2++;
                cout<<Num<<" "<<R[i]<<i<<" "<<k2<<" t[0]"<<t[0]<<" t2[0]"<<t2[0]<<endl;
            }
        }
        if (i > t[0] && t[0] != 0 && i < t2[0])
        {
            D1[i] = 0.25;//0.25; // 0.15 can be used to simulate OS
            D2[i] = D1[i];
            //cout << "D=" << D2[i] << " " << "t0= " << t[0] << endl;
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
        }
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            break;
        }
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << D2[i] << " " << UST << " " << UST2 << endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}

/* CASE X simulate the virtual patient from CASE 1, with therapy */
void RunTherapyNoStoStop(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.50;
    double sigma_p = 0.07;//0.065;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    double  t[2048] = {0.0}; // save the beginning time of chemo
    double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    int    k2 = 0;
    
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "CombTherapy-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);

    uniform_real_distribution<double> unistoTherapy(0.01,0.30);
    double UST = unistoTherapy(generator);
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << D2[0] << " "  << UST << endl;
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        if (abs(R[i] - UST)<1e-5)
        {
                t[k] = i;
                k++;
                //cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        /*if (abs(R[i] - 0.001)<1e-5 && i > t[0] && t[0] != 0)
        {
            t2[k2] = i;
            k2++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k2<<" t[0]"<<t[0]<<" t2[0]"<<t2[0]<<endl;
        }*/
        double beta_L0New = 0.0;
        double beta_L2New = 0.0;
        double beta_B2New = 0.0;
        //if (i > t[0] && t[0] != 0 && i < t2[0])
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = 0.25;
            D2[i] = D1[i];
            beta_L0New = beta_L0*0.50;//0.75;
            beta_L2New = beta_L2*0.50;//0.75;
            beta_B2New = beta_B2*1.00;
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
            beta_L0New = beta_L0;
            beta_L2New = beta_L2;
            beta_B2New = beta_B2;
        }
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2New/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0New/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2New*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            //cout << "D=" << D2[i+1] << " " << "beta_L0New= " << beta_L0New << endl;
            break;
        }
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << D2[i] << " " << UST << endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}

/* CASE X simulate the virtual patient from CASE 1, with therapy */
void RunTherapyNoStoStopOS(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500/1.0;//160;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.5692850210234907;//0.40717185181383;//0.60;//0.50;
    double sigma_p = 0.0789246091076159;//0.041554444168216645;//0.085;//0.065;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    double  t[81920] = {0.0}; // save the beginning time of chemo
    //double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    //int    k2 = 0;
    
    /* open files and write into initial value */
    ofstream fout;
    fout.open("CombTherapyOSB.txt",ios::app);

    uniform_real_distribution<double> unistoTherapy(0.05,0.20);
    double UST = unistoTherapy(generator);
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        if (abs(R[i] - UST)<1e-3)
        {
                t[k] = i;
                k++;
                //cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        /*if (abs(R[i] - 0.001)<1e-5 && i > t[0] && t[0] != 0)
        {
            t2[k2] = i;
            k2++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k2<<" t[0]"<<t[0]<<" t2[0]"<<t2[0]<<endl;
        }*/
        double beta_B2New = 0.0;
        double beta_L1New = 0.0;
        double beta_L2New = 0.0;
        //if (i > t[0] && t[0] != 0 && i < t2[0])
        if (i > t[0] && t[0] != 0)
        {
            // No therapy
            /*D1[i] = 0.00;
            D2[i] = 0.00;//D1[i];
            beta_B2New = beta_B2*1.0;
            beta_L1New = beta_L0*1.0;
            beta_L2New = beta_L2*1.0;*/
            // chemotherapy
            /*D1[i] = 0.3610983400073908;
            D2[i] = D1[i];
            beta_B2New = beta_B2*1.0;
            beta_L1New = beta_L0*1.0;
            beta_L2New = beta_L2*1.0;*/
            // targeted therapy
            D1[i] = 0.0;
            D2[i] = D1[i];
            beta_B2New = beta_B2*1.00;
            beta_L1New = beta_L0*0.10;
            beta_L2New = beta_L2*0.10;
            // chemo + target
            /*D1[i] = 0.24210781936136816/2.0;
            D2[i] = D1[i];
            beta_B2New = beta_B2*0.1;*/
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
            beta_B2New = beta_B2;
            beta_L1New = beta_L0;
            beta_L2New = beta_L2;
        }
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2New/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L1New/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2New*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            //cout << "D=" << D2[i+1] << " " << "beta_L0New= " << beta_L0New << endl;
            fout << Num << " " << (i+1)*dt << " " << B2[i+1] << " " << L1[i+1] << " " << L2[i+1] << " " << R[i+1] << endl;
            break;
        }
        
        /* write data into patient-Num.txt*/
        if(i == N-1)
        {
            fout << Num << " " << (i+1)*dt << " " << B2[i+1] << " " << L1[i+1] << " " << L2[i+1] << " " << R[i+1] << endl;
        }
    }
    fout.close();
}

/* CASE X simulate the virtual patient from CASE 1, with combination therapy */
void RunCombTherapy(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500/1.0;//160;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.5692850210234907;//0.40717185181383;//0.60;//0.50;
    double sigma_p = 0.0789246091076159;//0.041554444168216645;//0.085;//0.065;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    double  t[40960] = {0.0}; // save the beginning time of chemo
    double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    int    k2 = 0;
    
    /* open files and write into initial value */
    ofstream fout;
    fout.open("CombTherapyOSComb.txt",ios::app);

    uniform_real_distribution<double> unistoTherapy(0.05,0.20);
    double UST = unistoTherapy(generator);
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        if (abs(R[i] - UST)<1e-3)
        {
                t[k] = i;
                k++;
                //cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        //if (abs(R[i] - 0.001)<1e-5 && i > t[0] && t[0] != 0)
        if (abs(R[i] - 0.0001)<1e-6 && i > t[0] && t[0] != 0)
        {
            t2[k2] = i;
            k2++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k2<<" t[0]"<<t[0]<<" t2[0]"<<t2[0]<<endl;
        }
        double beta_B2New = 0.0;
        double beta_L1New = 0.0;
        double beta_L2New = 0.0;
        //if (i > t[0] && t[0] != 0 && i < t2[0])
        if (i > t[0] && t[0] != 0)
        {
            if(i < t2[0])
            {
                // chemo + target
                D1[i] = 0.3610983400073908/2.0;
                D2[i] = D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*0.1;
                beta_L2New = beta_L2*0.1;
                // chemo  1st + target 2nd
                /*D1[i] = 0.3610983400073908;
                D2[i] = D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*1.0;
                beta_L2New = beta_L2*1.0;*/
                // target 1st + chemo  2nd
                /*D1[i] = 0.0;
                D2[i] = 0.0;//D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*0.1;
                beta_L2New = beta_L2*0.1;*/
            }
            else
            {
                // chemo + target
                D1[i] = 0.3610983400073908/2.0;
                D2[i] = D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*0.1;
                beta_L2New = beta_L2*0.1;
                // chemo  1st + target 2nd
                /*D1[i] = 0.0;
                D2[i] = 0.0;//D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*0.1;
                beta_L2New = beta_L2*0.1;*/
                // target 1st + chemo  2nd
                /*D1[i] = 0.3610983400073908/2.0;
                D2[i] = D1[i];
                beta_B2New = beta_B2*1.0;
                beta_L1New = beta_L0*1.0;
                beta_L2New = beta_L2*1.0; */           }
            
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
            beta_B2New = beta_B2;
            beta_L1New = beta_L0;
            beta_L2New = beta_L2;
        }
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2New/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L1New/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2New*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            //cout << "D=" << D2[i+1] << " " << "beta_L0New= " << beta_L0New << endl;
            fout << Num << " " << (i+1)*dt << " " << B2[i+1] << " " << L1[i+1] << " " << L2[i+1] << " " << R[i+1] << endl;
            break;
        }
        
        /* write data into patient-Num.txt*/
        if(i == N-1)
        {
            fout << Num << " " << (i+1)*dt << " " << B2[i+1] << " " << L1[i+1] << " " << L2[i+1] << " " << R[i+1] << endl;
        }
    }
    fout.close();
}


/* CASE X simulate the virtual patient from CASE 1, with therapy, without death*/
void RunTherapyNoDeath(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability*/
    double    mu_p = 0.5700135324070382;//0.40717185181383;//0.60;//0.50;
    double sigma_p = 0.07906152604773091;//0.041554444168216645;//0.085;//0.065;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    
    /* parameters of chemotherapy */
    //double  t[2048] = {0.0}; // save the beginning time of chemo
    //double  t[8192] = {0.0};
    double  t[409600] = {0.0};
    //double t2[8192] = {365000*5+1}; // save the beginning time of chemo
    int    k  = 0;
    //int    k2 = 0;
    
    /* open files and write into initial value */
    char fnc[20];
    snprintf(fnc, 20, "VTherapy-%d.txt", Num);
    ofstream fout;
    fout.open(fnc,ios::out);

    uniform_real_distribution<double> unistoTherapy(0.05,0.20);
    //uniform_real_distribution<double> unistoTherapy(0.01,0.30);
    double UST = unistoTherapy(generator);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << D2[0] << " "  << UST << " " << t[0] <<endl;
    
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        //if (abs(R[i] - UST)<1e-6)
        if (abs(R[i] - UST)<1e-3)
        {
            t[k] = i;
            k++;
            //cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        /*if (abs(R[i] - 0.001)<1e-5 && i > t[0] && t[0] != 0)
        {
            t2[k2] = i;
            k2++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k2<<" t[0]"<<t[0]<<" t2[0]"<<t2[0]<<endl;
        }*/
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = 0.3652124412257087;//0.24210781936136816;
            D2[i] = D1[i];
            //cout << "D=" << D2[i] << " " << "t0= " << t[0] << endl;
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
        }
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+ epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
        R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
        P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            break;
        }
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << D2[i] << " " << UST << " " << t[0] <<endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}

/* CASE X more patients with death with death and without therapy */
void Run(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.50;
    double sigma_p = 0.073;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    
    /* parameters of chemotherapy */
    /*double t[N+1] = {0}; // save the beginning time of chemo
    int    k = 0;*/
    
    /* open files and write into initial value */
    char fnc[20];
    //sprintf(fnc,"patient-%d.txt",Num);
    snprintf(fnc, 20, "Vpatient-%d.txt", Num);
    //cout << fnc << endl;
    ofstream fout;
    fout.open(fnc,ios::out);
    //fout.open("Patient.txt",ios::out);
    
    fout << 0 << " " << H1[0] << " " << H2[0] << " " << H3[0] << " " << B1[0] << " " << B2[0] << " " << B3[0] << " " << L1[0] << " "  << L2[0] << " " << R[0] << " " << P[0] << " "  << endl;

    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        /*if (abs(R[i] - 0.2)<1e-7)
        {
            t[k] = i;
            k++;
        }
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = 1.0;
            D2[i] = 1.0;
            //cout << "D=" << D1[i] << " " << "t0= " << t[0] << endl;
        }*/
        
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-0.1*0) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+0.1*0)*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            //cout << "death probability" << P << endl;
            break;
        }
        
        /* write data into patient-Num.txt */
        if(i%100 == 0)
        {
            fout << (i+1)*dt << " " << H1[i+1] << " " << H2[i+1] << " " << H3[i+1] << " " << B1[i+1] << " " << B2[i+1] << " " << B3[i+1] << " "  << L1[i+1] << " "  << L2[i+1] << " " << R[i+1] << " " << P[i+1] << " " << endl;
            //cout <<"D1[i]=" << D1[i] << endl;
        }
    }
    fout.close();
}


/* CASE X simulation more patients with stochastic therapy */
 void RunStoTherapy(int Num, double beta_L0,double beta_L2, double d_L2)
{
    /* parameters of ODE */
    double   beta_H0 = 0.5;//0.0455;
    double   beta_H2 = 0.2450;
    double   beta_B1 = 0.4800;
    double   beta_B2 = 0.1852;
    //double   beta_L0 = 0.5000;
    //double   beta_L2 = 0.1240;

    double  kappa_H0 = 0.0455;
    double  kappa_H2 = 0.5;
    double  kappa_B1 = 0.1000;
    double  kappa_B2 = 0.1000;
    double  kappa_L1 = 0.2500;

    double      d_H1 = 0.0005*0;
    double      d_H2 = 0.0213*0;
    double      d_H3 = 0.0080;
    double      d_B1 = 0.0600;
    double      d_B2 = 0.0625;
    double      d_B3 = 0.1000;
    double      d_L1 = 0.0688;
    //double      d_L2 = 0.1330;
    
    /* parameters of sigma(B2) */
    double   delta_L = 4.5e-3*0.005;
    double epsilon_0 = 7.0;
    double      m_B2 = 5.0;
    double  theta_B2 = 1e8;
    
    /* parameters of death probability */
    double    mu_p = 0.50;
    double sigma_p = 0.05;//0.084375;
    double death;
    
    /* Calculate initial value */
    double   B1star = (beta_B1/(kappa_B1+d_B1) - 1)/theta_B1B1;
    double    first = B1star*kappa_B1*theta_B2B2 + beta_B2 - d_B2 - kappa_B2;
    double   B2star = (sqrt(pow(first,2) + 4*B1star*kappa_B1*theta_B2B2*(d_B2+kappa_B2)) + first)/2/theta_B2B2/(d_B2+kappa_B2);
    double   B3star = kappa_B2*B2star/d_B3;
    double     B_H1 = beta_H0*(1+B1star/(B1star+K_H1B1));
    double     K_H2 = kappa_H2*(1+B3star/(B3star+K_H2B3));
    double    FstH1 = -K_H2*(pow(B_H1,2)*theta_H1H3+d_H3*theta_H1H1*(B_H1-2*kappa_H0))+beta_H2*d_H3*theta_H1H1*(B_H1-2*kappa_H0);
    double    SecH1 = pow((B_H1*K_H2*theta_H1H3+d_H3*theta_H1H1*(beta_H2-K_H2)),2) - 4*d_H3*K_H2*theta_H1H1*theta_H1H3*(B_H1-kappa_H0)*(beta_H2-K_H2);
    double   H1star = (FstH1 - B_H1*sqrt(SecH1))/(2*d_H3*kappa_H0*theta_H1H1*theta_H1H1*(beta_H2-K_H2));
    double    FstH2 = B_H1*K_H2*theta_H1H3 + beta_H2*d_H3*theta_H1H1 - d_H3*K_H2*theta_H1H1;
    double   H2star = (-FstH2 -sqrt(SecH1))/(2*K_H2*theta_H1H1*theta_H1H3*(beta_H2-K_H2));
    double   H3star = K_H2*H2star/d_H3;
    H1[0] = H1star;
    H2[0] = H2star;
    H3[0] = H3star;
    B1[0] = B1star;
    B2[0] = B2star;
    B3[0] = B3star;
    L1[0] = 1e2;
    L2[0] = 1e2;
     R[0] = L2[0]/(H3[0]+L2[0]);
     P[0] = 0.0;
    D1[0] = 0.0;
    D2[0] = 0.0;
    //double D1[N+1] = {0.0};
    //double D2[N+1] = {0.0};
    
    /* parameters of chemotherapy */
    double t[2048] = {0.0}; // save the beginning time of chemo
    int    k = 0;
    
    /* open files and write into initial value */
    ofstream fout;
    fout.open("UniStoTherapy.txt",ios::app);

    uniform_real_distribution<double> unistoTherapy(0.05,0.20);
    double UST = unistoTherapy(generator);
    /* Euler mathod */
    for (int i = 0; i < N; i++)
    {
        //--------------------------------------------------------------------------
        //                                (B2(i)/theta_B2)^m_B2
        // sigma = delta_L*[epsilon_0 + -------------------------------]
        //                               1   +   (B2(i)/theta_B2)^m_B2
        //--------------------------------------------------------------------------
       lambda[i] = R[i]*sigma[i]*H1[i];
        sigma[i] = delta_L*(epsilon_0 + pow(B2[i]/theta_B2,m_B2))/(1+pow(B2[i]/theta_B2,m_B2));
        
        /* Dosage function */
        if (abs(R[i] - UST)<1e-6)
        {
            t[k] = i;
            k++;
            cout<<Num<<" "<<R[i]<<i<<" "<<k<<" "<<t[0]*dt<<endl;
        }
        if (i > t[0] && t[0] != 0)
        {
            D1[i] = 0.0;
            D2[i] = D1[i];
        }
        else
        {
            D1[i] = 0.0;
            D2[i] = 0.0;
        }
        /* Generate poisson distribution random number */
        poisson_distribution<int> distribution(lambda[i]*dt);
        uniform_real_distribution<double> unirand(0.0,1.0);
        int mutation = distribution(generator);
        
        /* Euler method */
        H1[i+1] = H1[i] + dt*H1[i]*(beta_H0/(1+theta_H1H1*H1[i]+theta_H1L1*L1[i])*(1 + epsilon_1*B1[i]/(B1[i]+K_H1B1))) -  dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) - dt*d_H1*H1[i] - mutation;
        //--------------------------------------------------------------------------
        //  dH1              beta_H0                             B1                      1
        //  --- = [-----------------------------(1 + epsilon_1---------) - kappa_H1---------------]*H1 - dN(t)
        //  dt     1+theta_H1H1*H1+theta_H1L1*L1              B1+K_H1B1            1+theta_H1H3*H3
        //--------------------------------------------------------------------------
        
        H2[i+1] = H2[i] + dt*H1[i]*kappa_H0/(1+theta_H1H3*H3[i]) + dt*H2[i]*beta_H2 - dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*d_H2*H2[i];
        //--------------------------------------------------------------------------
        //  dH2                  1                                                       B3
        //  --- = kappa_H1---------------*H1 + beta_H2*H2 - kappa_H2*H2*(1 + epsilon_2---------)
        //  dt            1+theta_H1H3*H3                                             B3+K_H2B3
        //--------------------------------------------------------------------------
        
        H3[i+1] = H3[i] + dt*H2[i]*kappa_H2*(1+epsilon_2*B3[i]/(B3[i]+K_H2B3)) - dt*H3[i]*d_H3;
        //--------------------------------------------------------------------------
        //  dH3                                B3
        //  --- = kappa_H2*H2*(1 + epsilon_2---------) - d_H3*H3
        //  dt                              B3+K_H2B3
        //--------------------------------------------------------------------------
        
        B1[i+1] = B1[i] + dt*B1[i]*beta_B1/(1+theta_B1B1*B1[i]) - dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) - dt*d_B1*B1[i];
        //--------------------------------------------------------------------------
        //  dB1                  1                                       L2
        //  --- = beta_B1*B1--------------- - kappa_B1*B1*(1 + epsilon_3---------) - d_B1*B1
        //  dt              1+theta_B1B1*B1                             L2+K_B1L2
        //--------------------------------------------------------------------------
        
        B2[i+1] = B2[i] + dt*B1[i]*kappa_B1*(1+epsilon_3*L2[i]/(K_B1L2+L2[i])) + dt*B2[i]*beta_B2/(1+theta_B2B2*B2[i]) - dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B2*B2[i];
        //--------------------------------------------------------------------------
        //   dB2                               L2                       1                              1
        //   --- = kappa_B1*B1*(1 + epsilon_3---------) + beta_B2*B2--------------- - kappa_B2*B2*--------------- - d_B2*B2
        //   dt                              L2+K_B1L2              1+theta_B2B2*B2               1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        B3[i+1] = B3[i] + dt*B2[i]*kappa_B2/(1+theta_B2L2*L2[i]) - dt*d_B3*B3[i];
        //--------------------------------------------------------------------------
        //  dB3                   1
        //  --- = kappa_B2*B2*--------------- - d_B3*B3
        //  dt                1+theta_B2L2*L2
        //--------------------------------------------------------------------------
        
        L1[i+1] = L1[i] + dt*L1[i]*(beta_L0/(1+theta_L1H1*H1[i]+theta_L1L1*L1[i])-kappa_L1-d_L1-0.1*0-D1[i]) + mutation;
        L2[i+1] = L2[i] + dt*L1[i]*kappa_L1 + dt*L2[i]*beta_L2*(1+epsilon_4*B2[i]/(K_L2B2+B2[i])) - dt*(d_L2+0.1*0+D2[i])*L2[i];
        //--------------------------------------------------------------------------
        //  dL1                    1
        //  --- = (beta_L1----------------------------- - kappa_L1 - d_L1)*L1 + dN(t) - D_1(t)*L1
        //  dt            1+theta_L1H1*H1+theta_L1L1*L1
        //  dL2                                             B2
        //  --- = kappa_L1*L_1 + beta_L2*L2*(1+epsilon_4*---------) - d_L2*L2 - D_2(t)*L1
        //  dt                                           K_L2B2+B2
        //--------------------------------------------------------------------------
         R[i+1] = L2[i+1]/(H3[i+1]+L2[i+1]);
         P[i+1] = unirand(generator);
          death = (1.0 + exp(-(1 - mu_p)/sigma_p))/(1.0 + exp(-(R[i+1] - mu_p)/sigma_p));
        if( P[i+1] < death*dt)
        {
            cout << Num <<": " << "death probability=" << P[i+1] << " " << death*dt << " R=" <<R[i+1] << " "<< endl;
            /* write data into UniStoTherapy.txt */
            fout << Num << " " << UST << " " << (i+1)*dt << endl;
            break;
        }
    }
    fout.close();
}
