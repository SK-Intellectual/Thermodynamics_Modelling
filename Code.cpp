#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

long double V0(long double a){
    // input in celsius
    return (1+18.1597*pow(10,-3)*a)/(0.9998+(18.2249*pow(10,-3)*a) -(7.9222*pow(10,-6)*a*a) -(55.4485*(pow(10,-9))*a*a*a)+(149.7562*pow(10,-12)*pow(a,4))-(393.2952*pow(10,-15)*pow(a,5)));
}
long double B(long double b){
    // input in celsius
    return ((19654.32) + (147.037*b) - (2.2155*pow(b,2)) + (1.0478*pow(10,-2)*pow(b,3)) - (2.2789*pow(10,-5)*pow(b,4)));
}
long double A1(long double a){
    // input in celsius
    return ((3.2891) - (2.391*pow(10,-3)*a) + (2.8446*pow(10,-4)*pow(a,2)) - (2.82*pow(10,-6)*pow(a,3)) + (8.477*pow(10,-9)*pow(a,4)));
}
long double A2(long double a){
    // input in celsius
    return ((6.245*pow(10,-5)) - (3.913*pow(10,-6)*a) -(3.499*pow(10,-8)*pow(a,2)) +( 7.942*pow(10,-10)*pow(a,3)) - (3.299*pow(10,-12)*pow(a,4)));
}
long double rho(long double v, long double b, long double a1, long double a2, long double p){
    //rho in g/cm^3
    return pow( (v - ((v*p)/(b+a1*p+a2*p*p))) ,-1);
}
long double f0H2O(long double Ps, long double P, long double rho, long double T){
    return Ps*exp(18.0152*(P-Ps)/(rho*83.14*(T+273.15)));
}
long double Ps(long double T, long double Tc, long double Pc){
    //Ps in bar
    T+=273.15;
    Tc+=273.15;
    long double c= (1-(T/Tc));
    return Pc*exp(Tc*( -7.8595178*c + 1.8440825*pow(c,1.5) -11.786649*pow(c,3) + 22.68071*pow(c,3.5) -15.9618719*pow(c,4) + 1.8012250*pow(c,7.5) )/T);
}
long double ai(long double Tc, long double Pc, long double T, long double mi){
    Tc+=273.15;
    T+=273.15;
    return 0.457236*pow((83.14*Tc),2)*pow( (1+mi*(1-pow((T/Tc),0.5))) , 2)/Pc;
}
long double bi(long double Tc, long double Pc){
    Tc+=273.15;
    return 0.077796*83.14*Tc/Pc;
}
long double mi(long double wi){
    return (0.37464 + 1.54226*wi - 0.26992*wi*wi);
}
long double f(long double z, long double A, long double B){
    return ( pow(z,3) - (1-B)*z*z +(A-2*B-3*B*B)*z -(A*B-B*B-B*B*B) );
}
long double fprime(long double z,long double A,long double B){
    return ( 3*z*z - 2*(1-B)*z + (A-2*B-3*B*B) );
}
long double gamma(long double mi, long double lambda, long double Epsilon){
    return exp((2*mi*lambda) + (2*mi*mi*Epsilon));
}
long double Par( long double T, long double P, long double c1, long double c2, long double c3, long double c4, long double c5, long double c6, long double c7, long double c8, long double c9, long double c10){
    T=T+273.15;
    return ( c1+ c2*T + c3/T + c4*P + c5/P + ((c6*P)/T)+ (c7*T)/(P*P) + (c8*P)/(630-T) + c9*T*log(P)+ ((c10*P)/(T*T)) );
}
long double h(long double f0, long double T, long double ph20, long double deltaB){
    T+=273.15;
    long double n=-0.114535;
    return exp( (1-n)*log(f0) +n*log(83.14*T*ph20/18.015) +2*ph20*deltaB );
}
long double deltaB(long double T){
    long double tau = -5.279063;
    long double beta = 6.187967;
    long double n=-0.114535;
    T=T+273.15;
    return ( tau + beta*pow((1000/T) , 0.5) );
}
long double K0H2O(long double T){
    // input in celsius
    return pow(10, (-2.209 + 3.097*0.01*T - 1.098*pow(10,-4)*T*T + 2.048*pow(10,-7)*pow(T,3)) );
}
int main(){
    ifstream inputFile("Input.txt");
    ifstream input("Input1.txt");
    ofstream y_H2Ovalues("y_H2Ovalues.txt");
    ofstream x_co2values("x_co2values.txt");
    ofstream phivalues("phivalues.txt");
    ofstream Hvalues("Hvalues.txt");
    ofstream Kh2Ovalues("kh2Ovalues.txt");
    ofstream ERRORvalues("ERRORvalues.txt");
    vector <long double> Temp, Press, Molality;
    if (inputFile.is_open()) {
        long double temp, press, molality;
        while (inputFile >> temp >> press >> molality) {
            Temp.push_back(temp);
            Press.push_back(press);
            Molality.push_back(molality);
        }
        inputFile.close();
    } else {
        cout << "Unable to open file." << endl;
    }
    long double Tc, Pc, TcH2O, PcH2O;
    input >> Pc >> Tc >> TcH2O >> PcH2O;
    vector <long double> X_CO2_actual(101);
    vector <long double> RHO(101);
    vector <long double> PressS(101);
    vector <long double> FugacityH20(101);
    vector <long double> AParameter(101);
    vector <long double> BParameter(101);
    for( int i =0;i<101;i++){
        input >> X_CO2_actual[i];
        Temp[i]-=273.15;
    } 
    for( int i =0;i<101;i++){
        RHO[i]=rho(V0(Temp[i]),B(Temp[i]),A1(Temp[i]),A2(Temp[i]),Press[i]);

        PressS[i]=Ps(Temp[i],TcH2O,PcH2O);

        FugacityH20[i]=f0H2O(PressS[i],Press[i],RHO[i],Temp[i]);

        AParameter[i]=ai(Tc, Pc, Temp[i], mi(0.224))*Press[i]/pow((83.14*(Temp[i]+273.15)) , 2);

        BParameter[i]=bi(Tc, Pc)*Press[i]/(83.14*(Temp[i]+273.15));
    }
    
    
    // CALCULATION OF Z
    long double d2= 1-pow(2,0.5);
    long double d1= 1+pow(2,0.5);
    vector <long double> Z(101);
    for (int i=0;i<101;i++){
        long double xi=2, xi1=5;
        long double temp=1 ,Zl,Zg;
        while(abs(xi1-xi)>0.0001){
            xi=temp;
            xi1= xi-(f(xi, AParameter[i] ,BParameter[i])/fprime(xi, AParameter[i] ,BParameter[i]));
            temp=xi1;
        }
        Zg=xi1;
        xi=BParameter[i];
        temp=1000000;
        xi1= 100;
        while(abs(xi1-xi)>0.0001){
            xi=temp;
            xi1= xi-(f(xi, AParameter[i] ,BParameter[i])/fprime(xi, AParameter[i] ,BParameter[i]));
            temp=xi1;
        }
        Zl = xi1;
        long double q=Zg-Zl + log((Zl-BParameter[i])/(Zg-BParameter[i])) -AParameter[i]/(BParameter[i]*(d2-d1))*log( ((Zl+d1*BParameter[i])/(Zl+d2*BParameter[i]))*((Zg+d2*BParameter[i])/(Zg+d1*BParameter[i])) ) ;
        if(q>=0){
            Z[i]=Zl;
        }else {
            Z[i]=Zg;
        }
    }

    vector <long double> Phi(101);
    vector <long double> H(101);
    vector <long double> Gamma(101);
    vector <long double> K_CO2(101);
    vector <long double> K_H2O(101);
    vector <long double> y_H2O(101);
    vector <long double> y_CO2(101);
    vector <long double> x_CO2(101);
    vector <long double> ERRor(101);
    for (int i =0;i<101;i++){

        Phi[i] = exp( (Z[i]-1) - log(Z[i]-BParameter[i]) - AParameter[i]/(BParameter[i]*(d2-d1))*log( (Z[i]+d2*BParameter[i])/(Z[i]+d1*BParameter[i]) ));
        phivalues << Phi[i]<< endl;

        H[i] = h(FugacityH20[i],Temp[i],RHO[i], deltaB(Temp[i]));
        Hvalues << H[i] << endl;

        Gamma[i]=gamma(Molality[i], Par( Temp[i], Press[i], -0.0652869, 1.6790636*pow(10,-4), 40.838951, 0, 0, -3.9266518*0.01, 0, 2.1157167*0.01, 6.5486487*pow(10,-6), 0), Par( Temp[i], Press[i], -1.144624*0.01, 2.8274958*pow(10,-5), 0, 0, 0, 1.3980876*0.01, 0, -1.4349005*0.01, 0 ,0));
        
        K_CO2[i]=H[i]*Gamma[i]/(Press[i]*Phi[i]);

        K_H2O[i]=(K0H2O(Temp[i])*exp(((Press[i]-1)*18.18)/(83.14*(Temp[i]+273.15))))/(FugacityH20[i]*Press[i]);
        
        Kh2Ovalues << K_H2O[i] << endl;
        
        y_H2O[i]= (1-1/K_CO2[i])/(1/K_H2O[i]-1/K_CO2[i]);
        y_H2Ovalues << y_H2O[i] << endl;
        
        y_CO2[i]= 1/(1+y_H2O[i]);
        
        x_CO2[i]= y_CO2[i]/K_CO2[i];
        x_co2values << x_CO2[i]*100 << endl;
        
        ERRor[i] = abs(1-X_CO2_actual[i]/(100*x_CO2[i]))*100;
        ERRORvalues << ERRor[i] << endl;
    }
    return 0;
}