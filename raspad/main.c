//
//  main.c
//  raspad
//
//  Created by Smirnov Ivan on 04.04.17.
//  Copyright © 2017 bmstu sm3. All rights reserved.
//

#include <stdio.h>      
#include <math.h>

double _ro(double p, double R, double T){
    
    return p/(R*T);
}

double _a(double k,double R,double T){
    
    return sqrt(k*R*T);
    }


double _psi(double s, double k){
    
    if (s<1) {
        
        return (2.0/(k-1.0))*(pow(s,(k-1.0)/(2.0*k))-1.0);
        
    }
    
    if (s>=1) {
        
        return sqrt(2.0/k)*(s-1.0)/sqrt(s*(k+1.0)+(k-1.0));
    }
    return 0.0;
}



double _phi(double s, double k){
    
    
    if (s<1) {
        
        return pow(s,1/k);
    }
    
    if (s>=1) {
        
        return ((k+1.0)*s+(k-1.0))/((k-1.0)*s+(k+1.0));
    }
        return 0.0;
}



double solve(double pObsh,double u[],double k[], double R[],double T[], double p[]){
    
    
    return u[0]+_a(k[0], R[0], T[0])*_psi(pObsh/p[0], k[0])
    -u[1]+_a(k[1], R[1], T[1])*_psi(pObsh/p[1], k[1]);
}


double hi(double pOb,double u[],double k[], double R[],double T[], double p[]){
    

    return _a(k[0], R[0], T[0])*_psi(pOb/p[0], k[0])+_a(k[1], R[1], T[1])*_psi(pOb/p[1], k[1]);
}


double resh(double u[],double k[], double R[],double T[], double p[]){
    
    double x0 = 2*pow(10, 3);
    double x1 = 1*pow(10, 3);
    double eps = 0.0000001;
    double tmp = 0.0;
    
    while (fabs(x1-x0)>eps) {
        tmp=x1;
        x1=x1-(x1-x0)*solve(x1, u, k, R, T, p)/
        (solve(x1, u, k, R, T, p)-solve(x0, u, k, R, T, p));
        x0=tmp;
    }
    
    return x0;
    
}

void paramVR(int r, double Dz0,double Dzobsh,double Uobsh,double u,double a,double k,double Ro,double p,double R,double aObshPr){
    
  //  double Dzobsh=Uobsh+aObshPr;
    double point=4;
    double Kolich=fabs(Dz0-Dzobsh)/point;
    printf("Dz aSech RoPrav1 PPrav TPrav M1 point\n");
    double DzSech=0;
    for (double i=0; i<point+1; i++) {
       DzSech=Dz0+i*Kolich;
       double aSech=0;
        double M1=0;
        if (!(r==1)) {
             aSech=a+((k-1)/(k+1))*fabs(DzSech-fabs(Dzobsh));
             M1=fabs(u+aSech-a)/aSech;
        }else{
             aSech=a-((k-1)/(k+1))*(DzSech-fabs(Dzobsh));
             M1=fabs(-u+aSech-a)/aSech;
        }
       double RoPrav1=Ro*pow((aSech/a),(2/(k-1)));
       double PPrav1=p*pow((aSech/a),(2*k/(k-1)));
       double TPrav1=PPrav1/(RoPrav1*R);
      // point=point+2;
        printf("%f %f %f %f %f %f %f\n",DzSech,aSech,RoPrav1,PPrav1,TPrav1,M1,i);
        
    }

    
}

double* paramObl(double uObsh,double pObsh ,double ro[],double phi[],double R[],double k[],double T[],double u[],double a[],double p[]){

    
    double roObsh[2];
    double TObsh[2];
    static double aObsh[2];
    double Mobsh[2];
    double D[2];
    char info[2] = {'L', 'R'};
    printf("volna pObsh roObsh TObsh aObsh Mobsh D uObsh\n");
    for (int i=0; i<2; i++) {
        
        roObsh[i]=ro[i]*phi[i];
        TObsh[i]=pObsh/(R[i]*roObsh[i]);
        aObsh[i]=_a(k[i], R[i], TObsh[i]);
        Mobsh[i]=fabs(uObsh/aObsh[i]);
        D[i]=(roObsh[i]*uObsh-pow(ro[i], u[i]))/(roObsh[i]-ro[i]);
        printf("%c %f %f %f %f %f %f %f\n",info[i],pObsh,roObsh[i],TObsh[i],aObsh[i],Mobsh[i],D[i],uObsh);
        
        
    }
   
    return aObsh;
}

double* paramObl2(double uObsh[],double pObsh ,double ro[],double phi[],double R[],double k[],double T[],double u[],double a[],double p[]){
    
    
    double roObsh[2];
    double TObsh[2];
    static double aObsh[2];
    double Mobsh[2];
    double D[2];
    char info[2] = {'L', 'R'};
    printf("volna pObsh roObsh TObsh aObsh Mobsh D uObsh\n");
    for (int i=0; i<2; i++) {
        
        roObsh[i]=ro[i]*phi[i];
        TObsh[i]=pObsh/(R[i]*roObsh[i]);
        aObsh[i]=_a(k[i], R[i], T[i]);
        Mobsh[i]=fabs(uObsh[i]/aObsh[i]);
        D[i]=(roObsh[i]*uObsh[i]-pow(ro[i], u[i]))/(roObsh[i]-ro[i]);
        printf("%c %f %f %f %f %f %f %f\n",info[i],pObsh,roObsh[i],TObsh[i],aObsh[i],Mobsh[i],D[i],uObsh[i]);
        
    }
    return aObsh;
}







int main(int argc, const char * argv[]) {
    // insert code here...
    double p[2]={0.1*pow(10, 6),1*pow(10, 6)};
    double u[2]={0,0};
    double T[2]={273,1000};
    double R[2]={4125,287};
    double k[2]={1.41,1.4};
    double hiP[2];
    double a[2];
    double ro[2];
    double psi[2];
    double phi[2];
    double M[2];
    
    double pObsh=resh(u, k, R, T, p);
    double hi0=hi(0, u, k, R, T, p);
    double rasV=u[1]-u[0];
    
    printf("obl p T ro u a M\n");
    for (int i=0; i<2; i++) {
        hiP[i]=hi(p[i], u, k, R, T,p);
        a[i]=_a(k[i], R[i], T[i]);
        M[i]=u[i]/a[i];
        ro[i]=_ro(p[i], R[i], T[i]);
        psi[i]=_psi(pObsh/p[i], k[i]);
        phi[i]=_phi(pObsh/p[i], k[i]);
        printf("%d %f %f %f %f %f %f\n",i,p[i],T[i],ro[i],u[i],a[i],M[i]);
    }

    
    
    
    
    if (p[0]<p[1]) {
        
        if (hi0>rasV) {

            printf("Структура a  Ошибка\n");
            
        }
        if ((hi0<rasV)&&(hiP[0]>rasV)) {
            
            printf("Структура b \n");
            double *aObsh;
            double uObsh=u[0]+a[0]*psi[0];
            aObsh=paramObl(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            paramVR(0,uObsh-aObsh[1],u[1]-a[1], uObsh, u[1], a[1], k[1], ro[1], p[1], R[1], aObsh[1]);
            paramVR(0,uObsh+aObsh[0],u[0]+a[0], uObsh, u[0], a[0], k[0], ro[0], p[0], R[0], aObsh[0]);
            
        }
 
        if ((hiP[0]<rasV)&&(hiP[1]>rasV)) {
            
            printf("Структура v \n");
            double *aObsh;
            double uObsh=u[1]-a[1]*psi[1];
            aObsh=paramObl(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            double Dz00P=u[0]+a[0];
            double Dzobsh=uObsh+aObsh[0];
            paramVR(0,Dz00P,Dzobsh, uObsh, u[0], a[0], k[0], ro[0], p[0], R[0], aObsh[0]);
        }

        if ((hiP[1]<rasV)) {
            
            printf("Структура g \n");
            double uObsh[2];
            uObsh[0]=u[0]+a[0]*psi[0];
            uObsh[1]=u[1]-a[1]*psi[1];
            paramObl2(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            
        }
        
    }else{
        
        if (hi0>rasV) {
            
            printf("Структура a  Ошибка\n");
            
        }
        if ((hi0<rasV)&&(hiP[1]>rasV)) {
            
            printf("Структура b \n");
            double *aObsh;
            double uObsh=u[0]+a[0]*psi[0];
            aObsh=paramObl(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            paramVR(0,uObsh-aObsh[1],u[1]-a[1], uObsh, u[1], a[1], k[1], ro[1], p[1], R[1], aObsh[1]);
            paramVR(0,uObsh+aObsh[0],u[0]+a[0], uObsh, u[0], a[0], k[0], ro[0], p[0], R[0], aObsh[0]);
            
        }
        
        if ((hiP[0]>rasV)&&(hiP[1]<rasV)) {
            
            printf("Структура v \n");
            double *aObsh;
            double uObsh=u[1]-a[1]*psi[1];
            aObsh=paramObl(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            double Dz00P=u[0]+a[0];
            double Dzobsh=uObsh+aObsh[0];
            paramVR(0,Dzobsh,Dz00P, uObsh, u[0], a[0], k[0], ro[0], p[0], R[0], aObsh[0]);
            
        }
        
        if ((hiP[0]<rasV)) {
            
            printf("Структура g \n");
            double uObsh[2];
            uObsh[0]=u[0]+a[0]*psi[0];
            uObsh[1]=u[1]-a[1]*psi[1];
            paramObl2(uObsh, pObsh, ro, phi, R, k, T,u,a,p);
            
        }
        
        
    }
    
    
    return 0;
}
