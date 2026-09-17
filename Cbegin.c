#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#define N 20
#define M 300000
#define therm 50000
#define Nt 40

    int jx = 1;
    int jy = 1;

    float tin = 2.25;
    float tfin = 2.28;

    int S[N][N];
    float T[Nt] = {0};

    double beta[Nt];
    float Ct[Nt] = {0};
    float Mt[Nt] = {0};

    float Eavg[M] = {0};
    float E2avg[M] = {0};
    float Mavg[M] = {0};

float siteE(int sx, int sy){
    double energy; 
    energy = S[sx][sy] * ( -jx * ( S[(sx+1)%N][sy] + S[(sx+N-1)%N][sy])  - jy * ( S[sx][(sy+1)%N] + S[sx][(sy+N-1)%N]) );
    return energy;
}

float Etot() {
    double sum = 0; 
    for(int sx=0; sx<N; sx++){
        for(int sy=0; sy<N; sy++){
    sum += siteE(sx, sy);
        }}
    return sum/2;
}



void flip(int sx, int sy, int p, double *Ecurrent)
{
    double diffE = -2 * siteE(sx, sy);
    if(diffE < 0){
        S[sx][sy] *= -1;
        *Ecurrent += diffE;
    }
    else{
        double r = (double)rand() / RAND_MAX;
        if(expf(-diffE*beta[p]) > r){
            S[sx][sy] *= -1;
            *Ecurrent += diffE;
        }
}}










float Etherm(int p){
    for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
        S[i][j] = 2 * (rand() % 2) - 1;
    }}

    double Ecurrent = Etot();
    double sum = 0;

    for(int t=0; t<therm; t++){
        for(int l=0; l<N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);
        }
    }

    for(int k =0; k<M; k++){

        for(int l=0; l<N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);}
    Eavg[k] = Ecurrent/(N*N);
    }
    
    for(int k=0; k<M; k++){
        sum += Eavg[k];
    }
return sum / M;
}




float Ctherm(int p){
    for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
        S[i][j] = 2 * (rand() % 2) - 1;
    }}

    double Ecurrent = Etot();

    double E = 0;
    double E2 = 0;

    for(int t = 0; t < therm; t++){
        for(int l=0; l < N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);}
    }

    for(int k = 0; k < M; k++){

        for(int l=0; l<N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);}
    Eavg[k] = Ecurrent/(N*N);
    E2avg[k] = Eavg[k]*Eavg[k];
    }
    
    for(int k=0; k<M; k++){
        E += Eavg[k]/M;
        E2 += E2avg[k]/M;
    }
return beta[p] * beta[p] * N * N * (E2 - E*E);
}

void C_t(){
    for(int p=0; p<Nt; p++){
        Ct[p] = Ctherm(p);
    }
}



double magn(){
    int mag = 0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            mag += S[i][j]; 
        }
    }
    return (double)abs(mag);
}


double Mtherm(int p){

    double Mmean = 0;

    for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
        S[i][j] = 2 * (rand() % 2) - 1;
    }}

    double Ecurrent = Etot();

    for(int t=0; t<therm; t++){
        for(int l=0; l<N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);}
    }

    for(int k =0; k<M; k++){

        for(int l=0; l<N*N; l++){
        int sx = rand() % N;
        int sy = rand() % N;
        flip(sx, sy, p, &Ecurrent);}
    Mavg[k] = magn()/ (double)(N*N);
    }
    
    for(int k=0; k<M; k++){
        Mmean += Mavg[k];
    }
return Mmean/M;
}


void M_t(){
    for(int p=0; p<Nt; p++){
        Mt[p] = Mtherm(p);
    }
}






int main()
{
    T[0] = tin;
    T[Nt] = tfin;
    for(int i=0; i<Nt; i++){
        T[i+1] = T[i] + fabs(tfin - tin)/Nt;
        printf("%f", T[i]);
        printf("\n");
    }


    for(int i = 0; i < Nt; i++){
        beta[i] = 1.0 / T[i];
    }


    srand(time(NULL));


printf("\n");



printf("T, C: \n");
for(int p=0; p<Nt; p++){
    Ct[p] = Ctherm(p);
    printf("%f, %f", T[p], Ct[p]);
    printf("\n");
}

FILE *file = fopen("C_t20.csv", "w");

fprintf(file, "temp,C\n");

for(int i = 0; i < Nt; i++){
    fprintf(file, "%f, %f \n", T[i], Ct[i]);
}
fclose(file);


    return 0;
}

