#include<math.h>
#include<iostream>
#include<fstream>

//Afisarea valorilor unei functii – L1

double f (double x) { 
double rezultat;
if(x<0) rezultat=2+3*x+x;
if ((unsigned)x<5) rezultat =1/(2+3*x);
if(x=5) rezultat=(sin(x)/2+2*cos(3*x)); 
return rezultat; 
}

int main()
{
x=2;
i=1;
    while(x<=10)
    {
        if(i%24==0)
        {
            cout<<"Apasati o tasta pentru a continua "<<endl;
            getch();
        }
        cout<<"f("<<x<<")="<<calcululfunctiei(x)<<endl;
        x+=0.5;
        i++;
    }
}

//Factorialul unui numar (recursiv si iterativ) - L1 

double factorial_recursiv(int n)
{
    if (n==0) return 1.0;
    else return (factorial_recursiv(n-1)*n);
}

double factorial_iterativ(int n)
{
    int i;
    double f=1.0;
    for(i=1;i<=n;i++) f*=i;
    return f;
}



//Ecuatia de gradul I – L1


int main()
{
cout << "Rezolvarea ecuatiei de gradul intai" << endl;
float a,b,x;
cout<<"a="; cin>>a;
cout<<"b="; cin>>b;

if(a==0)
    if(b==0)
       cout<<"Exista o infinitate de solutii!";
    else
       cout<<"Ecuatie imposibila!";
else
  {
	x=-b/a;
      cout<<"Solutia ecuatiei este "<<x;
  }
}



//Ecuatia de gradul II – L1


int main(){
cout << "Rezolvarea ecuatiei de gradul al doilea" << endl;
    float a,b,c,d,x1,x2;
    cout<<"a=";cin>>a;
    cout<<"b=";cin>>b;
    cout<<"c=";cin>>c;

    if(a!=0&&b!=0)
     {
        {
			d=b*b-4*a*c;
        }

    if(d>=0)

         {
			 x1=(-b+sqrt(d))/(2*a);
			 x2=(-b-sqrt(d))/(2*a);
			 cout<<"x1="<<x1<<endl;
			 cout<<"x2="<<x2<<endl;
         }
      else
    if(d<0)
          {cout<<"Ecuatia nu are solutii in R.";}

      }

    if (a==0&&b==0&&c==0)

          {cout<<"Ecuatia are o infinitate de solutii!";}
      else
    if(a==0&&b==0)

          {cout<<"Ecuatia nu are solutii!";}
      else
    if(a==0)
          {if(c!=0)
             {
				 x1=-b/c;
				 cout<<"x="<<x1;
             }
        else
             cout<<"Ecuatia nu are solutii!";
          }
      else
    if(b==0)
          {
			  if(-c/a>=0&&c!=0)

             {
				 x1=sqrt(-c/a);
				 cout<<"x="<<x1;
             }
        else
           if(c==0)
             {
				 cout<<"x=0";
             }
      else
             {cout<<"Ecuatia nu are solutii!";}

          }

}

//Ecuatia de gradul III (metoda lui Cardan) - L1 

// int main()
// {
//     double a,b,c,d;
//     double p,q,v,u,y1,y2r,y2i,x1,x2r,x2i,x3r,x3i,delta;
//     cout<<a<<"*x^3+ "<<b<<"*x^2 + "<<c<<"x +"<<d<<" = 0"<<endl;
//     cout<<"Solutiile ecuatiei sunt: "<<endl;
//     p=c/a-b*b/(3*a*a);
//     q=2*b*b/(27*a*a*a)-b*c/(3*a*a)+d/a;
//     delta=p*p*p/27+q*q/4;
//     u=exp(log(sqrt(delta) - q/2)/3);
//     v=-exp(log(q/2+sqrt(delta))/3);
//     y1=u+v;
//     y2r=y1/-2;
//     y2i=(u-v)*sqrt(3.0)/2;
//     x1=y1-b/3/a;
//     x2i=y2i;
//     x3r=x2r;
//     x3i=-x2r;
//     cout<<"x1: "<<x1<<endl;
//     cout<<"x2: "<<x2r<<" + "<<x2i<<" *i"<<endl;
//     cout<<"x3: "<<x3r<<" + "<<x3i<<" *i"<<endl;
// Permutari grad n, Combinari si Aranjamente de n luate cate k  - L1
// Main
//     int n, k, f_n, f_k, f_nk, nCk, nAk, nP, i;
// citire de la tastatura n si k;
// initializare factoriale cu 1;
//     for(i = 1; i <= n; i++)
//         f_n = f_n * i;
//     for(i = 1; i <= k; i++)
//         f_k = f_k * i;
//     for(i = 1; i <= n-k; i++)
//         f_nk = f_nk * i;
// nP  = f_n;
// nAk = f_n / f_nk;
// nCk = f_n / (f_k * f_nk);
// Conversii incompatibile – L2
// int nr; // cu semn
// unsigned int nr1; // fara semn
// float rez = -2.1;
// // folosirea Mixta a tipurilor de date - fara conversie implicita
//     // Numarul intreg este cu semn
//     nr = rez;
//     cout<<"\nFara conversie explicita -> Numarul intreg este "<<nr;

//     // Folosirea mixta a tipurilor de date -> cu conversie implicita
//     // Numarul intreg este cu semn
//     nr = (int)rez;
//     cout<<"\nCu conversie implicita -> Numarul intreg este: "<<nr;

//     //Folosirea mixta a tipurilor de date - fara conversie implicita
//     // Numarul intreg este fara semn
//     nr1=rez;
//     cout<<"\nFara conversie implicita -> Numarul intreg este: "<<nr1;

//     // Folosirea mixta a tipurilor de date - cu conversie implicita
//     // Numarul intreg este fara semn
//     nr1= (int)rez;
//     cout<<"\nCu conversie implicita -> numarul intreg este: "<<nr1;

// Rotunjiri limita – L2
// float nr1=4.3275,nr2=4.3285;

//     // Impunem reprezentarea cu 3 zecimale
//     printf("Rotunjire la limita -> Nr1: %.3f, Nr2. %.3f\n",nr1,nr2);



//Situatii in care apar numere speciale – L2



// union sitSpeciala
//     {
//         long int nrI;
//         float nrR;
//     }situatie;
//     float numerator = 0;
//     float denominator = 0;
//     float underSqrt = -1;

//     // Zero
//     situatie.nrI = 0x00000000; // zero cu plus
//     cout<<hex<<"\nCalcul eronat... zero 1... "<<situatie.nrR;
//     situatie.nrI = 0x80000000; // zero cu minus
//     cout<<"\nCalcul eronat... zero 2... "<<situatie.nrR;

//     // Infinit - cu plus si cu minus
//     cout<<"\nCalcul eronat... infinit 1... "<<1 / denominator;
//     cout<<"\nCalcul eronat... infinit 2... "<<1 / -denominator;

//     // NaN - Calcul eronat - Nedeterminari
//     cout<<"\nCalcul eronat... Nedeterminare 1... "<<numerator/denominator;
//     cout<<"\nCalcul eronat... Nedeterminare 2... "<<sqrt(underSqrt);

//     // Final
//     // cin>>numerator;
//     return 0;


//Tabelarea mai multor functii – L3


float f(float x)
{
    return(5/2*pow(x,4)-1/3*pow(x,3)+4/3*pow(x,2)-5*x+4);
}

float f1(float x)
{
    return(2*pow(x,4)-5*pow(x,2)+15*x-5);
}

float f2(float x)
{
    return(5*pow(x,3)-10*pow(x,2)+2*x-3);
}

float f3(float x)
{
    return(-3*pow(x,3)+5*pow(x,2)-10*x+8);
}

void TabelareaFunctiei(float functie(float x), float xmin, float xmax,
                       float n, float x[100], float y[100])
{
    int i;
    float d;
    d=(xmax-xmin)/(n-1);
    for(i=0;i<n;i++)
    {
        x[i] = xmin+i*d;
        y[i] = functie(x[i]);
    }
}


//Citirea fisierelor – L4

void citfis(int n_col, float a1[100],float a2[100],float a3[100], float a4[100], float a5[100])
{
    ifstream f("amestec.txt",ios::in);
    i=0;
    while(f.eof()==0)
    {
        i++;
        switch(n_col)
        {
            case 1:
                f>>a1[i];break;
            case 2:
                f>>a1[i]>>a2[i];break;
            case 3:
                f>>a1[i]>>a2[i]>>a3[i];break;
            case 4:
                f>>a1[i]>>a2[i]>>a3[i]>>a4[i];break;
            case 5:
                f>>a1[i]>>a2[i]>>a3[i]>>a4[i]>>a5[i];break;
            default:
                cout<<"S-a depasit numarul de coloane admise:";
        }
    }
    ndat=i;
    f.close();
}

int main()
{
    float y[100];
    citfis(2,M,G,ax,ax,ax);
    float S=0,Mmed=0;
    for(i=1;i<=ndat;i++)
    {
        cout<<i<<" "<<M[i]<<" "<<G[i]<<"\n";
        S+=G[i]/M[i];
    }
    for(i=1;i<=ndat;i++)
    {
        y[i] = (G[i] / M[i])/S;
        Mmed+=y[i]*M[i];
    }
    cout<<"Masa molara medie este "<<Mmed<<"\n";
}

// Tabelarea functiei in fisier – L4

float fct (float x)
{
    return (x*x*x*x*1/4-x*x*x*5/2+x*x*7/2-3*x+3);
}

void TabFunctie(float xmin, float xmax, float n, float z[300],float y[300])
{
    int i;
    float d;
    d=(xmax-xmin)/(n-1);
    for(i=0;i<n;i++)
    {
        z[i] = xmin + i * d;
        y[i] = f(z[i]);
    }
}

    ofstream f("tabelarefunctie.txt",ios::out);
    TabFunc(xmin,xmax,n,x,fx);
    for(i=0;i<n;i++)
    {
        f<<"\n"<<i+1<<"\t"<<x[i]<<"\t"<<fx[i];
    }


//Citirea matricei in fisier – L4


void citmat(int n,float a[10][10])
{
    int i,j;
    ifstream f("mat_a.dat",ios::in);
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            f>>a[i][j];
        }
    }
    f.close();
}
Citirea si scrierea in matrice – L5
void citirematrice(float a[10][10], int n, int m)
{
    ifstream f("mat_a.txt",ios::in);
    int i,j;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            f>>a[i][j];
        }
    }
    f.close();
}

void scrierematrice(int n,int m)
{
    ofstream f("mat_a.txt",ios::out);
    int i,j;
    float x;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            cout<<"a["<<i<<"]["<<j<<"]: ";
            cin>>x;
            f<<x<<"\t";
        }
        f<<"\n";
    }
    f.close();
}
Adunarea a doua matrici – L5
void citirematriceaa(float a[10][10],int n,int m)
{
    ifstream f("mat_a.txt",ios::in);
    int i,j;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=m;j++)
        {
            f>>a[i][j];
        }
    }
    f.close();
}

void citirematriceab(float b[10][10],int p,int q)
{
    ifstream f("mat_b.txt",ios::in);
    int i,j;
    for(i=1;i<=p;i++)
    {
        for(j=1;j<=q;j++)
        {
            f>>b[i][j];
        }
    }
    f.close();
}

void citirematriceas(float s[10][10],int n,int m)
{
    ifstream f("suma_mat.txt",ios::in);
    int i,j;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=m;j++)
        {
            f>>s[i][j];
        }
    }
    f.close();
}

void scrierematriceaa(int n,int m)
{
    ofstream f("mat_a.txt",ios::out);
    int i,j;
    float x;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=m;j++)
        {
            cout<<"a["<<i<<"]["<<j<<"]: ";
            cin>>x;
            f<<x<<"\t";
        }
        f<<"\n";
    }
    f.close();
}

void scrierematriceab(int p,int q)
{
    ofstream f("mat_b.txt",ios::out);
    int i,j;
    float x;
    for(i=1;i<=p;i++)
    {
        for(j=1;j<=q;j++)
        {
            cout<<"b["<<i<<"]["<<j<<"]: ";
            cin>>x;
            f<<x<<"\t";
        }
        f<<"\n";
    }
    f.close();
}

void suma(int n, int m, int p, int q, float a[10][10],float b[10][10])
{
    int i,j;
    float s[10][10];
    if((n==p) && (m==q))
    {
        for(i=1;i<=n;i++)
        {
            for(j=1;j<=m;j++)
            {
                s[i][j] = a[i][j] + b[i][j];
            }
        }

        ofstream f("suma_mat.txt",ios::out);
        for(i=1;i<=n;i++)
        {
            for(j=1;j<=m;j++)
            {
                f<<s[i][j]<<"\t";
            }
            f<<"\n";
        }
        f.close();
    }
    else cout<<"\nMatricele nu se pot aduna! (Nu sunt egale)";
}

// Inmultirea a doua matrici – L5
void citire(int x[10][10],int n,int m,int k)
{
    if(k==1)
    {
        for(int i=1;i<=n;i++)
        {
            for(int j=1;j<=m;j++)
            {
                f>>x[i][j];
            }
        }
        f.close();
    }
    else
    {
        fstream f("matrice2.txt",ios::in);
        for(int i=1;i<=n;i++)
        {
            for(int j=1;j<=m;j++)
            {
                f>>x[i][j];
            }
        }
        f.close();
    }
}

void afisare(int x[10][10],int n,int m)
{
    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=m;j++)
        {
            cout<<"\t"<<x[i][j];
        }
        cout<<"\n";
    }
}

void produs(int c[10][10],int d[10][10],int e[10][10],int n1,int m1,int m2)
{
    for(int i=1;i<=n1;i++)
    {
        for(int j=1;j<=m1;j++)
        {
            int s=0;
            for(int k=1;k<=m2;k++) s=s+c[i][k]*d[k][j];
            e[i][j]=s;
        }
    }
}
main
    int a[10][10],b[10][10],rez[10][10],n1,m1,n2,m2,l;
    citire(a,n1,m1,1);
    afisare(a,n1,m1);
    citire(b,n2,m2,2);
    afisare(b,n2,m2);
    if(n2==m1)
    {
        produs(a,b,rez,n1,m1,n2);
        cout<<"Matricea rezultanta este:\n";afisare(rez,n1,m2);
    }
    else cout<<"Inmultirea nu se poate efectua";
}

//Transpunerea unei matrici – L5

void transpusa (int x[10][10], int y[10][10], int n,int m)
{
for (int i=1;i<=m; i++)
 for (int j=1;j<=n; j++)
 y[i][j]=x[j][i];
}

Inverseara matricelor (Gauss-Jordan) – L6
void citire(double x[10][10],int n,int m)
{
    fstream f("mat_a.txt",ios::in);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            f>>x[i][j];
    f.close();
}

void afisare(double x[10][10],int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            cout<<"\t"<<x[i][j];
        cout<<"\n";
    }
}
Inversarea matricei A in matricea B (Gauss-Jordan) – L6
void inversa(int dim, double A[10][10],double B[10][10])
{
    int l, lu, c, ind;
    double maxim,interm;
    fstream f("matrice.txt",ios::out);
    for(l=0;l<dim-1;l++)
    {
        maxim=A[l][l];
        ind=l;
        for(lu=l+1;lu<dim;lu++)
            if(A[lu][l]>maxim) ind=lu;
        if(ind!=l)
        {
            for(int i=0;i<dim;i++)
            {
                interm=A[l][i];
                A[l][i]=A[ind][i];
                A[ind][i]=interm;
                interm=B[l][i];
                B[l][i]=B[ind][i];
                B[ind][i]=interm;
            }
        }
        for(lu=l+1;lu<dim;lu++)
        {
            maxim=-(A[lu][l]/A[l][l]);
            for(c=0;c<dim;c++)
            {
                A[lu][c]+=(maxim)*A[l][c];
                B[lu][c]+=(maxim)*B[l][c];
            }
        }
        f<<"\nAm interschimbat doua linii \n";
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
                f<<"\t"<<B[i][j];
            f<<"\n";
        }
    }
    for(l=dim-1;l>0;l--)
    {
        for(lu=l-1;lu>=0;lu--)
        {
            maxim=-(A[lu][l]/A[l][l]);
            for(c=0;c<dim;c++)
            {
                A[lu][c]+=(maxim)*A[l][c];
                B[lu][c]+=(maxim)*B[l][c];
            }
        }
        f<<"\n";
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
                f<<"\t"<<B[i][j];f<<"\n";
        }
    }
    for(l=0;l<dim;l++)
    {
        maxim=(1/A[l][l])-1;
        for(c=0;c<dim;c++)
        {
            A[l][c]+=(maxim)*A[l][c];
            B[l][c]+=(maxim)*B[l][c];
        }
    }
    f<<"\n";
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
            f<<"\t"<<B[i][j];
        f<<"\n";
    }
    f.close();
}

Solutia sistemului prin algoritmul Gauss – L6
void citire (float x[10][10],int n,int m)
{
    fstream f("mat_a.txt",ios::in);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=m;j++)
            f>>x[i][j];
    f.close();
}

void afisare(float x[10][10],int n,int m)
{
    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=n;j++) cout<<"\t"<<x[i][j];
        cout<<"\n";
    }
}

void liber(float y[10],int n)
{
    fstream f("sir_b.txt",ios::in);
    for(int i=1;i<=n;i++) f>>y[i];
    f.close();
}

void gauss1(float a[10][10], float b[10],int n, float x[10], int print)
{
    int l,k,i,j;
    float p, mult;
    char car;
    ofstream f("afisare.txt",ios::out);
    for(k=1;k<=n-1;k++)
    {
        l=k;
        for(i=k+1;i<=n;i++)
            if(abs(a[i][k])<=abs(a[l][k]))
                l=i;
        if(l!=k)
        {
            for(j=k;j<=n;j++)
            {
                p=a[k][j];
                a[k][j]=a[l][j];
                a[l][j]=p;
            }
            p=b[k];
            b[k]=b[l];
            b[l]=p;
        }

        for(i=k+1;i<=n;i++)
        {
            mult=a[i][k]/a[k][k];
            for(j=1;j<=n;j++) a[i][j]=a[i][j]-mult*a[k][j];
            b[i]=b[i]-mult*b[k];
        }
        if(print==1)
        {
            f<<"Iteraria "<<k<<endl;
            for(i=1;i<=n;i++)
            {
                for(j=1;j<=n;j++) f<<"\t"<<a[i][j];
                f<<"*\t"<<b[i]<<"*\t";
                f<<"\n";
            }
            f<<"\n";
        }
        x[n]=b[n]/a[n][n];
        for(i=n-1;i>=1;i--)
        {
            p=0;
            for(j=i+1;j<=n;j++) p=p+a[i][j]*x[j];
            x[i]=(b[i]-p)/a[i][i];
        }
    }
    if(print==1)
    {
        f<<"\nSolutia sistemului prin algoritmul Gauss este: ";
        for(i=1;i<=n;i++)
            f<<"\nx("<<i<<"): "<<x[i];
    }
}

Algoritmul bisectiei succesive – L7
float f(float x)
{
    return(4*x*x*x-8*x*x+12*x-8);
}

void bisuc(float x1,float x2,float epsif, float epsix,int nmax, float &xr, int &kod)
{
    int iter=0,contor,linie,nlin;
    float fx1,fx2,fx3,x3;
    ofstream fwrite("fisierbisectie.txt",ios::out);
    kod=0;
    fx1=f(x1);
    fx2=f(x2);
    fwrite<<"Procedura bisectie succesiva\n\n";
    if(fx1*fx2>0)
    {
        fwrite<<"Intervalul dat nu contine solutia\n";
        kod=2;
    }
    else
    {
        fwrite<<"Iter\t\tx1/fx1\t\tx2/fx2\t\tx3/fx3\n";
        iter++;
        x3=(x1+x2)/2;
        fx3=f(x3);
        while((fabs(fx3)>epsif) && (fabs(x1-x2)>epsix))
        {
            fwrite<<iter<<"\t\t"<<x1<<"\t\t"<<x3<<"\t\t"<<x2<<"\n";
            fwrite<<"\t\t"<<fx1<<"\t\t"<<fx2<<"\t\t"<<fx3<<"\n\n";
            if(fx1*fx3<0)
            {
                x2=x3;
                fx2=fx3;
            }
            else
            {
                x1=x3;
                fx1=fx3;
            }
            x3=(x1+x2)/2;
            fx3=f(x3);
            iter=iter+1;
            if(iter==nmax)
            {
                kod=1;
                xr=x3;
                goto a1;
            }
            if(fabs((x1-x2)/x3)<1e-11)
            {
                kod=0;
                xr=x3;
                goto a1;
            }
        }
        a1:xr=x3;
        if(iter<10) fwrite<<iter<<"\t\t"<<x1<<"\t\t"<<x3<<"\t\t"<<x2<<"\n";
        else fwrite<<iter<<"\t\t"<<x1<<"\t\t"<<x2<<"\t\t"<<x3<<"\n";
        fwrite<<"\t\t"<<fx1<<"\t\t"<<fx2<<"\t\t"<<fx3<<"\n";
        fwrite<<"Solutia ecuatiei: "<<xr;
        fwrite.close();
    }
}

Rezolvarea sistemului de ecuatii neliniare (Jacobianul) – L7
double func1(double x,double y)
{
    double k;
    return 4*x*x+2*x-3*y*y+3;
}

double func2(double x,double y)
{
    return x*x+x*y+y*y-61;
}

void jacobi(double x,double y,double j[10][10])
{
    j[1][1]=8*x+2;
    j[1][2]=-6*y;
    j[2][1]=2*x+y;
    j[2][2]=x+2*y;
}

void newton2(double &x1,double &y1,double eps)
{
    double f1,g1,fx,fy,gx,gy,del,x2,y2,dx,dy,J[10][10];
    fstream f("f_newton.txt",ios::out);
    int iter=99;
    int i=0;
    do
    {
        i++;
        f1=func1(x1,y1);
        g1=func2(x1,y1);
        jacobi(x1,y1,J);
        fx=J[1][1];fy=J[1][2];
        gx=J[2][1];gy=J[2][2];
        del=fx*gy-fy*gx;
        dx=(fy*g1-f1*gy)/del;
        dy=(f1*gx-fx*g1)/del;
        f<<"\n\n iteratia = "<<i<<"\nx = "<<x1<<"n y= "<<y1<<"\n f= "<<f1<<"\n g = "<<g1;
        x2=x1+dx;
        y2=y1+dy;
        x1=x2;y1=y2;
        if(i>=iter) break;
    }while(fabs(dx) >= eps && fabs(dy) >= eps);
    i++;
    f<<"\n\n Solutia sistemului de ecuatii: ";
    f<<"\n x = "<<x1<<"\n y = "<<y1<<"\n f = "<<f1<<"\n g = "<<g1;
    f.close();
}

Algoritmul Newton-Raphson – L7
void N_Raph(float x1,float epsif,float epsix,int nmax,float func(float),float dfunc(float), char* file, float &xr, int &kod)
{
    int iter;
    float fx1,fx2,dfx1;
    float x2;
    ofstream fwrite(file,ios::out);
    fwrite<<"Procedura Newton - Raphson\n\nIter \t\t x1/fx1 \t\t x2/fx2\n";
    kod=0;
    iter=1;
    fx1=func(x1);
    dfx1=dfunc(x1);
    x2=x1-fx1/dfx1;
    fx2=func(x2);

    while((fabs (fx2)> epsif) && (fabs (x1-x2)> epsix))
    {
        fwrite<<"\n"<<iter<<"\t\t"<<x1<<"\t\t"<<x2<<"\n"<<"\t\t"<<fx1<<"\t\t"<<fx2;
        x1=x2;
        fx1=fx2;
        dfx1=dfunc(x1);
        x2=x1-fx1/dfx1;
        fx2=func(x2);
        iter=iter+1;
        if(iter==nmax)
        {
            kod=1;
            fwrite<<"\nDepasirea numarului de iteratii";
            goto a1;
        }
        if(fabs((x1-x2)/x1)<1e-9)
        {
            kod=1;
            fwrite<<"\nS-a atins limita masinii de calcul";
            goto a1;
        }
    }
    xr=x2;
    fwrite<<"\n\nSolutia ecuatiei = "<<xr;
    fwrite.close();
    a1:
        cout;
}

float f(float x)
{
    return 4*pow(x,3)-8*pow(x,2)+12*x-8;
}

float df(float x)
{
    return 12*pow(x,2)-16*x+12;
}

Polinoame de interpolare (Cebisev) – L8
double t0,t1,tk,x;
    int n,k;
    cout<<"Introduceti datele polinomului: \n";
    cout<<"x: ";cin>>x;
    cout<<"Dati gradul polinomului Cebisev - n: ";cin>>n;
    t0=1.0;
    t1=x;
    for(k=2;k<=n;k++)
    {
        tk=1.0/k*(2*x*t1-t0);
        t0=t1;
        t1=tk;
    }
    cout<<"Polinomul Cebisev este: \n";
    cout<<"\t\t\t"<<tk;

Polinoame de interpolare (Lagrange) – L8
float lagrange(float xi[],float yi[],int ni, float x)
{
    int i,j;
    float p,y;
    for(i=1;i<=ni;i++)
    {
        p=1;
        for(j=1;j<=ni;j++)
            if(j!=i)
                p=p*(x-xi[j])/(xi[i]-xi[j]);
        y=y+p*yi[i];
    }
    return y;
}

    cout<<"n: ";cin>>ni;
    cout<<"Vector x: "<<endl;
    for(int i=1;i<=ni;i++)
    {
        cout<<"x["<<i<<"]: ";
        cin>>xi[i];
    }
    cout<<"Vector y: "<<endl;
    for(int i=1;i<=ni;i++)
    {
        cout<<"y["<<i<<"]: ";
        cin>>yi[i];
    }
    cout<<"Dati valoarea punctului pt. care se calculeaza polinomul: ";cin>>val;
    rezultat=lagrange(xi,yi,ni,val);
    cout<<"Valoarea pentru "<<val<<" este "<<rezultat;

Polinoame de interpolare (Laguerre – mai multe valori x) – L8
double x[100],y[100];
    double p0,p1,pk;
    int m,n,k,i;
    cout<<"Introduceti datele polinomului : \n";
    cout<<"Introduceti gradul polinomului Laguerre - n: ";
    cin>>n;
    cout<<"Introduceti numarul de puncte in care se face calculul m: ";
    cin>>m;
    for(int i=1;i<=m;i++)
    {
        cout<<"Punctul x_"<<i<<": ";cin>>x[i-1];
    }
    cout<<"-------------- Rezultate --------------\n";
    for(int i=1;i<=m;i++)
    {
        p0=1.0;
        p1=1-x[i-1];
        for(k=2;k<=n;k++)
        {
            pk=1.0/k*((k*2-1)*p1-(k-1)*p0);
            p0=p1;
            p1=pk;
        }
        y[i-1]=pk;
        cout<<x[i-1]<<"\t\t"<<y[i-1]<<"\n";
    }

Polinoame de interpolare (Laguerre – o singura valoare x) – L8
double p0,p1,pk,x;
    int n,k;
    cout<<"Introduceti datele polinomului: \n";
    cout<<"x: ";cin>>x;
    cout<<"Dati gradul polinomului n: ";cin>>n;
    p0=1.0;p1=1-x;
    for(k=2;k<=n;k++)
    {
        pk=1.0/k*((k*2-1)*p1-(k-1)*p0);
        p0=p1;
        p1=pk;
    }
    cout<<"Valoarea polinomului Laguerre este: \n";
    cout<<"\t\t\t\t"<<pk;

Factorizare QR – Probleme Suplimentare
double a[10][10], q[10][10], r[10][10], h[10][10], sigma, beta, v[10], s;
int i, j, k, n, m;

void initializare_mat(double a[10][10])
{
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			a[i][j] = 0;

}

void mat_unitate(double a[10][10])
{
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			if (i == j)
				a[i][j] = 1;
			else
				a[i][j] = 0;
}

void produs(double a[10][10], double b[10][10], double c[10][10])
{
	double ct[10][10];

	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
		{
			ct[i][j] = 0;
			for (int k = 1; k <= n; k++)
				ct[i][j] += a[i][k] * b[k][j];
		}
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
			c[i][j] = ct[i][j];

}

int main()
{
	cout << "\nDimensiune mat. A: "; cin >> n;

	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
		{
			cout << "a[" << i << "][" << j << "]= ";
			cin >> a[i][j];
		}
	mat_unitate(q);
	initializare_mat(r);
	for (k = 1; k <= n - 1; k++)
	{
		s = 0;
		for (i = k; i <= n; i++)
			s += a[i][k] * a[i][k];
		sigma = sqrt(s);
		if (sigma != 0)
		{
			if (a[k][k] < 0)
				sigma = -sigma;
			v[k] = a[k][k] + sigma;
			for (i = 1; i <= k - 1; i++)
				v[i] = 0;
			for (i = k + 1; i <= n; i++)
				v[i] = a[i][k];
			s = 0;
			for (i = k; i <= n; i++)
				s += v[i] * v[i];
			beta = s / 2;
			for (i = 1; i <= n; i++)
				for (j = 1; j <= n; j++)
					if (i != j)
						h[i][j] = -v[i] * v[j] / beta;
					else  h[i][j] = 1 - v[i] * v[j] / beta;

			produs(h, a, a);
			produs(q, h, q);
		}
	}

	cout << "\nMatricea Q:\n\n";
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
			cout << setprecision(3) << fixed << "\t" << q[i][j] << " ";
		cout << "\n";
	}

	cout << "\nMatricea R:\n\n";
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
			cout << setprecision(3) << fixed << "\t" << a[i][j] << " ";
		cout << "\n";
	}

	return 0;
}
Factorizare LR – Probleme Suplimentare
int main()
{
	int i, j, n, k;
	double a[10][10], L[10][10], R[10][10];
	cout << "\nDati dimensiunea matricei A:";
	cout << "\nn="; cin >> n;

	cout << "\nDati matricea A:\n";
		for (i = 1; i <= n; i++)
			for (j = 1; j <= n; j++) 
			{
				cout << "a[" << i << "][" << j << "]="; cin >> a[i][j];
			}

		for (k = 1; k <= n - 1; k++) 
		{
			if (a[k][k] != 0) 
			{
				for (i = k + 1; i <= n; i++) 
				{
					for (j = k + 1; j <= n; j++)
						a[i][j] -= a[i][k] * a[k][j] / a[k][k];
					a[i][k] = a[i][k] / a[k][k];
				}
			}
			else
			{
				cout << "Matricea A nu poate fi factorizata LR cu acest algoritm.";
				exit(1);
			}
		}

		cout << "\nMatricea L este:\n";
		for (i = 1; i <= n; i++)
			for (j = 1; j <= n; j++)
				L[i][j] = 0;

		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= i; j++)
				L[i][j] = a[i][j];
			L[i][i] = 1;
		}
		cout << "\n";
		for (i = 1; i <= n; i++) 
		{
			for (j = 1; j <= n; j++)
				cout << L[i][j] << "  ";
			cout << "\n";
		}

		cout << "\nMatricea R este:\n";
		for (i = 1; i <= n; i++)
			for (j = 1; j <= n; j++)
				R[i][j] = 0;

		for (i = 1; i <= n; i++)
			for (j = i; j <= n; j++)
				R[i][j] = a[i][j];
		cout << "\n";
		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
				cout << R[i][j] << "  ";
			cout << "\n";
		}


		return 0;
}