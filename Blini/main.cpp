//
//  main.cpp
//  Blini
//
//  Created by ian.moss on 06/07/2023.
//
//Remember to set paths to includes and library files
//in build settings and bulid phases (sift command G)
//e.g /opt/homebrew/Cellar/armadillo/12.4.1/lib
//

#include "ConfigFile.hpp"
#include "Constants.hpp"
#include "Field.hpp"
//#define ARMA_USE_HDF5
#include "armadillo"
#include <thread>

using namespace Blini;
using namespace arma;

double V(double r) {
    return (a/(chi*r))*exp(-r/TF);
}
void JThread(int i,List& W,cx_mat& vec,Tensor& J,int N) {
    
    for(int j=0;j<i;j++)
        for(int k=0;k<N;k++)
            for(int l=0;l<k;l++) {
                Complex Jc = 0.0;
                for(int p=0;p<W.NSize();p++)
                    for(int q=0;q<p;q++) {
                        Complex Omega1 = vec(p,i) * vec(q,j) - vec(p,j) * vec(q,i);
                        Complex Omega2 = vec(p,k) * vec(q,l) - vec(p,l) * vec(q,k);
                        double r = W.distance(p,q);
                        Jc += 0.25 * conj(Omega1) * Omega2 * V(r);
                    }
                J[i][j][k][l] = Jc;
            }
}

void MakeWafer(int Num,double f, int top, int hole,List& Wafer) {
    Dice dice;
    
    Wafer.Grow(Point(0,0),Num+top*hole,f);
    for(int i=0;i<top;i++) {
        List Hole;
        Hole.Grow(Wafer.Random(dice),hole,f);
        Wafer.Subtract(Hole);
    }
    Wafer.MakeLattice();
    std::cout << "Wafer generated with " << Wafer.HSize() << " hexes and " << Wafer.NSize() << " atoms" << std::endl;
}
void MakeHamiltonian(List& Wafer,sp_cx_mat& H,double B) {
    int ASize = Wafer.ASize();
    int BSize = Wafer.BSize();

    for(int i=0;i<ASize;i++) {
        double x =  Wafer.xcoord(i);
        double y =  Wafer.ycoord(i);
        for(int j=0;j<BSize;j++) {
            direction d = Wafer.Link(i,j);
            if(d == N) H(i,ASize+j)=EScale * exp(-0.5*I*b*FluxUnit*B*x);
            if(d == SE) H(i,ASize+j)=EScale * exp(-0.5*I*b*FluxUnit*B*(-s3*y-c3*x));
            if(d == SW) H(i,ASize+j)=EScale * exp(-0.5*I*b*FluxUnit*B*(s3*y-c3*x));
            if(d != C) H(ASize+j,i)=conj((Complex)H(i,ASize+j));
        }
    }
}
void MakeJTensor(List& Wafer,cx_mat& zerovec,Tensor& J) {
    int Flux = (int)zerovec.n_cols;
    std::vector<std::thread> threads;
    
    for(int i=0;i<Flux;i++) threads.push_back(std::thread(JThread,i,std::ref(Wafer),std::ref(zerovec),std::ref(J),Flux));
    
    for(auto& th : threads) th.join();
    
    //for(int i=0;i<Flux;i++) JThread(i,Wafer,zerovec,J,Flux);
                                    
    std::cout << "J coefficients calculated " << std::endl;
    
    for(int i=0;i<Flux;i++)
        for(int j=0;j<i;j++)
            for(int k=0;k<Flux;k++)
                for(int l=0;l<k;l++) {
                    J[j][i][k][l] = -J[i][j][k][l];
                    J[i][j][l][k] = -J[i][j][k][l];
                    J[j][i][l][k] = J[i][j][k][l];
                }
    std::cout << "J coefficients completed " << std::endl;
}


void Levels(int Num, double f, double Bmax, int maxlevel, int top, int hole, std::ofstream& os) {
    
    List Wafer;
    //Dice dice;
    
    MakeWafer(Num,f,top,hole,Wafer);
    Wafer.PrintPoints(os);

    std::cout << " Flux=" << FluxUnit * Bmax * Wafer.HSize() * area / twopi << std::endl;

    int nvals = 10;
    double dB = Bmax/nvals;
    double B = 0.0;
    
    eigs_opts  opts;
    opts.subdim=2*maxlevel+40;
    opts.maxiter=50000;
    mat realval(maxlevel,nvals);
    int NAtom = Wafer.NSize();
    
    for(int i=0;i<nvals;i++) {
        
        sp_cx_mat H(NAtom,NAtom);
        
        MakeHamiltonian(Wafer,H,B);
        
        cx_vec eigval;
        if(!eigs_gen(eigval, H, maxlevel,1.0e-6,opts))
            throw std::runtime_error("Diagonalisation incomplete");
        vec rval = sort(real(eigval));
        for(int j=0;j<maxlevel;j++) realval(j,i)=rval(j);
        B+=dB;
    }
    
    realval.save("levels.txt",raw_ascii);
    
}
void States(int Num, double f, double B, double tol, int maxlevel, int top, int hole, std::ofstream& os) {
    
    List Wafer;
    MakeWafer(Num,f,top,hole,Wafer);

    std::cout << " Flux=" << FluxUnit * B * Wafer.HSize() * area / twopi << std::endl;
    
    int NAtom = Wafer.NSize();
    int Flux = (int)(FluxUnit * B * Wafer.HSize() * area / twopi);
    
    std::cout << "Wafer generated with " << Wafer.HSize() << " hexes and " << NAtom << " atoms";
    std::cout << " flux=" << Flux << std::endl;
    
    sp_cx_mat H(NAtom,NAtom);
    
    MakeHamiltonian(Wafer,H,B);
    
    eigs_opts  opts;
    opts.subdim=2*maxlevel+40;
    opts.maxiter=50000;
    
    cx_vec eigval;
    cx_mat eigvec;
    if(!eigs_gen(eigval, eigvec, H, maxlevel,1.0e-6,opts))
        throw std::runtime_error("Diagonalisation incomplete");
    
    uvec indices = sort_index(real(eigval));        //extract the zero levels
    cx_mat zerovec(NAtom,Flux);
    vec    zeroval(Flux);
    int izero = maxlevel/2;
    for(int j=0;j<Flux;j++) {
        uword k = indices(j+izero-Flux/2);
        zerovec.col(j) = eigvec.col(k);
        zeroval(j) = real(eigval(k));
    }
    zeroval.save("levels.txt",raw_ascii);
    std::cout << "zero levels extracted" << std::endl;
  
    mat state(NAtom,2+Flux*2);                      //save the zero levels to file
    for(int i=0;i<NAtom;i++) {
        state(i,0) = Wafer.xcoord(i);
        state(i,1) = Wafer.ycoord(i);
        for(int j=0;j<Flux;j++) {
            Complex w = zerovec(i,j);
            state(i,2*j+2) = real(w);
            state(i,2*j+3) = imag(w);
        }
    }
    state.save("state.txt",raw_ascii);
    std::cout << "zero levels saved" << std::endl;

    
    Tensor J(Flux,Cubic(Flux,Matrix(Flux,CVector(Flux,0.0))));
    
    std::cout << "J coefficients initialised" << std::endl;
    
    MakeJTensor(Wafer,zerovec,J);
    
     for(int i=0;i<Flux;i++)
         for(int j=0;j<Flux;j++)
             for(int k=0;k<Flux;k++)
                 for(int l=0;l<Flux;l++)
                     if( i!=j && k!=l ) {
                     os.precision(5);
                     os << std::setw(12) << std::left << real(J[i][j][k][l]);
                     os << std::setw(12) << std::left << imag(J[i][j][k][l]);
                     os << std::endl;
                 }
}

int main(int argc, char * const argv[]) {
    
    string runtypes[]={"atoms","states"};
    char outdefault[]="data.txt";
    char indefault[]="params.txt";
    char *infilename=(argc==2?argv[1]:indefault);
    char *outfilename=(argc==3?argv[2]:outdefault);
    std::ofstream outfile;
    clock_t ticks = clock();
    time_t then = time(NULL);
    srand((int)time(NULL));
    
    try {
        ConfigFile params(infilename);
        outfile.open(outfilename);
        if(!outfile.good()) throw std::runtime_error("cannot open output file");
        string   type=params.read<string>("runtype","atoms");
        int      NHex=params.read<int>("area",64);
        double   fuzzy=params.read<double>("fuzzy",1.0);
        double   tol=params.read<double>("tol",1.0e-4);
        int   maxlevel=params.read<int>("maxlevel",1);
        double   B=params.read<double>("field",0.0);
        int      topology=params.read<int>("topology",0);
        int      holesize=params.read<int>("hole",64);
        
        std::cout << "running " << type << " B=" << B << " NHex=" << NHex;
        std::cout << " area=" << NHex * area << std::endl;
        
        if(type=="levels") Levels(NHex,fuzzy,B,maxlevel,topology,holesize,outfile);
        if(type=="states") States(NHex,fuzzy,B,tol,maxlevel,topology,holesize,outfile);
        std::cout << "Execution time/CPU time " << difftime(time(NULL),then) << "/";
        std::cout << (clock()-ticks)/CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    catch(std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    catch( ConfigFile::file_not_found& e ) {
        std::cout << "File '" << e.filename << "' not found." << std::endl;
    }
    catch( ConfigFile::key_not_found& e ) {
        std::cout << "Parameter '" << e.key << "' not found." << std::endl;
    }
    return 0;
}
