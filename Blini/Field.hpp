//
//  Field.hpp
//  Blini
//
//  Created by ian.moss on 06/07/2023.
//

#ifndef Field_hpp
#define Field_hpp

#include <stdio.h>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <time.h>

namespace Blini
{
typedef std::complex<double> Complex;
typedef std::vector<Complex> CVector;
typedef std::vector<CVector> Matrix;
typedef std::vector<Matrix> Cubic;
typedef std::vector<Cubic> Tensor;
Complex const I=Complex(0.0,1.0);
enum run_t {atoms,levels,states};
enum direction {C,N,SE,SW};
int const max_attempts = 1000000;

struct Dice {
public:
    double Shake(void) {return (double)rand()/RAND_MAX;}
    int Shake(int n) {return rand() % n;}
};

struct Point{
private:
    int n;
    int m;
public:
    Point() : n(0), m(0) {};
    Point(int a, int b) :  n(a), m(b) {};
    Point NE(void) { return Point(n,m+1);}
    Point NW(void) { return Point(n-1,m+1);}
    Point W(void) { return Point(n-1,m);}
    Point SW(void) { return Point(n,m-1);}
    Point SE(void) { return Point(n+1,m-1);}
    Point E(void) { return Point(n+1,m);}
    int getn()  const { return n;}
    int getm()  const { return m;}
    Point Rand(int n) {return Point(rand() % n,rand() % n);}
    bool operator==(const Point &rhs) {return n==rhs.n && m==rhs.m;}
    void Print(std::ostream& os) {os << std::setw(12) << std::left << n << m;}
};

class List{
protected:
    std::vector<Point> hexList;         //Hexagons
    std::vector<Point> AList;           //type A atom lattice
    std::vector<Point> BList;           //type B atom lattive
    std::vector<double> xList;          //atom x positions
    std::vector<double> yList;          //atom y positions
public:
    void Grow(Point p0, int maxSize, double prob);      //Grow the hexagons
    void MakeLattice(void);                             //List the atoms
    unsigned long AddHex(Point p);
    void AddA(Point p);
    void AddB(Point p);
    Point& A(int i) {return AList.at(i);}
    Point& B(int i) {return BList.at(i);}
    Point const& A(int i) const {return AList.at(i);}
    Point const& B(int i) const {return BList.at(i);}
    Point& Random(Dice d);
    int  HSize(void) {return (int)hexList.size();}
    int  ASize(void) {return (int)AList.size();}
    int  BSize(void) {return (int)BList.size();}
    int  NSize(void) {return ASize() + BSize();}
    double xcoord(int i) {return xList.at(i);}
    double ycoord(int i) {return yList.at(i);};
    double distance(int i,int j);
    void Subtract(List &L);                                   //Remove a list of hexagons
    direction Link(int a, int b);                             //Determine how B is linked to A
    void PrintList(std::ostream& os);
    void PrintPoints(std::ostream& os);
    
};

}


#endif /* Field_hpp */
