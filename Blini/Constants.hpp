//
//  Constants.hpp
//  deS
//
//  Created by Ian Moss on 19/07/2016.
//  Copyright © 2016 Ian Moss. All rights reserved.
//

#ifndef Constants_hpp
#define Constants_hpp

#include <cmath>

double const  a = 0.246;                    //nm
double const  b = a/sqrt(3.0);              //nm
double const  area = 1.5 * a * b;           //nm^2
double const  xi = 25.6556418;              //nm for 1 tesla
double const  FluxUnit = 1.0/(xi*xi);       //per Tesla per nm^2
double const  EScale = 2.7;                 //energy scale in eV
double const  TF = 50;                      //Thomas Fermi in nm
double const  chi = 2;                      //relative permeability
double const  c3 = 0.5;
double const  s3 = 0.5 * sqrt(3.0);
double const  pi = 4.0 * atan(1.0);
double const  twopi = 2.0 * pi;


#endif /* Constants_hpp */
