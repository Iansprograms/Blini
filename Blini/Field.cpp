//
//  Field.cpp
//  Blini
//
//  Created by ian.moss on 06/07/2023.
//

#include "Field.hpp"
#include "Constants.hpp"

namespace Blini

{
unsigned long List::AddHex(Point p) {
    
    if (find(hexList.begin(), hexList.end(), p) == hexList.end())
        hexList.push_back(p);
    return hexList.size();
}
void List::AddA(Point p) {
    
    if (find(AList.begin(), AList.end(), p) == AList.end())
        AList.push_back(p);
}

void List::AddB(Point p) {
    
    if (find(BList.begin(), BList.end(), p) == BList.end())
        BList.push_back(p);
}
Point& List::Random(Dice d) {
    
    return hexList.at(d.Shake(HSize()));
    
}
void List::Grow(Point p0, int maxSize, double prob) {
    
    int attempt=0;
    Dice dice;
    
    hexList.clear();
    hexList.push_back(p0);
    
    while(hexList.size() < maxSize && attempt++<max_attempts) {
        
        unsigned long Size = hexList.size();
        for(int i = 0; i < Size; i++) {
            Point p = hexList.at(i);
            if(dice.Shake()<prob) if(AddHex(p.NE()) >= maxSize) return;
            if(dice.Shake()<prob) if(AddHex(p.E())  >= maxSize) return;
            if(dice.Shake()<prob) if(AddHex(p.NW()) >= maxSize) return;
            if(dice.Shake()<prob) if(AddHex(p.W())  >= maxSize) return;
            if(dice.Shake()<prob) if(AddHex(p.SW()) >= maxSize) return;
            if(dice.Shake()<prob) if(AddHex(p.SE()) >= maxSize) return;
            }
        }
}
void List::Subtract(List& L) {
    
    std::vector<Point>::iterator i;
    std::vector<Point>::iterator ip;
    
    for(i = L.hexList.begin(); i<L.hexList.end();i++) {
        
        ip = find(hexList.begin(), hexList.end(), *i);
        if(ip != hexList.end()) hexList.erase(ip);
               
    }
}
void List::MakeLattice(void) {
    
    std::vector<Point>::iterator ihex;
    
    for(ihex=hexList.begin(); ihex<hexList.end(); ihex++) {
        Point p=*ihex;
        AddA(p);
        AddA(p.NE());
        AddA(p.E());
        AddB(p);
        AddB(p.W());
        AddB(p.SW());
    }
    std::vector<Point>::iterator i;
    
    for(i = AList.begin(); i < AList.end(); i++) {
        
        int n=i->getn();
        int m=i->getm();
    
        xList.push_back( a * (n + c3 * m) );
        yList.push_back( a * s3 * m );
    }
    for(i = BList.begin(); i < BList.end(); i++) {
        
        int n=i->getn();
        int m=i->getm();
    
        xList.push_back( a * (n + c3 * m + 1.0) );
        yList.push_back( a * (s3 * m + c3/s3) );
    }
}
direction List::Link(int a,int b) {
    
    Point A = AList.at(a);
    Point B = BList.at(b);
    
    if(B.getn()==A.getn()-1 && B.getm()==A.getm()) return N;
    if(B.getn()==A.getn()-1 && B.getm()==A.getm()-1) return SW;
    if(B.getn()==A.getn() && B.getm()==A.getm()-1) return SE;
    
    return C;           //no link found
}
void List::PrintList(std::ostream& os) {
        
        std::vector<Point>::iterator i;
    
    os << "AB  n          " << "m" << std::endl;
        
    for(i = AList.begin(); i < AList.end(); i++) {
            
            os << "A   ";
            i->Print(os);
            os << std::endl;
        }
    for(i = BList.begin(); i < BList.end(); i++) {
            os << "B   ";
            i->Print(os);
            os << std::endl;
        }
    }
double List::distance(int i,int j) {
    return sqrt( pow(xList.at(i)-xList.at(j),2)+pow(yList.at(i)-yList.at(j),2) );
}
void List::PrintPoints(std::ostream& os) {
    
    std::vector<double>::iterator ix = xList.begin();
    std::vector<double>::iterator iy = yList.begin();
    
    os << "x          " << "y" << std::endl;
    
    for(; ix < xList.end(); ix++,iy++) {
        os << std::setw(12) << std::left << *ix;
        os << std::setw(12) << std::left << *iy;
        os << std::endl;
    }
 
}


}
