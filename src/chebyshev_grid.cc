#include "Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{
    
    
    StandardGrid :: StandardGrid(size_t p)
    {
    _p = p;
    // T_j = cos (j pi/p) J = 0, ..p
    for (size_t j=0; j<=p; j++){
    _tj.push_back(cos(j * M_PI / static_cast<double>(p) ) );
      } // Costructer has to build everything, static_cast to make it a double
    
    
    
    //beta_j = (-1)^j * (1 if J!=0, p otherwise 1/2)
    for (size_t j = 0; j <= p; j++){
    double sign = j % 2 == 0 ? +1 : -1;
    //after question mark is what happens if it is and the : = otherwise
    // more elaborate questions are difficult to write in one line use if
    double scaling = 1 ;
    if (j == 0 || j==p) scaling = 0.5; 
    // j equal to zero OR p
    _betaj.push_back(sign * scaling);
    }
    
    //we already know the dimensions so compiler please use it
    _Dij.resize(p+1, vector_d(p+1, 0.));
    // D is a vector of dim p+1 of vectors of dim p+1 of zeros
    //we dont need to push back anymore only say at that index I want that value
    _Dij[0][0] = (2. * p*p +1) / 6;
    _Dij[p][p] = - _Dij[0][0];
    for (size_t j=1; j<p; j++){
     _Dij[j][j] = -0.5 * _tj[j] / (1 - pow(_tj[j], 2));
     //DO NOT X^2 for power, it's xor binary operation
     // NO x**2
     // multiplication or pow
    }
    //off diagonal now

    for (size_t i=0; i <= p; i++){
       for (size_t j =0; j<= p; j++){
       if (j==i) continue;
       _Dij[i][j] = - (_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
    }
    }


}

} // namespace Chebyshev
} // namespace Interpolation