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


  } //end constructor

// in interpolation we use second baricentric formula and we call the denominators from other function
// so that it takes less time and they are already precompiled
double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
   if (t < -1 || t > 1 || (end - start) != _p) {
      throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
                              + " \\notin [-1, +1] OR view "
                                "into fj of wrong size: ["
                              + std::to_string(start) + ", " + std::to_string(end) + "]");
   }

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return fj[l + start];
      den += _betaj[l] / (t - _tj[l]);
   }

   double sum = 0;
   for (size_t i = start, j = 0; i <= end; i++, j++) {
      sum += poli_weight(t, j, den) * fj[i];
   }
   return sum;
}

double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{
   if (t < -1 || t > 1 || (end - start) != _p) {
      throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
                              + " \\notin [-1, +1] OR view "
                                "into fj of wrong size: ["
                              + std::to_string(start) + ", " + std::to_string(end) + "]");
   }

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) {
         double sum = 0;
         for (size_t i = start, j = 0; i <= end; i++, j++) {
            sum += fj[i] * _Dij[j][l];
         }
         return sum;
      }
      den += _betaj[l] / (t - _tj[l]);
   }

   double sum = 0;
   for (size_t i = start, j = 0; i <= end; i++, j++) {
      sum += poli_weight_der(t, j, den) * fj[i];
   }
   return sum;
}

double StandardGrid::poli_weight(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   if (fabs(t - _tj[j]) < 1.0e-15) return 1.0;

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return 0.0;

      den += _betaj[l] / (t - _tj[l]);
      // i have to compute the denominator which is in general not neccesary
      //but will be useful for later
      // waste of time if I want to build the complete polinomial
   }

   return _betaj[j] / den / (t - _tj[j]);
}

//in this case we have denominator precomputed
double StandardGrid::poli_weight(double t, size_t j, double den) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   if (fabs(t - _tj[j]) < 1.0e-15) return 1.0;

   return _betaj[j] / den / (t - _tj[j]);
}


double StandardGrid::poli_weight_der(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   
   double res=0;
   for (size_t i=0; i <= _p; i++) {
    
     if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
     // should never happen se I call the function from interpolate
     // ma se chiamo from outside potrebbe
     res += _Dij[j][i]* poli_weight(t, i);
   }
   return res;
}

double StandardGrid::poli_weight_der(double t, size_t j, double den) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   
   double res=0;
   for (size_t i=0; i <= _p; i++) {
    
     if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
     res += _Dij[j][i]* poli_weight(t, i, den);
   }
   return res;
}


//apply derivative matrix directly to the vector
// I get vector of derivatives, modified in place, returns void
 void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
 {
   
   if (end - start != _p) {
      throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply "
                                  "derivative matrix to partial vector.");
   }
   vector_d temp((end - start + 1), 0.);
   for (size_t i = 0; i <= _p; i++) {
      for (size_t j = 0, k = start; k <= end; k++, j++) {
         temp[i] += fj[k] * _Dij[j][i];
         //  fj[i] += fj[k] * _Dij[j][i]; do not do that because loop over all 
         // and it expects the fj untouched till end of computation
      }
   }

   for (size_t i = start; i <= end; i++) {
      fj[i] = temp[i - start];
   }
   // I just substitute the vector
 }

 //just return from function a vector of function evaluated on grid points
vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
   vector_d result(_p+1, 0.);
   for (size_t i = 0; i<= _p; i++){
      result[i]= fnc(_tj[i]);
   }
   return result;
}




} // namespace Chebyshev
} // namespace Interpolation