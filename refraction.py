# /* Copyright (C) 2002-2012 Patrick Eriksson <patrick.eriksson@chalmers.se>

#    This program is free software; you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the
#    Free Software Foundation; either version 2, or (at your option) any
#    later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#    USA. */

# /*===========================================================================
#   === File description
#   ===========================================================================*/

# /*!
#   \file   refraction.cc
#   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
#   \date   2003-01-17
  
#   \brief  Functions releated to calculation of refractive index.
# */

# /*===========================================================================
#   === External declarations
#   ===========================================================================*/

#include "refraction.h"
#include <cmath>
#include "auto_md.h"
#include "matpack_complex.h"
#include "geodetic.h"
#include "interpolation.h"
#include "special_interp.h"

# extern const Numeric DEG2RAD;
# extern const Numeric RAD2DEG;
# extern const Numeric TEMP_0_C;

# /*===========================================================================
#   === The functions (in alphabetical order)
#   ===========================================================================*/

# //! complex_n_water_liebe93
# /*! 
#   Complex refractive index of liquid water according to Liebe 1993.

#   The method treats liquid water without salt. Thus, not valid below 10 GHz.
#   Upper frequency limit not known, here set to 1000 GHz. Model parameters taken
#   from Atmlab function epswater93 (by C. Maetzler), which refer to Liebe 1993
#   without closer specifications.
 
#   Temperature must be between 0 and 100 degrees Celsius.

#   The output matrix has two columns, where column 0 is real part and column 1
#   is imaginary part. And rows matches f_grid.
   
#    \param   complex_n   Out: Complex refractive index.        
#    \param   f_grid      As the WSV with the same name.
#    \param   t           Temperature

#    \author Patrick Eriksson
#    \date   2003-08-15
# */
# void complex_n_water_liebe93(Matrix& complex_n,
#                              const Vector& f_grid,
#                              const Numeric& t) {
import numpy as np
def     complex_n_water_liebe93(freq,t):  
    #INPUT 
    # freq              frequency array in GHz
    # t                 temperature in K
    #OUTPUT
    # complex_n         complex refractive index

    nf=len(freq)
    complex_n=np.empty((nf,2))
    complex_n[:]=np.nan
    
    TEMP_0_C=273.15
    if (t<(TEMP_0_C - 40 )) or    (t>  (  TEMP_0_C + 100) ):        
        return complex_n
    if (np.min(freq)<(10)) or    (np.max(freq)>  (  1000) ):        
        return complex_n   


#  // Implementation following epswater93.m (by C. Mätzler), part of Atmlab,
#  // but numeric values strictly following the paper version (146, not 146.4)
    theta = 1 - 300 / t
    e0 = 77.66 - 103.3 * theta
    e1 = 0.0671 * e0
    f1 = 20.2 + 146 * theta + 316 * theta * theta
    e2 = 3.52
    f2 = 39.8 * f1

    for iv in range(nf):
        ifGHz= (0.0+1j* freq[iv])
        n= np.sqrt(e2 + (e1 - e2) / ((1.0) - ifGHz / f2) + (e0 - e1) / ((1.0) - ifGHz / f1))

        complex_n[iv, 0] = np.real(n)
        complex_n[iv, 1] = np.imag(n)
    return complex_n


  

# //! complex_n_ice_matzler06
# /*! 
#   Complex refractive index of water ice according to Matzler 2006 (equivalent to
#   Warren 2008).

#   The method treats pure water ice (no impurities like salt). Valid from 10 MHz
#   up to 3 THz. Thus, not valid below 10 GHz. Follows the atmlab implementation,
#   including some relaxation of upper temperature limit to 280K.

#   The output matrix has two columns, where column 0 is real part and column 1
#   is imaginary part; rows match f_grid.
   
#    \param   complex_n   Out: Complex refractive index.        
#    \param   f_grid      As the WSV with the same name.
#    \param   t           Temperature

#    \author Jana Mendrok
#    \date   2016-03-21
# */
# void complex_n_ice_matzler06(Matrix& complex_n,
#                              const Vector& f_grid,
#                              const Numeric& t) {
import numpy as np
def     complex_n_ice_matzler06(freq,t):  
    #INPUT 
    # freq              frequency array in GHz
    # t                 temperature in K
    #OUTPUT
    # complex_n         complex refractive index
    nf=len(freq)
    complex_n=np.empty((nf,2))
    complex_n[:]=np.nan
    TEMP_0_C=273.15
    if (t<(20.)) or    (t>  (  280.) ):        
        return complex_n
    if (np.min(freq)<(0.01)) or    (np.max(freq)>  (  3000) ):        
        return complex_n   
    
#   chk_if_in_range("t", t, 20., 280.);
#   chk_if_in_range("min of f_grid", min(f_grid), 10e6, 3000e9);
#   chk_if_in_range("max of f_grid", max(f_grid), 10e6, 3000e9);

#   const Index nf = f_grid.nelem();

#   complex_n.resize(nf, 2);

#   // some parametrization constants
    B1 = 0.0207
    B2 = 1.16e-11
    b = 335.

    deltabeta = np.exp(-9.963 + 0.0372 * (t - 273))
    ebdt = np.exp(b / t)
    betam = (B1 / t) * ebdt / ((ebdt - 1.) * (ebdt - 1.))

    theta = 300. / t - 1
    alfa = (0.00504 + 0.0062 * theta) * np.exp(-22.1 * theta)
    reps = 3.1884 + 9.1e-4 * (t - 273)

    for iv in range( nf):
        f = freq[iv] 
        beta = betam + B2 * f * f + deltabeta
        ieps = alfa / f + beta * f

#     Complex eps(reps, ieps);
        eps=reps+1j*ieps
        n = np.sqrt(eps);
        complex_n[iv, 0] = np.real(n)
        complex_n[iv, 1] = np.imag(n)
    return complex_n