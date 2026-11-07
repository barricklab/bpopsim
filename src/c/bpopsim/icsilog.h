/*
    ICSI header file of the logarithm function based on a lookup table
    icsi_log source code and fill_icsi_log_table header

    Version 0.6 beta
    Build date: November 13th, 2007

    Copyright (C) 2007 International Computer Science Institute
    1947 Center Street. Suite 600
    Berkeley, CA 94704
    
    Contact information: 
         Oriol Vinyals	vinyals@icsi.berkeley.edu
         Gerald Friedland 	fractor@icsi.berkeley.edu

    Acknowledgements:
    Thanks to Harrison Ainsworth (hxa7241@gmail.com) for his idea that
    doubled the accuracy.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef __ICSI_LOG__
#define __ICSI_LOG__

/* ICSIlog v2.0 */
inline double icsi_log_v2(const float val, float* const pTable, const unsigned precision)
{
     /* get access to float bits */
     const int* const pVal = (const int*)(&val);

     /* extract exponent and mantissa (quantized) */
     const int exp = ((*pVal >> 23) & 255) - 127;
     const int man = (*pVal & 0x7FFFFF) >> (23 - precision);

     /* exponent plus lookup refinement */
     return ((double)(exp) + pTable[man]) * 0.69314718055995f;
}

#endif	

