# Created by and Copyright (C)  Lamin Juwara (McGill)(2017/18)
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# sample size calculation from incprecision

ss_calc_precision <- function(I = 0.017, RSE_I = 'out', PrevH = 0.1,
                    CR = 1, MDRI = 211, RSE_MDRI = 0.05, FRR = 0.009, RSE_FRR = 0.2,
                    BigT = 720, n = 5000, step = 1){
  
  temp<-incprecision(I = I, RSE_I = RSE_I, PrevH = PrevH,
                     CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, 
                     RSE_FRR = RSE_FRR,
                     BigT = BigT, n = n, step = step)
  rse<-temp$RSE_I
  return(rse)
}
#ss_calc_precision()

