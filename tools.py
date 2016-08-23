# This Python file uses the following encoding: utf-8

# Copyright (C) 2016 Dumur Étienne
# etienne.dumur@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



class tools(object):


    def __init__(self):
        pass



    def db2ma(self, x):
        """
            Transform dB in magnitude
        """

        return 10.**(x/20.)



    def ma2db(self, x):
        """
            Transform magnitude in dB
        """

        return 20.*np.log10(x)



    def frequency_factor(self, frequency_unit):
        """
        Give the factor of the frequency and set the attribut factor

        Parameters
        ----------
        style : str {'string', 'float'}
            Determined the style of the output

        Return
        ----------
        factor : float
            factor of the frequency
                1.  for  Hz
                1e3 for KHz
                1e6 for MHz
                1e9 for GHz
        """

        if frequency_unit == 'ghz' :
            return 1e9
        elif frequency_unit == 'mhz' :
            return 1e6
        elif frequency_unit == 'khz' :
            return 1e3
        else :
            return 1.
