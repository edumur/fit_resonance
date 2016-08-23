# This Python file uses the following encoding: utf-8

# s2p readout file

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

import numpy as np
import pandas as pd

from tools import tools

class s2p(tools):



    def __init__(self, full_name):
        """
        Initialize the class by reading all the information in the s2p file.

        Parameters
        ----------
        full_name : str
            Full path of the file.
        """

        self._full_name = full_name

        self._frequency_unit = None
        self._unit_mag   = None
        self._data     = None
        self._R        = None

        #We get parameters of the file and we store them into private attributes
        self._get_parameter()

        #We get data and store them into private attribute
        self._get_data()



    def __str__(self):

        print 's2p class'



    def __repr__(self):

        a  = 'File name: '+str(self.get_name())+'\n'
        a += '\n'
        a += '    Format of data:    '+str(self.get_format())+'\n'
        a += '    Impedance:    '+str(self.get_impedance())+'\n'

        return a



    def get_name(self):
        """
        Return the file name

        Return
        ----------
        filename : str
            Name of the file
        """

        return self._full_name.split('/')[-1].split('.')[0]


    def get_frequency_unit(self):
        """
        Return the unit of frequency
        """

        return self._frequency_unit



    def get_number_point(self):
        """
        Return the number of point

        Return
        ----------
        nb_point : float
            Number of data points
        """

        return len(self._data[0])




    def _get_parameter(self):
        """
        Extract parameters from the s2p file and order them in the class
        """

        f = open(self._full_name, 'r')

        # These variables will contain the last description lines of the file
        last_line     = ''
        last_last_line = ''

        # A counter to know how many rows at the beginning og the file doesn't
        # contain data
        self._skiprows = 0

        # Browse the file
        for i in f.readlines():
            self._skiprows += 1
            # The # character indicate the last description line
            # We then break the reading loop
            if i[0] == '#':
                last_line = i[1:].lower().split()
                break
            else:
                last_last_line = i[1:]

        # Contain a list of measured quantities
        self._measured = last_last_line.lower().replace(',','').replace(':','').split()[3:]

        self._frequency_unit = last_line[0]
        self._parameter      = last_line[1]
        self._unit_mag         = last_line[2]
        self._R              = float(last_line[4])




    def get_format(self):
        """
        Give the format of the data and set the attribut format

        Return
        ----------
        format : str
            Data format
        """

        if self._unit_mag == 'db' :
            return 'Decibel-Angle'

        elif self._unit_mag == 'ma' :
            return 'Magnitude-Angle'

        elif self._unit_mag == 'ri' :
            return 'Real-Imaginary'
        else:
            raise ValueError('Impossible to find what is the format of the s2p file (not "DB" "MA" or "RI")')




    def get_impedance(self):
        """
        Return the matching impedance

        Return
        ----------
        impedance : float
            Value of the matching impedance
        """

        return self._R




    def _get_data(self):
        """
        Extract data from the s2p file and format them in a ndarray
        """

        self._data = pd.read_csv(self._full_name,
                           delimiter=' ',
                           skiprows=self._skiprows-1).get_values().transpose()

        if max(self._data[4]) > 4:
            self._unit_phase = 'deg'
        else:
            self._unit_phase = 'rad'



    def get_SParameters(self, s='S21', data_format='db'):
        """
        Return desired S parameter from data.

        Parameters
        ----------
        s : str {'s11', 's12', 's21', 's22'}
            Set which parameter we want to get
        data_format : str {'ma', 'db', 'ri'}
            Set in which format we get data

        Return
        ----------
        (x, y, z) : tupple
            x: Frequency in Hertz
            y: if ma magnitude, if db attenuation in dB, if ri real part
            z: if ma or db angle in degrees, if ri imaginary part
        """

        s = s.lower()

        #We check the input parameters
        if s not in ('s11', 's12', 's21', 's22'):
            raise ValueError('The argument "s" should be in the form : "S11", "S12", "S21", "S22"')

        if data_format.lower() not in ('ma', 'db', 'ri'):
            raise ValueError('The argument "data_format" should be in the form : "ma", "ri", "db"')

        #For concision we create short variable names
        d = self._data
        f = self.frequency_factor(self._frequency_unit)
        m = (self._measured.index(s) + 1)*2 # Magic operation to easily access data

        # Depending on the file format and the asked format we return data
        if self._unit_mag == 'ma' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(d[m-1]),
                    d[m])

        elif self._unit_mag == 'db' and data_format.lower() == 'ma':

            return (d[0]*f,
                    self.db2ma(d[m-1]),
                    d[m])

        elif self._unit_mag == 'db' and data_format.lower() == 'ri':

            return (d[0]*f,
                    np.cos(np.radians(d[m]))*self.db2ma(d[m-1]),
                    np.sin(np.radians(d[m]))*self.db2ma(d[m-1]))

        elif self._unit_mag == 'ma' and data_format.lower() == 'ri':

            return (d[0]*f,
                    np.cos(np.radians(d[m]))*d[m-1],
                    np.sin(np.radians(d[m]))*d[m-1])

        elif self._unit_mag == 'ri' and data_format.lower() == 'ma':

            return (d[0]*f,
                    np.sqrt(d[m-1]**2 + d[m]**2),
                    np.degrees(np.angle(d[m-1] + 1j*d[m])))

        elif self._unit_mag == 'ri' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(np.sqrt(d[m-1]**2 + d[m]**2)),
                    np.degrees(np.angle(d[m-1] + 1j*d[m])))

        else:
            return (d[0]*f,
                    d[m-1],
                    d[m])
