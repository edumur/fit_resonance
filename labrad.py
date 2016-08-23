# This Python file uses the following encoding: utf-8

# labrad readout file

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
from ConfigParser import SafeConfigParser


from tools import tools

class labrad(tools):



    def __init__(self, full_name):
        """
        Initialize the class by reading all the information in the labrad file.

        Parameters
        ----------
        full_name : str
            Full path of the file.
        """

        self._full_name = full_name

        self._frequency_unit = None

        #We get parameters of the file and we store them into private attributes
        self._get_parameter()

        #We get data and store them into private attribute
        self._get_data()



    def get_name(self):
        """
        Return the file name

        Return
        ----------
        filename : str
            Name of the file
        """

        return self._full_name.split('/')[-1].split('.')[0]



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

        s = SafeConfigParser()

        s.read(self._full_name[:-3]+'ini')

        self._frequency_unit = s.get('Independent 1', 'units')

        self._parameters = {'s11', 's12', 's21', 's22'}
        self._units = []
        self._measured = []

        # We browse all sections and options of the file
        for section in s.sections():
            # We select only the Dependent sections
            if section.split()[0] == 'Dependent':
                # We browse all the options
                for option in s.options(section):

                    if option == 'units':
                        self._units.append(s.get(section, option).lower())
                    elif option == 'category':
                        self._measured.append(s.get(section, option).lower())

        self._unit_mag = self._units[0]
        self._unit_phase = self._units[-1]



    def _get_data(self):
        """
        Extract data from the s2p file and format them in a ndarray
        """

        self._data = pd.read_csv(self._full_name, header=None).get_values().transpose()



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
        m = self._measured.index(s) + 1 # Return position of the magnitude data

        # Depending on the file format and the asked format we return data
        if self._unit_mag == 'ma' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(d[m]),
                    d[m+4])

        elif self._unit_mag == 'db' and data_format.lower() == 'ma':

            return (d[0]*f,
                    self.db2ma(d[m]),
                    d[m+4])

        elif self._unit_mag == 'db' and data_format.lower() == 'ri':

            if self._unit_phase == 'rad':

                return (d[0]*f,
                        np.cos(d[m+4])*self.db2ma(d[m]),
                        np.sin(d[m+4])*self.db2ma(d[m]))
            else:

                return (d[0]*f,
                        np.cos(np.radians(d[m+4]))*self.db2ma(d[m]),
                        np.sin(np.radians(d[m+4]))*self.db2ma(d[m]))

        elif self._unit_mag == 'ma' and data_format.lower() == 'ri':

            if self._unit_phase == 'rad':

                return (d[0]*f,
                        np.cos(d[m+4])*d[m],
                        np.sin(d[m+4])*d[m])

            else:

                return (d[0]*f,
                        np.cos(np.radians(d[m+4]))*d[m],
                        np.sin(np.radians(d[m+4]))*d[m])

        elif self._unit_mag == 'ri' and data_format.lower() == 'ma':

            return (d[0]*f,
                    np.sqrt(d[m]**2 + d[m+4]**2),
                    np.degrees(np.angle(d[m] + 1j*d[m+4])))

        elif self._unit_mag == 'ri' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(np.sqrt(d[m]**2 + d[m+4]**2)),
                    np.degrees(np.angle(d[m] + 1j*d[m+4])))

        else:
            return (d[0]*f,
                    d[m],
                    d[m+4])
