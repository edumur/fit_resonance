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
import os
from ConfigParser import ConfigParser, SafeConfigParser


from tools import Tools

class Labrad(Tools):



    def __init__(self, file_number, folder=None):
        """
        Initialize the class by reading all the information in the labrad file.

        Parameters
        ----------
        file_number : int
            The number being at the beginning of the file.
        folder : str {None}
            If precised the folder in which the file is.
            If None the folder of the main script is used
        """

        self._folder = folder
        self._file_number = file_number

        # _folder and _file_number have to be defined before _get_full_name is
        # called
        self._full_name = self._get_full_name()

        self._frequency_unit = None

        #We get parameters of the file and we store them into private attributes
        self._get_parameter()

        #We get data and store them into private attribute
        self._get_data()


    def _get_full_name(self):

        if self._folder is None:
            folder = os.getcwd()
        else:
            folder = self._folder
        # Get files in the directory
        files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

        # Get files whose the extension is csv
        files = [f for f in files if 'csv' == f.split('.')[-1]]

        # Get file whose the name corresponds to the file number
        f = [f for f in files if self._file_number == int(f.split('-')[0])]

        if len(f) < 1:
            raise ValueError('Your file_number is incorrect, no file found.')

        return f[0]



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
        nb_point : int
            Number of data points
            If the number of dimensions in greater than one, return a tuple
        """

        if self._dimension == 2:
            return len(np.unique(self._data[0])), len(np.unique(self._data[1]))
        elif self._dimension == 1:
            return len(self._data[0])



    def get_number_dimension(self):
        """
        Return the number of sweeped parameters

        Return
        ----------
        nb_dimension : int
            Number of dimension
        """

        return self._dimension



    def _get_parameter(self):
        """
        Extract parameters from the s2p file and order them in the class
        """

        s = ConfigParser()

        if self._folder is None:
            folder = os.getcwd()
        else:
            folder = self._folder

        s.read(os.path.join(folder, self._full_name[:-3]+'ini'))

        # Contain information on the sweeped parameters
        self._sweeped_labels = []
        self._sweeped_units = []

        # Contain information about the measured parameters
        self._parameters = {'s11', 's12', 's21', 's22'}
        self._units = []
        self._measured = []

        # Attribut which will contain the number dimension of the measurement
        self._dimension = 0

        # We browse all sections and options of the file
        for section in s.sections():
            # We select only the Independent sections
            if section.split()[0] == 'Independent':
                self._sweeped_labels.append(s.get(section, 'label').lower())
                self._sweeped_units.append(s.get(section, 'units').lower())

                self._dimension += 1

            # We select only the Dependent sections
            if section.split()[0] == 'Dependent':
                # We browse all the options
                for option in s.options(section):

                    if option == 'units':
                        self._units.append(s.get(section, option).lower())
                    elif option == 'category':
                        self._measured.append(s.get(section, option).lower())

        self._frequency_unit = self._sweeped_units[self._sweeped_labels.index('frequency')]

        self._unit_mag = self._units[0]
        self._unit_phase = self._units[-1]



    def _get_data(self):
        """
        Extract data from the s2p file and format them in a ndarray
        """

        if self._folder is None:
            folder = os.getcwd()
        else:
            folder = self._folder

        os.path.join(folder, self._full_name[:-3]+'ini')

        self._data = pd.read_csv(os.path.join(folder, self._full_name),
                                 header=None).get_values().transpose()



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

        # Depending on the file format and the asked format we get data
        if self._unit_mag == 'ma' and data_format.lower() == 'db':

            a = self.ma2db(d[m])
            b = d[m+4]

        elif self._unit_mag == 'db' and data_format.lower() == 'ma':

            a = self.db2ma(d[m]),
            b = d[m+4]

        elif self._unit_mag == 'db' and data_format.lower() == 'ri':

            if self._unit_phase == 'rad':

                a = np.cos(d[m+4])*self.db2ma(d[m])
                b = np.sin(d[m+4])*self.db2ma(d[m])
            else:

                a = np.cos(np.radians(d[m+4]))*self.db2ma(d[m]),
                b = np.sin(np.radians(d[m+4]))*self.db2ma(d[m])

        elif self._unit_mag == 'ma' and data_format.lower() == 'ri':

            if self._unit_phase == 'rad':

                a = np.cos(d[m+4])*d[m]
                b = np.sin(d[m+4])*d[m]

            else:

                a = np.cos(np.radians(d[m+4]))*d[m]
                b = np.sin(np.radians(d[m+4]))*d[m]

        elif self._unit_mag == 'ri' and data_format.lower() == 'ma':

            a = np.sqrt(d[m]**2 + d[m+4]**2)
            b = np.degrees(np.angle(d[m] + 1j*d[m+4]))

        elif self._unit_mag == 'ri' and data_format.lower() == 'db':

            a = self.ma2db(np.sqrt(d[m]**2 + d[m+4]**2))
            b = np.degrees(np.angle(d[m] + 1j*d[m+4]))

        else:
            a = d[m]
            b = d[m+4]

        # Depending of the number of dimensions we return data
        if self._dimension == 2:

            # Find index for sweeped parameters
            l  = self._sweeped_labels.index('frequency')
            temp = np.array([0, 1])
            ll = temp[temp!=l][0]

            # Return 1D array for sweeped parameters and 2D array for measured
            # one.
            n = len(np.unique(d[l]))
            p = d[ll][::n]
            f = d[l][:n]*f
            a = np.reshape(a, (len(p), len(f)))
            b = np.reshape(b,  (len(p), len(f)))

            return p, f, a, b
        elif self._dimension == 1:

            return d[0]*f, a, b
