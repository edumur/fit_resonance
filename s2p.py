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
import matplotlib.pyplot as plt

from matplotlib.widgets import CheckButtons, RadioButtons

class s2p(object):



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
        self._format   = None
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
        self._format         = last_line[2]
        self._R              = float(last_line[4])



    def get_frequency_unit(self, style='String'):
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
               : str
                     Hertz
                Kilo Hertz
                Mega Hertz
                Giga Hertz
        """

        answer = ''
        factor = 0

        if self._frequency_unit == 'ghz' :
            answer = 'Giga Hertz'
            factor = 1e9

        elif self._frequency_unit == 'mhz' :
            answer = 'Mega Hertz'
            factor = 1e6

        elif self._frequency_unit == 'khz' :
            answer = 'Kilo Hertz'
            factor = 1e3

        else :
            answer = 'Hertz'
            factor = 1.

        if style.lower() == 'string':
            return answer
        elif style.lower() == 'float':
            return factor
        else :
            raise ValueError('The style input should be "String" or "Float"')




    def get_format(self):
        """
        Give the format of the data and set the attribut format

        Return
        ----------
        format : str
            Data format
        """

        if self._format == 'db' :
            return 'Decibel-Angle'

        elif self._format == 'ma' :
            return 'Magnitude-Angle'

        elif self._format == 'ri' :
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
        f = self.get_frequency_unit('float')
        m = (self._measured.index(s) + 1)*2 # Magic operation to easily access data

        # Depending on the file format and the asked format we return data
        if self._format == 'ma' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(d[m-1]),
                    d[m])

        elif self._format == 'db' and data_format.lower() == 'ma':

            return (d[0]*f,
                    self.db2ma(d[m-1]),
                    d[m])

        elif self._format == 'db' and data_format.lower() == 'ri':

            return (d[0]*f,
                    np.cos(np.radians(d[m]))*self.db2ma(d[m-1]),
                    np.sin(np.radians(d[m]))*self.db2ma(d[m-1]))

        elif self._format == 'ma' and data_format.lower() == 'ri':

            return (d[0]*f,
                    np.cos(np.radians(d[m]))*d[m-1],
                    np.sin(np.radians(d[m]))*d[m-1])

        elif self._format == 'ri' and data_format.lower() == 'ma':

            return (d[0]*f,
                    np.sqrt(d[m-1]**2 + d[m]**2),
                    np.degrees(np.angle(d[m-1] + 1j*d[m])))

        elif self._format == 'ri' and data_format.lower() == 'db':

            return (d[0]*f,
                    self.ma2db(np.sqrt(d[m-1]**2 + d[m]**2)),
                    np.degrees(np.angle(d[m-1] + 1j*d[m])))

        else:
            return (d[0]*f,
                    d[m-1],
                    d[m])



    def _parse_number(self, x):

        power_ten = int(np.log10(abs(x[-1])))/3*3

        prefix = {-24 : 'y',
                  -21 : 'z',
                  -18 : 'a',
                  -15 : 'p',
                  -12 : 'p',
                   -9 : 'n',
                   -6 : 'µ',
                   -3 : 'm',
                    0 : '',
                    3 : 'k',
                    6 : 'M',
                    9 : 'G',
                   12 : 'T',
                   15 : 'p',
                   18 : 'E'}

        return x/10.**power_ten, prefix[power_ten]



    def _data_range(self, x, y, z, data_range):
        """
        Return a zoom in the data.
        User specify a range in frequency (a, b).
        """

        if data_range[0] < x[0] or data_range[-1] > x[-1]:
            raise ValueError('Your data range is out of measured frequency.')

        a = abs(x - data_range[0]).argmin()
        b = abs(x - data_range[1]).argmin()

        return x[a:b], y[a:b], z[a:b]



    def plot_SParameter(self, s='s21', data_format='db', grid=False,
                        data_range=None, unwrap_phase=False, phase_unit='deg'):
        """
        Plot  S parameter from data

        Parameters
        ----------
        s : str {'s11', 's12', 's21', 's22'}
            Set which parameter you want to plot
        asked_format : str {'ma', 'db', 'ri'}
            Set in which format data will be displayed
        grid : booleen
            Grid displayed
        data_range : tupple
            Start and Stop point of plotted data in frequency unit of the s2p
            file (see get_frequency_unit method).

        Return
        ----------
        matplotlib 2d figure
        """

        x, y, z = self.get_SParameters(s, data_format)

        # If user plots part of the data
        if data_range is not None:
            x, y, z = self._data_range(x, y, z, data_range)




        if data_format.lower() == 'ri':

            fig, ax1 = plt.subplots(1, 1)

            ax1.plot(y, z)

            # We draw axes lines
            ax1.axhline(color='black')
            ax1.axvline(color='black')

            ax1.set_xlabel('Real part')
            ax1.set_ylabel('Imaginary part')

        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

            # Obtain the frequency data in a good format (usually GHz)
            f_data, f_unit = self._parse_number(x)

            # We format the phase
            if phase_unit == 'rad':
                if unwrap_phase:
                    z = np.unwrap(np.deg2rad(z))
                else:
                    z = np.deg2rad(z)
            if phase_unit == 'deg':
                if unwrap_phase:
                    z = np.rad2deg(np.unwrap(np.deg2rad(z)))

            ax1.plot(f_data, y)
            ax2.plot(f_data, z)

            ax2.set_xlabel('Frequency '+f_unit+'Hz')
            ax2.set_ylabel('Phase [deg]')

            if data_format.lower() == 'ma':
                ax1.set_ylabel('Magnitude')
            elif data_format.lower() == 'db':
                ax1.set_ylabel('Attenuation [dB]')

        # Display grid or not
        if grid:
            for ax in fig.get_axes():
                ax.grid()

        fig.suptitle(s.upper())
        plt.show()



    def plot_SParameters(self, data_format='db', grid=False, data_range=None):
        """
        Plot  S parameters from data

        Parameters
        ----------
        asked_format : str {'ma', 'db', 'ri'}
            Set in which format data will be displayed
        grid : booleen
            Grid displayed
        data_range : tupple
            Start and Stop point of plotted data in frequency unit of the s2p
            file (see get_frequency_unit method).

        Return
        ----------
        matplotlib 2d figure
        """

        if data_format.lower() == 'ri':

            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

            for s, ax in zip (('s11', 's21', 's12', 's22'), (ax1, ax2, ax3, ax4)):

                x, y, z = self.get_SParameters(s, data_format)

                # If user plots part of the data
                if data_range is not None:
                    x, y, z = self._data_range(x, y, z, data_range)

                ax.plot(y, z)
                ax.set_title(s.upper())

                # We draw axes lines
                ax.axhline(color='black')
                ax.axvline(color='black')

            ax1.set_ylabel('Imaginary part')
            ax3.set_ylabel('Imaginary part')
            ax3.set_xlabel('Real part')
            ax4.set_xlabel('Real part')

        else:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)

            for s, ax in zip (('s11', 's21', 's12', 's22'), (ax1, ax2, ax3, ax4)):

                x, y, z = self.get_SParameters(s, data_format)

                # If user plots part of the data
                if data_range is not None:
                    x, y, z = self._data_range(x, y, z, data_range)

                f_data, f_unit = self._parse_number(x)

                ax.plot(f_data, y)
                ax.set_title(s.upper())

            if data_format.lower() == 'ma':
                ax1.set_ylabel('Magnitude')
                ax3.set_ylabel('Magnitude')
            elif data_format.lower() == 'db':
                ax1.set_ylabel('Attenuation [dB]')
                ax3.set_ylabel('Attenuation [dB]')

            ax3.set_xlabel('Frequency '+f_unit+'Hz')
            ax4.set_xlabel('Frequency '+f_unit+'Hz')

        # Display grid or not
        if grid:
            for ax in fig.get_axes():
                ax.grid()

        plt.show()



    def plot_interactive_SParameter(self, s='s21'):

        def zoom_func(event):

            a, b = ax3.get_xlim()

            x, y, z = self.get_SParameters(s, data_format='ri')

            if x[0]<a*1e9 and x[-1]>b*1e9:
                x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))
                line_ri.set_data(y, z)


                phase_func('a')
                mag_func('a')

                fig.canvas.draw()

        def phase_func(label):

            z = line_phase.get_ydata()

            a, b = ax3.get_xlim()

            x, y, z = self.get_SParameters(s, data_format='db')

            if x[0]<a*1e9 and x[-1]>b*1e9:
                x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))

            # If deg
            if ax5_button.circles[0].get_facecolor()[0] == 0.:
                # If unwrap
                if ax4_button.circles[1].get_facecolor()[0] == 0.:
                    z = np.rad2deg(np.unwrap(np.deg2rad(z)))
                    ax3.set_ylabel('Phase (unwrap) [deg]')
                else:
                    ax3.set_ylabel('Phase [deg]')
            else:
                # If unwrap
                if ax4_button.circles[1].get_facecolor()[0] == 0.:
                    z = np.unwrap(np.deg2rad(z))
                    ax3.set_ylabel('Phase (unwrap) [rad]')
                else:
                    z = np.deg2rad(z)
                    ax3.set_ylabel('Phase [rad]')


            f_data, f_unit = self._parse_number(x)
            line_phase.set_data(f_data, z)

            ax3.relim()
            ax3.autoscale_view()

            fig.canvas.draw()

        def mag_func(label):

            z = line_phase.get_ydata()

            a, b = ax3.get_xlim()

            # If db
            if ax6_button.circles[0].get_facecolor()[0] == 0.:
                x, y, z = self.get_SParameters(s, data_format='db')

                if x[0]<a*1e9 and x[-1]>b*1e9:
                    x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))


                ax2.set_ylim(y.min()*1.05, y.max()*0.95)
                ax2.set_ylabel('Attenuation [dB]')
            else:
                x, y, z = self.get_SParameters(s, data_format='ma')

                if x[0]<a*1e9 and x[-1]>b*1e9:
                    x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))

                ax2.set_ylim(y.min()*0.9, y.max()*1.1)
                ax2.set_ylabel('Attenuation')


            f_data, f_unit = self._parse_number(x)
            line_ma.set_data(f_data, y)


            ax2.autoscale_view()

            fig.canvas.draw()


        fig = plt.figure()
        fig.suptitle(s.upper())

        ax1 = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=2)
        ax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2)
        ax3 = plt.subplot2grid((3, 4), (1, 2), colspan=2, sharex=ax2)
        ax4 = plt.subplot2grid((3, 4), (2, 0))
        ax5 = plt.subplot2grid((3, 4), (2, 1))
        ax6 = plt.subplot2grid((3, 4), (2, 2))



        ax4.set_title('Phase wrapping')
        ax4_button = RadioButtons(ax=ax4,
                                  labels=('wrap', 'unwrap'),
                                  active=0)

        ax5.set_title('Phase unit')
        ax5_button = RadioButtons(ax=ax5,
                                  labels=('deg', 'rad'),
                                  active=0)

        ax6.set_title('Magnitude unit')
        ax6_button = RadioButtons(ax=ax6,
                                  labels=('dB', 'mag'),
                                  active=0)


        x, y, z = self.get_SParameters(s, data_format='ri')
        line_ri, = ax1.plot(y, z)

        # We draw axes lines
        ax1.axhline(color='black')
        ax1.axvline(color='black')

        x, y, z = self.get_SParameters(s, data_format='db')
        f_data, f_unit = self._parse_number(x)

        line_ma,   = ax2.plot(f_data, y)
        line_phase, = ax3.plot(f_data, z)

        plt.setp(ax2.get_xticklabels(), visible=False)

        ax1.set_xlabel('Real part')
        ax1.set_ylabel('Imaginary part')
        ax2.set_ylabel('Attenuation [dB]')
        ax3.set_ylabel('Phase [deg]')
        ax3.set_xlabel('Frequency '+f_unit+'Hz')


        fig.canvas.mpl_connect('button_release_event', zoom_func)
        ax4_button.on_clicked(phase_func)
        ax5_button.on_clicked(phase_func)
        ax6_button.on_clicked(mag_func)


        # fig.tight_layout()

        plt.show()
