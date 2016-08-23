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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, RadioButtons

from s2p import s2p
from labrad import labrad



class plot(object):


    def __init__(self, full_name):

        if full_name[-3:] == 's2p':
            self.data = s2p(full_name)
        elif full_name[-3:] == 'csv':
            self.data = labrad(full_name)
        else:
            raise ValueError('File format unknown: '+full_name[:-3])



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

        x, y, z = self.data.get_SParameters(s, data_format)

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
                    if self.data._unit_phase == 'deg':
                        z = np.unwrap(np.deg2rad(z))
                    else:
                        z = np.unwrap(z)
                else:
                    if self.data._unit_phase == 'deg':
                        z = np.deg2rad(z)
            if phase_unit == 'deg':
                if unwrap_phase:
                    if self.data._unit_phase == 'deg':
                        z = np.rad2deg(np.unwrap(np.deg2rad(z)))
                    else:
                        z = np.rad2deg(np.unwrap(z))
                else:
                    if self.data._unit_phase == 'rad':
                        z = np.rad2deg(z)


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

                x, y, z = self.data.get_SParameters(s, data_format)

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

                x, y, z = self.data.get_SParameters(s, data_format)

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

            x, y, z = self.data.get_SParameters(s, data_format='ri')

            if x[0]<a*1e9 and x[-1]>b*1e9:
                x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))
                line_ri.set_data(y, z)


                phase_func('a')
                mag_func('a')

                fig.canvas.draw()

        def phase_func(label):

            z = line_phase.get_ydata()

            a, b = ax3.get_xlim()

            x, y, z = self.data.get_SParameters(s, data_format='db')
            if self.data._unit_phase == 'rad':
                z = np.rad2deg(z)

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
                x, y, z = self.data.get_SParameters(s, data_format='db')

                if x[0]<a*1e9 and x[-1]>b*1e9:
                    x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))


                a, b = ax2.get_ylim()
                if y.min() > b or y.max() < b:
                    ax2.set_ylim(y.min()*1.05, y.max()*0.95)

                ax2.set_ylabel('Attenuation [dB]')
            else:
                x, y, z = self.data.get_SParameters(s, data_format='ma')

                if x[0]<a*1e9 and x[-1]>b*1e9:
                    x, y, z = self._data_range(x, y, z, (a*1e9, b*1e9))

                a, b = ax2.get_ylim()
                if y.min() > b or y.max() < b:
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


        x, y, z = self.data.get_SParameters(s, data_format='ri')
        line_ri, = ax1.plot(y, z)

        # We draw axes lines
        ax1.axhline(color='black')
        ax1.axvline(color='black')

        x, y, z = self.data.get_SParameters(s, data_format='db')
        f_data, f_unit = self._parse_number(x)
        if self.data._unit_phase == 'rad':
            z = np.rad2deg(z)

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
