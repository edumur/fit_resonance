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
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, minimize

from plot import plot

class fit(plot):



    def __init__(self):

        # Fit parameters are saved as class attributes
        self.qi               = None
        self.qc               = None
        self.f0               = None
        self.phi              = None
        self.background_db    = None
        self.background_phase = None
        self.phase_delay      = None



    @property
    def frequency_range(self):

        return (self.x[0], self.x[-1])

    @frequency_range.setter
    def frequency_range(self, (a, b)):

        self.z = self.z[self.x<b]
        self.z = self.z[self.x[self.x<b]>a]

        self.y = self.y[self.x<b]
        self.y = self.y[self.x[self.x<b]>a]

        self.x = self.x[self.x[self.x<b]>a]



################################################################################
#
#
#                   Phase shift and electronic delay
#
#
################################################################################



    def model_phase_shift_electronic_delay(self, p, x):

        a, b = p

        return np.angle(np.exp(-1j*2.*np.pi*(a*x + b)))



    def fit_phase_shift_electronic_delay(self, p, x, y):

        def func(p, x, y):

            return self.model_phase_shift_electronic_delay(p, x) - y


        return leastsq(func, p, args=(x, y))[0]



    def print_phase_shift_electronic_delay(self, p):

        a, b = p

        result  = 'Results:\n'\
                  '    Electronic delay: '+str(a*1e9)+' ns\n'\
                  '    Phase shift: '+str(b)+' rad or '+str(np.rad2deg(b))+' deg'

        return result



    def plot_phase_shift_electronic_delay(self, p, x, y, grid=True):

        x_data = np.linspace(x[0], x[-1], 1e4)
        y_model = self.model_phase_shift_electronic_delay(p, x_data)

        fig, ax = plt.subplots(1, 1)

        # Obtain the frequency data in a good format (usually GHz)
        f_data, f_unit = self._parse_number(x)

        ax.plot(f_data, y, '.', label='data')
        ax.plot(x_data/1e9, y_model, '--', label='model')

        ax.set_xlabel('Frequency '+f_unit+'Hz')
        ax.set_ylabel('Phase [rad]')

        # Display grid or not
        if grid:
            ax.grid()

        fig.suptitle('Phase shift and electronic delay')
        plt.show()



################################################################################
#
#
#                    S21
#
#
################################################################################



    def model(self, p, x, style='db', phase='rad'):

        qi, qc, f0, phi, background_db, background_phase, phase_delay = p

        # Calculate S21
        dx = (x - f0)/f0
        y = 1./(1. + qi/qc*np.exp(1j*phi)/(1. + 2j*qi*dx))

        # Adding the background in db, the background phase and phase delay
        y = y*10.**(background_db/20.)\
            *np.exp(-1j*(2.*np.pi*x*phase_delay - background_phase))

        if style.lower() == 'db':

            a = 20.*np.log10(abs(y))

            if phase.lower() == 'rad':
                b = np.angle(y, deg=False)
            elif phase.lower() == 'deg':
                b = np.angle(y, deg=True)
            else:
                raise ValueError("phase argument must be: 'rad' or 'deg'.")
        elif style.lower() == 'mag':

            a = abs(y)

            if phase.lower() == 'rad':
                b = np.angle(y, deg=False)
            elif phase.lower() == 'deg':
                b = np.angle(y, deg=True)
            else:
                raise ValueError("phase argument must be: 'rad' or 'deg'.")
        elif style.lower() == 'ri':

            a = np.real(y)
            b = np.imag(y)
        elif style.lower() == 'inverse':

            a = np.real(1./y)
            b = np.imag(1./y)
        else:
            raise ValueError("style argument must be: 'db', 'mag' or 'ri'.")

        return a, b



    def fit_db_phase(self, p, x, y, z):

        def func(p, x, y, z):

            y_model, z_model = self.model(p, x, style='db')

            y_error = y - y_model
            z_error = z - z_model

            # return np.sum(np.concatenate((y_error, z_error))**2.)
            return np.concatenate((y_error, z_error))

        # return minimize(func, p, args=(self.x, self.y, self.z),
        #                 method='L-BFGS-B',
        #                 bounds=((0., None),
        #                         (0., None),
        #                         (0., None),
        #                         (None, None),
        #                         (None, None),
        #                         (None, None),
        #                         (None, None))).x
        p, cov_p, infodict, mesg, ier = leastsq(func, p, args=(x, y, z), full_output=True)

        p_err = [np.nan]*len(p)
        if cov_p is not None:
            dof = len(x) - len(p)
            chi_sq = np.sum(func(p, x, y, z)**2.)

            for i in range(len(p)):
                p_err[i] = np.sqrt(cov_p[i][i]) * np.sqrt(chi_sq / dof)

        return p, p_err



    def fit_inverse_circle(self, p, x, y, z):

        def func(p, x, y, z):

            y_model, z_model = self.model(p, x, style='inverse')

            y_error = y - y_model
            z_error = z - z_model

            # return np.sum(np.concatenate((y_error, z_error))**2.)
            return np.concatenate((y_error, z_error))

        # return minimize(func, p, args=(self.x, self.y, self.z),
        #                 method='L-BFGS-B',
        #                 bounds=((0., None),
        #                         (0., None),
        #                         (0., None),
        #                         (None, None),
        #                         (None, None),
        #                         (None, None),
        #                         (None, None))).x
        p, cov_p, infodict, mesg, ier = leastsq(func, p, args=(x, y, z), full_output=True)

        p_err = [np.nan]*len(p)
        if cov_p is not None:
            dof = len(x) - len(p)
            chi_sq = np.sum(func(p, x, y, z)**2.)

            for i in range(len(p)):
                p_err[i] = np.sqrt(cov_p[i][i]) * np.sqrt(chi_sq / dof)

        return p, p_err



    def print_s21(self, p, p_err=None):

        result  = 'Results:\n'
        result += '    Qi: {0:.2E} +- {1:.2E}\n'.format(p[0], p_err[0])
        result += '    Qc: {0:.2E} +- {1:.2E}\n'.format(p[1], p_err[1])
        result += '    f0: {0:.2E} +- {1:.2E} GHz\n'.format(p[2], p_err[2])
        result += '    phi: {0:.2E} +- {1:.2E} deg\n'.format(p[3], p_err[3])
        result += '    background_db: {0:.2E} +- {1:.2E} dB\n'.format(p[4], p_err[4])
        result += '    background_phase: {0:.2E} +- {1:.2E} deg\n'.format(p[5], p_err[5])
        result += '    phase_delay: {0:.2E} +- {1:.2E} deg'.format(p[6], p_err[6])

        return result



    def plot_db_fit(self, p, x, y, z, grid=False):

        x_data = np.linspace(x[0], x[-1], 1e4)
        y_model, z_model = self.model(p, x_data, style='db')


        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        # Obtain the frequency data in a good format (usually GHz)
        f_data, f_unit = self._parse_number(x)


        ax1.plot(f_data, y, '.')
        ax1.plot(x_data/1e9, y_model)

        ax2.plot(f_data, z, '.')
        ax2.plot(x_data/1e9, z_model)

        ax2.set_xlabel('Frequency '+f_unit+'Hz')
        ax2.set_ylabel('Phase [rad]')

        ax1.set_ylabel('Attenuation [dB]')

        # Display grid or not
        if grid:
            for ax in fig.get_axes():
                ax.grid()

        fig.suptitle('S21')
        plt.show()




    def plot_inverse_circle_fit(self, p, x, y, z, grid=False):

        x_data = np.linspace(x[0], x[-1], 1e4)
        y_model, z_model = self.model(p, x_data, style='inverse')


        fig, ax = plt.subplots(1, 1)


        ax.plot(y, z, '.')
        ax.plot(y_model, z_model, '.-')

        ax.set_xlabel('Re[1/S21]')
        ax.set_ylabel('Im[1/S21]')

        # Display grid or not
        if grid:
            ax.grid()

        fig.suptitle('S21')
        plt.show()
