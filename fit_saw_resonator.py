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
import scipy.constants as cst
import lmfit
from lmfit import Parameters
from scipy.stats import pearsonr

from tools import Tools

class FitSawResonator(Tools):



    def __init__(self, x=None, y=None, z=None):

        self.x = x
        self.y = y
        self.z = z


    @property
    def frequency_range(self):

        return (self.x[0], self.x[-1])

    @frequency_range.setter
    def frequency_range(self, (a, b)):

        if self.z is not None:
            self.z = self.z[self.x<b]
            self.z = self.z[self.x[self.x<b]>a]

        self.y = self.y[self.x<b]
        self.y = self.y[self.x[self.x<b]>a]

        self.x = self.x[self.x[self.x<b]>a]



################################################################################
#
#
#                   Model
#
#
################################################################################



    def model_phase_shift_electronic_delay(self, p, x):

        p = p.valuesdict()
        a = p['electronic_delay']
        b = p['phase_shift']

        return np.unwrap(np.angle(np.exp(-1j*(2.*np.pi*(a*x) + b))))



    def model_s21(self, p, x, style='db', phase='rad'):

        p   = p.valuesdict()
        qi  = p['qi']
        qc  = p['qc']
        qc2 = p['qc2']
        f0  = p['f0']

        d = (x - f0)/f0
        y = (-d + 1j/2./qi)/(1./4./qi/qc2 - d\
             + 1j*(d/2./qc2 + 1./2./qc + 1./2./qi))

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



################################################################################
#
#
#                   Residual
#
#
################################################################################



    def residual_phase_shift_electronic_delay(self, p, weight=None):

        if weight is None:
            weight = np.ones_like(self.x)

        residual = self.model_phase_shift_electronic_delay(p, self.x) - self.y

        return  residual/weight



    def residual_inverse_circle(self, p, weight=None):

        if weight is None:
            weight = np.ones(len(self.x)*2)

        y_model, z_model = self.model_s21(p, self.x, style='inverse')

        y_error = self.y - y_model
        z_error = self.z - z_model

        residual = np.concatenate((y_error, z_error))

        return  residual/weight



    def residual_db(self, p, weight=None):

        if weight is None:
            weight = np.ones(len(self.x)*2)

        y_model, z_model = self.model_s21(p, self.x, style='inverse')

        y_error = self.y - y_model
        z_error = self.z - z_model

        residual = np.concatenate((y_error, z_error))

        return  residual/weight



################################################################################
#
#
#                   Fit
#
#
################################################################################



    def fit_phase_shift_electronic_delay(self, p, weight=None):

        self.result = lmfit.minimize(self.residual_phase_shift_electronic_delay, p, args=(weight))

        return self.result




    def fit_inverse_circle(self, p, weight=None):

        self.result = lmfit.minimize(self.residual_inverse_circle, p, args=(weight))

        return self.result



################################################################################
#
#
#                   Plot
#
#
################################################################################



    def plot_phase_shift_electronic_delay(self, title, file_name,
                                          file_format='png', grid=True, show=False):

        # Obtain the frequency data in a good format (usually GHz)
        f_data, f_unit = self._parse_number(self.x)

        fig, ax = plt.subplots(1, 1)

        ax.plot(f_data, self.y, '.-')
        ax.plot(f_data, self.model_phase_shift_electronic_delay(self.result.params, self.x))

        ax.set_xlabel('Frequency '+f_unit+'Hz')
        ax.set_ylabel('Unwrap phase [rad]')

        if grid:
            ax.grid(which='both')

        textstr = u'phase shift={0:.2E}rad, std={1:.2E}rad\n'\
                  u'electronic delay={2:.2E}ns, std={3:.2E}ns\n'\
                  .format(self.result.params['phase_shift'].value,
                          self.result.params['phase_shift'].stderr,
                          self.result.params['electronic_delay'].value,
                          self.result.params['electronic_delay'].stderr)

        props = dict(boxstyle='round', facecolor='white', alpha=1.)
        ax.text(0.25, 1., textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        fig.suptitle(title)

        if show:
            plt.show()
        else:
            plt.savefig('{0}.{1}'.format(file_name, file_format))
            plt.close(fig)



    def plot_inverse_circle(self, title, file_name,
                            file_format='png', grid=True, show=False,
                            zoom_in=False):


        y_model, z_model = self.model_s21(self.result.params, self.x, style='inverse')

        fig, ax = plt.subplots(1, 1)

        ax.plot(self.y, self.z, '.-')


        ax.plot(y_model, z_model, '-')

        ax.set_xlabel('Re(1/s21)')
        ax.set_ylabel('Im(1/s21)')

        if grid:
            ax.grid(which='both')

        textstr = u'Qi={0:.2E}, std={1:.2E}\n'\
                  u'Qc={2:.2E}, std={3:.2E}\n'\
                  u'Qc2={4:.2E}, std={5:.2E}\n'\
                  u'f0={6:.2E}, std={7:.2E}\n'\
                  .format(self.result.params['qi'].value,
                          self.result.params['qi'].stderr,
                          self.result.params['qc'].value,
                          self.result.params['qc'].stderr,
                          self.result.params['qc2'].value,
                          self.result.params['qc2'].stderr,
                          self.result.params['f0'].value,
                          self.result.params['f0'].stderr)

        props = dict(boxstyle='round', facecolor='white', alpha=1.)
        ax.text(0.625, 1.1, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        fig.suptitle(title)

        if zoom_in:
            ax.set_xlim(0., 2.)
            ax.set_ylim(-1., 1.)

        if show:
            plt.show()
        else:
            plt.savefig('{0}.{1}'.format(file_name, file_format))
            plt.close(fig)



    def plot_s21_conf_interval2d(self, a, b, title, file_name,
                                 a_nb_point=50, b_nb_point=50,
                                 cmap=plt.cm.jet, file_format='png',
                                 grid=True, show=False):


        mini = lmfit.Minimizer(self.residual_inverse_circle, self.result.params)

        fig, ax = plt.subplots(1, 1)

        cx, cy, data_grid = lmfit.conf_interval2d(mini, self.result,
                                             a.lower(), b.lower(),
                                             a_nb_point, b_nb_point)

        cax = plt.imshow(data_grid,
                         interpolation='none',
                         origin='bottom',
                         extent=[cx[0], cx[-1], cy[0], cy[-1]],
                         aspect='auto',
                         cmap=cmap)

        cb = fig.colorbar(cax)
        cb.solids.set_rasterized(True)
        cb.solids.set_edgecolor('face')
        cb.set_label('Probability')

        ax.set_ylabel(b)
        ax.set_xlabel(a)

        ax.ticklabel_format(style='scientific', scilimits=(0,0))

        if grid:
            ax.grid(which='both', color='w')

        fig.suptitle(title)

        if show:
            plt.show()
        else:
            plt.savefig('{0}.{1}'.format(file_name, file_format))
            plt.close(fig)



################################################################################
#
#
#                   Print
#
#
################################################################################



    def print_results(self):

        lmfit.printfuncs.report_fit(self.result.params)



################################################################################
#
#
#                   Others
#
#
################################################################################


    def get_pearsonr(self, x, y):

        return pearsonr(x, y)[0]



    def power2photon_number(self, power, qi=None, qc=None, qc2=None, f0=None):
        """
        Calculate the average photon number in a microwave resonator.

        Parameters
        ----------
        power : float
            Input power at the entrance of the resonator in watt.
        qi : float
            Internal quality factor.
            If None the result of inverse circle fit is used.
        qc : float
            Coupling quality factor.
            If None the result of inverse circle fit is used.
        qc2 : float
            Second coupling quality factor.
            If None the result of inverse circle fit is used.
        f0 : float
            Resonance frequency in hertz.
            If None the result of inverse circle fit is used.

        Return
        ----------
        average photon number : float
        """

        if qi is None and qc is None and qc2 is None and f0 is None:
            qi = self.result.params['qi'].value
            qc = self.result.params['qc'].value
            qc2 = self.result.params['qc2'].value
            o0 = self.result.params['f0'].value*2.*np.pi

        return power/cst.hbar/o0**2.*qc*((1./2./qc2)**2. + (qi/qc)**2.)



    def get_L(self, qc=None, qc2=None, f0=None, z0=50.):

        if qc is None and qc2 is None and f0 is None:
            qc  = self.result.params['qc'].value
            qc2 = self.result.params['qc2'].value
            o0  = self.result.params['f0'].value*2.*np.pi

        return z0/2./o0*(qc + 2.*qc2)/qc2**2.



    def get_R(self, qi=None, qc=None, qc2=None, z0=50.):

        if qi is None and qc is None and qc2 is None:
            qi  = self.result.params['qi'].value
            qc  = self.result.params['qc'].value
            qc2 = self.result.params['qc2'].value

        return z0/2./qi/qc2**2.*np.sqrt(qc*(qc + 2.*qc2))



    def get_C(self, qc=None, qc2=None, f0=None, z0=50.):

        if qc is None and qc2 is None and f0 is None:
            qc  = self.result.params['qc'].value
            qc2 = self.result.params['qc2'].value
            o0  = self.result.params['f0'].value*2.*np.pi

        return 2.*qc2**2./z0/o0/qc



    def get_C0(self, qc2=None, f0=None, z0=50.):

        if qc2 is None and f0 is None:
            qc2 = self.result.params['qc2'].value
            o0  = self.result.params['f0'].value*2.*np.pi

        return qc2/z0/o0
