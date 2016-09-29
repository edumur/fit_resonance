import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('../../fit_resonance/')

from fit_superconducting_resonance import FitSuperconductingResonance as Fit
from fit_superconducting_resonance import Parameters

################
#    Data import
################

dat = pd.read_csv('dat.csv', header=None).get_values().transpose()
f = dat[0] #In GHz
s21 = 10.**(dat[3]/20.)*np.exp(1j*dat[7])

################
#    Fit electronic delay, phase shift
################

p = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
p.add_many(('phase_shift',      np.deg2rad(150),  True, -np.pi, np.pi,  None),
           ('electronic_delay', 32.105,           True,     0.,  None,  None))

# We select data out of resonance
x = np.concatenate((f[:9], f[-9:]))
y = np.concatenate((np.angle(s21)[:9], np.angle(s21)[-9:]))

# Algorithme work much better with unwrap data
y = np.unwrap(y)

# Instance fit class with the data
fm = Fit(x, y)

# Fit the data
result = fm.fit_phase_shift_electronic_delay(p)

# Print fit report
fm.print_results()

# Plot data and fit
fm.plot_phase_shift_electronic_delay(title='Fit electronic delay and phase shift',
                                     file_name='fit_phase_shift_electronic_delay',
                                     file_format='svg',
                                     grid=True)

phase_shift = result.params['phase_shift']
electronic_delay = result.params['electronic_delay']*1e-9

################
#    Find background
################

background = np.mean(20.*np.log10(abs(s21))[:9])

print 'Background: {0}dBm'.format(background)

fig, ax = plt.subplots(1, 1)

ax.plot(f, 20.*np.log10(abs(s21)) - background, '.-')

ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('S21 [dB]')

ax.grid(which='both')

textstr = u'background={0:.2E}dB'.format(background)

props = dict(boxstyle='round', facecolor='white', alpha=1.)
ax.text(0.25, 1., textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

fig.suptitle('Find background')

plt.savefig('find_background.svg')
plt.close(fig)

################
#    Fit inverse circle
################

p = Parameters()
#           (Name, Value,      Vary, Min,   Max,  Expr)
p.add_many(('qi',  3e5,        True,   0., None,  None),
           ('qc',  1e3,        True,   0., None,  None),
           ('f0',  7.58982e9,  True,   0., None,  None),
           ('phi', -0.5,       True, None, None,  None))

# Remove background, phase shift and electronic delay
s21 *= 10**(-background/20.)\
       *np.exp(1j*(2.*np.pi*f*1e9*electronic_delay + phase_shift))

re = np.real(1./s21)
im = np.imag(1./s21)

# Instance fit class with data
fm = Fit(f*1e9, re, im)

# Fit the data
result = fm.fit_inverse_circle(p)

# Print fit report
fm.print_results()

# Get pearson correlation coefficient
re_model, im_model = fm.model_s21(result.params, f*1e9, style='inverse')
pc = fm.get_pearsonr(np.concatenate((re_model, im_model)), np.concatenate((re, im)))
print 'Pearson correlation coefficient: {0}'.format(pc)

# Plot data and fit
fm.plot_inverse_circle(title='Nice inverse fit',
                       file_name='fit_inverse_circle',
                       file_format='svg',
                       grid=True)
