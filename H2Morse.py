import numpy as np
from scipy.constants import h, hbar, c, u
from scipy.special import factorial
from scipy.special import genlaguerre, gamma
from matplotlib import rc
import matplotlib.pyplot as plt
from simulation.shrodinger_1D.classes.class_secant import Secant_Method
from simulation.shrodinger_1D.classes.class_numerov import Numerov_Method
from simulation.shrodinger_1D.classes.class_simpson import Simpson_Method
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
rc('text', usetex=True)
# Factor for conversion from cm-1 to J
FAC = 100 * h * c
#Factor for conversion from W/cm2 to V/m
43416586692.18

class Morse:
    """A class representing the Morse oscillator model of a diatomic."""

    def __init__(self, mA, mB, we, wexe, re, Te):
        """Initialize the Morse model for a diatomic molecule.

        mA, mB are the atom masses (atomic mass units).
        we, wexe are the Morse parameters (cm-1).
        re is the equilibrium bond length (m).
        Te is the electronic energy (minimum of the potential well; origin
            of the vibrational state energies).

        """

        self.mA, self.mB = mA, mB
        self.mu = mA*mB/(mA+mB) * u
        self.we, self.wexe = we, wexe
        self.re = re
        self.Te = Te

        self.De = we**2 / 4 / wexe * FAC
        self.ke = (2 * np.pi * c * 100 * we)**2 * self.mu
        #  Morse parameters, a and lambda.
        self.a = self.calc_a()
        self.lam = np.sqrt(2 * self.mu * self.De) / self.a / hbar
        # Maximum vibrational quantum number.
        self.vmax = int(np.floor(self.lam - 0.5))

        self.make_rgrid()
        #self.V = self.Vmorse(self.r)

    def make_rgrid(self, n=1000, rmin=None, rmax=None, retstep=False):
        """Make a suitable grid of internuclear separations."""

        self.rmin, self.rmax = rmin, rmax
        if rmin is None:
            # minimum r where V(r)=De on repulsive edge
            self.rmin = self.re - np.log(2) / self.a
        if rmax is None:
            # maximum r where V(r)=f.De
            f = 0.999
            self.rmax = self.re - np.log(1-f)/self.a
        self.r, self.dr = np.linspace(self.rmin, self.rmax, n,
                                      retstep=True)
        if retstep:
            return self.r, self.dr
        return self.r

    def calc_a(self):
        """Calculate the Morse parameter, a.

        Returns the Morse parameter, a, from the equilibrium
        vibrational wavenumber, we in cm-1, and the dissociation
        energy, De in J.

        """

        return (self.we * np.sqrt(2 * self.mu/self.De) * np.pi *
                c * 100)

    def Vmorse(self, r, F):
        """Calculate the Morse potential, V(r).

        Returns the Morse potential at r (in m) for parameters De
        (in J), a (in m-1) and re (in m).

        """
        

        return self.De * (1 - np.exp(-self.a*(r - self.re)))**2-1.602176634E-19*5.143E11*F*(self.r-self.r[0])*6.24181E+18*FAC/0.00012

    def Emorse(self, v):
        """Calculate the energy of a Morse oscillator in state v.

        Returns the energy of a Morse oscillator parameterized by
        equilibrium vibrational frequency we and anharmonicity
        constant, wexe (both in cm-1).

        """
        vphalf = v + 0.5
        return (self.we * vphalf - self.wexe * vphalf**2) * FAC

    def calc_turning_pts(self, E):
        """Calculate the classical turning points at energy E.

        Returns rm and rp, the classical turning points of the Morse
        oscillator at energy E (provided in J). rm < rp.

        """

        b = np.sqrt(E / self.De)
        return (self.re - np.log(1+b) / self.a,
                self.re - np.log(1-b) / self.a)

    def calc_psi(self, v, r=None, normed=True, psi_max=1):
        """Calculates the Morse oscillator wavefunction, psi_v.

        Returns the Morse oscillator wavefunction at vibrational
        quantum number v. The returned function is "normalized" to
        give peak value psi_max.

        """

        if r is None:
            r = self.r
        z = 2 * self.lam * np.exp(-self.a*(r - self.re))
        alpha = 2*(self.lam - v) - 1
        psi = (z**(self.lam-v-0.5) * np.exp(-z/2) *
               genlaguerre(v, alpha)(z))
        psi *= psi_max / np.max(psi)
        return psi

    def calc_psi_z(self, v, z):
        alpha = 2*(self.lam - v) - 1
        psi = (z**(self.lam-v-0.5) * np.exp(-z/2) *
               genlaguerre(v, alpha)(z))
        Nv = np.sqrt(factorial(v) * (2*self.lam - 2*v - 1) /
                     gamma(2*self.lam - v))
        return Nv * psi

    def plot_V(self, ax, **kwargs):
        """Plot the Morse potential on Axes ax."""

        ax.plot(self.r*1.e10, self.V / FAC*0.00012 + self.Te*0.00012, **kwargs)

    def get_vmax(self):
        """Return the maximum vibrational quantum number."""

        return int(self.we / 2 / self.wexe - 0.5)

    def draw_Elines(self, vlist, ax, **kwargs):
        """Draw lines on Axes ax representing the energy level(s) in vlist."""

        if isinstance(vlist, int):
            vlist = [vlist]
        for v in vlist:
            E = self.Emorse(v)
            rm, rp = self.calc_turning_pts(E)
            ax.hlines(E / FAC*0.00012 + self.Te*0.00012, rm*1.e10, rp*1e10, **kwargs)

    def label_levels(self, vlist, ax):
        if isinstance(vlist, int):
            vlist = [vlist]

        for v in vlist:
            E = self.Emorse(v)
            rm, rp = self.calc_turning_pts(E)
            ax.text(s=r'$v={}$'.format(v), x=rp*1e10 + 0.6,
                    y=E / FAC + self.Te, va='center')

    def plot_psi(self, vlist, ax, r_plot=None, scaling=1, **kwargs):
        """Plot the Morse wavefunction(s) in vlist on Axes ax."""
        if isinstance(vlist, int):
            vlist = [vlist]
        for v in vlist:
            E = self.Emorse(v)
            if r_plot is None:
                rm, rp = self.calc_turning_pts(E)
                x = self.r[self.r<rp*1.2]
            else:
                x = r_plot
            psi = self.calc_psi(v, r=x, psi_max=self.we/2)
            psi_plot = psi*scaling*0.00012 + self.Emorse(v)/FAC*0.00012 + self.Te*0.00012
            ax.plot(x*1.e10, psi_plot, **kwargs)




COLOUR1 = (0.6196, 0.0039, 0.2588, 1.0)

# Atom masses and equilibrium bond length for (1H)(35Cl).
mA, mB = 1., 1.
X_re = 0.74144E-10
X_Te = 0
X_we, X_wexe = 4401.213, 121.336

X = Morse(mA, mB, X_we, X_wexe, X_re, X_Te)
X.make_rgrid()
X.V = X.Vmorse(X.r,0.0)

fig, ax = plt.subplots(figsize=(6, 8))
X.plot_V(ax, color='k')
#X.V = X.Vmorse(X.r,0.05)
#X.plot_V(ax, color='r')
#X.V = X.Vmorse(X.r,0.1)
#X.plot_V(ax, color='r')

X.draw_Elines(range(2), ax)
X.draw_Elines(X.get_vmax(), ax, linestyles='--', linewidths=1)
#X.plot_psi([0, 0], ax, scaling=0.5, color=COLOUR1)
#X.label_levels([0, 0], ax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.set_ylim([-1,7])

# Atom masses and equilibrium bond length for (1H)(35Cl).
mA, mB = 1., 1.
X_re = 1.052E-10
X_Te = 124418
X_we, X_wexe = 2321.7, 66.2

X = Morse(mA, mB, X_we, X_wexe, X_re, X_Te)
X.make_rgrid()
X.V = X.Vmorse(X.r,F=0)

X.plot_V(ax, color='k')

X.draw_Elines(range(11), ax)
X.draw_Elines(X.get_vmax(), ax, linestyles='--', linewidths=1)
#X.plot_psi([0, 10], ax, scaling=0.5, color=COLOUR1)
#X.label_levels([0, 10], ax)
#
ax.set_xlabel(r'$\mathrm{Internuclear\ distance \ R}\ (\mathrm{\AA}$)')
ax.set_ylabel(r'Energy (eV)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlim([0, 3])
plt.tight_layout()
plt.savefig('morse-psi.png',dpi=300,transparent=True)
plt.show()
