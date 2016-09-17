import cmath
import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol, simplify
from sympy.solvers import solve
from scipy.signal import sawtooth

# We set our unknowns: vo, vr, ir, ic and il
vo = Symbol("V_O")
vr = Symbol("V_R")
ir = Symbol("I_R")
ic = Symbol("I_C")
il = Symbol("I_L")   
r = Symbol("R")          # resistance
omega = Symbol("\omega") # angular frequency
c = Symbol("C")          # capacitance
l = Symbol("L")          # inductance

eq1 = (vr + vo - 1, 
       ir - ic - il, 
       vr - ir*r,
       # what does 1j mean? (1j represents the imaginary number i)
       vo - ic/(1j*omega*c),
       vo - 1j*omega*l*il)
# what does the following line do?
sol = solve(eq1, (vo, vr, ir, ic, il))
vos = simplify(sol[vo])
"""
The system of equations is solved for the 5 specified variables.
Subsequently, the variable vo is simplified.
"""
# compare the output of the following line if vos = sol[vo]
print vos
#vos = sol[vo]
print sol[vo]
"""
When 'vos=sol[vo]' instead of 'vos=simplify(sol[vo])', the output is still
the same because the expression is already simplified after solving.
"""

numvalue = {c: 10**-6, l: 10**-3}
# what does subs() do? is vos.subs(c=10**-6, l=10**-3) allowed? Try it.
vosnum = vos.subs(numvalue)
flist = [vosnum.subs({r: 100.0*3**s}) for s in range(0, 4)]
omega_axis = np.linspace(20000, 43246, 100)
"""
'subs()' substitutes a variable expression with an actual number.
Using 'numvalue = {c: 10**-6, l: 10**-3}', the result of the substitution
is i*omega/1000.
"""
# what does 121 in the following line mean?
# what are the other possible parameters of subplot()?
plt.subplot(121)
"""
'121' indicates the position of the subplot. The three numbers represent
numrows, numcols, fignum. fignum ranges from 1 to numrows*numcols. 
"""
# describe (python type, dimensions, etc) of the input parameter/s of zip() below
# what does zip(*a) do if a is a 2-D list or numpy array?
plt.plot(omega_axis, zip(*[[abs(f.subs({omega: o})) for o in omega_axis] 
                                                    for f in flist]))
"""
The input parameter is a 2-D list. It is a list of absolute omega values
within another list.
When there is an asterisk before the arguments (e.g. zip(*a) where a is 2-D),
it outputs a list of tuples.
"""
plt.xlim(20000, 43246)
plt.ylim(0, 1)
plt.xlabel('$\omega$')
plt.ylabel('$V_O$')
plt.xticks([20000, 30000, 40000])
# Replicate Fig. 2.6, right pane following the code for Fig. 2.6, left pane
"""Code shown below."""
plt.subplot(122)
plt.plot(omega_axis, zip(*[[cmath.phase(f.subs({omega: o})) for o in omega_axis] 
                                                            for f in flist]))
plt.xlim(20000, 43246)
plt.ylim(-1.5, 1.5)
plt.xlabel('$\omega$')
plt.ylabel('$\phi$')
plt.xticks([20000, 30000, 40000])
plt.tight_layout()
plt.savefig('fig2.6.png', dpi=300)
plt.show()

def vsaw(t, T=1.0): 
    """Output a sawtooth wave over a given time array t.
    
    In order to create a sawtooth wave, I utilized 'scipy.signal.sawtooth()'.
    The instructions call for a period T = 1.0. However, 'sawtooth()' has a
    period T = 2pi.
    To get around this, the argument for the 'sawtooth()' function must be
    multiplied by 2pi.
    """
    return sawtooth(t*2*np.pi)
    
omegares = 1./np.sqrt(np.prod(numvalue.values()))
alist = (1/np.sqrt(256)) * vsaw(np.arange(256)/256.0)
blist = np.sqrt(256) * np.fft.fft(alist)

def plot3(fac, w):
    # add a docstring for this function
    """Plots the output voltage of a sawtooth voltage input.
    
    Parameters
    ----------
    fac : Ratio of input fundamental frequency to the resonance frequency
    w : Resistor value
    
    Returns
    -------
    Output voltage V_out vs. t/T plot showing the filtered sawtooth input.
    """
    omegai = fac * omegares
    # How were the limits of arange() in the following line chosen?
    volist = np.concatenate(([complex(vosnum.subs({omega: omegai*s, r:
                                                   w}).evalf()) 
                                 for s in np.arange(1, 129)],
                             [0.0],
                             [complex(vosnum.subs({omega: omegai*s, r:
                                                   w}).evalf()) 
                                 for s in np.arange(-127, 0)]))
    vtrans = np.fft.ifft(blist * volist)
    plotlist = np.array([[(k+1)/256., vtrans[k%256]] for k in range(768)])
    plt.plot(plotlist[:,0], plotlist[:,1])
    # what does the following line do?
    plt.axhline(0, color='black')
    """'axhline()' creates a horizontal axis line."""
    # add labels
    """Labels shown below."""
    plt.xlabel('$t/T$')
    plt.ylabel('$V_O(t)$')
    fname = 'fig2.7and8_f' + str(fac) + 'r' + str(w) + '.png'
    plt.savefig(fname, dpi=300)
    plt.show()

plot3(1, 2700.0)
plot3(1/3., 200.0)
plot3(3.0, 5.0)

#eq2 = (ir * (r + 1/(1j*omega*c) + 1j*omega*l) + vo - 1,
#       ir - (1j*omega*c + 1/(1j*omega*l)) * vo)
#sol2 = # complete this line
#vos2 = simplify(sol2[vo])
#irs = simplify(sol2[ir])
## why should irs be passed to sympy.abs() before squaring?
#power = (r**2) *( sympy.abs(irs)**2)
#flist3 = [sympy.abs(vos2.subs(numvalue).subs({r: 10.0*3**s})) 
#            for s in range(0, 3)]
#omega_axis = np.linspace(10000, 70000, 1000)
#lines = # ...
## what does plt.setp() do?
#plt.setp(lines[0], lw=2)
#plt.setp(lines[1], ls='--'
## add labels and ticks
#plt.minorticks_on()
#plt.show()
#
## replicate fig. 2.10