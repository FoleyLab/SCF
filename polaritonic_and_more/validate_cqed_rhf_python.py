from __future__ import print_function

"""
A reference implementation of cavity quantum electrodynamics 
configuration interactions singles.
"""

__authors__   = ["Jon McTague", "Jonathan Foley"]
__credits__   = ["Jon McTague", "Jonathan Foley"]

__copyright_amp__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__   = "BSD-3-Clause"
__date__      = "2021-01-15"

# ==> Import Psi4, NumPy, & SciPy <==
import psi4
import numpy as np
import scipy.linalg as la
import time
from helper_cqed_rhf import *

# Set Psi4 & NumPy Memory Options
psi4.set_memory('2 GB')
#psi4.core.set_output_file('output.dat', False)

numpy_memory = 2

# basis set etc
psi4.set_options({'basis':        'def2-tzvppd',
                  'scf_type':     'pk',
                  'reference':    'rhf',
                  'mp2_type':     'conv',
                  'save_jk': True,
                  'e_convergence': 1e-10,
                  'd_convergence': 1e-10})

NaF_string = """

0 1
    NA           0.000000000000     0.000000000000    -0.875819904077
    F            0.000000000000     0.000000000000     1.059820520433
no_reorient
#nocom
symmetry c1
"""

NaCl_string = """

0 1
    NA           0.000000000000     0.000000000000    -1.429419641344
    CL           0.000000000000     0.000000000000     0.939751385626
no_reorient
#nocom
symmetry c1
"""

expected_NaF =  -261.371070718358
expected_NaCl = -621.438985539266

# electric field
Ex = 0.
Ey = 0.
Ez = 0.01

lam = np.array([Ex, Ey, Ez])
# run cqed_rhf on NaF and compare to expected answer
cqed_rhf_dict = cqed_rhf(lam, NaF_string)
em_cqed_rhf_e = cqed_rhf_dict['cqed_rhf_energy']
em_rhf_e = cqed_rhf_dict['rhf_energy']
assert np.isclose(em_cqed_rhf_e, expected_NaF,5e-5)
print("def2-tzvppd RHF energy of NaF:               ", em_rhf_e)
print("def2-tzvppd CQED-RHF energy of NaF:          ", em_cqed_rhf_e)
print("reference def2-tzvppd CQED-RHF energy of NaF:", expected_NaF)

cqed_rhf_dict = cqed_rhf(lam, NaCl_string)
em_rhf_e = cqed_rhf_dict['rhf_energy']
em_cqed_rhf_e = cqed_rhf_dict['cqed_rhf_energy']
assert np.isclose(em_cqed_rhf_e, expected_NaCl,5e-5)
print("def2-tzvppd RHF energy of NaCl:               ", em_rhf_e)
print("def2-tzvppd CQED-RHF energy of NaCl:          ", em_cqed_rhf_e)
print("reference def2-tzvppd CQED-RHF energy of NaCl:", expected_NaCl)
# change to the cc-pVDZ basis set
psi4.set_options({'basis':        'cc-pVDZ'})
# template for the z-matrix for MgH+
mol_tmpl = """Mg
H 1 **R**
symmetry c1
1 1"""
# array of bondlengths for MgH+
#r_array = np.array([1.0, 1.1, 1.2, 1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3])
r_array = np.array([1.0, 1.1, 1.2, 1.3,1.4,1.5,1.6])
# lambda arrays
l_big = np.array([0.0, 0.0, 0.1])
l_med = np.array([0.0, 0.0, 0.01])

# empty arrays for the energies 
# for "big lambda" -> (0.1, 0.1, 0.1)
cqed_energy_array_l_big = np.zeros_like(r_array)
# for "medium lambda" -> (0.01, 0.01, 0.01)
cqed_energy_array_l_med = np.zeros_like(r_array)
rhf_energy_array = np.zeros_like(r_array)
cqed_one_array_med = np.zeros_like(r_array)
cqed_two_array_med = np.zeros_like(r_array)
cqed_dipole_energy_array_med = np.zeros_like(r_array)
cqed_1e_dipole_energy_array_med = np.zeros_like(r_array)
cqed_quadrupole_energy_array_med = np.zeros_like(r_array)

cqed_one_array_big = np.zeros_like(r_array)
cqed_two_array_big = np.zeros_like(r_array)
cqed_dipole_energy_array_big = np.zeros_like(r_array)
cqed_1e_dipole_energy_array_big = np.zeros_like(r_array)
cqed_quadrupole_energy_array_big = np.zeros_like(r_array)

# loop over the different bond-lengths, create different instances
# of HF molecule
ctr = 0
for r in r_array:
    molstr = mol_tmpl.replace("**R**", str(r))
    print(molstr)
    med_dict = cqed_rhf(l_med, molstr)
    big_dict = cqed_rhf(l_big, molstr)
    rhf_energy_array[ctr] = med_dict['rhf_energy']
    cqed_energy_array_l_med[ctr] = med_dict['cqed_rhf_energy']
    cqed_energy_array_l_big[ctr] = big_dict['cqed_rhf_energy']
    cqed_one_array_med[ctr] = med_dict['One Electron Energy Contribution']
    cqed_two_array_med[ctr] = med_dict['Two Electron Energy Contribution']
    cqed_dipole_energy_array_med[ctr] = med_dict['Nuclear Dipolar Energy']
    cqed_1e_dipole_energy_array_med[ctr] = med_dict['1 e- Dipole Energy Contribution']
    cqed_quadrupole_energy_array_med[ctr] = med_dict['Quadrupole Energy Contribution']

    cqed_one_array_big[ctr] = big_dict['One Electron Energy Contribution']
    cqed_two_array_big[ctr] = big_dict['Two Electron Energy Contribution']
    cqed_dipole_energy_array_big[ctr] = big_dict['Nuclear Dipolar Energy']
    cqed_1e_dipole_energy_array_big[ctr] = big_dict['1 e- Dipole Energy Contribution']
    cqed_quadrupole_energy_array_big[ctr] = big_dict['Quadrupole Energy Contribution']
    ctr+=1

from matplotlib import pyplot as plt
plt.plot(r_array, rhf_energy_array, 'r-o', label="$\lambda$=(0,0,0) a.u.")
plt.plot(r_array, cqed_energy_array_l_med, 'b-*', label="$\lambda$=(0, 0, 0.01) a.u.")
plt.plot(r_array, cqed_energy_array_l_big, 'p--', label="$\lambda$=(0, 0, 0.1) a.u.")
plt.xlabel("Bondlength (Angstroms)")
plt.ylabel("Energy (Hartrees)")
plt.legend()
plt.show()

print(cqed_one_array_med,'\n',cqed_two_array_med,'\n', cqed_dipole_energy_array_med, '\n', cqed_1e_dipole_energy_array_med, '\n', cqed_quadrupole_energy_array_med)
print('med energy array', cqed_energy_array_l_med)
print('Summed med energies: ', (cqed_one_array_med[-1] +  cqed_two_array_med[-1] + cqed_dipole_energy_array_med[-1] + cqed_1e_dipole_energy_array_med + cqed_quadrupole_energy_array_med))
print('big energy array', cqed_energy_array_l_big)
print('Summed big energies: ', (cqed_one_array_big[-1] +  cqed_two_array_big[-1] + cqed_dipole_energy_array_big[-1] + cqed_1e_dipole_energy_array_big + cqed_quadrupole_energy_array_big))

