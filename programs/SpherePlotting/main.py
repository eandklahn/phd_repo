import argparse
import crystallographyClasses as cc
from polarized_neutrons import _plot_peak

# CIF-navn og navn på den _data-blok i CIF-filen, som skal bruges
_cifname = 'CoClTu.cif'
_blockname = 'shelx'
_structure = cc.crystalStructure(_cifname, blockname=_blockname)

# Bragg-peak, feltstyrke og polariseringsværdi
hkl = [0,0,3]
_mag_field = 0.5
_polarization = 0.94

# Tilføj information om magnetiske atomer i strukturen.
# For eksemplet DBM er det

Dy1 = _structure.atoms[_structure.atomdict['Co1']]
Dy1.ion = 'Co2'
Dy1._type_magnetic_form = 'j0'
Dy1._angular_L = 1
Dy1._angular_S = 3/2
Dy1._angular_J = 5/2

_plot_peak(_structure, *hkl, _mag_field, _polarization)