"""
"""
import os
import numpy as np
from numpy.testing import assert_array_equal
import pygamma

profile_filename = 'example.dat'
profile_path = os.path.join(pygamma.datadir, profile_filename)
with open(profile_path, 'rb') as profile:
    slon, slat, sn = profile.readline().split()
    lon, lat, n = float(slon), float(slat), int(sn)
    S, T, P = np.loadtxt(profile, unpack=True)

gamma, dg_lo, dg_hi = pygamma.gamma_n(S, T, P, lon, lat)

glevels = np.array([26.8, 27.9, 28.1], dtype=np.float64)

surfaces = pygamma.neutral_surfaces(glevels, S, T, P, gamma)
# Test alternative signature:
def test_alternative_signature():
    surfaces2 = pygamma.neutral_surfaces(glevels, S, T, P, lon=lon, lat=lat)
    # print("Test signatures: outputs match: ", np.all(surfaces == surfaces2))
    assert_array_equal(surfaces, surfaces2)

# Test 2-d inputs:
SS = np.row_stack((S,S))
TT = np.row_stack((T,T+0.1)) # slight offset
PP = np.row_stack((P,P))

surfaces3 = pygamma.neutral_surfaces(glevels, SS, TT, PP,
                    lon=[lon]*2, lat=[lat]*2)

#print("Test 1-D versus 2-D, outputs match:", np.all(surfaces == surfaces3[0]))
def test_1D_vs_2D():
    assert_array_equal(surfaces, surfaces3[0])

# Test masked input:
Sm = np.ma.array(S)
Sm[-3:] = np.ma.masked
Sm[0] = np.ma.masked

gamma_masked = pygamma.gamma_n(Sm, T, P, lon, lat)[0]

SSm = np.ma.array(SS)
SSm[1, -3:] = np.ma.masked
SSm[0,:2] = np.ma.masked

gamma_masked2 = pygamma.gamma_n(SSm, TT, PP, [lon]*2, [lat]*2)[0]
def test_masking():
    assert_array_equal(SSm.mask, gamma_masked2.mask)

# Now compare to the example.out from a slight modification of the original
# package example.  (The fortran print doesn't seem to be able to write
# a string without putting a space at the beginning of the line, so we
# put a space in front of each label here, too.)

outlines = []
outlines.extend(['', ' location', '', '%12.4f%12.4f' % (lon, lat), ''])
outlines.extend([' labels', ''])

P, gamma, dg_lo, dg_hi

for rec in zip(P, gamma, dg_lo, dg_hi):
    outlines.append('%8.2f%12.6f%12.6f%12.6f' % rec)

outlines.extend(['', ' surfaces', ''])

for i, lev in enumerate(glevels):
    surf = list(surfaces[i])
    rec = tuple([lev] + surf[:3] + surf[-1:])
    outlines.append('%8.2f%12.6f%12.6f%14.6f%6.1f' % rec)

outlines.extend(['',''])
with open('example_py.out', 'w') as f:
    f.write('\n'.join(outlines))

example = os.path.join(pygamma.datadir, 'example.out')

def test_example_output():
    with open(example) as target:
        with open('example_py.out') as new:
            found = new.readline()
            expected = target.readline()
            assert found == expected
