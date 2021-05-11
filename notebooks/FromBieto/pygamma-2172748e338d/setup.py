
import sys
import os, shlex  # needed for hack; see below
from setuptools import setup
from numpy.distutils.core import Extension

try:
    from pycurrents.setup_helper import write_hg_status, use_currents_prefix
    write_hg_status()
    use_currents_prefix()
except ImportError:
    pass

print("Possibly modified setup argument list is:\n  ", ' '.join(sys.argv[1:]))


SOURCES = ['atg.f',
          'depth-ns.f',
          'depth-scv.f',
          'derthe.f',
          'eosall.f',
          'alphabeta.f',
          'eos8d.f',
          'e-solve.f',
          'gamma-errors.f',
          'gamma-n.f',
          'goor.f',
          'gamma-qdr.f',
          'get-lunit.f',
          'goor-solve.f',
          'index.f',
          'neutral-surfaces.f',
          'ocean-test.f',
          'scv-solve.f',
          'sig-vals.f',
          'stp-interp.f',
          'svan.f',
          'theta.f',
          'read-nc.f']

sources = ["src/_gamma.pyf"]
sources.extend(["src/%s" % f for f in SOURCES])

# Hack from
# https://github.com/DougBurke/sherpa/commit/334c23494737555046712f5afcfa449063fb39fb

f77_ldflags = os.environ.get('LDFLAGS', '')
if f77_ldflags == '':
    f77_extra_link_args = None
else:
    f77_extra_link_args = ["-shared"] + shlex.split(f77_ldflags)

# Another hack: homebrew doesn't seem to put the libnetcdff location into a
# standard search path.

library_dirs = None
if sys.platform == 'darwin':
    f77 = os.environ.get('F77', None)
    if f77 and 'darwin' in f77:
        # Looks like conda
        pass
    else:
        # Maybe it's homebrew and needs the following
        library_dirs = ['/usr/local/lib']

# Originally I thought (perhaps incorrectly) that hdf5_hl needed to be
# included in the library list, since netcdff uses it.  Apparently this
# is not true now, and including it can cause linking to fail.
_gamma = Extension(name='pygamma._gamma',
                 sources=sources,
                 extra_f77_compile_args=['-ffixed-line-length-none', '-fPIC'],
                 libraries=['netcdff',],
                 extra_link_args=f77_extra_link_args,  # part of hack, above
                 library_dirs = library_dirs, # part of second hack, above
                 )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'pygamma',
          packages=['pygamma'],
          package_dir={'pygamma':''},
          description       = "Interface to Jackett and McDougall gamma-n",
          author            = "Eric Firing",
          author_email      = "efiring@hawaii.edu",
          ext_modules = [_gamma],
          package_data={'pygamma': ['data/gamma.nc', 'data/example.dat',
                                    'data/example.out']},
          )


