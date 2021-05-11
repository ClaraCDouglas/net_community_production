"""
Simple python wrapper for the fortran gamma-n code of
Jackett and McDougall, 1997.

'A neutral density variable for the world's oceans' by David Jackett
and Trevor McDougall, Journal of Physical Oceanography, 1997, Vol.27(2),
237-263.


See the src/README and src/ConditionsOfUse.

The fortran code has been modified to make it more suitable
for wrapping, for example, by returning an error code instead
of executing 'stop'.

The module provides two functions: gamma_n and neutral_surfaces.

"""

import os
import numpy as np
import numpy.ma as ma
from functools import reduce

_thisdir = os.path.split(os.path.realpath(__file__))[0]

# If necessary, provision can be made later for alternative
# data locations.
datadir = os.path.join(_thisdir, 'data')
ncfile = os.path.join(datadir, 'gamma.nc')
if not os.path.exists(ncfile):
    raise ImportError("Cannot find data file, %s" % ncfile)


from pygamma import _gamma
from pygamma._gamma import neutral_surfaces as _neutral_surfaces

def _check(ierr):
    """
    Handle errors returned by fortran code.
    """
    if ierr == 1:
        # Should not happen; we check for it at the python level.
        raise ValueError("Latitude or Longitude is out of range")
    elif ierr == 2:
        raise ValueError("ERROR 1 in neutral-surfaces.f : missing gamma value")
    elif ierr == 3:
        raise RuntimeError("ERROR 2 in neutral-surfaces.f")
    elif ierr == 4:
        raise RuntimeError("ERROR 3 in neutral-surfaces.f")
    elif ierr == 5:
        raise RuntimeError('ERROR 1 in gamma-out-of-range.f')
    elif ierr == 6:
        raise RuntimeError('ERROR 1 in gamma-errors: negative scv error')
    elif ierr == 7:
        raise RuntimeError('n > nmax (100) in depth-ns.f')
    elif ierr == 8:
        raise RuntimeError('ERROR 1 in depth-ns.f')
    elif ierr == 9:
        raise RuntimeError('n > nmax (100) in depth-scv.f')
    elif ierr == 10:
        raise RuntimeError('ERROR 1 in depth-scv.f')
    elif ierr == 11:
        raise RuntimeError('ERROR 2 in depth-scv.f')
    elif ierr == 12:
        raise IOError('cannot open file in get-lunit.f')
    elif ierr == 13:
        raise IOError('cannot open data file in read-nc.f')



def _gamma_n(*args):
    ret = _gamma.gamma_n(*args)
    _check(ret[-1])
    return ret[:-1]

def _neutral_surfaces(*args):
    ret = _gamma.neutral_surfaces(*args)
    _check(ret[-1])
    return ret[:-1]

# Copied from matplotlib.cbook:
def _safe_masked_invalid(x, copy=False):
    x = np.array(x, subok=True, copy=copy)
    if not x.dtype.isnative:
        # Note that the argument to `byteswap` is 'inplace',
        # thus if we have already made a copy, do the byteswap in
        # place, else make a copy with the byte order swapped.
        # Be explicit that we are swapping the byte order of the dtype
        x = x.byteswap(copy).newbyteorder('S')

    try:
        xm = np.ma.masked_invalid(x, copy=False)
        xm.shrink_mask()
    except TypeError:
        return x
    return xm

def broadcast(*args):
    def _mask_or(a, b):
        return ma.mask_or(a, b, shrink=True)
    args = [_safe_masked_invalid(arg) for arg in args]
    if any([ma.isMA(arg) for arg in args]):
        vars = [ma.getdata(var) for var in args]
        mvars = [ma.getmaskarray(var) for var in args]
        outargs = list(map(np.array, np.broadcast_arrays(*vars)))
        masks = list(map(np.array, np.broadcast_arrays(*mvars)))
        mask = reduce(_mask_or, masks)
    else:
        mask = ma.nomask
        # Using map(np.array, ...) to get contiguous copies.
        outargs = list(map(np.array, np.broadcast_arrays(*args)))
    if outargs[0].ndim == 0:
        scalar = True
        for arg in outargs:
            arg.shape = (1,)
        if mask is not ma.nomask:
            mask.shape = (1,)
    else:
        scalar = False
    return scalar, mask, outargs


def gamma_n(s, t, p, lon, lat):
    """
    Arguments:
        salinity
        temperature
        pressure
        lon
        lat

    s, t, p must be broadcastable against each other,
    and 0-D, 1-D, or 2-D.  They may be masked, or nan
    may be used for invalid values.
    lon, lat must be at most 1-D, and not masked.

    If lon and/or lat is 1-D, its length must match the *first*
    dimension of s, t, p; in that case, if the latter are 2-D, their
    *second* dimension must be the depth dimension.

    Returns:
        gamma (generally a masked array or a scalar)
        dg_lo
        dg_hi

    """
    scalar, mask, (s, t, p) = broadcast(s, t, p)
    masked = (mask is not ma.nomask)
    if masked:
        gmask = ~mask

    lon, lat = list(map(np.array, np.broadcast_arrays(lon, lat)))

    if lon.ndim == 0:
        if s.ndim != 1:
            raise ValueError("lat, lon dimensions must match number of profiles")
        if lat < -80 or lat > 64:
            g = np.ma.masked_all(s.shape, dtype=float)
            dl = g.copy()
            dh = g.copy()
            return g, dl, dh
        if masked:
            g = ma.zeros(s.shape, dtype=float)
            g.mask = mask
            dl = g.copy()
            dh = g.copy()
            if gmask.any():
                g[gmask], dl[gmask], dh[gmask] = _gamma_n(ncfile,
                                        s[gmask], t[gmask], p[gmask],
                                        float(lon), float(lat))
        else:
            g, dl, dh =  _gamma_n(ncfile, s, t, p, float(lon), float(lat))
    elif lon.ndim == 1:
        if s.ndim != 2:
            raise ValueError("incompatible dimensions")
        if len(lon) != s.shape[0]:
            raise ValueError("len(lat), len(lon) must match s.shape[0]")
        if masked:
            g = ma.zeros(s.shape, dtype=float)
            g.mask = mask
        else:
            g = np.zeros(s.shape, dtype=float)
        dl = g.copy()
        dh = g.copy()
        if masked:
            for i in range(s.shape[0]):
                if lat[i] < -80 or lat[i] > 64:
                    g[i] = np.ma.masked
                    dl[i] = np.ma.masked
                    dh[i] = np.ma.masked
                    continue
                gm = gmask[i]
                if gm.any():
                    g[i,gm], dl[i,gm], dh[i,gm] = _gamma_n(ncfile,
                                                  s[i,gm], t[i,gm], p[i,gm],
                                                         lon[i], lat[i])

        else:
            badlat = False
            for i in range(s.shape[0]):
                if lat[i] < -80 or lat[i] > 64:
                    g[i] = np.nan
                    dl[i] = np.nan
                    dh[i] = np.nan
                    badlat = True
                    continue
                g[i], dl[i], dh[i] = _gamma_n(ncfile, s[i], t[i],
                                                         p[i], lon[i], lat[i])
            if badlat:
                g = np.ma.masked_invalid(g)
                dl = np.ma.masked_invalid(dl)
                dh = np.ma.masked_invalid(dh)
    else:
        raise ValueError("Unsupported combination of input dimensions")


    bad = (g < -98)
    if bad.any():
        g = np.ma.masked_where(bad, g)
        masked = True
    if scalar:
        g.shape = ()
        dl.shape = ()
        dh.shape = ()
        if not masked:
            g, dl, dh = [float(a) for a in (g, dl, dh)]
    return g, dl, dh

def neutral_surfaces(glevels, s, t, p, gamma=None, lon=None, lat=None):
    """
    Arguments:
        glevels, length ng (1-D)
        s, t, p, gamma, all the same shape, (n,) or (nprofs, n);
            they may be masked, or use nan as a bad value.
        gamma is optional; if omitted, it will be calculated,
            using the values of lon and lat, which must be specified
            only if gamma is omitted.

    Returns:
        structured array, shape (ng,) or (nprofs, ng),
            field names s, t, p, ds, dt, dp


    """
    glevels = np.asarray(glevels, dtype=float).flatten()
    ng = len(glevels)

    if gamma is None:
        gamma = gamma_n(s, t, p, lon, lat)[0]

    scalar, mask, (s, t, p, gamma) = broadcast(s, t, p, gamma)
    masked = (mask is not ma.nomask)
    if masked:
        gmask = ~mask

    # We *cannot* use T here if we want to use a recarray as the
    # output, because it will be interpreted as "transpose" when
    # attribute lookup is done.

    dtype = np.dtype(dict(names=['s', 't', 'p', 'ds', 'dt', 'dp'],
                          formats=[np.float_]*6))

    if s.ndim > 2 or s.ndim < 1:
        raise ValueError("Only 1-D or 2-D inputs are supported")

    inshape = s.shape

    if s.ndim == 1:
        if masked:
            out = ma.zeros((ng,6), dtype=float)
            out[mask] = ma.masked
            vars = _neutral_surfaces(s[gmask], t[gmask],
                                p[gmask], gamma[gmask], glevels)
            for i, var in enumerate(vars):
                out[gmask,i] = var
        else:
            out = np.empty((ng,6), dtype=float)
            vars = _neutral_surfaces(s, t, p, gamma, glevels)
            for i, var in enumerate(vars):
                out[:,i] = var
    else:
        nprofs = s.shape[0]
        if masked:
            out = np.ma.zeros((nprofs, ng, 6), dtype=float)
            out[:] = ma.masked
            for iprof in range(nprofs):
                gm = gmask[iprof]
                if gm.any():
                    vars = _neutral_surfaces(s[iprof,gm], t[iprof,gm],
                                                p[iprof,gm],
                                                gamma[iprof,gm], glevels)
                    for i, var in enumerate(vars):
                        out[iprof, :, i] = var
            out[out < -98] = ma.masked
        else:
            out = np.empty((nprofs, ng, 6), dtype=float)
            for iprof in range(nprofs):
                vars = _neutral_surfaces(s[iprof], t[iprof], p[iprof],
                                                gamma[iprof], glevels)
                for i, var in enumerate(vars):
                    out[iprof, :,i] = var

    out = out.view(dtype=dtype)
    out.shape = out.shape[:-1]
    return out                            #.view(np.recarray)



