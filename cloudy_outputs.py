"""
Script for processing a grid of cloudy ouputs. In active development!
"""
from astropy.io import ascii
from astropy.table import Table, join
import re
import numpy as np
from astropy import units as u
from scipy import constants


# Function to read a spectrum in cloudy and return the correct spectrum
def get_spectrum(f):
    c_angs = (constants.c*u.m/u.s).to(u.Angstrom/u.s) 
    # Split into lines
    lines = f.splitlines()

    # Load the various spectra
    wl = [float(line.split()[0]) for line in lines if line.strip()] # Energy in rydbergs
    wl = np.array(wl) * u.Ry
    ev = wl.to(u.eV)
    joules = ev.to(u.J)
    wl_ang = (constants.h * c_angs) / joules # Final wavelength in angstroms

    source = [float(line.split()[1]) for line in lines if line.strip()] 
    total = [float(line.split()[4]) for line in lines if line.strip()] 
    nebular = [float(line.split()[3]) for line in lines if line.strip()] 

    total= np.array(total) / wl_ang.value # in erg/s/A
    nebular = np.array(nebular) / wl_ang.value # in erg/s/A
    source = np.array(source) / wl_ang.value # in erg/s/A

    return wl_ang.value[::-1], source[::-1], total[::-1], nebular[::-1]

# Parse the cloudy GRID file into a machine readable format
def read_cloudy_grid(path):
    with open(path) as f:
        header = f.readline().lstrip('#').rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
 
    param_names = [h.split('=')[0].strip() for h in header[6:-1]]
    cols = list(zip(*rows))
 
    t = Table()
    t['index'] = [int(v) for v in cols[0]]
    t['failure'] = [v.strip() == 'T' for v in cols[1]]
    t['warnings'] = [v.strip() == 'T' for v in cols[2]]
    t['exit_code'] = [v.strip() for v in cols[3]]
    t['rank'] = [int(v) for v in cols[4]]
    t['seq'] = [int(v) for v in cols[5]]
    for i, name in enumerate(param_names):
        t[name] = [float(v) for v in cols[6 + i]]
    return t

def read_cloudy_spectrum(path):
    with open(path, "r") as f:
        data = f.read()
    # We don't need the header here since there's spaces in each file
    models = re.split(r"#+\s*GRID_DELIMIT\s*--\s*grid\d+", data)
    # Remove any empty chunks
    models = [m.strip() for m in models if m.strip()]
    wl_all, inc_all, comb_all, neb_all = [],[],[],[]
    ind = 0
    inds = []
    for m in models:
        #print("Processing new continuum model: %i"%count)
        test = m
        lines = test.splitlines()
        new_cont = "\n".join(lines[1:])
        tf = new_cont
        tf = re.sub(r'^#+\s*$', '', tf, flags=re.M).strip()
        wl, inc, comb, neb = get_spectrum(tf)
        wl_all.append(wl)
        inc_all.append(inc)
        comb_all.append(comb)
        neb_all.append(neb)
        inds.append(ind)
        ind+=1
    t = Table()
    t['index'] = inds
    t['wave'] = wl_all
    t['flux_incident'] = inc_all
    t['flux_nebular'] = neb_all
    t['flux_total'] = comb_all
    return t


def _clean_label(label):
    parts = label.split()
    m = re.match(r'([\d.]+)(\D*)$', parts[-1])
    if m:
        wave = np.format_float_positional(float(m.group(1)), precision=4,
                                          unique=False, fractional=False, trim='-')
        parts[-1] = wave + m.group(2)
    return ''.join(parts[:-1]) + '_' + parts[-1]
 
def read_cloudy_lines(path):
    with open(path) as f:
        labels = f.readline().rstrip('\n').split('\t')[1:]
        models, current = [], None
        for line in f:
            if line.startswith('#'):
                if 'GRID_DELIMIT' in line and current is not None:
                    models.append(current)
                    current = None
                continue
            if line.strip():
                current = [float(v) for v in line.rstrip('\n').split('\t')[1:]]
        if current is not None:                 # file without trailing delimiter
            models.append(current)
 
    t = Table()
    t['index'] = list(range(len(models)))
    names = [_clean_label(lab) for lab in labels]
    for i, (name, label) in enumerate(zip(names, labels)):
        # make duplicate labels unique so no column is overwritten
        if name in t.colnames:
            n = 2
            while f'{name}_{n}' in t.colnames:
                n += 1
            name = f'{name}_{n}'
        t[name] = [m[i] for m in models]
    return t


def get_cloudy_table(root):
    grid = read_cloudy_grid(root+'.grd')
    lines = read_cloudy_lines(root+'.lin')
    spec = read_cloudy_spectrum(root+'.con')
    j1 = join(grid,lines, keys='index')
    return join(j1,spec, keys='index')
