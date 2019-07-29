#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table

from Chandra.Time import DateTime
from Chandra.Time import Time
from kadi import events

import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime



plt.style.use('ggplot')
labelsizes = 15
plt.rcParams['font.size'] = labelsizes
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = labelsizes
plt.rcParams['xtick.labelsize'] = labelsizes
plt.rcParams['ytick.labelsize'] = labelsizes

hk_msidset = ['2SMTRATM',
        '2TLEV1RT',
        '2VLEV1RT', 
        '2SHEV1RT', 
        '224PCAST', 
        '215PCAST',
        '215NCAST',
        '2SPTPAST',
        '2SPBPAST',
        '2IMTPAST',
        '2IMBPAST',
        '2NYMTAST',
        '2PYMTAST',
        '2CLMTAST',
        '2DRMTAST',
        '2ALMTAST',
        '2MSMDARS',
        '2MDIRAST',
        '2MSNBAMD',
        '2MSNAAMD', '2MSLBAMD', '2MSLAAMD', '2MSPRAMD', '2MSDRAMD', '2MCMDARS',
        '2MCNBAMD', '2MCNAAMD', '2MCLBAMD', '2MCLAAMD', '2MCPRAMD', '2MDRVAST', 
        '2OSLSAST', '2OPLSAST', '2CSLSAST', '2CPLSAST', '2OSLSADT', '2OSLSAAC',
        '2OPLSAAC', '2CSLSADT', '2CSLSAAC', '2CPLSAAC', '2FCPUAST', '2FCPVAST',
        '2CBHUAST', '2CBLUAST']



temp_msidset = [
        "2FE00ATM",  # Front-end Temperature (c)
        "2LVPLATM",  # LVPS Plate Temperature (c)
        "2IMHVATM",  # Imaging Det HVPS Temperature (c)
        "2IMINATM",  # Imaging Det Temperature (c)
        "2SPHVATM",  # Spectroscopy Det HVPS Temperature (c)
        "2SPINATM",  # Spectroscopy Det Temperature (c)
        "2FEPRATM"
        # "2PMT1T"  ,  # PMT 1 EED Temperature (c)
        # "2PMT2T"  ,  # PMT 2 EED Temperature (c)
        # "2DCENTRT",  # Outdet2 EED Temperature (c)
        # "2FHTRMZT",  # FEABox EED Temperature (c)
        # "2CHTRPZT",  # CEABox EED Temperature (c)
        # "2FRADPYT",  # +Y EED Temperature (c)
        # "2CEAHVPT",  # -Y EED Temperature (c)
        # "2CONDMXT",  # Conduit Temperature (c)
        # "2UVLSPXT",  # Snout Temperature (c)
        # #"2CE00ATM",  # CEA Temperature 1 (c) THESE HAVE FEWER POINTS AS THEY WERE RECENTLY ADDED BY TOM
        # #"2CE01ATM",  # CEA Temperature 2 (c) THESE HAVE FEWER POINTS AS THEY WERE RECENTLY ADDED BY TOM
        # "2FEPRATM",  # FEA PreAmp (c)
        # # "2SMTRATM",  # Selected Motor Temperature (c) THIS IS ALWAYS 5 DEGREES THROUGHOUT ENTIRE MISSION
        # "2DTSTATT"   # OutDet1 Temperature (c)
]


jupiter_ssc_start = '2019:149:10:00:00.000'
jupiter_ssc_stop =  '2019:149:17:00:00.000'


hk_data = fetch.MSIDset(hk_msidset, jupiter_ssc_start, jupiter_ssc_stop, filter_bad=True)
temp_data = fetch.MSIDset(temp_msidset, jupiter_ssc_start, jupiter_ssc_stop, filter_bad=True)


ave_table = Table()
all_names = []
all_means = []

for msid in temp_data:
    mean = np.mean(temp_data[msid].vals)

    all_names.append(msid)
    all_means.append(mean)

ave_table["MSIDname"] = all_names
ave_table["Average"] = all_means

ave_table.sort("Average")
ordered_temp_msidlist = ave_table["MSIDname"]


fig, ax = plt.subplots(figsize=(12,12))
for item in hk_msidset:
    # ax.plot(data[item].times, data[item].vals, label=f'{item}')
    # print(hk_data[item].raw_vals)

    if hk_data[item].raw_vals is not None:
        ax.plot(hk_data[item].times, hk_data[item].raw_vals, label=f'{item}')

ax.legend()

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(12,18))

n_lines = len(ordered_temp_msidlist)
color_idx = np.linspace(0, 1, n_lines)

for i, msid in zip(color_idx,ordered_temp_msidlist):
    ax1.plot(temp_data[msid].times, temp_data[msid].vals, label=f'{msid}', color=plt.cm.Oranges(i))

for msid in hk_msidset[:4]:
    ax2.plot(hk_data[msid].times, hk_data[msid].vals, label=f'{msid}')

for item in hk_msidset[4:]:
    if hk_data[item].raw_vals is not None:
        ax3.plot(hk_data[item].times, hk_data[item].raw_vals, label=f'{item}')


ax1.legend()
ax2.legend()
ax3.legend()

ax2.set_xlabel("Chandra Time")
plt.tight_layout()
plt.show()