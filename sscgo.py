#!/usr/bin/env python

from Chandra.Time import DateTime
from Chandra.Time import Time
from kadi import events

import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime

import matplotlib.pyplot as plt

plt.style.use('ggplot')
labelsizes = 15
plt.rcParams['font.size'] = labelsizes
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = labelsizes
plt.rcParams['xtick.labelsize'] = labelsizes
plt.rcParams['ytick.labelsize'] = labelsizes

msidset = ['2SMTRATM',
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


jupiter_ssc_start = '2019:149:10:00:00.000'
jupiter_ssc_stop =  '2019:149:17:00:00.000'

data = fetch.MSIDset(msidset, jupiter_ssc_start, jupiter_ssc_stop)

fig, ax = plt.subplots(figsize=(12,12))

for item in msidset[:3]:
    ax.plot(data[item].times, data[item].vals, label=f'{item}')

ax.legend()

plt.show()