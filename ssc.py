import os

import datetime as dt

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib import dates

import numpy as np
from scipy.interpolate import spline
from scipy.interpolate import interp1d
from scipy.signal import hilbert

try:
    from hrcsentinel import hrccore as hrc
except ImportError:
    raise ImportError(
        "hrcsentinel required. Download here: \
        https://github.com/granttremblay/HRCsentinel")

hrc.styleplots()

home_directory = os.path.expanduser("~")
msid_directory = home_directory + "/Dropbox/Work/HRCOps/MSIDCloud/"
figure_save_directory = home_directory + "/HRCOps/Projects/HRCSSC/figures/"

# Make that directory if it doesn't exist
if not os.path.exists(figure_save_directory):
    os.mkdir(figure_save_directory)
    print("Made directory {}".format(figure_save_directory))

rasterized = True

temperature_msids = [   "2FE00ATM",  # Front-end Temperature (c)
                		"2LVPLATM",  # LVPS Plate Temperature (c)
                		"2IMHVATM",  # Imaging Det HVPS Temperature (c)
                		"2IMINATM",  # Imaging Det Temperature (c)
                		"2SPHVATM",  # Spectroscopy Det HVPS Temperature (c)
                		"2SPINATM",  # Spectroscopy Det Temperature (c)
                		"2PMT1T",    # PMT 1 EED Temperature (c)
                		"2PMT2T",    # PMT 2 EED Temperature (c)
                		"2DCENTRT",  # Outdet2 EED Temperature (c)
                		"2FHTRMZT",  # FEABox EED Temperature (c)
                		"2CHTRPZT",  # CEABox EED Temperature (c)
                		"2FRADPYT",  # +Y EED Temperature (c)
                		"2CEAHVPT",  # -Y EED Temperature (c)
                		"2CONDMXT",  # Conduit Temperature (c)
                		"2UVLSPXT",  # Snout Temperature (c)
                		"2CE00ATM",  # CEA Temperature 1 (c)
                		"2CE01ATM", # CEA Temperature 2 (c)
                        "2FEPRATM", # FEA PreAmp (c)
                        "2SMTRATM", # Selected Motor Temperature (c)
                        "2DTSTATT" # OutDet1 Temperature (c)
    	            ]

voltage_msids = [       "2PRBSVL",   # Primary Bus Voltage (V)
                		"2PRBSCR",   # Primary Bus Current (amps)
                        "2C05PALV",  # +5V Bus Monitor
                        "2C15NALV",  # -15V Bus Monitor
                        "2C15PALV",  # +15V Bus Monitor
                        "2C24PALV",  # +24V Bus Monitor
    	            ]

temperature_data = {}
voltage_data = {}
binning = "full"
valtype = "vals"
timeperiod = "pastyear"

print("\nReading Temperature MSIDs")
for msidname in temperature_msids:
    print("Reading {} {} for {}".format(binning, valtype.upper(), msidname))
    times, values = hrc.parse_generic_msid(msid_directory + "{}".format(msidname) + "_{}_{}.csv".format(binning, timeperiod), valtype)
    temperature_data["{}_times".format(msidname)] = times
    temperature_data["{}_values".format(msidname)] = values

print("\nReading Voltage MSIDs")
for msidname in voltage_msids:
    print("Reading {} {} for {}".format(binning, valtype.upper(), msidname))
    times, values = hrc.parse_generic_msid(msid_directory + "{}".format(msidname) + "_{}_{}.csv".format(binning, timeperiod), valtype)
    voltage_data["{}_times".format(msidname)] = times
    voltage_data["{}_values".format(msidname)] = values

print("\nCollecting SSC Events flagged by bit 10 in 2SMTRATM")
ssc_events = {}
ssc_flags_file = msid_directory + "HRC_SS_HK_BAD_full_pastfiveyears.csv"
ssc_times, ssc_flags = hrc.parse_generic_msid(ssc_flags_file, valtype="vals")

ssc_events['SSC_times'] = ssc_times
ssc_events['SSC_flags'] = ssc_flags

print("\nCreating SSC Byte Shift Anomaly Table")
ssc_events_table = Table(ssc_events)

byte_shift_mask = ssc_events_table['SSC_flags'] == 1024
print("There are {} byte shift anomaly events over the past three years.".format(sum(byte_shift_mask)))
byte_shift_times = ssc_events_table[byte_shift_mask]['SSC_times']


temperature_figure_savename = figure_save_directory + "temperature.pdf"

fig, ax = plt.subplots(figsize=(12, 8))

for msidname in temperature_msids:
    ax.plot_date(temperature_data["{}_times".format(msidname)],
                 temperature_data["{}_values".format(msidname)], '.', alpha=1.0, markersize=2.5, label='{}'.format(msidname), rasterized=rasterized)
ax.set_ylabel('Temperature (C)')
ax.set_xlabel('Date')
ax.set_ylim(0, 40)

ssccolor = 'gray'
for byte_shift_time in byte_shift_times[-10000:]:
    ax.axvline(byte_shift_time, lw=0.5, color=ssccolor, alpha=0.3, rasterized=rasterized)


    ax.set_title("HRC Thermistor Temperatures over the past year")

ax.legend(loc=2, prop={'size': 13})
plt.show()
fig.savefig(temperature_figure_savename, dpi=300)




voltage_figure_savename = figure_save_directory + "voltages.pdf"

fig, ax = plt.subplots(figsize=(12, 8))

for msidname in voltage_msids:
    ax.plot_date(voltage_data["{}_times".format(msidname)],
                 voltage_data["{}_values".format(msidname)], '.', alpha=1.0, markersize=2.5, label='{}'.format(msidname), rasterized=rasterized)
ax.set_ylabel('Voltage (V)')
ax.set_xlabel('Date')

ssccolor = 'gray'
for byte_shift_time in byte_shift_times[-10000:]:
    ax.axvline(byte_shift_time, lw=0.5, color=ssccolor, alpha=0.3, rasterized=rasterized)


    ax.set_title("HRC Voltages over the past year")

ax.legend(loc=2, prop={'size': 13})
plt.show()
fig.savefig(voltage_figure_savename, dpi=300)
