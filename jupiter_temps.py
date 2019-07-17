#!/usr/bin/env python

# Compare distributions of an MSID (temperature) for normal and corrupted
# secondary science.  The MSID is interpolated to get values during SS
# corruption.

import sys
import math
import matplotlib.pyplot as plt
from Ska.engarchive import fetch
import numpy as np
import argparse
import datetime as dt
from scipy.stats import norm
import Chandra.Time

dtype = np.dtype('f8')


plt.style.use('ggplot')
labelsizes = 15
plt.rcParams['font.size'] = labelsizes
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = labelsizes
plt.rcParams['xtick.labelsize'] = labelsizes
plt.rcParams['ytick.labelsize'] = labelsizes


def getStartTime(tend, duration):
    """Determine start time from end time string the duration
    as a timedelta.
    """
    # tend format is YYYY:DDD or YYYY:DD:HH:MM:SS
    t = len(tend.split(':'))
    if t == 2:
        finishtime = dt.datetime.strptime(tend, "%Y:%j")
    elif t == 5:
        finishtime = dt.datetime.strptime(tend, "%Y:%j:%H:%M:%S")
    else:
        print("Unknown end time format:", tend)
        sys.exit(1)

    begintime = finishtime - duration
    tstart = begintime.strftime("%Y:%j:18:00:00")
    print("Time interval:", tstart, "to", tend)
    return tstart


def getBadTemps(Temp, Tmin, Tmax):
    """Interpolated temperatures for times of corrupted secondary science.
    """
    print("Number initially flagged as corrupt:", sum(Temp.bads))
    Tbad = []
    # Temperature estimates for bad data
    nTemp = len(Temp.vals)
    for i in range(nTemp):
        # MSID values outside [Tmin, Tmax] are also assumed to be bad
        if Temp.vals[i] < Tmin or Temp.vals[i] > Tmax:
            Temp.bads[i] = True
        if Temp.bads[i]:
            # Interpolate the temperature from good bracketing values
            ibelow = i - 1
            while ibelow >= 0 and Temp.bads[ibelow]:
                ibelow -= 1
            iabove = i + 1
            while iabove < nTemp and Temp.bads[iabove]:
                iabove += 1
            if ibelow < 0:
                Tinterp = Temp.vals[iabove]
            elif iabove >= nTemp:
                Tinterp = Temp.vals[ibelo]
            else:
                Tinterp = ((iabove - i) * Temp.vals[ibelow] + (i - ibelow)
                           * Temp.vals[iabove]) / (iabove - ibelow)
            Tbad.append(Tinterp)

    # Sanity checks
    nbad = sum(Temp.bads)
    print("Total number of samples:", nTemp)
    print("Number of temperatures flagged as bad:", nbad)
    print("len (Tbad) =", len(Tbad))
    print("Overall fraction of bad temperatures:", nbad / float(nTemp))

    return Tbad


def goodHist(Temp):
    """Make a histogram of good temperatures and plot it as the
    distribution of fractions.
    """
    # Good temperatures
    Tgood = [x for x, y in zip(Temp.vals, Temp.bads) if not y]
    # Counts of good temperatures
    Tgcount = dict()
    for x in Tgood:
        if x in Tgcount:
            Tgcount[x] += 1
        else:
            Tgcount[x] = 1

    # Make normalized histogram of the good temperatures.
    # The ordered list of discrete temperature values
    Tval = sorted(Tgcount.keys())
    print("Number of distinct temperatures:", len(Tval))
    print("T values:", Tval)
    # For plotting
    Tgpl = np.array(Tval, dtype)
    norm = 1.0 / float(len(Tgood))
    # Counts for good temperatures as a list ordered by temperature
    Tgcl = [Tgcount[T] for T in Tval]
    fgood = np.array([c * norm for c in Tgcl], dtype)
    # These errors are generally negligible
    fgerr = np.array([math.sqrt(c) * norm for c in Tgcl], dtype)
    plt.errorbar(Tgpl, fgood, yerr=fgerr, fmt='o')
    return Tval, Tgpl, Tgcl


def badHist(Tbad, Tval, Tgpl, Tgcl):
    """Make histogram of temperatures for corrupted secondary science
    and plot it.
    Tbad = the list of temperatures for corrupted SS
    Tval = ordered list of discrete (good) temperatures seen
    Tgpl = Tval as a numpy array
    Tgcl = counts of good temperatures corresponding to Tval
    """
    # Assemble counts of bad temperatures, binned to Tval
    nval = len(Tval)
    # Temperature cuts are not evenly spaced
    Tcuts = [0.5 * (Tval[i] + Tval[i + 1]) for i in range(nval - 1)]
    Tbcount = [0] * nval
    for T in Tbad:
        if T < Tcuts[0]:
            Tbcount[0] += 1
        elif T >= Tcuts[-1]:
            Tbcount[-1] += 1
        else:
            ilo = 0
            ihi = nval - 2
            while ihi > ilo + 1:
                i = (ilo + ihi) // 2
                if T < Tcuts[i]:
                    ihi = i
                else:
                    ilo = i
            Tbcount[ihi] += 1

    # Simple check that the extreme bins are not over-populated
    if Tbcount[0] > Tbcount[1]:
        print("****Lower MSID limit too high?****:", Tbcount[0], Tbcount[1])
    if Tbcount[-1] > Tbcount[-2]:
        print("****Upper MSID limit too low?****:", Tbcount[-1], Tbcount[-2])
    # Determine count fractions and probabilities, with errors
    bnorm = 1.0 / float(len(Tbad))
    # Gaussian 1 sigma confidence level
    q = 2.0 * norm.cdf(1.0) - 1.0
    fbl = []  # Bad fraction
    fbel = []  # Error on bad fraction
    prob = []  # Probability of being bad at this temperature
    pe = []  # Error in the probability
    for i in range(nval):
        bad = Tbcount[i]
        good = Tgcl[i]
        fbl.append(bad * bnorm)
        prob.append(bad / float(bad + good))
        if bad != 0:
            # Poisson error estimate for bad count > 0.
            err = math.sqrt(bad)
            fbel.append(bnorm * err)
            # The number of samples for this temperature is assumed fixed,
            # which hardly makes a difference when good >> bad
            pe.append(err / (bad + good))
        else:
            # When there are no bad temperatures, use the 1 sigma upper
            # limit on probability of being bad.  This is the binomial
            # result, but it is close to Poission unless the total number
            # of samples is small.
            pup = 1.0 - (1.0 - q)**(1.0 / float(good))
            fbel.append(pup * good * bnorm)
            pe.append(pup)
            print("Upper limit: T =", Tval[i], "; p =", pup, "; n =", good, \
                "; n p =", good * pup)

    fbad = np.array(fbl, dtype)
    fberr = np.array(fbel, dtype)
    plt.errorbar(Tgpl, fbad, yerr=fberr, fmt='o')
    # Force lower limit to zero
    ylims = plt.ylim()
    plt.ylim([0.0, ylims[1]])
    plt.xlabel(r'$T\ (\rm C)$')
    plt.ylabel("fraction")
    return prob, pe


############################################################
#
# For Chandra.Time formats see
# http://cxc.cfa.harvard.edu/mta/ASPECT/tool_doc/eng_archive/fetch_tutorial.html#date-and-time-formats

# Option parsing
parser = argparse.ArgumentParser(
    description='Compare good and bad MSID distributions.')
parser.add_argument('-m', '--msid',
                    help='MSID',
                    default='2FEPRATM')
parser.add_argument('-e', '--end',
                    help='ending time')
parser.add_argument('-d', '--duration',
                    help='time elapsed (days)',
                    default='3')
parser.add_argument('-l', '--low-threshold',
                    help='lower limit on good MSID values',
                    default=15.0)
parser.add_argument('-u', '--upper-threshold',
                    help='upper limit on good MSID values',
                    default=35.0)

argdata = parser.parse_args()
msid = argdata.msid
tend = argdata.end
if not tend:
    # Use 6 pm UTC to ensure getting today's updates
    tend = dt.datetime.utcnow().strftime("2019:149:17:00:00")
duration = dt.timedelta(float(argdata.duration))
Tmin = float(argdata.low_threshold)
Tmax = float(argdata.upper_threshold)

tstart = getStartTime(tend, duration)

# Collect MSID data
Temp = fetch.MSID(msid, tstart, tend)
tlast = Chandra.Time.DateTime(Temp.times[-1])
print("Final time in data: ", tlast.date)

# Get interpolated MSID values during secondary science corruption.
# NB: Updates temp.bads, so must precede goodHist
Tbad = getBadTemps(Temp, Tmin, Tmax)
nbad = len(Tbad)

# Plot the histogram of good MSID values
plt.figure(1)
Tval, Tgpl, Tgcl = goodHist(Temp)
legstr = ["Good " + msid]

# Plot the histogram of MSID values for corrupted secondary science
prob, pe = badHist(Tbad, Tval, Tgpl, Tgcl)
legstr.append("Bad " + msid)
plt.legend(legstr, loc='upper left')
plt.show()

# Plot probability of corrupt data vs temperature
plt.figure(2)
plt.errorbar(Tgpl, prob, yerr=pe, fmt='o')
lims = plt.ylim()
plt.ylim([0.0, lims[1]])
plt.xlabel(r'$T\ (\rm C)$')
plt.ylabel(r'$P_{\rm bad}$')
plt.legend([msid], loc='upper left')
plt.show()
