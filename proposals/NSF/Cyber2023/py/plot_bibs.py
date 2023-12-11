import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2

fn = "pypeit_sum.csv"

t = ascii.read(fn)
t['yearfrac'] = t['months'] / 12.
t['fullyear'] = t['years'] + t['yearfrac']

refereed = t['refereed'] == 1

cites = np.ones_like(t['fullyear'])

if False:
    plt.plot(t['fullyear'][refereed][::-1], np.cumsum(cites[refereed]), label='Refereed')
    plt.plot(t['fullyear'][::-1], np.cumsum(cites),label='All')
    plt.xlabel('Year')
    plt.ylabel('Cumulative Citations')
    plt.legend()
    plt.semilogy()
    #plt.savefig('cumu_cites.jpg')
    plt.show()

refereed = t['refereed'][::-1] == 1

tdates = t['fullyear'][::-1]
months = t['months'] + 1
pdates = [ f"{t['years'][i]}-{months[i]:02d}" for i in range(0,len(tdates))]
pdates = np.asarray(pdates)
pdates = pdates[::-1]

fig, ax = plt.subplots(1, 1, figsize=(8, 3), layout='constrained', linewidth=2)
fig.set_linewidth(2)
ax.plot(mdates.date2num(pdates[refereed]),  np.cumsum(cites[refereed]), label='Refereed')
ax.plot(mdates.date2num(pdates), np.cumsum(cites),label='All')
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
ax.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))

for label in ax.xaxis.get_ticklabels():
#    label.set_rotation(45)
    label.set_fontsize(12)
for label in ax.yaxis.get_ticklabels():
#    label.set_rotation(45)
    label.set_fontsize(12)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    
plt.legend(fontsize=12)
#plt.xlabel('Date')
ax.xaxis.set_label_text('Date', fontsize=12)
ax.yaxis.set_label_text('Cumulative Citations', fontsize=12)

plt.savefig('Figures/cumu_cites.jpg')

plt.show()



