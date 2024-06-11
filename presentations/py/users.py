from matplotlib import pyplot, rc, dates

import numpy
from IPython import embed

user_dates = ["2021-03-11", "2022-04-29", "2022-11-07", "2022-12-06", "2023-06-08", "2023-06-29", "2023-07-11", "2023-09-03", "2023-10-13", "2023-12-01", "2023-12-15", "2024-02-22", "2024-03-21", "2024-04-09", "2024-05-02", "2024-05-19", "2024-06-06", "2024-06-10"]
user_dates = numpy.array([numpy.datetime64(date) for date in user_dates])
user_number = numpy.array([125, 293, 390, 394, 477, 487, 506, 518, 531, 544, 551, 568, 579, 588, 596, 603, 616, 620])

user_pred_dates = numpy.array([numpy.datetime64(date)
                            for date in ["2024-06-10", "2024-12-31", "2025-12-31", "2026-12-31",
                                         "2027-12-31"]])
user_pred_num = numpy.array([620, 620+0.5*160, 620+1.5*160, 620+2.5*160, 620+3.5*160])


cite_dates = ["2020-12-31", "2021-12-31", "2022-12-31", "2023-12-31", "2024-06-10"]
cite_dates = numpy.array([numpy.datetime64(date) for date in cite_dates])
cite_ref = numpy.cumsum([7, 24, 31, 56, 53])
cite_all = numpy.cumsum([8, 25, 33, 68, 81])

cite_pred_dates = numpy.array([numpy.datetime64(date)
                            for date in ["2024-06-10", "2024-12-31", "2025-12-31", "2026-12-31", "2027-12-31"]])
cite_pred_all = numpy.cumsum([8, 25, 33, 68, 81, numpy.sqrt(1.5)*81, 1.5*68, 1.5**2*68, 1.5**3*68])
cite_pred_all = cite_pred_all[cite_all.size-1:]

rc('font', size=14)

w,h = pyplot.figaspect(1)
fig = pyplot.figure(figsize=(1.5*w,1.5*h))

ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
ax.plot(user_dates, user_number, ls='-', color='k', label='Slack Users')
ax.plot(user_pred_dates, user_pred_num, ls=':', color='k')
ax.set_ylim([0,1400])
ax.set_ylabel("Cumulative Usage Metric")
ax.set_xlabel("Date")
ax.xaxis.set_major_locator(dates.MonthLocator(bymonth=[1,7]))
fig.autofmt_xdate()

#axt = ax.twinx()
#ax.plot(cite_dates, cite_ref, ls='--', color='0.5', label='All Cite')
ax.plot(cite_dates, cite_all, ls='-', color='C0', label='Citations')
ax.plot(cite_pred_dates, cite_pred_all, ls=':', color='C0')
#axt.set_ylim([0,250])

#ax.scatter(user_dates, user_number, marker='.', lw=0, s=200, color='k')
#fig.canvas.print_figure('pypeit_users.pdf', bbox_inches='tight')

ax.legend()
pyplot.show()
fig.clear()
pyplot.close(fig)

