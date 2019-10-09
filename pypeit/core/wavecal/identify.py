""" Module for finding patterns in arc line spectra
"""
import os
import sys
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
import matplotlib.transforms as mtransforms

matplotlib.use('Qt5Agg')

#from pypeit import utils
#from pypeit.core.wavecal import autoid
#from pypeit import msgs


class Identify(object):
    """
    Manually identify arc lines.

    """

    def __init__(self, canvas, ax, spec):
        """
        Description to go here
        """
        self.ax = ax
        self.spec = spec
        self.fb = None
        self.lineidx = 0
        self._addsub = 0  # Adding a region (1) or removing (0)
        self._changes = False
        self.annlines = []
        self.anntexts = []

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas = canvas

        # Sort the lines by decreasing probability of detection
        self.sort_lines()
        self.get_current_line()

        # Draw the first line
        self.next_spectrum()

    def draw_lines(self):
        # annotations = [child for child in self.ax.get_children() if isinstance(child, matplotlib.text.Annotation)]
        for i in self.annlines: i.remove()
        for i in self.anntexts: i.remove()
        self.annlines = []
        self.anntexts = []
        default = False
        molecules = False
        # Plot the lines
        if default:
            xmn, xmx = self.ax.get_xlim()
            ymn, ymx = self.ax.get_ylim()
            xmn /= (1.0 + self.prop._zabs)
            xmx /= (1.0 + self.prop._zabs)
            w = np.where((self.atom._atom_wvl > xmn) & (self.atom._atom_wvl < xmx))[0]
            for i in range(w.size):
                dif = i % 5
                self.annlines.append(self.ax.axvline(self.atom._atom_wvl[w[i]] * (1.0 + self.prop._zabs), color='r'))
                txt = "{0:s} {1:s} {2:.1f}".format(self.atom._atom_atm[w[i]], self.atom._atom_ion[w[i]],
                                                   self.atom._atom_wvl[w[i]])
                ylbl = ymn + (ymx - ymn) * (dif + 1.5) / 8.0
                self.anntexts.append(
                    self.ax.annotate(txt, (self.atom._atom_wvl[w[i]] * (1.0 + self.prop._zabs), ylbl), rotation=90.0,
                                     color='b', ha='center', va='bottom'))
            if molecules:
                # Plot molecules
                labls, lines = np.loadtxt("/Users/rcooke/Software/ALIS_dataprep/molecule.dat",
                                          dtype={'names': ('ion', 'wvl'), 'formats': ('S6', 'f8')}, unpack=True,
                                          usecols=(0, 1))
                w = np.where((lines > xmn) & (lines < xmx))[0]
                for i in range(w.size):
                    dif = i % 5
                    self.annlines.append(self.ax.axvline(lines[w[i]] * (1.0 + self.prop._zabs), color='g'))
                    txt = "{0:s} {1:.1f}".format(labls[w[i]], lines[w[i]])
                    ylbl = ymn + (ymx - ymn) * (dif + 1.5) / 8.0
                    self.anntexts.append(
                        self.ax.annotate(txt, (lines[w[i]] * (1.0 + self.prop._zabs), ylbl), rotation=90.0, color='b',
                                         ha='center', va='bottom'))
        else:
            zabs = np.array(
                [0.0, 0.72219, 0.86431, 1.15729, 1.724087, 1.84506, 2.16797, 2.30812, 2.315, 2.316, 2.37995, 2.433,
                 2.439, 2.44059, 2.57813, 2.57884, 2.7125])
            zabs = np.append(zabs, self.prop._zabs)
            labls, lines = np.loadtxt("/Users/rcooke/Software/ALIS_dataprep/atomic_UV.dat",
                                      dtype={'names': ('ion', 'wvl'), 'formats': ('S6', 'f8')}, unpack=True)
            for z in range(zabs.size):
                xmn, xmx = self.ax.get_xlim()
                ymn, ymx = self.ax.get_ylim()
                xmn /= (1.0 + zabs[z])
                xmx /= (1.0 + zabs[z])
                w = np.where((self.atom._atom_wvl > xmn) & (self.atom._atom_wvl < xmx))[0]
                for i in range(w.size):
                    dif = i % 5
                    self.annlines.append(self.ax.axvline(self.atom._atom_wvl[w[i]] * (1.0 + zabs[z]), color='r'))
                    txt = "{0:s} {1:s} {2:.1f}".format(self.atom._atom_atm[w[i]], self.atom._atom_ion[w[i]],
                                                       self.atom._atom_wvl[w[i]])
                    ylbl = ymn + (ymx - ymn) * (dif + 1.5) / 8.0
                    self.anntexts.append(
                        self.ax.annotate(txt, (self.atom._atom_wvl[w[i]] * (1.0 + zabs[z]), ylbl), rotation=90.0,
                                         color='b', ha='center', va='bottom'))
                # Do far UV list
                w = np.where((lines > xmn) & (lines < xmx))[0]
                for i in range(w.size):
                    dif = i % 5
                    self.annlines.append(self.ax.axvline(lines[w[i]] * (1.0 + zabs[z]), color='r'))
                    txt = "{0:s} {1:.1f}".format(labls[w[i]], lines[w[i]])
                    ylbl = ymn + (ymx - ymn) * (dif + 1.5) / 8.0
                    self.anntexts.append(
                        self.ax.annotate(txt, (lines[w[i]] * (1.0 + zabs[z]), ylbl), rotation=90.0, color='b',
                                         ha='center', va='bottom'))
        return

    def draw_callback(self, event):
        trans = mtransforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        if self.fb is not None:
            self.fb.remove()
        # Find all regions
        regwhr = np.copy(self.prop._regions == 1)
        # Fudge to get the leftmost pixel shaded in too
        regwhr[np.where((self.prop._regions[:-1] == 0) & (self.prop._regions[1:] == 1))] = True
        self.fb = self.ax.fill_between(self.prop._wave, 0, 1, where=regwhr, facecolor='green', alpha=0.5,
                                       transform=trans)
        self.ax.draw_artist(self.spec)
        self.draw_lines()

    def get_ind_under_point(self, event):
        """
        Get the index of the spectrum closest to the cursor
        """
        ind = np.argmin(np.abs(self.prop._wave - event.xdata))
        return ind

    def sort_lines(self, method="ion"):
        """
        Sort lines by decreasing probability of detection
        """
        if method == "sig":
            coldens = 10.0 ** (self.atom.solar - 12.0)
            ew = coldens * (self.atom._atom_wvl ** 2 * self.atom._atom_fvl)
            # snr = 1.0/self.prop._flue
            # indices = np.abs(np.subtract.outer(self.prop._wave, self.atom._atom_wvl*(1.0+self.prop._zabs))).argmin(0)
            sigdet = ew  # *snr[indices]
            self.sortlines = np.argsort(sigdet)[::-1]
        elif method == "ion":
            self.sortlines = np.arange(self.atom._atom_wvl.size)
        elif method == "wave":
            self.sortlines = np.argsort(self.atom._atom_wvl)
        return

    def get_current_line(self):
        if self.lineidx < 0:
            self.lineidx += self.sortlines.size
        if self.lineidx >= self.sortlines.size:
            self.lineidx -= self.sortlines.size
        self.linecur = self.sortlines[self.lineidx]

    def button_press_callback(self, event):
        """
        whenever a mouse button is pressed
        """
        if event.inaxes is None:
            return
        if self.canvas.toolbar.mode != "":
            return
        if event.button == 1:
            self._addsub = 1
        elif event.button == 3:
            self._addsub = 0
        self._start = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """
        whenever a mouse button is released
        """
        if event.inaxes is None:
            return
        if self.canvas.toolbar.mode != "":
            return
        self._end = self.get_ind_under_point(event)
        if self._end != self._start:
            if self._start > self._end:
                tmp = self._start
                self._start = self._end
                self._end = tmp
            self.update_regions()
        trans = mtransforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.canvas.restore_region(self.background)
        if self.fb is not None:
            self.fb.remove()
        # Find all regions
        regwhr = np.copy(self.prop._regions == 1)
        # Fudge to get the leftmost pixel shaded in too
        regwhr[np.where((self.prop._regions[:-1] == 0) & (self.prop._regions[1:] == 1))] = True
        self.fb = self.ax.fill_between(self.prop._wave, 0, 1, where=regwhr, facecolor='green', alpha=0.5,
                                       transform=trans)
        self.canvas.draw()

    def key_press_callback(self, event):
        """
        whenever a key is pressed
        """
        if not event.inaxes:
            return
        if event.key == '?':
            print("============================================================")
            print("       MAIN OPERATIONS")
            print("p       : toggle pan/zoom with the cursor")
            print("w       : write the spectrum with the associated region")
            print("q       : exit")
            print("------------------------------------------------------------")
            print("       SHORTCUTS TO MOVE BETWEEN LINES")
            print("[ / ]   : go to previous/next element")
            print(", / .   : go to previous/next ion")
            print("b / n   : go to previous/next line")
            print("------------------------------------------------------------")
            print("       ATOMIC DATA OF THE CURRENT LINE")
            print("{0:s} {1:s}  {2:f}".format(self.atom._atom_atm[self.linecur].strip(),
                                              self.atom._atom_ion[self.linecur].strip(),
                                              self.atom._atom_wvl[self.linecur]))
            print("Observed wavelength = {0:f}".format(self.atom._atom_wvl[self.linecur] * (1.0 + self.prop._zabs)))
            print("f-value = {0:f}".format(self.atom._atom_fvl[self.linecur]))
            print("------------------------------------------------------------")
        elif event.key == 'w':
            self.write_data()
        elif event.key == 'n':
            self.lineidx += 1
            self.next_spectrum()
        elif event.key == 'b':
            self.lineidx -= 1
            self.next_spectrum()
        elif event.key == ']':
            self.next_element(1)
            self.next_spectrum()
        elif event.key == '[':
            self.next_element(-1)
            self.next_spectrum()
        elif event.key == '.':
            self.next_element(1, ion=True)
            self.next_spectrum()
        elif event.key == ',':
            self.next_element(-1, ion=True)
            self.next_spectrum()
        elif event.key == 'q':
            if self._changes:
                print("WARNING: There are unsaved changes!!")
                print("Press q again to exit")
                self._changes = False
            else:
                sys.exit()
        self.canvas.draw()

    def next_element(self, pm, ion=False):
        if ion == True:
            arrsrt = np.core.defchararray.add(self.atom._atom_atm, self.atom._atom_ion)
        else:
            arrsrt = self.atom._atom_atm
        unq, idx = np.unique(arrsrt, return_index=True)
        unq = unq[idx.argsort()]
        nxt = np.where(unq == arrsrt[self.linecur])[0][0] + pm
        if nxt >= unq.size:
            nxt = 0
        ww = np.where(arrsrt == unq[nxt])[0]
        self.lineidx = ww[0]
        return

    def next_spectrum(self):
        self.get_current_line()
        # Update the wavelength range of the spectrum being plot
        self.update_waverange()
        # See if any regions need to be loaded
        self.prop._regions[:] = 0
        idtxt = "{0:s}_{1:s}_{2:.1f}".format(self.atom._atom_atm[self.linecur].strip(),
                                             self.atom._atom_ion[self.linecur].strip(),
                                             self.atom._atom_wvl[self.linecur])
        tstnm = self.prop._outp + "_" + idtxt + "_reg.dat"
        if os.path.exists(tstnm):
            wv, reg = np.loadtxt(tstnm, unpack=True, usecols=(0, 3))
            mtch = np.in1d(self.prop._wave, wv, assume_unique=True)
            self.prop._regions[np.where(mtch)] = reg.copy()
            self.ax.set_xlim([np.min(wv), np.max(wv)])
        # Other stuff
        self.canvas.draw()
        return

    def update_regions(self):
        self.prop._regions[self._start:self._end] = self._addsub
        std = np.std(self.prop._flux[self._start:self._end])
        med = np.median(self.prop._flue[self._start:self._end])
        mad = 1.4826 * np.median(
            np.abs(self.prop._flux[self._start:self._end] - np.median(self.prop._flux[self._start:self._end])))
        print(np.mean(self.prop._wave[self._start:self._end]), std / med, mad / med)

    def update_waverange(self):
        self.get_current_line()
        wcen = self.atom._atom_wvl[self.linecur] * (1.0 + self.prop._zabs)
        xmn = wcen * (1.0 - self.veld / 299792.458)
        xmx = wcen * (1.0 + self.veld / 299792.458)
        self.ax.set_xlim([xmn, xmx])
        medval = np.median(self.prop._flux[~np.isnan(self.prop._flux)])
        self.ax.set_ylim([-0.1 * medval, 1.1 * medval])
        # print("Updated wavelength range:", xmn, xmx)
        self.canvas.draw()

    def write_data(self):
        # Plot the lines
        xmn, xmx = self.ax.get_xlim()
        wsv = np.where((self.prop._wave > xmn) & (self.prop._wave < xmx))
        idtxt = "{0:s}_{1:s}_{2:.1f}".format(self.atom._atom_atm[self.linecur].strip(),
                                             self.atom._atom_ion[self.linecur].strip(),
                                             self.atom._atom_wvl[self.linecur])
        outnm = self.prop._outp + "_" + idtxt + "_reg.dat"
        sclfct = 1.0
        np.savetxt(outnm, np.transpose((self.prop._wave[wsv], sclfct * self.prop._flux[wsv] * self.prop._cont[wsv],
                                        sclfct * self.prop._flue[wsv] * self.prop._cont[wsv], self.prop._regions[wsv])))
        print("Saved file:")
        print(outnm)


def initialise(spec, detns, lines):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    # Ignore lines outside of wavelength range
    wmin, wmax = np.min(prop._wave) / (1.0 + prop._zabs), np.max(prop._wave) / (1.0 + prop._zabs)

    spec = Line2D(np.arange(), prop._flux, linewidth=1, linestyle='solid', color='k', drawstyle='steps', animated=True)

    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    ax.add_line(spec)
    reg = SelectRegions(fig.canvas, ax, spec, prop, atom)

    ax.set_title("Press '?' to list the available options")
    # ax.set_xlim((prop._wave.min(), prop._wave.max()))
    ax.set_ylim((0.0, prop._flux.max()))
    plt.show()


if __name__ == '__main__':
    # Start with a linelist and detections
    dir = "/Users/rcooke/Work/Research/vmp_DLAs/observing/WHT_ISIS_2019B/N1/wht_isis_blue_U/"
    spec = np.loadtxt(dir+"spec.dat")
    detns = np.loadtxt(dir+"detns.dat")
    lines = np.loadtxt(dir+"lines.dat")
    initialise(spec, detns, lines)
