import sys
import getopt
# Import PYPEIT routines
import arload
import armsgs as msgs

last_updated = "Last updated 4th June 2015"

def usage(prognm):
	print "\n#####################################################################"
	print msgs.pypeitheader()
	print "##  -----------------------------------------------------------------"
	print "##  Options: (default values in brackets)"
	print "##   -c or --cpus      : (all) Number of cpu cores to use"
	print "##   -h or --help      : Print this message"
	print "##   -v or --verbose   : (2) Level of verbosity (0-2)"
	print "##  -----------------------------------------------------------------"
	print "##  %s" % last_updated
	print "#####################################################################\n"
	sys.exit()

class ClassMain:

	def __init__(self, argflag, quick=False):
		"""
		argflag :: A list of arguments and flags for the reduction
		quick   :: Results in a quick reduction, if a quicker (but less accurate) reduction exists
		---------------------------------------------------
		
		"""

		# Set parameters
		self._argflag = argflag

		# First send all signals to messages to be dealt
		# with (i.e. someone hits ctrl+c)
		signal.signal(signal.SIGINT, msgs.signal_handler)

		# Ignore all warnings given by python
		warnings.resetwarnings()
		warnings.simplefilter("ignore")

		# Record the starting time
		self._tstart=time.time()

		# Load the Input file
		self._parlines, self._datlines, self._spclines = arload.load_input(self)

		# Determine the type of data that is being reduced
		msgs.work("TO BE DONE")

		# If a quick reduction has been requested, make sure the requested pipeline
		# is the quick implementation (if it exists), otherwise run the standard pipeline.
		if quick:
			# Change to a "quick" settings file
			msgs.work("TO BE DONE")

		# Load the Spectrograph settings
		self._spect = arload.load_spect(self)

		# Load any changes to the spectrograph settings
		self._spect = arload.load_spect(self, lines=self._spclines)

		# Load the important information from the fits headers
		self._fitsdict = arload.load_headers(self)

		# Load the list of standard stars
		self._standardStars = arload.load_standards(self)

		# Reduce the data!
		msgs.work("Send the data away to a definition of the type of reduction needed")
		status = 0
		if quick:
			msgs.work("define what is needed here")
		else:
			success = self.armed()
		if status==0:
			msgs.info("Reduction complete")

	def armed(self):
		success = False
		# Insert series of reduction steps here
		return success

if __name__ == "__main__":
	prognm = sys.argv[0]
	debug = True
	quick = False

	# Load options from command line
	try:
		opt,arg=getopt.getopt(argv[1:],'hqc:v:', ['help', 
                                                 'quick'
                                                ])
	except getopt.GetoptError, err:
		msgs.error(err.msg)
		usage(prognm)
	for o,a in opt:
		if   o in ('-h', '--help')      : usage(argflag)
		elif o in ('-q', '--quick')     : quick = True
#		elif o in ('-c', '--cpus')      : argflag['run']['ncpus']     = a
#		elif o in ('-v', '--verbose')   : argflag['out']['verbose']   = int(a)

	if debug:
		argflag = arload.optarg(sys.argv, last_updated)
		ClassMain(argflag, quick=quick)
	else:
		try:
			argflag = arload.optarg(sys.argv, last_updated)
			ClassMain(argflag, quick=quick)
		except Exception:
			# There is a bug in the code, print the file and line number of the error.
			et, ev, tb = sys.exc_info()
			while tb:    
				co = tb.tb_frame.f_code
				filename = str(co.co_filename)
				line_no =  str(traceback.tb_lineno(tb))
				tb = tb.tb_next
			filename=filename.split('/')[-1]
			msgs.bug("There appears to be a bug on Line "+line_no+" of "+filename+" with error:"+msgs.newline()+str(ev)+msgs.newline()+"---> please contact the author")
	