import sys

class colors:
	"""
	Create coloured text for messages printed to screen.

	For further details on colours see the following example:
	http://ascii-table.com/ansi-escape-sequences.php
	"""
	# Start and end coloured text
	start = "\x1B["
	end   = "\x1B[" + "0m"
	# Clear Backgrounds
	black_CL  = "1;30m"
	yellow_CL = "1;33m"
	blue_CL   = "1;34m"
	green_CL  = "1;32m"
	red_CL    = "1;31m"
	# Coloured Backgrounds
	white_RD  = "1;37;41m"
	white_GR  = "1;37;42m"
	white_BK  = "1;37;40m"
	white_BL  = "1;37;44m"
	black_YL  = "1;37;43m"
	yellow_BK = "1;33;40m"

	def disable(self):
		self.black_CL = ''
		self.red_CL   = ''
		self.green_CL = ''
		self.black_RD = ''
		self.black_GR = ''

def armedheader(prognm):
	header = "##  "
	header += colors.start + colors.white_GR + "ARMED : "
	header += "Automated Reduction and Modelling of Echelle Data v1.0" + colors.end + "\n"
	header += "##  "
	header += "Usage : "
	header += "python %s [options] filelist" % (prognm)
	return header

def pypitheader(prognm):
	header = "##  "
	header += colors.start + colors.white_GR + "PYPIT : "
	header += "The Python Spectroscopic Data Reduction Pipeline v1.0" + colors.end + "\n"
	header += "##  "
	header += "Usage : "
	header += "python %s [options] filelist" % (prognm)
	return header

def signal_handler(signalnum, handler):
	if signalnum == 2:
		info("Ctrl+C was pressed. Ending processes...")
		sys.exit()

def error(msg):
	premsg="\n"+colors.start + colors.white_RD + "[ERROR]   ::" + colors.end + " "
	print >>sys.stderr,premsg+msg
	sys.exit()

def info(msg):
	premsg=colors.start + colors.green_CL + "[INFO]    ::" + colors.end + " "
	print >>sys.stderr,premsg+msg

def info_update(msg,last=False):
	premsg="\r" + colors.start + colors.green_CL + "[INFO]    ::" + colors.end + " "
	if last:
		print >>sys.stderr,premsg+msg
	else:
		print >>sys.stderr,premsg+msg,

def test(msg):
	premsg=colors.start + colors.white_BL   + "[TEST]    ::" + colors.end + " "
	print >>sys.stderr,premsg+msg

def warn(msg):
	premsg=colors.start + colors.red_CL   + "[WARNING] ::" + colors.end + " "
	print >>sys.stderr,premsg+msg

def bug(msg):
	premsg=colors.start + colors.white_BK   + "[BUG]     ::" + colors.end + " "
	print >>sys.stderr,premsg+msg

def work(msg):
	premsgP=colors.start + colors.black_CL   + "[WORK IN ]::" + colors.end + "\n"
	premsgS=colors.start + colors.yellow_CL   + "[PROGRESS]::" + colors.end + " "
	print >>sys.stderr,premsgP+premsgS+msg

def prindent(msg):
	premsg = "             "
	print >>sys.stderr,premsg+msg

def input():
	premsg=colors.start + colors.blue_CL  + "[INPUT]   ::" + colors.end + " "
	return premsg

def newline():
	return "\n             "

def indent():
	return "             "
