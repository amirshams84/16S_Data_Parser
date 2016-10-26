# ################################### IMPORTS ################################ #
import sys
import argparse
import time
import os
import multiprocessing
import itertools
import errno
import signal
import datetime
import logging as log
import subprocess
import traceback
import shutil
import collections
import math
import platform


# ################################## GLOBAL DESC ############################# #

fasta_extensions = ['.fasta', '.fa', '.fna']
fastq_extensions = ['.fastq', '.fq']
CHECK_MARK = "OK"
FAILED_MARK = ":("
# ################################### OBJECTS ################################ #


class mothur_process:
	def __init__(self, mothur_input_dictionary):
		for var_name, var_value in mothur_input_dictionary.items():
			setattr(self, var_name, var_value)

	def build_mothur_command(self):
		space = " "
		string = ''
		if hasattr(self, 'nohup_in'):
			string += self.nohup_in + space
		string += self.mothur_exec_path + space + '"#set.dir(output=' + self.outputdir + ');' + space
		string += 'set.logfile(name=mothur.' + self.command + '.logfile);' + space
		string += 'set.current(processors=' + self.processors + ');' + space
		string += self.command + '('
		if hasattr(self, 'parameters') and len(self.parameters) > 0:
			for each_element in self.parameters:
				string += each_element
		string += ');"'
		if hasattr(self, 'nohup_out'):
			string += space + self.nohup_out
		if hasattr(self, 'pid_file'):
			string += ' echo $! > ' + self.pid_file
		report(string)
		print string
		return string

	def execute_mothur_command(self):
		exec_dict = {}
		exec_dict = self.execute([self.build_mothur_command()])
		if exec_dict['exitCode'] != 0:
			print "ERROR occurred!!!!"
			return (False, exec_dict)
		else:
			return(True, exec_dict)

	def execute(self, command, ** kwargs):
		assert isinstance(command, list), "Expected 'command' parameter to be a list containing the process/arguments to execute. Got %s of type %s instead" % (command, type(command))
		assert len(command) > 0, "Received empty list of parameters"
		retval = {
					"exitCode": -1,
					"stderr": u"",
					"stdout": u"",
					"execTime": datetime.timedelta(0),
					"command": None,
					"pid": None
				}
		retval["command"] = command
		log.info("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(command)))
		cwd = kwargs.get("cwd", os.getcwd())
		#user = kwargs.get("user", os.getuid)
		sheel = kwargs.get("shell", True)
		startDatetime = datetime.datetime.now()
		myPopy = subprocess.Popen(command, cwd=cwd, preexec_fn=os.seteuid(os.getuid()), shell=sheel, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		retval["pid"] = myPopy.pid
		log.debug("::singleProcessExecuter > Command \"%s\" got pid %s" % (" ".join(command), myPopy.pid))
		try:
			retval["stdout"], retval["stderr"] = myPopy.communicate()
			myPopy.wait()
		except OSError, osErr:
			log.debug("::singleProcessExecuter > Got %s %s in myPopy.communicate() when trying get output of command %s. It is probably a bug (more info: http://bugs.python.org/issue1731717)" % (osErr, type(osErr), command[0]))
		except Exception, e:
			log.warn("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s" % (type(e), e, " ".join(command)))
			log.debug("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s. Showing traceback:\n%s" % (type(e), e, " ".join(command), traceback.format_exc()))
			raise
		retval["exitCode"] = myPopy.returncode
		retval["execTime"] = datetime.datetime.now() - startDatetime
		return retval


class usearch_process:

	def __init__(self, usearch_input_dictionary):
		for var_name, var_value in usearch_input_dictionary.items():
			setattr(self, var_name, var_value)

	def build_usearch_command(self):
		space = " "
		string = ''
		if hasattr(self, 'nohup_in'):
			string += self.nohup_in + space
		string += self.usearch_exec_path + space
		for each_element in self.commandline:
			string += each_element
		if hasattr(self, 'nohup_out'):
			string += space + self.nohup_out
		if hasattr(self, 'pid_file'):
			string += ' echo $! > ' + self.pid_file
		report(string)
		print string
		return string

	def execute_usearch_command(self):
		exec_dict = {}
		exec_dict = self.execute([self.build_usearch_command()])
		if exec_dict['exitCode'] != 0:
			print "ERROR occurred!!!"
			return (False, exec_dict)
		else:
			#print "Execution started!!!"
			return(True, exec_dict)

	def execute(self, command, ** kwargs):
		assert isinstance(command, list), "Expected 'command' parameter to be a list containing the process/arguments to execute. Got %s of type %s instead" % (command, type(command))
		assert len(command) > 0, "Received empty list of parameters"
		retval = {
				"exitCode": -1,
				"stderr": u"",
				"stdout": u"",
				"execTime": datetime.timedelta(0),
				"command": None,
				"pid": None
				}
		retval["command"] = command
		log.info("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(command)))
		cwd = kwargs.get("cwd", os.getcwd())
		#user = kwargs.get("user", os.getuid)
		sheel = kwargs.get("shell", True)
		startDatetime = datetime.datetime.now()
		myPopy = subprocess.Popen(command, cwd=cwd, preexec_fn=os.seteuid(os.getuid()), shell=sheel, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		retval["pid"] = myPopy.pid
		log.debug("::singleProcessExecuter > Command \"%s\" got pid %s" % (" ".join(command), myPopy.pid))
		try:
			retval["stdout"], retval["stderr"] = myPopy.communicate()
			myPopy.wait()
		except OSError, osErr:
			log.debug("::singleProcessExecuter > Got %s %s in myPopy.communicate() when trying get output of command %s. It is probably a bug (more info: http://bugs.python.org/issue1731717)" % (osErr, type(osErr), command[0]))
		except Exception, e:
			log.warn("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s" % (type(e), e, " ".join(command)))
			log.debug("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s. Showing traceback:\n%s" % (type(e), e, " ".join(command), traceback.format_exc()))
			raise
		retval["exitCode"] = myPopy.returncode
		retval["execTime"] = datetime.datetime.now() - startDatetime
		return retval

# ################################### MAIN ################################### #


def main(argv):
	report_string = ''
	# ############################# PARSE INPUT ARGUMENTS
	parser = argparse.ArgumentParser()
	
	# ########### INPUT/OUTPUT GROUP
	input_output = parser.add_argument_group('Input/Output parameters')
	input_output.add_argument("--inputdir", help="Input file directory", required=True)
	input_output.add_argument("--outputdir", help="Output file directory", required=True)
	input_output.add_argument("--execdir", help="Executables file directory", required=True)
	
	# ########### TAXONOMY GROUP
	taxonomy = parser.add_argument_group('taxonomy parameters')
	taxonomy.add_argument("--reference", help="16S reference fasta file", required=True)
	taxonomy.add_argument("--taxonomy", help="16S reference taxonomy file", required=True)
	taxonomy.add_argument("--chimera_reference", help="16S chimera_removal reference file", required=True)
	
	# ########### SEQUENCING MODE GROUP
	sequencing = parser.add_argument_group('Sequencing parameters')
	sequencing.add_argument("--filetype", help="Fastq or fasta file type", required=True)
	sequencing.add_argument("--filesystem", help="System of input file['paired', 'single', 'mixed']", required=True)
	sequencing.add_argument("--ngs", help="['illumina', 'sanger', 'pacbio', 'iontorrent', 'auto']", required=True)
	
	# ########### GENERAL PARAMETERS GROUP
	general = parser.add_argument_group('general parameters')
	general.add_argument("--relabel", help="in mixed filetype relabel each headers based on ; in other filetypes relabel each headers based on filename", required=True)
	general.add_argument("--separator", help="in mixed filetype separator is on each headers ; in other filetypes separators is applied on each filename", required=True)
	general.add_argument("--name", help="name of output files", required=True)
	general.add_argument("--processors", help="number of processors assigned", required=True)
	general.add_argument("--design", help="Design file for grouping the data", action='store')

	# ########### FILTERING PARAMETERS GROUP
	filtering = parser.add_argument_group('Filtering parameteres')
	filtering.add_argument("--expected_error_rate", help="filtering by expected_error rate(0.25-0.5-0.75-1.0-2.0)(false for disable)", required=True)
	filtering.add_argument("--maxambig", help="Filtering by ambiguous base limit(false for disable)", required=True)
	filtering.add_argument("--maxhomop", help="Filtering by homopolymer limit(false for disable)", required=True)
	filtering.add_argument("--minqscore", help="Filtering by phred score limit(false for disable)", required=True)
	filtering.add_argument("--maxlength", help="Filtering by length limit(digit-auto-false)", required=True)
	filtering.add_argument("--precluster", help="Filtering by preclustering method [Huse et all](true for active and False for disable)", required=True)
	
	# ########### PIPELINE PARAMETERS GROUP
	pipeline = parser.add_argument_group('Pipeline parameters')
	pipeline.add_argument("--pipeline", help="Pipeline methods(uparse-uclust-upgma-supervised-unifrac)", required=True)
	pipeline.add_argument("--min_abund_size", help="Minimum abundance of cluster size", required=True)
	pipeline.add_argument("--identity", help="Identity percentage", required=True)
	pipeline.add_argument("--maxhits", help="Number of reporting hits", required=True)
	
	args = parser.parse_args()
	# #############################BEURACRATICS PROCEDURES
	report_string += "######################################################################################################################################\n"
	print "######################################################################################################################################"
	report_string += "16S SIMPLE PARSER EXECUTION HAS INITIATED" + '\n'
	print "16S SIMPLE PARSER EXECUTION HAS INITIATED"
	report_string += "Initiation time: " + time.strftime("%Y-%m-%d %H:%M:%S") + '\n'
	print "Initiation time: ", time.strftime("%Y-%m-%d %H:%M:%S")
	report_string += "###################################################################" + '\n'
	print "###################################################################"
	report_string += "INFORMATION ABOUT THE ENVIRONMENT, EXECUTABLES AND PROVIDED DATA" + '\n'
	print "INFORMATION ABOUT THE ENVIRONMENT, EXECUTABLES AND PROVIDED DATA"
	report_string += "###################################################################" + '\n'
	print "###################################################################"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "COMMAND LINE:"
	report_string += "COMMAND LINE:\n"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	commandline_string = 'python ' + ' '.join(sys.argv) + '\n'
	print commandline_string
	report_string += commandline_string + '\n'
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report_string += "ARGUMENTS:" + '\n'
	print "ARGUMENTS:"
	report_string += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" + '\n'
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	
	# ########### OUTPUT DIRECTORIES
	if args.outputdir[-1] != '/':
		args.outputdir += '/'
	if isPathExist(args.outputdir) is False:
		print "[--outputdir]: Output file directory has Access/Exist issue!!!"
		report_string += "[--outputdir]: Output file directory has Access/Exist issue!!!" + '\n'
		print "ABORTING!!!"
		sys.exit(2)
	else:
		print "[--outputdir]: Output file directory is: ", args.outputdir
		report_string += "[--outputdir]: Output file directory is: " + args.outputdir
	global report_file
	report_file = args.outputdir + "16S_simple_parser_report.txt"
	check_it_and_remove_it(report_file, True)
	report(report_string)

	# ########### INPUT DIRECTORIES
	if isPathExist(args.inputdir) is False:
		error("[--inputdir]: Input file directory has Access/Exist issue!!!")
		print "[--inputdir]: Input file directory has Access/Exist issue!!!"
		error("ABORTING!!!")
		print "ABORTING!!!"
		sys.exit(2)
	else:
		report("[--inputdir]: Input file directory is: " + args.inputdir)
		print "[--inputdir]: Input file directory is: ", args.inputdir

	# ########### REFERENCE FASTA FILE
	if isFileExist(args.reference) is False:
		print "[--reference]: 16S reference_fasta file has Access/Exist issue!!!"
		error("[--reference]: 16S reference_fasta file has Access/Exist issue!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--reference]: 16S reference_fasta file is: " + args.reference)
		print "[--reference]: 16S reference_fasta file is: ", args.reference

	# ########### REFERENCE TAXONOMY FILE
	if isFileExist(args.taxonomy) is False:
		print "[--taxonomy]: 16S reference_taxonomy file has Access/Exist issue!!!"
		error("[--taxonomy]: 16S reference_taxonomy file has Access/Exist issue!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--taxonomy]: 16S reference_taxonomy file is: " + args.taxonomy)
		print "[--taxonomy]: 16S reference_taxonomy file is: ", args.taxonomy

	# ########### CHIMERA REFERENCE FASTA FILE
	if isFileExist(args.chimera_reference) is False:
		print "[--chimera_reference]: 16S chimera reference_fasta file has Access/Exist issue!!!"
		error("[--chimera_reference]: 16S chimera reference_fasta file has Access/Exist issue!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--chimera_reference]: 16S chimera reference_fasta file is: " + args.chimera_reference)
		print "[--chimera_reference]: 16S chimera reference_fasta file is: ", args.chimera_reference

	# ########### FILE TYPE CHECKING
	if args.filetype.lower() not in ['fasta', 'fastq']:
		error("[--filetype]: The only supported file types is (fastq or fasta)!!!")
		print "[--filetype]: The only supported file types is (fastq or fasta)!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--filetype]: File type is: " + args.filetype.lower())
		print "[--filetype]: File type is: ", args.filetype.lower()
		args.filetype = args.filetype.lower()

	# ########### filesystem MODE CHECKING
	if args.filesystem.lower() not in ['paired', 'single', 'mixed']:
		error("[--filesystem]: The supported file system is ('paired', 'single', 'mixed')!!!")
		print "[--filesystem]: The supported file system is ('paired', 'single', 'mixed')!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--filesystem]: The file system is: " + args.filesystem)
		print "[--filesystem]: The file system is: ", args.filesystem
		args.filesystem = args.filesystem.lower()

	# ########### NGS MODE CHECKING
	if args.ngs.lower() not in ['illumina', 'sanger', 'pacbio', 'iontorrent', 'auto']:
		error("[--ngs]: The supported ngs modes is ('illumina', 'sanger', 'pacbio', 'iontorrent', 'auto')")
		print "[--ngs]: The supported ngs modes is ('illumina', 'sanger', 'pacbio', 'iontorrent', 'auto')"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--ngs]: The ngs mode is: " + args.ngs.lower())
		print "[--ngs]: The ngs mode is: ", args.ngs.lower()
		if args.ngs.lower() == 'auto':
			args.ngs = 'illumina'
		ngs_qscore = set_qscore(args.ngs.lower())

	# ########### PROCESSORS CHECKING
	if isfloat(args.processors) is False:
		error("[--processors]: Number of processors should be a digit!!!")
		print "[--processors]: Number of processors should be a digit!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif float(args.processors) < 1.0:
		error("[--processors]: Atleast one processor is needed!!!")
		print "[--processors]: Atleast one processor is needed!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif float(args.processors) > float(multiprocessing.cpu_count()):
		error("[--processors]: Number of assigned processors is higher than available CPUs: " + str(multiprocessing.cpu_count()))
		print "[--processors]: Number of assigned processors is higher than available CPUs: ", str(multiprocessing.cpu_count())
		error("Please select a value lower than " + str(multiprocessing.cpu_count()))
		print "Please select a value lower than ", str(multiprocessing.cpu_count())
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--processors]: The Processors number is: " + args.processors)
		print "[--processors]: The Processors number is: ", args.processors

	# ########### SEPARATOR CHECKING
	if args.separator.lower() == 'false':
		report("[--separator]: Read headers separator is disabled")
		print "[--separator]: Read headers separator is disabled"
		args.separator = False
	else:
		report("[--separator]: Read headers separator is: " + args.separator)
		print "[--separator]: Read headers separator is: ", args.separator

	# ########### relabel CHECKING
	if args.relabel.lower() not in ['true', 'false']:
		report("[--relabel]: the only supported value for relabel is (true, false)!!!")
		print "[--relabel]: the only supported value for relabel is (true, false)!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
	elif args.relabel.lower() == 'false':
		report("[--relabel]: relabel headers is disabled")
		print "[--relabel]: relabel headers is disabled"
		args.relabel = False
	elif args.relabel.lower() == 'true':
		report("[--relabel]: Relabel headers is enabled.")
		print "[--relabel]: Relabel headers is enabled."
		args.relabel = True

	# ########### NAME CHECKING
	if args.name and args.name.strip():
		report("[--name]: Specified name for output files is: " + args.name)
		print "[--name]: Specified name for output files is: ", args.name
	else:
		error("[--name]: Specified name for output files can not be empty!!!")
		print "[--name]: Specified name for output files can not be empty!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)

	# ########### DESIGN FILE CHECKING
	if not args.design:
		report("[--design]: Design file is not specified")
		print "[--design]: Design file is not specified"
		args.design = False
	elif isFileExist(args.design) is False:
		print "[--design]: Design file has Access/Exist issue!!!"
		error("[--design]: Design file has Access/Exist issue!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--design]: Design file is: " + args.design)
		print "[--design]: Design file is: ", args.design

	# ########### EXPECTED ERROR RATE CHECKING
	if args.expected_error_rate.lower() == 'false':
		report("[--expected_error_rate]: Expected error rate is disabled")
		print "[--expected_error_rate]: Expected error rate is disabled"
		args.expected_error_rate = False
	elif args.expected_error_rate and isfloat(args.expected_error_rate) is False:
		print "[--expected_error_rate]: please provide a float number for expected error rate(0.25, 0.5, 1.0, 2.0)!!!"
		error("[--expected_error_rate]: please provide a float number for expected error rate(0.25, 0.5, 1.0, 2.0)!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif args.expected_error_rate and args.filetype == 'fasta':
		print "[--expected_error_rate]: Expected error rate is only applicable on fastq files, your filetype is fasta"
		report("[--expected_error_rate]: Expected error rate is only applicable on fastq files, your filetype is fasta")
		report("[--expected_error_rate]: Expected error rate is disabled")
		print "[--expected_error_rate]: Expected error rate is disabled"
		args.expected_error_rate = False
	elif args.expected_error_rate not in ['0.25', '0.5', '1.0', '2.0']:
		print "[--expected_error_rate]: Invalid value for expected error rate"
		print "[--expected_error_rate]: The only acceptable expected error rates is 0.25, 0.5, 1.0, 2.0"
		error("[--expected_error_rate]: Invalid value for expected error rate\nThe only acceptable expected error rates is 0.25, 0.5, 1.0, 2.0")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--expected_error_rate]: Expected error rate for fastq filtering is activated")
		print "[--expected_error_rate]: Expected error rate for fastq filtering is activated"
		report("[--expected_error_rate]: Expected error rate is: " + args.expected_error_rate)
		print "[--expected_error_rate]: Expected error rate is: ", args.expected_error_rate

	# ########## AMBIGUOUS BASE FILTERING
	if args.maxambig.lower() == 'false':
		report("[--maxambig]: Ambiguous base filtering mode is disabled")
		print "[--maxambig]: Ambiguous base filtering mode is disabled"
		args.maxambig = False

	elif args.maxambig and isfloat(args.maxambig) is False:
		error("[--maxambig]: Maximum Number of ambiguous base should be a digit")
		print "[--maxambig]: Maximum Number of ambiguous base should be a digit"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--maxambig]: Maximum Number of ambiguous base in each read is: " + args.maxambig)
		print "[--maxambig]: Maximum Number of ambiguous base in each read is: ", args.maxambig

	# ########### HOMOPOLYMER BASE FILTERING
	if args.maxhomop.lower() == 'false':
		report("[--maxhomop]: Homopolymer filtering mode is disabled")
		print "[--maxhomop]: Homopolymer filtering mode is disabled"
		args.maxhomop = False

	elif args.maxhomop and isfloat(args.maxhomop) is False:
		error("[--maxhomop]: Maximum Number of Homopolymer base should be a digit")
		print "[--maxhomop]: Maximum Number of Homopolymer base should be a digit"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--maxhomop]: Maximum Number of Homopolymer base in each read is: " + args.maxhomop)
		print "[--maxhomop]: Maximum Number of Homopolymer base in each read is: ", args.maxhomop

	# ########### QUALITY SCORE FILTERING
	if args.minqscore.lower() == 'false':
		report("[--minqscore]: fastq quality score filtering is disabled")
		print "[--minqscore]: fastq quality score filtering is disabled"
		args.minqscore = False
	elif args.minqscore and isfloat(args.minqscore) is False:
		print "[--minqscore]: please provide fastq quality score filtering as an integer between 0 and 41"
		error("[--minqscore]: please provide fastq quality score filtering as an integer between 0 and 41")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif args.minqscore and args.filetype == 'fasta':
		print "[--minqscore]: Fastq quality score filtering is only applicable on fastq files, your filetype(--filetype) is fasta"
		report("[--minqscore]: Fastq quality score filtering is only applicable on fastq files, your filetype(--filetype) is fasta")
		report("[--minqscore]: fastq quality score filtering is disabled")
		print "[--minqscore]: fastq quality score filtering is disabled"
		args.minqscore = False
	elif float(args.minqscore) < 0 or float(args.minqscore) > 41:
		print "[--minqscore]: Invalid value for fastq quality score filtering"
		print "[--minqscore]: The only acceptable fastq quality score filtering is between 0 and 41"
		error("[--minqscore]: Invalid value for fastq quality score filtering\nThe only acceptable fastq quality score filtering(--minqscore) is between 0 and 41")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--minqscore]: fastq quality score filtering is: " + args.minqscore)
		print "[--minqscore]: fastq quality score filtering is: ", args.minqscore

	# ########### GLOBAL TRIMMING LENGTH
	if args.maxlength.lower() == 'false':
		report("[--maxlength]: Global trimming length is disabled")
		print "[--maxlength]: Global trimming length is disabled"
		args.maxlength = False
	elif args.maxlength.lower() == 'auto':
		report("[--maxlength]: Global trimming length will calculated automatically")
		print "[--maxlength]: Global trimming length will calculated automatically"
		args.maxlength = 'auto'
	elif args.maxlength.lower() not in ['auto', 'false'] and isfloat(args.maxlength.lower()) is False and int(args.maxlength.lower()) < 0:
		report("[--maxlength]: Global trimming length should be a digit greater than zero")
		print "[--maxlength]: Global trimming length should be a digit greater than zero"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--maxlength]: Global trimming length is: " + args.maxlength)
		print "[--maxlength]: Global trimming length is: ", args.maxlength

	# ########## PRECLUSTER PROCESS
	if args.precluster.lower() == 'false':
		report("[--precluster]: Preclustering mode is disabled")
		print "[--precluster]: Preclustering mode is disabled"
		args.precluster = False
	elif args.precluster.lower() == 'true':
		report("[--precluster]: Preclustering mode enabled")
		print "[--precluster]: Preclustering mode enabled"
		args.precluster = True
	else:
		report("[--precluster]: Preclustering mode should be true or false")
		print "[--precluster]: Preclustering mode should be true or false"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)

	# ########## PIPELINE PROCESS
	if args.pipeline.lower() not in ['uparse', 'uclust', 'upgma', 'supervised', 'unifrac']:
		error("[--pipeline]: The only supported pipeline is ('uparse', 'uclust', 'upgma', 'supervised', 'unifrac')!!!")
		print "[--pipeline]: The only supported pipeline is ('uparse', 'uclust', 'upgma', 'supervised', 'unifrac')!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--pipeline]: The pipeline is: " + args.pipeline.lower())
		print "[--pipeline]: The pipeline is: ", args.pipeline.lower()
		args.pipeline = args.pipeline.lower()

	# ########################
	if isfloat(args.identity) is False:
		error("[--identity]: The identity percentage of alignment should be integer less than 100 and greather than 1")
		print "[--identity]: The identity percentage of alignment should be integer less than 100 and greather than 1"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif float(args.identity) < 1.0 or float(args.identity) > 101.0:
		error("[--identity]: The identity percentage of alignment should be integer less than 100 and greather than 1")
		print "[--identity]: The identity percentage of alignment should be integer less than 100 and greather than 1"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--identity]: The identity percentage of alignment is: " + args.identity)
		print "[--identity]: The identity percentage of alignment is: ", args.identity
		args.identity = float(args.identity) / 100.0

	# ########################
	if isint(args.maxhits) is False:
		error("[--maxhits]: The number of max hits to be included in abundance calculations should be integer greather than zero")
		print "[--maxhits]: The number of max hits to be included in abundance calculations should be integer greather than zero"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif int(args.maxhits) < 0:
		error("[--maxhits]: The number of max hits to be included in abundance calculations should be integer greather than zero")
		print "[--maxhits]: The number of max hits to be included in abundance calculations should be integer greather than zero"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--maxhits]: The number of best hits is: " + args.maxhits)
		print "[--maxhits]: The number of best hits is: ", args.maxhits

	# ########################
	if isint(args.min_abund_size) is False:
		error("[--min_abund_size]: The Minimum abundance of cluster size should be an integer greater than zero")
		print "[--min_abund_size]: The Minimum abundance of cluster size should be an integer greater than zero"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif int(args.min_abund_size) < 0:
		error("[--min_abund_size]: The Minimum abundance of cluster size should be an integer greater than zero")
		print "[--min_abund_size]: The Minimum abundance of cluster size should be an integer greater than zero"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("[--min_abund_size]: The Minimum abundance of cluster size is: " + args.min_abund_size)
		print "[--min_abund_size]: The Minimum abundance of cluster size is: ", args.min_abund_size

	# ########### EXECECUTIVE DIRECTORY CHECKING
	if args.execdir[-1] != '/':
		args.execdir += '/'
	if isPathExist(args.execdir) is False:
		error("[--execdir]: executables directory has Access/Exist issue!!!")
		print "[--execdir]: executables directory has Access/Exist issue!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)

	print "[--execdir]: The excutives directory path is: ", args.execdir
	report("[--execdir]: the excutives directory path is: " + args.execdir)
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY"
	report("VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY")
	print "###################################################################\n"
	report("###################################################################\n")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("0: ENVIRONMENT")
	print "0: ENVIRONMENT"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "Operating System version is: ", platform.platform()
	report("Operating System version is: " + platform.platform())
	python_version = sys.version.split(' (')[0]
	print "Python version is: ", python_version
	report("Python version is: " + python_version)
	if float(python_version[0:3]) < 2.7:
		error("python version is older than 2.7")
		print "python version is older than 2.7"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("1: MOTHUR")
	print "1: MOTHUR"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	mothur_exec_path = args.execdir + 'mothur'
	if isFileExist(mothur_exec_path) is False:
		error("Your mothur path has Access/Exist issue")
		print "Your mothur path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("mothur execution file is: " + mothur_exec_path)
		print "mothur execution file is: ", mothur_exec_path
		report("Testing mothur executables: ")
		print "Testing mothur executables: "
		flag, stderr = execute_functions(test_mothur, args.processors, args.outputdir, 'multi', 'mothur', mothur_exec_path)
		if flag is False:
			report("[" + FAILED_MARK + "]")
			print "[" + FAILED_MARK + "]"
			print "Execution of mothur failed!!!"
			error("Execution of mothur failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			target_lines = []
			flag = parse_mothur_logfile(stderr, 'mothur v.', target_lines)
			if flag is False:
				print "This keyword is not avalaible: mothur v."
				error("This keyword is not avalaible: mothur v.")
				print "ABORTING!!!"
				error("ABORTING!!!")
				sys.exit(2)
			
			report("mothur executables responding successfully!!!")
			print "mothur executables responding successfully!!!"
			report("version of mothur executables: " + target_lines[0])
			print "version of mothur executables:", target_lines[0]
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("2: USEARCH")
	print "2: USEARCH"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	usearch_exec_path = args.execdir + 'usearch'
	if isFileExist(usearch_exec_path) is False:
		error("Your usearch path has Access/Exist issue")
		print "Your usearch path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("usearch execution file is: " + usearch_exec_path)
		print "usearch execution file is: ", usearch_exec_path
		report("Testing usearch executables: ")
		print "Testing usearch executables:"
		flag, stderr = execute_functions(test_usearch, args.processors, args.outputdir, 'multi', 'usearch', usearch_exec_path)
		if flag is False:
			
			print "Execution of usearch failed!!!"
			error("Execution of usearch failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			report("usearch executables responding successfully!!!")
			print "usearch executables responding successfully!!!"
			report("version of usearch executables: " + stderr.split('\n')[0].split('_')[0])
			print "version of usearch executables:", stderr.split('\n')[0].split('_')[0]
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("3: MAFFT")
	print "3: MAFFT"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	mafft_exec_path = args.execdir + 'mafft/mafft.bat'
	if isFileExist(mafft_exec_path) is False:
		error("Your mafft path has Access/Exist issue")
		print "Your mafft path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("mafft execution file is: " + mafft_exec_path)
		print "mafft execution file is: ", mafft_exec_path
		report("Testing mafft executables: ")
		print "Testing mafft executables:"
		flag, stderr = execute_functions(test_mafft, args.processors, args.outputdir, 'multi', 'mafft', mafft_exec_path)
		if flag is False:
			print "Execution of mafft failed!!!"
			error("Execution of mafft failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			report("mafft executables responding successfully!!!")
			print "mafft executables responding successfully!!!"
			report("version of mafft executables: " + stderr.split('\n')[2])
			print "version of mafft executables:", stderr.split('\n')[2]
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("4: FASTTREE")
	print "4: FASTTREE"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	fasttree_exec_path = args.execdir + 'FastTreeMP'
	if isFileExist(fasttree_exec_path) is False:
		error("Your fasttree_exec_path path has Access/Exist issue")
		print "Your fasttree_exec_path path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("fasttree execution file is: " + fasttree_exec_path)
		print "fasttree execution file is: ", fasttree_exec_path
		report("Testing fasttree executables: ")
		print "Testing fasttree executables:"
		flag, stderr = execute_functions(test_fasttree, args.processors, args.outputdir, 'multi', 'mafft', fasttree_exec_path)
		if flag is False:
			print "Execution of fasttree failed!!!"
			error("Execution of fasttree failed!!!")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			report("fasttree executables responding successfully!!!")
			print "fasttree executables responding successfully!!!"
			report("version of fasttree executables: " + stderr.split('\n')[0].split(',')[0])
			print "version of fasttree executables:", stderr.split('\n')[0].split(',')[0]

	# #####################################################################
	# END OF BEURACRATICS PROCEDURES
	# #####################################################################
	# #####################################################################
	# INPUT DATA PREPARATION
	# #####################################################################
	print "#################################################################"
	report("#################################################################")
	report("STEP1: INPUT DATA PREPARATION")
	print "STEP1: INPUT DATA PREPARATION"
	report("################################################################\n")
	print "#################################################################\n"
	files_path = args.outputdir + args.name + '_files_path_STEP1.txt'
	flag = data_prep(args.inputdir, args.filetype, args.filesystem, files_path, usearch_exec_path, mothur_exec_path, args.processors, args.outputdir)
	if flag is True:
		print "INPUT DATA PREPARATION STEP: PASSED!!!"
		report("INPUT DATA PREPARATION STEP: PASSED!!!")
	else:
		error("INPUT DATA PREPARATION STEP: FAILED!!!")
		print "INPUT DATA PREPARATION STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	# #####################################################################
	# DATA MERGE AND RELABELING
	# #####################################################################
	report("###############################################################")
	print "################################################################"
	report("STEP2: DATA MERGE AND RELABELING")
	print "STEP2: DATA MERGE AND RELABELING"
	report("###############################################################\n")
	print "################################################################\n"
	total_merged_data = []
	flag, total_merged_data = data_merge_relabel(args.filetype, args.filesystem, files_path, args.relabel, args.separator, args.name, ngs_qscore, usearch_exec_path, args.processors, args.outputdir)
	if flag is True:
		if len(total_merged_data) == 1:
			total_merged_fasta = total_merged_data[0]
		else:
			total_merged_fastq = total_merged_data[0]
			total_merged_fasta = total_merged_data[1]
		read_count = fasta_read_count(total_merged_fasta)
		print "TOTAL NUMBER OF READ IS: ", read_count
		report("TOTAL NUMBER OF READ IS: " + read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "DATA MERGE AND RELABELING STEP: PASSED!!!"
		report("DATA MERGE AND RELABELING STEP: PASSED!!!")
	else:
		error("DATA MERGE AND RELABELING STEP: FAILED!!!")
		print "DATA MERGE AND RELABELING STEP: FAILED!!!"
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	# ###########################################################
	# DATA FILTERING
	# ###########################################################
	report("################################################################")
	print "#################################################################"
	report("STEP3: DATA FILTERING")
	print "STEP3: DATA FILTERING"
	report("################################################################\n")
	print "#################################################################\n"
	total_fasta_filtered_file = args.outputdir + args.name + '_total_filtered_STEP3.fasta'
	total_discarded_fasta = args.outputdir + args.name + '_total_discarded_STEP3.fasta'
	# ##########################################
	if args.filetype == 'fastq':
		# ##################
		trunc_len = False
		min_len = False
		maxambig = args.maxambig
		minqscore = args.minqscore
		maxhomop = args.maxhomop
		expected_error_rate = args.expected_error_rate
		fastq_filtering_process(total_merged_fastq, total_merged_fasta, total_fasta_filtered_file, total_discarded_fasta, trunc_len, min_len, maxambig, maxhomop, minqscore, ngs_qscore, expected_error_rate, usearch_exec_path, mothur_exec_path, args.name, args.processors, args.outputdir)
	# ##########################################
	elif args.filetype == 'fasta':
		# ##############
		trunc_len = False
		min_len = False
		maxambig = args.maxambig
		minqscore = args.minqscore
		maxhomop = args.maxhomop
		fasta_filtering_process(total_merged_fasta, total_fasta_filtered_file, total_discarded_fasta, trunc_len, min_len, maxambig, maxhomop, usearch_exec_path, mothur_exec_path, args.name, args.processors, args.outputdir)
	report("DATA FILTERING STEP: PASSED!!!")
	print "DATA FILTERING STEP: PASSED!!!"
	# ###########################################################
	# MEDIAN LENGTH AND GLOBAL TRIMMING
	# ###########################################################
	report("###############################################################")
	print "################################################################"
	report("STEP4: MEDIAN LENGTH AND GLOBAL TRIMMING")
	print "STEP4: MEDIAN LENGTH AND GLOBAL TRIMMING"
	report("###############################################################\n")
	print "################################################################\n"
	median_length = calculate_fasta_median_length(total_fasta_filtered_file, args.name, mothur_exec_path, args.processors, args.outputdir)

	if args.pipeline in ['supervised', 'uclust']:
		print "Global trimming for SUPERVISED and UCLUST will be skipped."
		report("Global trimming for SUPERVISED and UCLUST will be skipped.")
		report("GLOBAL TRIMMING: " + FAILED_MARK)
		print "GLOBAL TRIMMING: ", FAILED_MARK
		fasta_trimmed_file = args.outputdir + args.name + '_global_trimmed_STEP4.fasta'
		copy_file(total_fasta_filtered_file, fasta_trimmed_file)
	else:
		fasta_trimmed_file = args.outputdir + args.name + '_global_trimmed_STEP4.fasta'
		flag = global_trimming_fasta(total_fasta_filtered_file, fasta_trimmed_file, median_length, args.maxlength, args.name, usearch_exec_path, mothur_exec_path, args.processors, args.outputdir)
		if flag is False:
			print "Something is not right with the global_trimming_data."
		else:
			report("THE DATA GLOBAL-TRIMMING STEP: PASSED!!!")
			print "THE DATA GLOBAL-TRIMMING STEP: PASSED!!!"
	
	# ###########################################################
	# ERROR RATE CALCULATION
	# ###########################################################
	report("###############################################################")
	print "################################################################"
	report("STEP4: EXPECTED ERROR RATE/ AVERAGE PROBABILITY/ AVERAGE QUALITY CALCULATION")
	print "STEP4: EXPECTED ERROR RATE/ AVERAGE PROBABILITY/ AVERAGE QUALITY CALCULATION"
	report("###############################################################\n")
	print "################################################################\n"
	total_error_rate = args.outputdir + args.name + '_phred_error_rate_STEP5.txt'
	error_rate_calculation(total_merged_fasta, total_error_rate, mothur_exec_path, args.processors, args.outputdir)
	report("EXPECTED ERROR RATE/ AVERAGE PROBABILITY/ AVERAGE QUALITY CALCULATION STEP: PASSED!!!")
	print "EXPECTED ERROR RATE/ AVERAGE PROBABILITY/ AVERAGE QUALITY CALCULATION STEP: PASSED!!!"
	report("###############################################################")
	print "###############################################################"
	# ########################################################
	# CLUSTERING
	# ########################################################
	ready_to_align_fasta = args.outputdir + args.name + '_ready_for_abundance_STEP7.fasta'
	copy_file(total_merged_fasta, ready_to_align_fasta)

	ready_to_cluster_fasta = args.outputdir + args.name + '_ready_for_cluster_STEP7.fasta'
	copy_file(fasta_trimmed_file, ready_to_cluster_fasta)

	read_count = fasta_read_count(ready_to_align_fasta)
	print "READS FOR ABUNDANCE CALCULATION: ", read_count
	report("READS FOR ABUNDANCE CALCULATION: " + read_count)
	read_count = fasta_read_count(ready_to_cluster_fasta)
	print "READS FOR CENTROIDS DETECTION: ", read_count
	report("READS FOR CENTROIDS DETECTION: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	
	# ########################################################
	# UPARSE
	# ########################################################
	if args.pipeline == 'uparse':
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "UPARSE METHOD:"
		report("UPARSE METHOD:")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		flag, qiimeout, mothurout, biomout, sasout, phylotype_file = uparse_pipeline(ready_to_align_fasta, ready_to_cluster_fasta, args.reference, args.taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, args.name, args.chimera_reference, args.filetype, total_error_rate, args.identity, args.maxhits, args.precluster, median_length, args.min_abund_size, args.processors, args.outputdir)
		
		path, pipabsname, ext = split_file_name(args.pipeline)
		path, refabsname, ext = split_file_name(args.reference)
		
		sas_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_sas.txt'
		check_it_and_remove_it(sas_report)
		os.rename(sasout, sas_report)
		
		biom_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_biom.txt'
		check_it_and_remove_it(biom_report)
		os.rename(biomout, biom_report)
		
		qiime_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_qiime.txt'
		check_it_and_remove_it(qiime_report)
		os.rename(qiimeout, qiime_report)
		
		mothur_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_mothur.txt'
		check_it_and_remove_it(mothur_report)
		os.rename(mothurout, mothur_report)

		phylotype_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_phylotype.txt'
		check_it_and_remove_it(phylotype_report)
		os.rename(phylotype_file, phylotype_report)
		report("THE CLUSTERING PIPELINE STEP: STEP7 PASSED!!!")
		print "THE CLUSTERING PIPELINE STEP: STEP7 PASSED!!!"
		print "DeNOVO ClSUTERING USING UPARSE EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("DeNOVO ClSUTERING USING UPARSE EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))
	
	# ########################################################
	# UPGMA
	# ########################################################
	elif args.pipeline == 'upgma':
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "UPGMA METHOD:"
		report("UPGMA METHOD:")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		flag, qiimeout, mothurout, biomout, sasout, phylotype_file = upgma_pipeline(ready_to_align_fasta, ready_to_cluster_fasta, args.reference, args.taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, args.name, args.chimera_reference, args.filetype, total_error_rate, args.identity, args.maxhits, args.precluster, median_length, args.min_abund_size, args.processors, args.outputdir)
		
		path, pipabsname, ext = split_file_name(args.pipeline)
		path, refabsname, ext = split_file_name(args.reference)
		
		sas_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_sas.txt'
		check_it_and_remove_it(sas_report)
		os.rename(sasout, sas_report)
		
		biom_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_biom.txt'
		check_it_and_remove_it(biom_report)
		os.rename(biomout, biom_report)
		
		qiime_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_qiime.txt'
		check_it_and_remove_it(qiime_report)
		os.rename(qiimeout, qiime_report)
		
		mothur_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_mothur.txt'
		check_it_and_remove_it(mothur_report)
		os.rename(mothurout, mothur_report)

		phylotype_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_phylotype.txt'
		check_it_and_remove_it(phylotype_report)
		os.rename(phylotype_file, phylotype_report)

		print "DeNOVO ClSUTERING USING UPGMA EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("DeNOVO ClSUTERING USING UPGMA EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))
	
	# ########################################################
	# UCLUST
	# ########################################################
	elif args.pipeline == 'uclust':
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "UCLUST METHOD:"
		report("UCLUST METHOD:")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		flag, qiimeout, mothurout, biomout, sasout, phylotype_file = uclust_pipeline(ready_to_align_fasta, ready_to_cluster_fasta, args.reference, args.taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, args.name, args.chimera_reference, args.filetype, total_error_rate, args.identity, args.maxhits, args.precluster, median_length, args.min_abund_size, args.processors, args.outputdir)
		
		path, pipabsname, ext = split_file_name(args.pipeline)
		path, refabsname, ext = split_file_name(args.reference)
		
		sas_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_sas.txt'
		check_it_and_remove_it(sas_report)
		os.rename(sasout, sas_report)
		
		biom_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_biom.txt'
		check_it_and_remove_it(biom_report)
		os.rename(biomout, biom_report)
		
		qiime_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_qiime.txt'
		check_it_and_remove_it(qiime_report)
		os.rename(qiimeout, qiime_report)
		
		mothur_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_mothur.txt'
		check_it_and_remove_it(mothur_report)
		os.rename(mothurout, mothur_report)

		phylotype_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_phylotype.txt'
		check_it_and_remove_it(phylotype_report)
		os.rename(phylotype_file, phylotype_report)

		print "DeNOVO ClSUTERING USING UCLUST EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("DeNOVO ClSUTERING USING UCLUST EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))

	# ########################################################
	# SUPERVISED
	# ########################################################
	elif args.pipeline == 'supervised':
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "SUPERVISED METHOD:"
		report("SUPERVISED METHOD:")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		flag, qiimeout, mothurout, biomout, sasout, phylotype_file = supervised_pipeline(ready_to_align_fasta, ready_to_cluster_fasta, args.reference, args.taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, args.name, args.chimera_reference, args.filetype, total_error_rate, args.identity, args.maxhits, args.precluster, median_length, args.min_abund_size, args.processors, args.outputdir)
		
		path, pipabsname, ext = split_file_name(args.pipeline)
		path, refabsname, ext = split_file_name(args.reference)
		
		sas_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_sas.txt'
		check_it_and_remove_it(sas_report)
		os.rename(sasout, sas_report)
		
		biom_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_biom.txt'
		check_it_and_remove_it(biom_report)
		os.rename(biomout, biom_report)
		
		qiime_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_qiime.txt'
		check_it_and_remove_it(qiime_report)
		os.rename(qiimeout, qiime_report)
		
		mothur_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_mothur.txt'
		check_it_and_remove_it(mothur_report)
		os.rename(mothurout, mothur_report)

		phylotype_report = args.outputdir + args.name + '_' + pipabsname + '_' + refabsname + '_phylotype.txt'
		check_it_and_remove_it(phylotype_report)
		os.rename(phylotype_file, phylotype_report)

		print "SUPERVISED PIPELINE EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("SUPERVISED PIPELINE EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))

	# ########################################################
	# UNIFRAC
	# ########################################################
	elif args.pipeline == 'unifrac':
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "UNIFRAC METHOD:"
		report("UNIFRAC METHOD:")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		flag, unifrac_unweighted_distance_matrix, unifrac_weighted_distance_matrix = unifrac_pipeline(ready_to_align_fasta, ready_to_cluster_fasta, args.reference, args.taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, fasttree_exec_path, args.name, args.chimera_reference, args.filetype, total_error_rate, args.identity, args.maxhits, args.precluster, median_length, args.min_abund_size, args.processors, args.outputdir)
		
		path, pipabsname, ext = split_file_name(args.pipeline)
		path, refabsname, ext = split_file_name(args.reference)
		unifrac_unweighted_report = args.outputdir + args.name + '_' + pipabsname + '_unweighted_distance_matrix.txt'
		check_it_and_remove_it(unifrac_unweighted_report)
		os.rename(unifrac_unweighted_distance_matrix, unifrac_unweighted_report)

		print "UNIFRAC EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("UNIFRAC EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))

# ################################# PIPELINE_FUNCTIONS ####################### #


def uparse_pipeline(ready_to_align_fasta_file, ready_to_cluster_fasta_file, reference, taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, sample_name, chimera_reference, file_type, total_error_rate, similarity, maxhits, precluster, median_length, min_abund_size, processors, outputdir):
	
	# ########################################################
	# REMOVING DUPLICATES
	# ########################################################
	print "DUPLICATES REMOVAL:"
	report("DUPLICATES REMOVAL:")
	derep_fasta = outputdir + sample_name + '_dereplicated_STEP7.fasta'
	check_it_and_remove_it(derep_fasta)
	flag, stderr = execute_functions(usearch_dereplicate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_cluster_fasta_file, derep_fasta)
	if flag is False:
		print "Execution of usearch_dereplicate failed!!!"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	derep_read_count = fasta_read_count(derep_fasta)
	print "The number of reads after removing duplicates reads is: ", derep_read_count
	report("The number of reads after removing duplicates reads is: " + derep_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "DUPLICATES REMOVAL: PASSED!!"
	report("DUPLICATES REMOVAL: PASSED!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# PRECLUSTER
	# ########################################################
	if precluster is True:
		print "PRECLUSTER: "
		report("PRECLUSTER: ")
		maxdiffs = str(int(float(median_length) / 100))
		preclustered_fasta = outputdir + sample_name + '_preclustered_STEP7.fasta'
		check_it_and_remove_it(preclustered_fasta)
		#maxdiffs = '3'
		flag, stderr = execute_functions(usearch_precluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, derep_fasta, preclustered_fasta, maxdiffs)
		if flag is False:
			print "Execution of usearch_precluster failed!!!"
			sys.exit(2)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		preclustered_read_count = fasta_read_count(preclustered_fasta)
		print "The number of reads after precluster is: ", preclustered_read_count
		report("The number of reads after precluster reads is: " + preclustered_read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "PRECLUSTER: PASSED!!!"
		report("PRECLUSTER: PASSED!!!")
	else:
		preclustered_fasta = derep_fasta
		preclustered_read_count = derep_read_count
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# LOW ABUNDANCE REMOVAL
	# ########################################################
	print "LOW ABUNDANCE REMOVAL:"
	report("LOW ABUNDANCE REMOVAL:")
	sorted_fasta = outputdir + sample_name + '_sorted_STEP7.fasta'
	flag, stderr = execute_functions(usearch_sort_by_size, processors, outputdir, 'multi', 'usearch', usearch_exec_path, preclustered_fasta, sorted_fasta)
	if flag is False:
		print "Execution of usearch_precluster failed!!!"
		sys.exit(2)
	high_abundant_fasta = outputdir + sample_name + '_high_abundant_STEP7.fasta'
	check_it_and_remove_it(high_abundant_fasta)
	flag, stderr = execute_functions(usearch_sort_abundance, processors, outputdir, 'multi', 'usearch', usearch_exec_path, sorted_fasta, high_abundant_fasta, min_abund_size)
	if flag is False:
		print "Execution of usearch_sort_abundance failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	hq_read_count = fasta_read_count(high_abundant_fasta)
	print "The number of reads after removing low abundance reads is: ", hq_read_count
	report("The number of reads after removing low abundance reads is: " + hq_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = check_and_verify(preclustered_read_count, hq_read_count, 3)
	if flag is False:
		error("Too little high abundance reads remains for clustering processes, removing low abundance reads disabled")
		print "Too little high abundance reads remains for clustering processes, removing low abundance reads disabled"
		high_abundant_fasta = preclustered_fasta
	else:
		print "LOW ABUNDANCE REMOVAL: PASSED!!!"
		report("LOW ABUNDANCE REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CLUSTERING USING UPARSE
	# ########################################################
	print "CLUSTERING USING UPARSE: "
	report("CLUSTERING USING UPARSE: ")
	centroids_fasta = outputdir + sample_name + '_uparse_centroids_STEP7.fasta'
	check_it_and_remove_it(centroids_fasta)
	uparse_file = outputdir + sample_name + '_uparse_STEP7.txt'
	check_it_and_remove_it(uparse_file)
	flag, stderr = execute_functions(usearch_uparse_cluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, high_abundant_fasta, centroids_fasta, uparse_file)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(centroids_fasta)
	print "The number of centroids reads is: ", read_count
	report("The number of centroids reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "CLUSTERING USING UPARSE: PASSED!!!"
	report("CLUSTERING USING UPARSE: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CHIMERA REMOVAL
	# ########################################################
	print "CHIMERA REMOVAL: "
	report("CHIMERA REMOVAL: ")
	
	dechim_fasta = outputdir + sample_name + '_chimeric_free_centroids_STEP7.fasta'
	check_it_and_remove_it(dechim_fasta)
	
	chim_fasta = outputdir + sample_name + '_chimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(chim_fasta)
	
	chim_report = outputdir + sample_name + '_chimera_report_STEP7.txt'
	check_it_and_remove_it(chim_report)
	
	flag, stderr = execute_functions(usearch_chimera_removal, processors, outputdir, 'multi', 'usearch', usearch_exec_path, centroids_fasta, chimera_reference, chim_fasta, dechim_fasta, chim_report)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(dechim_fasta)
	print "The number of chimeric_free centroids reads is: ", read_count
	report("The number of chimeric_free centroids reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "CHIMERA REMOVAL: PASSED!!!"
	report("CHIMERA REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CLASSIFICATION OF CENTROIDS
	# ########################################################
	print "CLASSIFICATION OF CENTROIDS: "
	report("CLASSIFICATION OF CENTROIDS: ")
	centroids_classification_report = outputdir + sample_name + '_centroids_classification_report_STEP7.txt'
	check_it_and_remove_it(centroids_classification_report)
	flag = execute_functions(mothur_classify_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, dechim_fasta, reference, taxonomy)
	if flag is False:
		print "Execution of mothur_classify_seqs failed!!!"
		sys.exit(2)
	classify_seqs_container = []
	extension_list = ['.knn.taxonomy']
	flag = scandirs(outputdir, classify_seqs_container, extension_list, mode='ex_partial')
	if flag is False:
		print "This extension is not availble: ", list_to_string(extension_list, ',')
		sys.exit(2)
	else:
		os.rename(classify_seqs_container[0], centroids_classification_report)
	centroids_label = outputdir + sample_name + '_centroids_classification_label_STEP7.txt'
	check_it_and_remove_it(centroids_label)
	
	relabeled_nonchimera_centroids_fasta = outputdir + sample_name + '_relabeled_nonchimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(relabeled_nonchimera_centroids_fasta)
	
	flag = classify_centroids(dechim_fasta, relabeled_nonchimera_centroids_fasta, centroids_classification_report, taxonomy, centroids_label)
	if flag is False:
		print "Execution of classify_centroids failed!!!"
		sys.exit(2)
	print "CLASSIFICATION OF CENTROIDS:  PASSED!!!"
	report("CLASSIFICATION OF CENTROIDS:  PASSED!!!")
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# ABUNDANCE CALCULATION AND TABLE FORMATION
	# ########################################################
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: "
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: ")
	
	alignment_report = outputdir + sample_name + '_alignment_report_STEP7.txt'
	check_it_and_remove_it(alignment_report)
	
	qiimeout = outputdir + sample_name + '_qiime_report_STEP7.txt'
	check_it_and_remove_it(qiimeout)
	abundance_count = outputdir + sample_name + '_abundance_count_report_STEP7.txt'
	check_it_and_remove_it(abundance_count)

	mothurout = outputdir + sample_name + '_mothur_report_STEP7.txt'
	check_it_and_remove_it(mothurout)

	mothurout_temp = outputdir + sample_name + '_mothur_report_temp.txt'
	check_it_and_remove_it(mothurout_temp)
	
	biomout_temp = outputdir + sample_name + '_biom_report_temp.txt'
	check_it_and_remove_it(biomout_temp)

	biomout = outputdir + sample_name + '_biom_report_STEP7.txt'
	check_it_and_remove_it(biomout)

	unaligned_reads = outputdir + sample_name + '_unaligned_reads_STEP7.fasta'
	check_it_and_remove_it(unaligned_reads)
	aligned_reads = outputdir + sample_name + '_aligned_reads_STEP7.fasta'
	check_it_and_remove_it(aligned_reads)
	
	identity = str(similarity)
	refcount = fasta_read_count(dechim_fasta)
	tsegout = False
	aligned_reads = False
	unaligned_reads = False
	flag, stderr = execute_functions(usearch_usearch_global, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_align_fasta_file, relabeled_nonchimera_centroids_fasta, alignment_report, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout_temp, biomout_temp, identity, maxhits, refcount)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	
	flag = mothur_shared_correct(centroids_label, mothurout_temp, mothurout)
	if flag is False:
		print "Execution of mothur correct failed!!!"
		sys.exit(2)
	print "MOTHUR abundance table generated."
	report("MOTHUR abundance table generated.")
	flag = biom_report_correct(centroids_label, biomout_temp, biomout)
	if flag is False:
		print "Execution of biom_report_correct failed!!!"
		sys.exit(2)
	print "BIOM abundance table generated."
	report("BIOM abundance table generated.")

	print "QIIME abundance table generated."
	report("QIIME abundance table generated.")

	phylotype_file = outputdir + sample_name + '_mothur_phylotype_STEP7.txt'
	check_it_and_remove_it(phylotype_file)
	flag = phylotype(alignment_report, centroids_label, phylotype_file, mothur_exec_path, sample_name, processors, outputdir)
	if flag is False:
		print "Execution of phylotype failed!!!"
		sys.exit(2)
	print "Phylotype abundance table generated."
	report("Phylotype abundance table generated.")

	sasout = outputdir + sample_name + '_sas_report_STEP7.txt'
	check_it_and_remove_it(sasout)
	create_sas_table(alignment_report, centroids_label, sasout, sample_name, file_type, total_error_rate, processors, outputdir)
	print "SAS abundance table generated."
	report("SAS abundance table generated.")

	print "ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!"
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	return (flag, qiimeout, mothurout, biomout, sasout, phylotype_file)


def upgma_pipeline(ready_to_align_fasta_file, ready_to_cluster_fasta_file, reference, taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, sample_name, chimera_reference, file_type, total_error_rate, similarity, maxhits, precluster, median_length, min_abund_size, processors, outputdir):
	
	# ########################################################
	# REMOVING DUPLICATES
	# ########################################################
	print "DUPLICATES REMOVAL:"
	report("DUPLICATES REMOVAL:")
	derep_fasta = outputdir + sample_name + '_dereplicated_STEP7.fasta'
	check_it_and_remove_it(derep_fasta)
	flag, stderr = execute_functions(usearch_dereplicate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_cluster_fasta_file, derep_fasta)
	if flag is False:
		print "Execution of usearch_dereplicate failed!!!"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	derep_read_count = fasta_read_count(derep_fasta)
	print "The number of reads after removing duplicates reads is: ", derep_read_count
	report("The number of reads after removing duplicates reads is: " + derep_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "DUPLICATES REMOVAL: PASSED!!"
	report("DUPLICATES REMOVAL: PASSED!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# PRECLUSTER
	# ########################################################
	if precluster is True:
		print "PRECLUSTER: "
		report("PRECLUSTER: ")
		maxdiffs = str(int(float(median_length) / 100))
		preclustered_fasta = outputdir + sample_name + '_preclustered_STEP7.fasta'
		check_it_and_remove_it(preclustered_fasta)
		#maxdiffs = '3'
		flag, stderr = execute_functions(usearch_precluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, derep_fasta, preclustered_fasta, maxdiffs)
		if flag is False:
			print "Execution of usearch_precluster failed!!!"
			sys.exit(2)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		preclustered_read_count = fasta_read_count(preclustered_fasta)
		print "The number of reads after precluster is: ", preclustered_read_count
		report("The number of reads after precluster reads is: " + preclustered_read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "PRECLUSTER: PASSED!!!"
		report("PRECLUSTER: PASSED!!!")
	else:
		preclustered_fasta = derep_fasta
		preclustered_read_count = derep_read_count
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# LOW ABUNDANCE REMOVAL
	# ########################################################
	print "LOW ABUNDANCE REMOVAL:"
	report("LOW ABUNDANCE REMOVAL:")
	sorted_fasta = outputdir + sample_name + '_sorted_STEP7.fasta'
	flag, stderr = execute_functions(usearch_sort_by_size, processors, outputdir, 'multi', 'usearch', usearch_exec_path, preclustered_fasta, sorted_fasta)
	if flag is False:
		print "Execution of usearch_precluster failed!!!"
		sys.exit(2)
	high_abundant_fasta = outputdir + sample_name + '_high_abundant_STEP7.fasta'
	check_it_and_remove_it(high_abundant_fasta)
	flag, stderr = execute_functions(usearch_sort_abundance, processors, outputdir, 'multi', 'usearch', usearch_exec_path, sorted_fasta, high_abundant_fasta, min_abund_size)
	if flag is False:
		print "Execution of usearch_sort_abundance failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	hq_read_count = fasta_read_count(high_abundant_fasta)
	print "The number of reads after removing low abundance reads is: ", hq_read_count
	report("The number of reads after removing low abundance reads is: " + hq_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = check_and_verify(preclustered_read_count, hq_read_count, 3)
	if flag is False:
		error("Too little high abundance reads remains for clustering processes, removing low abundance reads disabled")
		print "Too little high abundance reads remains for clustering processes, removing low abundance reads disabled"
		high_abundant_fasta = preclustered_fasta
		hq_read_count = preclustered_read_count
	else:
		print "LOW ABUNDANCE REMOVAL: PASSED!!!"
		report("LOW ABUNDANCE REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# MULTIPLE ALIGNMENT USING MAFFT
	# ########################################################
	print "MULTIPLE ALIGNMENT USING MAFFT:"
	report("MULTIPLE ALIGNMENT USING MAFFT:")
	space = ' '
	mafft_parameter = ''
	if int(hq_read_count) > 10000:
		mafft_parameter = '--parttree' + space
	else:
		mafft_parameter = '--auto' + space
	multiple_aligned_fasta = outputdir + sample_name + '_mafft_multiple_aligned_STEP7.fasta'

	flag, stderr = execute_functions(mafft_align, processors, outputdir, 'multi', 'mafft', mafft_exec_path, high_abundant_fasta, multiple_aligned_fasta, mafft_parameter)
	if flag is False:
		print "Execution of usearch_dereplicate failed!!!"
		sys.exit(2)
	report("MULTIPLE ALIGNMENT USING MAFFT: PASSED!!!")
	print "MULTIPLE ALIGNMENT USING MAFFT: PASSED!!!"
	report("####################################################################")
	print "####################################################################"
	
	# ########################################################
	# DISTANCE MATRIX CALCULATION
	# ########################################################
	report("DISTANCE MATRIX CALCULATION: ")
	print "DISTANCE MATRIX CALCULATION: "
	mothur_distance_matrix = outputdir + sample_name + '_distance_matrix_STEP7.dist'
	check_it_and_remove_it(mothur_distance_matrix)
	cutoff_limit = '0.03'
	flag, stderr = execute_functions(mothur_dist_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, multiple_aligned_fasta, cutoff_limit)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
		sys.exit(2)
	else:
		cluster_distance_container = []
		extension_list = ['.dist']
		flag = scandirs(outputdir, cluster_distance_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			os.rename(cluster_distance_container[0], mothur_distance_matrix)

	report("DISTANCE MATRIX CALCULATION: PASSED!!!")
	print "DISTANCE MATRIX CALCULATION: PASSED!!!"
	report("####################################################################")
	print "####################################################################"
	
	# ########################################################
	# CLUSTERING USING UPGMA
	# ########################################################
	print "CLUSTERING USING UPGMA: "
	report("CLUSTERING USING UPGMA: ")
	
	count_table = outputdir + sample_name + '_multiple_aligned_count_STEP7.txt'
	group_table = outputdir + sample_name + '_multiple_aligned_group_STEP7.txt'
	flag = mothur_count_group_table(high_abundant_fasta, count_table, group_table)
	if flag is not True:
		error("can not create group and count_table")
		print "Can not create group and count_table"
		sys.exit(2)
	
	cluster_list_container = []
	mothur_cluster_list = outputdir + sample_name + '_cluster_list_STEP7.list'
	check_it_and_remove_it(mothur_cluster_list)

	flag, stderr = execute_functions(mothur_cluster_split, processors, outputdir, 'multi', 'mothur', mothur_exec_path, mothur_distance_matrix, count_table)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
	else:
		cluster_list_container = []
		extension_list = ['.an.unique_list.list']
		flag = scandirs(outputdir, cluster_list_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(cluster_list_container[0], mothur_cluster_list)

	flag, stderr = execute_functions(mothur_get_oturep, processors, outputdir, 'multi', 'mothur', mothur_exec_path, mothur_distance_matrix, mothur_cluster_list, high_abundant_fasta, count_table)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
	else:
		otu_rep_container = []
		extension_list = ['.0.03.rep.fasta', '.0.04.rep.fasta', '.0.02.rep.fasta', '.0.01.rep.fasta', '.unique.rep.fasta']
		for rate in extension_list:
			rate_list = [rate]
			flag = scandirs(outputdir, otu_rep_container, rate_list)
			if flag is False:
				print "This extension is not availble: ", rate_list
			else:
				break
		if len(otu_rep_container) < 1:
			print "Centroids file retreival failed!!!!"
			sys.exit(2)

	centroids_fasta = outputdir + sample_name + '_UPGMA_centroids_STEP7.fasta'
	check_it_and_remove_it(centroids_fasta)
	tag = 'OTU'
	suffix = ''
	flag = relable_fasta(otu_rep_container[0], tag, suffix, centroids_fasta)
	if flag is False:
		print "Execution of relabel fasta failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(centroids_fasta)
	print "The number of centroids reads is: ", read_count
	report("The number of centroids reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("CLUSTERING USING UPGMA:: PASSED!!!")
	print "CLUSTERING USING UPGMA:: PASSED!!!"
	report("####################################################################")
	print "####################################################################"

	# ########################################################
	# CHIMERA REMOVAL
	# ########################################################
	print "CHIMERA REMOVAL: "
	report("CHIMERA REMOVAL: ")
	
	dechim_fasta = outputdir + sample_name + '_chimeric_free_centroids_STEP7.fasta'
	check_it_and_remove_it(dechim_fasta)
	
	chim_fasta = outputdir + sample_name + '_chimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(chim_fasta)
	
	chim_report = outputdir + sample_name + '_chimera_report_STEP7.txt'
	check_it_and_remove_it(chim_report)
	
	flag, stderr = execute_functions(usearch_chimera_removal, processors, outputdir, 'multi', 'usearch', usearch_exec_path, centroids_fasta, chimera_reference, chim_fasta, dechim_fasta, chim_report)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(dechim_fasta)
	print "The number of chimeric_free centroids reads is: ", read_count
	report("The number of chimeric_free centroids reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "CHIMERA REMOVAL: PASSED!!!"
	report("CHIMERA REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CLASSIFICATION OF CENTROIDS
	# ########################################################
	print "CLASSIFICATION OF CENTROIDS: "
	report("CLASSIFICATION OF CENTROIDS: ")
	centroids_classification_report = outputdir + sample_name + '_centroids_classification_report_STEP7.txt'
	check_it_and_remove_it(centroids_classification_report)
	flag = execute_functions(mothur_classify_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, dechim_fasta, reference, taxonomy)
	if flag is False:
		print "Execution of mothur_classify_seqs failed!!!"
		sys.exit(2)
	classify_seqs_container = []
	extension_list = ['.knn.taxonomy']
	flag = scandirs(outputdir, classify_seqs_container, extension_list, mode='ex_partial')
	if flag is False:
		print "This extension is not availble: ", list_to_string(extension_list, ',')
		sys.exit(2)
	else:
		os.rename(classify_seqs_container[0], centroids_classification_report)
	centroids_label = outputdir + sample_name + '_centroids_classification_label_STEP7.txt'
	check_it_and_remove_it(centroids_label)
	relabeled_nonchimera_centroids_fasta = outputdir + sample_name + '_relabeled_nonchimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(relabeled_nonchimera_centroids_fasta)
	flag = classify_centroids(dechim_fasta, relabeled_nonchimera_centroids_fasta, centroids_classification_report, taxonomy, centroids_label)
	if flag is False:
		print "Execution of classify_centroids failed!!!"
		sys.exit(2)
	print "CLASSIFICATION OF CENTROIDS:  PASSED!!!"
	report("CLASSIFICATION OF CENTROIDS:  PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# ABUNDANCE CALCULATION AND TABLE FORMATION
	# ########################################################
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: "
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: ")
	
	alignment_report = outputdir + sample_name + '_alignment_report_STEP7.txt'
	check_it_and_remove_it(alignment_report)
	
	qiimeout = outputdir + sample_name + '_qiime_report_STEP7.txt'
	check_it_and_remove_it(qiimeout)
	
	mothurout = outputdir + sample_name + '_mothur_report_STEP7.txt'
	check_it_and_remove_it(mothurout)

	mothurout_temp = outputdir + sample_name + '_mothur_report_temp.txt'
	check_it_and_remove_it(mothurout_temp)
	
	biomout_temp = outputdir + sample_name + '_biom_report_temp.txt'
	check_it_and_remove_it(biomout_temp)

	biomout = outputdir + sample_name + '_biom_report_STEP7.txt'
	check_it_and_remove_it(biomout)

	unaligned_reads = outputdir + sample_name + '_unaligned_reads_STEP7.fasta'
	check_it_and_remove_it(unaligned_reads)

	aligned_reads = outputdir + sample_name + '_aligned_reads_STEP7.fasta'
	check_it_and_remove_it(aligned_reads)
	
	identity = str(similarity)
	refcount = fasta_read_count(dechim_fasta)
	tsegout = False
	aligned_reads = False
	unaligned_reads = False
	flag, stderr = execute_functions(usearch_usearch_global, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_align_fasta_file, relabeled_nonchimera_centroids_fasta, alignment_report, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout_temp, biomout_temp, identity, maxhits, refcount)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	
	flag = mothur_shared_correct(centroids_label, mothurout_temp, mothurout)
	if flag is False:
		print "Execution of mothur correct failed!!!"
		sys.exit(2)
	print "MOTHUR abundance table generated."
	report("MOTHUR abundance table generated.")
	flag = biom_report_correct(centroids_label, biomout_temp, biomout)
	if flag is False:
		print "Execution of biom_report_correct failed!!!"
		sys.exit(2)
	print "BIOM abundance table generated."
	report("BIOM abundance table generated.")

	print "QIIME abundance table generated."
	report("QIIME abundance table generated.")

	phylotype_file = outputdir + sample_name + '_mothur_phylotype_STEP7.txt'
	check_it_and_remove_it(phylotype_file)
	
	flag = phylotype(alignment_report, centroids_label, phylotype_file, mothur_exec_path, sample_name, processors, outputdir)
	print "Phylotype abundance table generated."
	report("Phylotype abundance table generated.")
	
	sasout = outputdir + sample_name + '_sas_report_STEP7.txt'
	check_it_and_remove_it(sasout)
	create_sas_table(alignment_report, centroids_label, sasout, sample_name, file_type, total_error_rate, processors, outputdir)
	print "SAS abundance table generated."
	report("SAS abundance table generated.")
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!"
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	return (flag, qiimeout, mothurout, biomout, sasout, phylotype_file)


def uclust_pipeline(ready_to_align_fasta_file, ready_to_cluster_fasta_file, reference, taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, sample_name, chimera_reference, file_type, total_error_rate, similarity, maxhits, precluster, median_length, min_abund_size, processors, outputdir):
	
	# ########################################################
	# REMOVING DUPLICATES
	# ########################################################
	print "DUPLICATES REMOVAL:"
	report("DUPLICATES REMOVAL:")
	derep_fasta = outputdir + sample_name + '_dereplicated_STEP7.fasta'
	check_it_and_remove_it(derep_fasta)
	flag, stderr = execute_functions(usearch_dereplicate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_cluster_fasta_file, derep_fasta)
	if flag is False:
		print "Execution of usearch_dereplicate failed!!!"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	derep_read_count = fasta_read_count(derep_fasta)
	print "The number of reads after removing duplicates reads is: ", derep_read_count
	report("The number of reads after removing duplicates reads is: " + derep_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "DUPLICATES REMOVAL: PASSED!!"
	report("DUPLICATES REMOVAL: PASSED!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# PRECLUSTER
	# ########################################################
	if precluster is True:
		print "PRECLUSTER: "
		report("PRECLUSTER: ")
		maxdiffs = str(int(float(median_length) / 100))
		preclustered_fasta = outputdir + sample_name + '_preclustered_STEP7.fasta'
		check_it_and_remove_it(preclustered_fasta)
		#maxdiffs = '3'
		flag, stderr = execute_functions(usearch_precluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, derep_fasta, preclustered_fasta, maxdiffs)
		if flag is False:
			print "Execution of usearch_precluster failed!!!"
			sys.exit(2)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		preclustered_read_count = fasta_read_count(preclustered_fasta)
		print "The number of reads after precluster is: ", preclustered_read_count
		report("The number of reads after precluster reads is: " + preclustered_read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "PRECLUSTER: PASSED!!!"
		report("PRECLUSTER: PASSED!!!")
	else:
		preclustered_fasta = derep_fasta
		preclustered_read_count = derep_read_count
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# LOW ABUNDANCE REMOVAL
	# ########################################################
	print "LOW ABUNDANCE REMOVAL:"
	report("LOW ABUNDANCE REMOVAL:")
	sorted_fasta = outputdir + sample_name + '_sorted_STEP7.fasta'
	flag, stderr = execute_functions(usearch_sort_by_size, processors, outputdir, 'multi', 'usearch', usearch_exec_path, preclustered_fasta, sorted_fasta)
	if flag is False:
		print "Execution of usearch_precluster failed!!!"
		sys.exit(2)
	high_abundant_fasta = outputdir + sample_name + '_high_abundant_STEP7.fasta'
	check_it_and_remove_it(high_abundant_fasta)
	flag, stderr = execute_functions(usearch_sort_abundance, processors, outputdir, 'multi', 'usearch', usearch_exec_path, sorted_fasta, high_abundant_fasta, min_abund_size)
	if flag is False:
		print "Execution of usearch_sort_abundance failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	hq_read_count = fasta_read_count(high_abundant_fasta)
	print "The number of reads after removing low abundance reads is: ", hq_read_count
	report("The number of reads after removing low abundance reads is: " + hq_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = check_and_verify(preclustered_read_count, hq_read_count, 3)
	if flag is False:
		error("Too little high abundance reads remains for clustering processes, removing low abundance reads disabled")
		print "Too little high abundance reads remains for clustering processes, removing low abundance reads disabled"
		high_abundant_fasta = preclustered_fasta
	else:
		print "LOW ABUNDANCE REMOVAL: PASSED!!!"
		report("LOW ABUNDANCE REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# CLUSTERING USING UCLUST
	# ########################################################
	print "CLUSTERING USING UCLUST: "
	report("CLUSTERING USING UCLUST: ")
	uclust_centroids_fasta = outputdir + sample_name + '_uclust_centroids_temp_STEP7.fasta'
	check_it_and_remove_it(uclust_centroids_fasta)
	uclust_file = outputdir + sample_name + '_uclust_STEP7.txt'
	check_it_and_remove_it(uclust_file)

	flag, stderr = execute_functions(usearch_uclust_cluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, high_abundant_fasta, uclust_centroids_fasta, uclust_file)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(uclust_centroids_fasta)
	print "The number of reads after clustering using uclust reads is: ", read_count
	report("The number of reads after clustering using uclust reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	
	print "CLUSTERING USING UCLUST: PASSED!!!"
	report("CLUSTERING USING UCLUST: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CHIMERA REMOVAL
	# ########################################################
	print "CHIMERA REMOVAL: "
	report("CHIMERA REMOVAL: ")
	
	dechim_fasta = outputdir + sample_name + '_chimeric_free_centroids_STEP7.fasta'
	check_it_and_remove_it(dechim_fasta)
	
	chim_fasta = outputdir + sample_name + '_chimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(chim_fasta)
	
	chim_report = outputdir + sample_name + '_chimera_report_STEP7.txt'
	check_it_and_remove_it(chim_report)
	
	flag, stderr = execute_functions(usearch_chimera_removal, processors, outputdir, 'multi', 'usearch', usearch_exec_path, uclust_centroids_fasta, chimera_reference, chim_fasta, dechim_fasta, chim_report)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	read_count = fasta_read_count(dechim_fasta)
	print "The number of chimeric_free centroids reads is: ", read_count
	report("The number of chimeric_free centroids reads is: " + read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "CHIMERA REMOVAL: PASSED!!!"
	report("CHIMERA REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# CLASSIFICATION OF CENTROIDS
	# ########################################################
	print "CLASSIFICATION OF CENTROIDS: "
	report("CLASSIFICATION OF CENTROIDS: ")
	centroids_classification_report = outputdir + sample_name + '_centroids_classification_report_STEP7.txt'
	check_it_and_remove_it(centroids_classification_report)
	flag = execute_functions(mothur_classify_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, dechim_fasta, reference, taxonomy)
	if flag is False:
		print "Execution of mothur_classify_seqs failed!!!"
		sys.exit(2)
	classify_seqs_container = []
	extension_list = ['.knn.taxonomy']
	flag = scandirs(outputdir, classify_seqs_container, extension_list, mode='ex_partial')
	if flag is False:
		print "This extension is not availble: ", list_to_string(extension_list, ',')
		sys.exit(2)
	else:
		os.rename(classify_seqs_container[0], centroids_classification_report)
	centroids_label = outputdir + sample_name + '_centroids_classification_label_STEP7.txt'
	check_it_and_remove_it(centroids_label)
	relabeled_nonchimera_centroids_fasta = outputdir + sample_name + '_relabeled_nonchimeric_centroids_STEP7.fasta'
	check_it_and_remove_it(relabeled_nonchimera_centroids_fasta)
	flag = classify_centroids(dechim_fasta, relabeled_nonchimera_centroids_fasta, centroids_classification_report, taxonomy, centroids_label)
	if flag is False:
		print "Execution of classify_centroids failed!!!"
		sys.exit(2)
	print "CLASSIFICATION OF CENTROIDS:  PASSED!!!"
	report("CLASSIFICATION OF CENTROIDS:  PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# ABUNDANCE CALCULATION AND TABLE FORMATION
	# ########################################################
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: "
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: ")
	
	alignment_report = outputdir + sample_name + '_alignment_report_STEP7.txt'
	check_it_and_remove_it(alignment_report)
	
	qiimeout = outputdir + sample_name + '_qiime_report_STEP7.txt'
	check_it_and_remove_it(qiimeout)
	
	mothurout = outputdir + sample_name + '_mothur_report_STEP7.txt'
	check_it_and_remove_it(mothurout)

	mothurout_temp = outputdir + sample_name + '_mothur_report_temp.txt'
	check_it_and_remove_it(mothurout_temp)
	
	biomout_temp = outputdir + sample_name + '_biom_report_temp.txt'
	check_it_and_remove_it(biomout_temp)

	biomout = outputdir + sample_name + '_biom_report_STEP7.txt'
	check_it_and_remove_it(biomout)

	unaligned_reads = outputdir + sample_name + '_unaligned_reads_STEP7.fasta'
	check_it_and_remove_it(unaligned_reads)

	aligned_reads = outputdir + sample_name + '_aligned_reads_STEP7.fasta'
	check_it_and_remove_it(aligned_reads)
	
	identity = str(similarity)
	refcount = fasta_read_count(dechim_fasta)
	tsegout = False
	aligned_reads = False
	unaligned_reads = False
	flag, stderr = execute_functions(usearch_usearch_global, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_align_fasta_file, relabeled_nonchimera_centroids_fasta, alignment_report, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout_temp, biomout_temp, identity, maxhits, refcount)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	
	flag = mothur_shared_correct(centroids_label, mothurout_temp, mothurout)
	if flag is False:
		print "Execution of mothur correct failed!!!"
		sys.exit(2)
	print "MOTHUR abundance table generated."
	report("MOTHUR abundance table generated.")
	flag = biom_report_correct(centroids_label, biomout_temp, biomout)
	if flag is False:
		print "Execution of biom_report_correct failed!!!"
		sys.exit(2)
	print "BIOM abundance table generated."
	report("BIOM abundance table generated.")
	print "QIIME abundance table generated."
	report("QIIME abundance table generated.")

	phylotype_file = outputdir + sample_name + '_mothur_phylotype_STEP7.txt'
	check_it_and_remove_it(phylotype_file)
	
	flag = phylotype(alignment_report, centroids_label, phylotype_file, mothur_exec_path, sample_name, processors, outputdir)
	print "Phylotype abundance table generated."
	report("Phylotype abundance table generated.")
	
	sasout = outputdir + sample_name + '_sas_report_STEP7.txt'
	check_it_and_remove_it(sasout)
	create_sas_table(alignment_report, centroids_label, sasout, sample_name, file_type, total_error_rate, processors, outputdir)
	print "SAS abundance table generated."
	report("SAS abundance table generated.")
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!"
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	return (flag, qiimeout, mothurout, biomout, sasout, phylotype_file)


def supervised_pipeline(ready_to_align_fasta_file, ready_to_cluster_fasta_file, reference, taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, sample_name, chimera_reference, file_type, total_error_rate, similarity, maxhits, precluster, median_length, min_abund_size, processors, outputdir):

	# ########################################################
	# ABUNDANCE CALCULATION AND TABLE FORMATION
	# ########################################################
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: "
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: ")
	
	alignment_report = outputdir + sample_name + '_alignment_report_STEP7.txt'
	check_it_and_remove_it(alignment_report)
	
	qiimeout = outputdir + sample_name + '_qiime_report_STEP7.txt'
	check_it_and_remove_it(qiimeout)
	
	mothurout = outputdir + sample_name + '_mothur_report_STEP7.txt'
	check_it_and_remove_it(mothurout)

	mothurout_temp = outputdir + sample_name + '_mothur_report_temp.txt'
	check_it_and_remove_it(mothurout_temp)
	
	biomout_temp = outputdir + sample_name + '_biom_report_temp.txt'
	check_it_and_remove_it(biomout_temp)

	biomout = outputdir + sample_name + '_biom_report_STEP7.txt'
	check_it_and_remove_it(biomout)

	unaligned_reads = outputdir + sample_name + '_unaligned_reads_STEP7.fasta'
	check_it_and_remove_it(unaligned_reads)

	aligned_reads = outputdir + sample_name + '_aligned_reads_STEP7.fasta'
	check_it_and_remove_it(aligned_reads)
	
	identity = str(similarity)
	tsegout = False
	unaligned_reads = False
	aligned_reads = False
	flag, stderr = execute_functions(usearch_usearch_global_supervised, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_align_fasta_file, reference, alignment_report, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout_temp, biomout_temp, identity, maxhits)
	if flag is False:
		print "Execution of usearch_uparse_cluster failed!!!"
		sys.exit(2)
	
	flag = mothur_shared_correct_supervised(taxonomy, mothurout_temp, mothurout)
	if flag is False:
		print "Execution of mothur correct failed!!!"
		sys.exit(2)
	print "MOTHUR abundance table generated."
	report("MOTHUR abundance table generated.")
	flag = biom_report_correct_supervised(taxonomy, biomout_temp, biomout)
	if flag is False:
		print "Execution of biom_report_correct failed!!!"
		sys.exit(2)
	print "BIOM abundance table generated."
	report("BIOM abundance table generated.")
	print "QIIME abundance table generated."
	report("QIIME abundance table generated.")
	phylotype_file = outputdir + sample_name + '_mothur_phylotype_STEP7.txt'
	check_it_and_remove_it(phylotype_file)
	
	flag = phylotype_supervised(alignment_report, taxonomy, phylotype_file, mothur_exec_path, sample_name, processors, outputdir)
	print "Phylotype abundance table generated."
	report("Phylotype abundance table generated.")
	sasout = outputdir + sample_name + '_sas_report_STEP7.txt'
	check_it_and_remove_it(sasout)
	create_sas_table_supervised(alignment_report, taxonomy, sasout, sample_name, file_type, total_error_rate, processors, outputdir)
	print "SAS abundance table generated."
	report("SAS abundance table generated.")
	print "ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!"
	report("ABUNDANCE CALCULATION AND TABLE FORMATION: PASSED!!!")
	report("###############################################################")
	print "###############################################################"
	return (flag, qiimeout, mothurout, biomout, sasout, phylotype_file)


def unifrac_pipeline(ready_to_align_fasta_file, ready_to_cluster_fasta_file, reference, taxonomy, usearch_exec_path, mafft_exec_path, mothur_exec_path, fasttree_exec_path, sample_name, chimera_reference, file_type, total_error_rate, similarity, maxhits, precluster, median_length, min_abund_size, processors, outputdir):
	
	# ########################################################
	# REMOVING DUPLICATES
	# ########################################################
	print "DUPLICATES REMOVAL:"
	report("DUPLICATES REMOVAL:")
	derep_fasta = outputdir + sample_name + '_dereplicated_STEP7.fasta'
	check_it_and_remove_it(derep_fasta)
	flag, stderr = execute_functions(usearch_dereplicate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, ready_to_cluster_fasta_file, derep_fasta)
	if flag is False:
		print "Execution of usearch_dereplicate failed!!!"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	derep_read_count = fasta_read_count(derep_fasta)
	print "The number of reads after removing duplicates reads is: ", derep_read_count
	report("The number of reads after removing duplicates reads is: " + derep_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "DUPLICATES REMOVAL: PASSED!!"
	report("DUPLICATES REMOVAL: PASSED!!")
	report("###############################################################")
	print "###############################################################"
	
	# ########################################################
	# PRECLUSTER
	# ########################################################
	if precluster is True:
		print "PRECLUSTER: "
		report("PRECLUSTER: ")
		maxdiffs = str(int(float(median_length) / 100))
		preclustered_fasta = outputdir + sample_name + '_preclustered_STEP7.fasta'
		check_it_and_remove_it(preclustered_fasta)
		#maxdiffs = '3'
		flag, stderr = execute_functions(usearch_precluster, processors, outputdir, 'multi', 'usearch', usearch_exec_path, derep_fasta, preclustered_fasta, maxdiffs)
		if flag is False:
			print "Execution of usearch_precluster failed!!!"
			sys.exit(2)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		preclustered_read_count = fasta_read_count(preclustered_fasta)
		print "The number of reads after precluster is: ", preclustered_read_count
		report("The number of reads after precluster reads is: " + preclustered_read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "PRECLUSTER: PASSED!!!"
		report("PRECLUSTER: PASSED!!!")
	else:
		preclustered_fasta = derep_fasta
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# LOW ABUNDANCE REMOVAL
	# ########################################################
	print "LOW ABUNDANCE REMOVAL:"
	report("LOW ABUNDANCE REMOVAL:")
	sorted_fasta = outputdir + sample_name + '_sorted_STEP7.fasta'
	flag, stderr = execute_functions(usearch_sort_by_size, processors, outputdir, 'multi', 'usearch', usearch_exec_path, preclustered_fasta, sorted_fasta)
	if flag is False:
		print "Execution of usearch_precluster failed!!!"
		sys.exit(2)
	high_abundant_fasta = outputdir + sample_name + '_high_abundant_STEP7.fasta'
	check_it_and_remove_it(high_abundant_fasta)
	flag, stderr = execute_functions(usearch_sort_abundance, processors, outputdir, 'multi', 'usearch', usearch_exec_path, sorted_fasta, high_abundant_fasta, min_abund_size)
	if flag is False:
		print "Execution of usearch_sort_abundance failed!!!"
		sys.exit(2)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	hq_read_count = fasta_read_count(high_abundant_fasta)
	print "The number of reads after removing low abundance reads is: ", hq_read_count
	report("The number of reads after removing low abundance reads is: " + hq_read_count)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = check_and_verify(preclustered_read_count, hq_read_count, 3)
	if flag is False:
		error("Too little high abundance reads remains for clustering processes, removing low abundance reads disabled")
		print "Too little high abundance reads remains for clustering processes, removing low abundance reads disabled"
		high_abundant_fasta = preclustered_fasta
		hq_read_count = preclustered_read_count
	else:
		print "LOW ABUNDANCE REMOVAL: PASSED!!!"
		report("LOW ABUNDANCE REMOVAL: PASSED!!!")
	report("###############################################################")
	print "###############################################################"

	# ########################################################
	# MULTIPLE ALIGNMENT USING MAFFT
	# ########################################################
	print "MULTIPLE ALIGNMENT USING MAFFT:"
	report("MULTIPLE ALIGNMENT USING MAFFT:")
	space = ' '
	mafft_parameter = ''
	if int(hq_read_count) > 10000:
		mafft_parameter = '--parttree' + space
	else:
		mafft_parameter = '--auto' + space
	multiple_aligned_fasta = outputdir + sample_name + '_mafft_multiple_aligned_STEP7.fasta'

	flag, stderr = execute_functions(mafft_align, processors, outputdir, 'multi', 'mafft', mafft_exec_path, high_abundant_fasta, multiple_aligned_fasta, mafft_parameter)
	if flag is False:
		print "Execution of mafft_align failed!!!"
		sys.exit(2)
	report("MULTIPLE ALIGNMENT USING MAFFT: PASSED!!!")
	print "MULTIPLE ALIGNMENT USING MAFFT: PASSED!!!"
	report("####################################################################")
	print "####################################################################"
	
	# ########################################################
	# PHYLOGENETIC TREE CALCULATION
	# ########################################################
	report("PHYLOGENETIC TREE USING FASTTREE: ")
	print "PHYLOGENETIC TREE USING FASTTREE: "
	phylogenetic_tree = outputdir + sample_name + '_phylogenetic_tree_STEP7.tre'
	check_it_and_remove_it(phylogenetic_tree)
	flag, stderr = execute_functions(fasttree_build, processors, outputdir, 'multi', 'mothur', fasttree_exec_path, multiple_aligned_fasta, phylogenetic_tree)
	if flag is False:
		print "Execution of fasttree_build failed!!!"
		sys.exit(2)
	report("PHYLOGENETIC TREE USING FASTTREE: PASSED!!!")
	print "PHYLOGENETIC TREE USING FASTTREE: PASSED!!!"
	report("####################################################################")
	print "####################################################################"
	
	# ########################################################
	# UNIFRAC WIGHTED CALCULATION
	# ########################################################
	report("WEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: ")
	print "WEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: "
	count_table = outputdir + sample_name + '_multiple_aligned_count_STEP7.txt'
	group_table = outputdir + sample_name + '_multiple_aligned_group_STEP7.txt'
	flag = mothur_count_group_table(high_abundant_fasta, count_table, group_table)
	if flag is not True:
		error("can not create group and count_table")
		print "Can not create group and count_table"
		sys.exit(2)
	unifrac_weighted_distance_matrix = outputdir + sample_name + '_unifrac_weighted_distance_matrix_STEP7.dist'
	check_it_and_remove_it(unifrac_weighted_distance_matrix)
	flag, stderr = execute_functions(mothur_unifrac_weighted, processors, outputdir, 'multi', 'mothur', mothur_exec_path, phylogenetic_tree, group_table)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
		sys.exit(2)
	else:
		unifrac_weighted_distance_container = []
		extension_list = ['.tre1.weighted.phylip.dist']
		flag = scandirs(outputdir, unifrac_weighted_distance_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			os.rename(unifrac_weighted_distance_container[0], unifrac_weighted_distance_matrix)
	report("WEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: PASSED!!!")
	print "WEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: PASSED!!!"
	report("####################################################################")
	print "####################################################################"

	# ########################################################
	# UNIFRAC UNWIGHTED CALCULATION
	# ########################################################
	report("UNWEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: ")
	print "UNWEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: "
	unifrac_unweighted_distance_matrix = outputdir + sample_name + '_unifrac_unweighted_distance_matrix_STEP7.dist'
	check_it_and_remove_it(unifrac_unweighted_distance_matrix)
	flag, stderr = execute_functions(mothur_unifrac_unweighted, processors, outputdir, 'multi', 'mothur', mothur_exec_path, phylogenetic_tree, group_table)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
		sys.exit(2)
	else:
		unifrac_unweighted_distance_container = []
		extension_list = ['.tre1.unweighted.phylip.dist']
		flag = scandirs(outputdir, unifrac_unweighted_distance_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			os.rename(unifrac_unweighted_distance_container[0], unifrac_unweighted_distance_matrix)
	report("UNWEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: PASSED!!!")
	print "UNWEIGHTED UNIFRAC DISTANCE MATRIX CALCULATION: PASSED!!!"
	report("####################################################################")
	print "####################################################################"

	return (True, unifrac_unweighted_distance_matrix, unifrac_weighted_distance_matrix)


# ################################# MERGE AND RELABEL FUNCTIONS ############## #

def data_merge_relabel(filetype, filesystem, files_path, relabel, separator, name, ngs_qscore, usearch_exec_path, processors, outputdir):
	total_data = []
	# ###############################################################################
	if filetype == 'fastq' and filesystem == 'paired':
		total_fastq = outputdir + name + '_total_merged_relabeled_STEP2.fastq'
		total_fasta = outputdir + name + '_total_merged_relabeled_STEP2.fasta'
		print "Fastq Paired parser started"
		report("Fastq Paired parser started")
		flag = fastq_paired_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir)
		if flag is True:
			pass
		else:
			error("DATA MERGE AND RELABELING: FAILED!!!")
			print "DATA MERGE AND RELABELING: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		total_data.append(total_fastq)
		total_data.append(total_fasta)
		print "Fastq Paired parser: Done!!!"
		report("Fastq Paired parser: Done!!!")
	# ###############################################################################
	elif filetype == 'fastq' and filesystem == 'single':
		total_fastq = outputdir + name + '_total_merged_relabeled_STEP2.fastq'
		total_fasta = outputdir + name + '_total_merged_relabeled_STEP2.fasta'
		print "Fastq Single parser started"
		report("Fastq Single parser started")
		flag = fastq_single_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir)
		if flag is True:
			pass
		else:
			error("DATA MERGE AND RELABELING: FAILED!!!")
			print "DATA MERGE AND RELABELING: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		total_data.append(total_fastq)
		total_data.append(total_fasta)
		print "Fastq Single parser: Done!!!"
		report("Fastq Single parser: Done!!!")
	# ###############################################################################
	elif filetype == 'fastq' and filesystem == 'mixed':
		total_fastq = outputdir + name + '_total_merged_relabeled_STEP2.fastq'
		total_fasta = outputdir + name + '_total_merged_relabeled_STEP2.fasta'
		print "Fastq Mixed parser started"
		report("Fastq Mixed parser started")
		flag = fastq_mixed_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir)
		if flag is True:
			pass
		else:
			error("DATA MERGE AND RELABELING: FAILED!!!")
			print "DATA MERGE AND RELABELING: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		total_data.append(total_fastq)
		total_data.append(total_fasta)
		print "Fastq Mixed parser: Done!!!"
		report("Fastq Mixed parser: Done!!!")
	# ###############################################################################
	elif filetype == 'fasta' and filesystem == 'single':
		total_fasta = outputdir + name + '_total_merged_relabeled_STEP2.fasta'
		flag = fasta_single_parser(files_path, total_fasta, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir)
		if flag is True:
			pass
		else:
			error("DATA MERGE AND RELABELING: FAILED!!!")
			print "DATA MERGE AND RELABELING: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		total_data.append(total_fasta)
	# ###############################################################################
	elif filetype == 'fasta' and filesystem == 'mixed':
		total_fasta = outputdir + name + '_total_merged_relabeled_STEP2.fasta'
		flag = fasta_mixed_parser(files_path, total_fasta, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir)
		if flag is True:
			pass
		else:
			error("DATA MERGE AND RELABELING: FAILED!!!")
			print "DATA MERGE AND RELABELING: FAILED!!!"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		total_data.append(total_fasta)
	
	# ###############################################################################
	return (True, total_data)


def fastq_paired_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	# we merge paires of fastq file, rename them, concat them into one large fastq file
	fastq_total_list = []
	fasta_total_list = []
	
	file_paths_handle = open(files_path, 'rU')
	for path_line in file_paths_handle:
		path_line = path_line.rstrip()
		forward_fastq = path_line.split('\t')[1]
		reverse_fastq = path_line.split('\t')[2]
		path, absname, ext = split_file_name(forward_fastq)
		if separator is False:
			sample_name = absname
		else:
			sample_name = absname.split(separator)[0]
		sample_name = slugify(sample_name)
		sample_outputdir = outputdir + sample_name + '/'
		flag = isPathExist(sample_outputdir)
		if flag is True:
			error("BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator/outputdir.")
			print "BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator/outputdir"
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		flag = create_folder(sample_outputdir)
		if flag is False:
			error("Can not create folder at specified path: " + sample_outputdir)
			print "Can not create folder at specified path: ", sample_outputdir
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		merge_and_relabeled_fasta = sample_outputdir + sample_name + '_merged_and_relabeled.fasta'
		merge_and_relabeled_fastq = sample_outputdir + sample_name + '_merged_and_relabeled.fastq'
		
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "MERGING PAIRS:"
		report("MERGING PAIRS:")
		print forward_fastq
		report(forward_fastq)
		print reverse_fastq, '\n'
		report(reverse_fastq + '\n')
		flag = fastq_merge_and_relable(forward_fastq, reverse_fastq, merge_and_relabeled_fastq, merge_and_relabeled_fasta, sample_name, ngs_qscore, relabel, usearch_exec_path, processors, sample_outputdir)
		if flag is False:
			error("fastq_merge_and_relable generated error while processing these files: " + forward_fastq + '\t' + reverse_fastq + '\n')
			print "fastq_merge_and_relable generated error while processing these files: ", forward_fastq, reverse_fastq
			print "ABORTING!!!"
			sys.exit(2)
		fastq_total_list.append(merge_and_relabeled_fastq)
		fasta_total_list.append(merge_and_relabeled_fasta)
		
		print "MERGING PAIRS: ", CHECK_MARK
		report("MERGING PAIRS: " + CHECK_MARK)
		read_count = fasta_read_count(merge_and_relabeled_fasta)
		print sample_name, " read count is: ", read_count
		report(sample_name + " read count is: " + read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = merge_files(fastq_total_list, total_fastq)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		sys.exit(2)
	flag = merge_files(fasta_total_list, total_fasta)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		sys.exit(2)
	return True


def fastq_single_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	# relabel and convert fastq files
	fastq_total_list = []
	fasta_total_list = []
	file_paths_handle = open(files_path, 'rU')
	for path_line in file_paths_handle:
		path_line = path_line.rstrip()
		forward_fastq = path_line.split('\t')[0]
		path, absname, ext = split_file_name(forward_fastq)
		if separator is False:
			sample_name = absname
		else:
			sample_name = absname.split(separator)[0]
		sample_name = slugify(sample_name)
		sample_outputdir = outputdir + sample_name + '/'
		flag = isPathExist(sample_outputdir)
		if flag is True:
			error("BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator.")
			print "BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator."
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		flag = create_folder(sample_outputdir)
		if flag is False:
			error("Can not create folder at specified path: " + sample_outputdir)
			print "Can not create folder at specified path: ", sample_outputdir
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		# ##############
		relabeled_fastq = sample_outputdir + sample_name + '_relabeled.fastq'
		relabeled_fasta = sample_outputdir + sample_name + '_relabeled.fasta'
		tag = sample_name + ';S_'
		flag, stderr = execute_functions(relable_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, forward_fastq, tag, relabeled_fastq, relabeled_fasta, ngs_qscore, relabel)
		if flag is False:
			print "Execution of relable_fastq failed!!!"
		# ##############
		fastq_total_list.append(relabeled_fastq)
		fasta_total_list.append(relabeled_fasta)
		read_count = fasta_read_count(relabeled_fasta)
		print sample_name, " read count is: ", read_count
		report(sample_name + " read count is: " + read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = merge_files(fastq_total_list, total_fastq)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	flag = merge_files(fasta_total_list, total_fasta)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	return True


def fastq_mixed_parser(files_path, total_fasta, total_fastq, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	# relabel and convert fastq files
	fastq_total_list = []
	fasta_total_list = []
	file_paths_handle = open(files_path, 'rU')
	for path_line in file_paths_handle:
		path_line = path_line.rstrip()
		forward_fastq = path_line.split('\t')[0]
		path, absname, ext = split_file_name(forward_fastq)
		sample_name = absname
		"""if separator is False:
			sample_name = absname
		else:
			sample_name = absname.split(separator)[0]"""
		sample_name = slugify(sample_name)
		sample_outputdir = outputdir + sample_name + '/'
		flag = isPathExist(sample_outputdir)
		if flag is True:
			error("BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator.")
			print "BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator."
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		flag = create_folder(sample_outputdir)
		if flag is False:
			error("Can not create folder at specified path: " + sample_outputdir)
			print "Can not create folder at specified path: ", sample_outputdir
			print "ABORTING!!!"
			report("ABORTING!!!")
			sys.exit(2)
		# ##############
		if relabel is False:
			relabeled_fastq = sample_outputdir + sample_name + '_relabeled.fastq'
			relabeled_fasta = sample_outputdir + sample_name + '_relabeled.fasta'
			copy_file(forward_fastq, relabeled_fastq)
			flag = execute_functions(fastq_to_fasta_converter, processors, outputdir, 'multi', 'usearch', usearch_exec_path, relabeled_fastq, relabeled_fasta, ngs_qscore)
			if flag is False:
				print "Execution of relable_fastq failed!!!"
			# ##############
			fastq_total_list.append(relabeled_fastq)
			fasta_total_list.append(relabeled_fasta)
			read_count = fasta_read_count(relabeled_fasta)
			print sample_name, " read count is: ", read_count
			report(sample_name + " read count is: " + read_count)
		else:
			# Grab headers of fastq file
			exploded_fastq_list = []
			flag = explode_fastq(forward_fastq, separator, exploded_fastq_list, usearch_exec_path, processors, sample_outputdir)
			for each_exploded_fastq in exploded_fastq_list:
				path, absname, ext = split_file_name(each_exploded_fastq)
				sample_name = absname
				relabeled_fastq = sample_outputdir + sample_name + '_relabeled.fastq'
				relabeled_fasta = sample_outputdir + sample_name + '_relabeled.fasta'
				tag = sample_name + ';S_'
				flag, stderr = execute_functions(relable_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, each_exploded_fastq, tag, relabeled_fastq, relabeled_fasta, ngs_qscore, relabel)
				if flag is False:
					print "Execution of relable_fastq failed!!!"
				# ##############
				fastq_total_list.append(relabeled_fastq)
				fasta_total_list.append(relabeled_fasta)
				read_count = fasta_read_count(relabeled_fasta)
				print sample_name, " read count is: ", read_count
				report(sample_name + " read count is: " + read_count)
				report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
				print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = merge_files(fastq_total_list, total_fastq)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	flag = merge_files(fasta_total_list, total_fasta)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		report("ABORTING!!!")
		sys.exit(2)
	return True


def fasta_single_parser(files_path, total_fasta, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	# relabeling fasta file and concat them into single big one
	fasta_total_list = []
	file_paths_handle = open(files_path, 'rU')
	for path_line in file_paths_handle:
		path_line = path_line.rstrip()
		forward_fasta = path_line.split('\t')[0]
		path, absname, ext = split_file_name(forward_fasta)
		if separator is False:
			sample_name = absname
		else:
			sample_name = absname.split(separator)[0]
		sample_name = slugify(sample_name)
		sample_outputdir = outputdir + sample_name + '/'
		flag = isPathExist(sample_outputdir)
		if flag is True:
			error("BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator.")
			print "BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator"
			print "ABORTING!!!"
			sys.exit(2)
		flag = create_folder(sample_outputdir)
		if flag is False:
			error("Can not create folder at specified path: " + sample_outputdir)
			print "Can not create folder at specified path: ", sample_outputdir
			print "ABORTING!!!"
			sys.exit(2)
		relabeled_fasta = sample_outputdir + sample_name + '_relabeled.fasta'
		if relabel is False:
			copy_file(forward_fasta, relabeled_fasta)
		else:
			flag = relable_single_fasta_reads(forward_fasta, sample_name, relabeled_fasta)
			if flag is False:
				print "Execution of relable_fasta failed!!!"

		fasta_total_list.append(relabeled_fasta)
		read_count = fasta_read_count(relabeled_fasta)
		print sample_name, " read count is: ", read_count
		report(sample_name + " read count is: " + read_count)
	flag = merge_files(fasta_total_list, total_fasta)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		sys.exit(2)
	return True


def fasta_mixed_parser(files_path, total_fasta, separator, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	# relabeling fasta file and concat them into single big one
	fasta_total_list = []
	file_paths_handle = open(files_path, 'rU')
	for path_line in file_paths_handle:
		path_line = path_line.rstrip()
		forward_fasta = path_line.split('\t')[0]
		path, absname, ext = split_file_name(forward_fasta)
		sample_name = slugify(absname)
		sample_outputdir = outputdir + sample_name + '/'
		flag = isPathExist(sample_outputdir)
		if flag is True:
			error("BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator.")
			print "BAD SEPARATOR/OUTPUTDIR(--outputdir) Contain similar directory: Overwriting prohibitted, please choose another separator"
			print "ABORTING!!!"
			sys.exit(2)
		flag = create_folder(sample_outputdir)
		if flag is False:
			error("Can not create folder at specified path: " + sample_outputdir)
			print "Can not create folder at specified path: ", sample_outputdir
			print "ABORTING!!!"
			sys.exit(2)

		relabeled_fasta = sample_outputdir + sample_name + '_relabeled.fasta'
		if relabel is False:
			copy_file(forward_fasta, relabeled_fasta)
		else:

			flag = relable_mixed_fasta_reads(forward_fasta, separator, relabeled_fasta)
			if flag is False:
				print "Execution of relable_fasta_reads failed!!!"

		fasta_total_list.append(relabeled_fasta)
		read_count = fasta_read_count(relabeled_fasta)
		print sample_name, " read count is: ", read_count
		report(sample_name + " read count is: " + read_count)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	flag = merge_files(fasta_total_list, total_fasta)
	if flag is not True:
		print "something went wrong on writing files during merging"
		error("something went wrong on writing files during merging")
		print "ABORTING!!!"
		sys.exit(2)
	return True


def fastq_merge_and_relable(R1_file, R2_file, total_fastq, total_fasta, sample_name, ngs_qscore, relabel, usearch_exec_path, processors, outputdir):
	#error_log = "something went wrong on writing files during merging/relabeling"
	merged_fastq = outputdir + sample_name + '_merged_temp.fastq'
	not_merged_fastq_fwd = outputdir + sample_name + '_not_merged_fwd_temp.fastq'
	not_merged_fastq_rev = outputdir + sample_name + '_not_merged_rev_temp.fastq'
	# ##############
	flag, stderr = execute_functions(usearch_fastq_mergepairs, processors, outputdir, 'multi', 'usearch', usearch_exec_path, R1_file, R2_file, merged_fastq, not_merged_fastq_fwd, not_merged_fastq_rev)
	if flag is False:
		print "Execution of usearch_fastq_mergepairs failed!!!"
		report("ABORTING!!!")
		print "ABORTING!!!"
		sys.exit(2)
	# ##############
	tag = sample_name + ';M_'
	merged_fastq_relabled = outputdir + sample_name + '_merged.fastq'
	merged_fasta_relabled = outputdir + sample_name + '_merged.fasta'
	#merged_fastq_discarded = outputdir + sample_name + '_merged_discarded.fastq'
	#merged_fasta_discarded = outputdir + sample_name + '_merged_discarded.fasta'
	flag, stderr = execute_functions(relable_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, merged_fastq, tag, merged_fastq_relabled, merged_fasta_relabled, ngs_qscore, relabel)
	if flag is False:
		print "Execution of relable_fastq failed!!!"
	# ##############
	tag = sample_name + ';F_'
	not_merged_fastq_fwd_relabled = outputdir + sample_name + '_not_merged_fwd.fastq'
	not_merged_fasta_fwd_relabled = outputdir + sample_name + '_not_merged_fwd.fasta'
	#not_merged_fastq_fwd_discarded = outputdir + sample_name + '_fwd_discarded.fastq'
	#not_merged_fasta_fwd_discarded = outputdir + sample_name + '_fwd_discarded.fasta'
	flag, stderr = execute_functions(relable_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, not_merged_fastq_fwd, tag, not_merged_fastq_fwd_relabled, not_merged_fasta_fwd_relabled, ngs_qscore, relabel)
	if flag is False:
		print "Execution of relable_fastq failed!!!"
	# ##############
	tag = sample_name + ';R_'
	not_merged_fastq_rev_relabled = outputdir + sample_name + '_not_merged_rev.fastq'
	not_merged_fasta_rev_relabled = outputdir + sample_name + '_not_merged_rev.fasta'
	#not_merged_fastq_rev_discarded = outputdir + sample_name + '_rev_discarded.fastq'
	#not_merged_fasta_rev_discarded = outputdir + sample_name + '_rev_discarded.fasta'
	flag, stderr = execute_functions(relable_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, not_merged_fastq_rev, tag, not_merged_fastq_rev_relabled, not_merged_fasta_rev_relabled, ngs_qscore, relabel)
	if flag is False:
		print "Execution of relable_fastq failed!!!"
	# ##############
	flag = merge_files([merged_fastq_relabled, not_merged_fastq_fwd_relabled, not_merged_fastq_rev_relabled], total_fastq)
	if flag is False:
		print "Execution of merge_files failed!!!"
		
	flag = merge_files([merged_fasta_relabled, not_merged_fasta_fwd_relabled, not_merged_fasta_rev_relabled], total_fasta)
	if flag is False:
		print "Execution of merge_files failed!!!"

	check_it_and_remove_it(merged_fastq)
	check_it_and_remove_it(not_merged_fastq_rev)
	check_it_and_remove_it(not_merged_fastq_fwd)
	check_it_and_remove_it(merged_fastq_relabled)
	check_it_and_remove_it(not_merged_fastq_fwd_relabled)
	check_it_and_remove_it(not_merged_fastq_rev_relabled)
	check_it_and_remove_it(merged_fasta_relabled)
	check_it_and_remove_it(not_merged_fasta_fwd_relabled)
	check_it_and_remove_it(not_merged_fasta_rev_relabled)

	return True


# ################################# USEARCH_FUNCTIONS ######################## #


def usearch_uclust_cluster(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, centroids_fasta, uparse_file):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-cluster_fast' + space
	usearch_string += fasta_file + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-id 0.97' + space
	usearch_string += '-sort size' + space
	usearch_string += '-strand both' + space
	usearch_string += '-sizein' + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-centroids' + space + centroids_fasta + space
	usearch_string += '-uc' + space + uparse_file + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_usearch_global(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, reference, ucout, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout, biomout, identity, maxhits, refcount):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-usearch_global' + space
	usearch_string += fasta_file + space
	usearch_string += '-db' + space + reference + space
	usearch_string += '-id' + space + identity + space
	usearch_string += '-strand both' + space
	usearch_string += '-maxhits' + space + maxhits + space
	usearch_string += '-maxaccepts' + space + refcount + space
	usearch_string += '-maxrejects 0' + space
	usearch_string += '-threads' + space + processors + space
	if ucout is not False:
		usearch_string += '-uc' + space + ucout + space
	if biomout is not False:
		usearch_string += '-biomout' + space + biomout + space
	if mothurout is not False:
		usearch_string += '-mothur_shared_out' + space + mothurout + space
	if qiimeout is not False:
		usearch_string += '-otutabout' + space + qiimeout + space
	if tsegout is not False:
		usearch_string += '-tsegout' + space + tsegout + space
	if unaligned_reads is not False:
		usearch_string += '-notmatched' + space + unaligned_reads + space
	if aligned_reads is not False:
		usearch_string += '-matched' + space + aligned_reads + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_usearch_global_supervised(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, reference, ucout, tsegout, aligned_reads, unaligned_reads, qiimeout, mothurout, biomout, identity, maxhits):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-usearch_global' + space
	usearch_string += fasta_file + space
	usearch_string += '-db' + space + reference + space
	usearch_string += '-id' + space + identity + space
	usearch_string += '-strand both' + space
	usearch_string += '-maxhits' + space + maxhits + space
	usearch_string += '-maxaccepts 32' + space
	usearch_string += '-maxrejects 32' + space
	usearch_string += '-threads' + space + processors + space
	if ucout is not False:
		usearch_string += '-uc' + space + ucout + space
	if biomout is not False:
		usearch_string += '-biomout' + space + biomout + space
	if mothurout is not False:
		usearch_string += '-mothur_shared_out' + space + mothurout + space
	if qiimeout is not False:
		usearch_string += '-otutabout' + space + qiimeout + space
	if tsegout is not False:
		usearch_string += '-tsegout' + space + tsegout + space
	if unaligned_reads is not False:
		usearch_string += '-notmatched' + space + unaligned_reads + space
	if aligned_reads is not False:
		usearch_string += '-matched' + space + aligned_reads + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_uparse_cluster(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, centroids_fasta, uparse_file):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-cluster_otus' + space
	usearch_string += fasta_file + space
	#usearch_string += '-relabel OTU' + space
	usearch_string += '-sizein' + space
	usearch_string += '-otu_radius_pct 3.0' + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-otus' + space + centroids_fasta + space
	usearch_string += '-uparseout' + space + uparse_file + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_sort_by_size(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, sorted_fasta):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-sortbysize' + space
	usearch_string += fasta_file + space
	usearch_string += '-fastaout' + space + sorted_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_sort_by_length(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, sorted_fasta):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-sortbylength' + space
	usearch_string += fasta_file + space
	usearch_string += '-fastaout' + space + sorted_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_precluster(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, preclustered_fasta, maxdiffs):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-cluster_smallmem' + space
	usearch_string += fasta_file + space
	usearch_string += '-id 0.99' + space
	usearch_string += '-sortedby size' + space
	usearch_string += '-strand both' + space
	usearch_string += '-maxdiffs' + space + maxdiffs + space
	usearch_string += '-centroids' + space + preclustered_fasta
	
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_sort_abundance(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, sorted_fasta, min_abund_size):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-sortbysize' + space
	usearch_string += fasta_file + space
	usearch_string += '-minsize' + space + min_abund_size + space
	usearch_string += '-fastaout' + space + sorted_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_dereplicate(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, derep_fasta):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-derep_fulllength' + space
	usearch_string += fasta_file + space
	# usearch_string += '-relabel UNIQ' + space
	usearch_string += '-strand both' + space
	usearch_string += '-sizeout' + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-fastaout' + space + derep_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fasta_truncate(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, hq_fasta, trunc_len):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_truncate' + space
	usearch_string += fasta_file + space
	usearch_string += '-trunclen' + space + trunc_len + space
	usearch_string += '-fastaout' + space + hq_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fastq_truncate(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_file, hq_fastq, hq_fasta, trunc_len):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_truncate' + space
	usearch_string += fastq_file + space
	usearch_string += '-trunclen' + space + trunc_len + space
	usearch_string += '-fastaout' + space + hq_fasta + space
	usearch_string += '-fastqout' + space + hq_fastq
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fastq_filter(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_file, hq_fasta, lq_fasta, expected_error_rate, maxambig, minqscore, trunc_len, min_len, ngs_qscore):
	# filtering fastq file sbase on in put parameters
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastq_filter' + space
	usearch_string += fastq_file + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += ngs_qscore + space
	if maxambig is not False:
		usearch_string += '-fastq_maxns' + space + maxambig + space
	if expected_error_rate is not False:
		usearch_string += '-fastq_maxee' + space + expected_error_rate + space
	if trunc_len is not False:
		usearch_string += '-fastq_trunclen' + space + trunc_len + space
	if min_len is not False:
		usearch_string += '-fastq_minlen' + space + min_len + space
	if minqscore is not False:
		usearch_string += '-fastq_truncqual' + space + minqscore + space
	usearch_string += '-fastaout' + space + hq_fasta + space
	usearch_string += '-fastaout_discarded' + space + lq_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_chimera_removal(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, chimera_reference, chimeric_fasta, non_chimeric_fasta, chimera_report):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-uchime_ref' + space
	usearch_string += fasta_file + space
	usearch_string += '-db' + space + chimera_reference + space
	usearch_string += '-strand plus' + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-self -selfid' + space
	usearch_string += '-nonchimeras' + space + non_chimeric_fasta + space
	usearch_string += '-chimeras' + space + chimeric_fasta + space
	usearch_string += '-uchimeout' + space + chimera_report
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fasta_split(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, split_value):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_split' + space
	usearch_string += fasta_file + space
	usearch_string += '-splits' + space + str(split_value) + space
	usearch_string += '-outname ' + outputdir + 'split@.fasta'
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def fastq_to_fasta_converter(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_file, fasta_file, ngs_qscore):
	# conversion of fast to fasta using usearch
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastq_filter' + space
	usearch_string += fastq_file + space + '-threads' + space + processors + space
	#usearch_string += '-fastq_minlen 10' + space
	usearch_string += ngs_qscore + space
	usearch_string += '-fastaout' + space + fasta_file + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def relable_fastq(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_file, tag, new_name_fastq, new_name_fasta, ngs_qscore, relabel):
	# relabel fastq files one by one based on tag
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastq_filter' + space
	usearch_string += fastq_file + space + '-threads' + space + processors + space
	if relabel is not False:
		usearch_string += '-relabel ' + space + '"' + tag + '"' + space
		usearch_string += '-fastq_eeout' + space
	usearch_string += ngs_qscore + space
	usearch_string += '-fastqout' + space + new_name_fastq + space
	usearch_string += '-fastaout' + space + new_name_fasta + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def relable_fasta2(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_file, tag, new_name_fasta, ngs_qscore, relabel):
	# relabel fastq files one by one based on tag
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_truncate' + space
	usearch_string += fasta_file + space
	usearch_string += '-stripright 0' + space
	if relabel is not False:
		usearch_string += '-relabel ' + space + '"' + tag + '"' + space
		usearch_string += '-label_suffix ' + space + '";ee=0.0"' + space
	usearch_string += '-fastaout' + space + new_name_fasta + space
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fastq_mergepairs(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, R1_file, R2_file, merged_fastq, not_merged_fastq_fwd, not_merged_fastq_rev):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastq_mergepairs' + space
	usearch_string += R1_file + space
	usearch_string += '-reverse' + space + R2_file + space 
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-fastqout' + space + merged_fastq + space
	usearch_string += '-fastqout_notmerged_fwd' + space + not_merged_fastq_fwd + space
	usearch_string += '-fastqout_notmerged_rev' + space + not_merged_fastq_rev + space

	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_extract_fasta(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fasta_headers_file, main_fasta_file, extracted_fasta):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_getseqs' + space
	usearch_string += main_fasta_file + space
	usearch_string += '-labels' + space + fasta_headers_file + space
	usearch_string += '-fastaout' + space + extracted_fasta
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_extract_fastq(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_headers_file, main_fastq_file, extracted_fastq):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_getseqs' + space
	usearch_string += main_fastq_file + space
	usearch_string += '-labels' + space + fastq_headers_file + space
	usearch_string += '-fastqout' + space + extracted_fastq
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_extract_fastq_partial(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, sample_name, main_fastq_file, extracted_fastq):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastx_getseqs' + space
	usearch_string += main_fastq_file + space
	usearch_string += '-label_word' + space + '"' + sample_name + '"' + space
	usearch_string += '-label_substr_match' + space
	usearch_string += '-fastqout' + space + extracted_fastq
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def usearch_fastq_header_extract(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path, fastq_file, header_file):
	space = ' '
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	usearch_string = ''
	usearch_string += '-fastq_filter' + space
	usearch_string += fastq_file + space
	usearch_string += '-fastq_trunclen 0' + space
	usearch_string += '-threads' + space + processors + space
	usearch_string += '-fastaout' + space + header_file
	commandline_list = []
	commandline_list.append(usearch_string)
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def test_usearch(processors, outputdir, stderr, stdout, run_pid, usearch_exec_path):
	# test mothur to see if it is working and grab the version of mothur by scanning its output log
	usearch_input_dictionary = {}
	usearch_input_dictionary['usearch_exec_path'] = usearch_exec_path
	commandline_list = []
	commandline_list.append('-help')
	usearch_input_dictionary['commandline'] = commandline_list
	usearch_input_dictionary['nohup_in'] = 'nohup'
	usearch_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	make_file_object = usearch_process(usearch_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_usearch_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True

# ################################# MOTHUR_FUNCTIONS ######################### #


def mothur_get_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, header_list_file):
	mothur_input_dictionary = {}
	space = ' '
	command = 'get.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'accnos=' + header_list_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_unifrac_unweighted(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, tree_file, group_file):
	mothur_input_dictionary = {}
	space = ' '
	command = 'unifrac.unweighted'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('tree=' + tree_file)
	parameter_list.append(',' + space + 'group=' + group_file + ',' + space + 'distance=square, iters=1000')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_unifrac_weighted(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, tree_file, group_file):
	mothur_input_dictionary = {}
	space = ' '
	command = 'unifrac.weighted'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('tree=' + tree_file)
	parameter_list.append(',' + space + 'group=' + group_file + ',' + space + 'distance=square, iters=1000')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_classify_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, reference, taxonomy):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'classify.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'reference=' + reference)
	parameter_list.append(',' + space + 'taxonomy=' + taxonomy)
	parameter_list.append(',' + space + 'probs=F, method=knn, cutoff=80, iters=100, numwanted=1')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_oturep(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, distance_file, list_file, fasta_file, count_table):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.oturep'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('column=' + distance_file)
	parameter_list.append(',' + space + 'list=' + list_file)
	parameter_list.append(',' + space + 'fasta=' + fasta_file)
	parameter_list.append(',' + space + 'count=' + count_table)
	parameter_list.append(',' + space + 'label=0.03')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_cluster_split(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, distance_file, count_table):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'cluster.split'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('column=' + distance_file)
	parameter_list.append(',' + space + 'count=' + count_table + ',' + space + 'large=T, method=average, cutoff=0.03')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_dist_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, cutoff_limit):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'dist.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	parameter_list.append(',' + space + 'cutoff=' + cutoff_limit)
	parameter_list.append(',' + space + 'output=column')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_summary_tax(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, taxonomy_file, group_table):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'summary.tax'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('taxonomy=' + taxonomy_file + ',' + space)
	parameter_list.append('group=' + group_table)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_fastq_info(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fastq_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'fastq.info'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fastq=' + fastq_file)
	parameter_list.append(', fasta=F, qfile=T')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_summary_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	command = 'summary.seqs'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_screen_seqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, fasta_file, maxambig, maxhomop, min_len):
	mothur_input_dictionary = {}
	command = 'screen.seqs'
	space = ' '
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('fasta=' + fasta_file)
	if maxambig is not False:
		parameter_list.append(',' + space + 'maxn=' + maxambig)
	if min_len is not False:
		parameter_list.append(',' + space + 'minlength=' + min_len)
	if maxhomop is not False:
		parameter_list.append(',' + space + 'maxhomop=' + maxhomop)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_make_file(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, input_file_directory):
	# grab path of paired fastq files and save it into container
	mothur_input_dictionary = {}
	command = 'make.file'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('inputdir=' + input_file_directory)
	mothur_input_dictionary['parameters'] = parameter_list
	#here we construct the object
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def test_mothur(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path):
	# test mothur to see if it is working and grab the version of mothur by scanning its output log
	mothur_input_dictionary = {}
	command = 'get.current'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	#parameter_list = []
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def parse_mothur_logfile(logfile, keyword, output_container, mode=None):
	if mode == 'file':
		f = open(logfile, 'rU')
		data = f.read().split('\n')
	else:
		data = logfile
	for line in data.split('\n'):
		if keyword in line:
			output_container.append(line)
	if len(output_container) < 1:
		return False
	return True


# ################################# FASTTREE FUNCTIONS ####################### #

def test_fasttree(processors, outputdir, stderr, stdout, run_pid, fasttree_exec_path):
	space = ' '
	fasttree_string = 'nohup' + space + fasttree_exec_path + space + '> ' + stderr + ' 2> ' + stdout + ' & echo $! > ' + run_pid
	print "EXECUTING: ", fasttree_string
	report(fasttree_string)
	exec_dict = {}
	exec_dict = execute([fasttree_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


def fasttree_build(processors, outputdir, stderr, stdout, run_pid, fasttree_exec_path, fasta_file, newick_tree):
	space = ' '
	fasttree_string = 'nohup' + space + fasttree_exec_path + space + '-noml -nome -fastest -out' + space + newick_tree + space + '-nt' + space + fasta_file + '> ' + stderr + ' 2> ' + stdout + ' & echo $! > ' + run_pid
	print fasttree_string
	report(fasttree_string)
	exec_dict = {}
	exec_dict = execute([fasttree_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


# ################################# MAFFT FUNCTIONS ########################## #

def mafft_align(processors, outputdir, stderr, stdout, run_pid, mafft_exec_path, fasta_file, mafft_out, mafft_parameter):
	space = ' '
	mafft_string = 'nohup' + space
	mafft_string += mafft_exec_path + space
	mafft_string += mafft_parameter + '--thread' + space + processors + space
	mafft_string += fasta_file + space
	#mafft_string += '> ' + mafft_out + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	mafft_string += '> ' + mafft_out + ' & echo $! > ' + run_pid
	exec_dict = {}
	print mafft_string
	report(mafft_string)
	exec_dict = execute([mafft_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


def test_mafft(processors, outputdir, stderr, stdout, run_pid, mafft_exec_path):
	space = ' '
	mafft_string = 'nohup' + space + mafft_exec_path + space + '--help' + space + '--thread' + space + processors + space + '> ' + stderr + ' 2> ' + stdout + ' & echo $! > ' + run_pid
	print "EXECUTING: ", mafft_string
	report(mafft_string)
	exec_dict = {}
	exec_dict = execute([mafft_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True

# ################################# SPECIFIC FUNCTIONS ####################### #


def qiime_to_count(qiime_table, count_table):
	string = ''
	f = open(qiime_table, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == '#OTU ID':
			string += 'Representative_Sequence' + '\t' + 'total' + '\t' + '\t'.join(line[1:]) + '\n'
		else:
			total_value = 0
			value_list = []
			for each_cell in line[1:]:
				total_value += float(format(each_cell))
				value_list.append(format(each_cell))
			string += line[0] + '\t' + str(total_value) + '\t' + '\t'.join(value_list) + '\n'
	o = open(count_table, 'w')
	o.write(string)
	o.close()
	return True


def biom_report_correct_supervised(taxonomy, biom_report_temp, biom_report):
	tax_dict = {}
	f = open(taxonomy, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		acc = line[0].split(';')[0]
		taxonom = line[1].split('_Strain')[0] + ';'
		tax_dict[acc] = taxonom
	f.close()
	f = open(biom_report_temp, 'rU')
	string = ''
	otu_flag = False
	for i in f:
		if 'metadata' in i:
			temp = i.split(',')[0]
			acc = temp.split(':')[1]
			acc = acc.replace('"', '').strip()
			if acc not in tax_dict:
				string += i
				continue
			tax_list = tax_dict[acc].split(';')
			tax_string = '"' + '","'.join(tax_list)[:-2]
			string += '\t\t{"id":"' + acc + '", ' + '"metadata"' + ':{"taxonomy":[' + tax_string + ']}},\n'
			otu_flag = True

		else:
			if otu_flag is True:
				string = string[:-2]
				string += '\n'
				otu_flag = False
			string += i
	o = open(biom_report, 'w')
	o.write(string)
	o.close()
	f.close()
	string = ''
	return True


def biom_report_correct(taxonomy, biom_report_temp, biom_report):
	tax_dict = {}
	f = open(taxonomy, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax_dict[line[1]] = line[3]
	f.close()
	f = open(biom_report_temp, 'rU')
	string = ''
	otu_flag = False
	for i in f:
		if 'metadata' in i:
			temp = i.split(',')[0]
			acc = temp.split(':')[1]
			acc = acc.replace('"', '').strip()
			if acc not in tax_dict:
				string += i
				continue
			tax_list = tax_dict[acc].split(';')
			tax_string = '"' + '","'.join(tax_list)[:-2]
			string += '\t\t{"id":"' + acc + '", ' + '"metadata"' + ':{"taxonomy":[' + tax_string + ']}},\n'
			otu_flag = True

		else:
			if otu_flag is True:
				string = string[:-2]
				string += '\n'
				otu_flag = False
			string += i
	o = open(biom_report, 'w')
	o.write(string)
	o.close()
	f.close()
	string = ''
	return True


def classify_centroids(centroids_fasta, relabeled_centroids_fasta, centroids_classification_report, taxonomy, centroids_label):
	tax_dict = {}
	f = open(taxonomy, 'r')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		acc = line[0].split(';')[0]
		length = line[0].split(';')[1]
		tax_dict[acc] = length
	f.close()

	classify_dict = {}
	f = open(centroids_classification_report)
	string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		acc = line[1].split('Strain_')[1]
		acc = acc.split(';')[0]
		acc = acc.split('_')[0]
		tax = line[1].split('_Strain_')[0]
		tax += ';'
		classify_dict[line[0]] = acc
		string += line[0] + '\t' + acc + '\t' + tax_dict[acc] + '\t' + tax + '\n'
	f.close()
	o = open(centroids_label, 'w')
	o.write(string)
	o.close()

	string = ''
	f = open(centroids_fasta, 'r')
	for i in f:
		if i[0] == '>':
			i = i.rstrip()
			i = i[1:]
			string += '>' + classify_dict[i] + '\n'
		else:
			string += i
	f.close()
	o = open(relabeled_centroids_fasta, 'w')
	o.write(string)
	o.close()
	return True


def classify_centroids2(centroids_fasta, centroids_taxonomy, centroids_classification_report, reference, taxonomy, usearch_exec_path, processors, outputdir):
	acc_labels = outputdir + 'acc_labels_STEP7.txt'
	f = open(centroids_classification_report, 'rU')
	acc_list = []
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		# centroids_label = line[0]
		tax = line[1]
		raw_acc = tax.split('_ACCESSION_')[1]
		acc = raw_acc.split('_END_ACCESSION')[0]
		acc = acc.replace('_', ';')
		acc_list.append(acc)
	f.close()
	acc_list = list(set(acc_list))
	acc_labels_string = ''
	for each_acc in acc_list:
		acc_labels_string += each_acc + '\n'
	o = open(acc_labels, 'w')
	o.write(acc_labels_string)
	o.close()
	# processing taxonomy
	tax_dict = {}
	f = open(taxonomy, 'r')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		acc = line[0]
		raw_tax = line[1]
		tax = raw_tax.split('_ACCESSION_')[0] + ';'
		species_name = acc.split(';')[0] + '_' + tax.split(';')[-2]
		tax_dict[acc] = (species_name, tax)
	f.close()
	
	# extract centroids from reference fasta file
	raw_centroids_fasta = outputdir + 'raw_centroids_STEP7.fasta'
	flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, acc_labels, reference, raw_centroids_fasta)
	if flag is False:
		print "Execution of usearch_extract_fasta failed!!!"
	#relabel fasta headers
	new_fasta_string = ''
	new_tax_string = ''

	f = open(raw_centroids_fasta, 'r')
	for i in f:
		if i[0] == '>':
			i = i.rstrip()
			i = i[1:]
			new_fasta_string += '>' + tax_dict[i][0] + '\n'
			new_tax_string += tax_dict[i][0] + '\t' + tax_dict[i][1] + '\n'
		else:
			new_fasta_string += i
	f.close()
	o = open(centroids_fasta, 'w')
	o.write(new_fasta_string)
	o.close()
	o = open(centroids_taxonomy, 'w')
	o.write(new_tax_string)
	o.close()

	return True


def mothur_count_group_table(fasta_file, count_table, group_table):
	count_string = 'Representative_Sequence\ttotal\n'
	group_string = ''
	f = open(fasta_file, 'rU')
	for i in f:
		i = i.rstrip()
		if i[0] == '>':
			count_string += i[1:] + '\t' + i.split('size=')[1][:-1] + '\n'
			group_string += i[1:] + '\t' + i[1:].split(';')[0] + '\n'
	f.close()
	f = open(count_table, 'w')
	f.write(count_string)
	f.close()
	f = open(group_table, 'w')
	f.write(group_string)
	f.close()
	group_string = ''
	count_string = ''
	return True


def phylotype(alignment_report, centroids_label, phylotype_file, mothur_exec_path, sample_name, processors, outputdir):
	# step 1: create taxonomy file
	tax_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax_dict[line[1]] = line[3]
	f.close()
	f = open(alignment_report, 'rU')
	tax_string = ''
	group_string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'H':
			tax_string += line[8] + '\t' + tax_dict[line[9]] + '\n'
			group_string += line[8] + '\t' + line[8].split(';')[0] + '\n'
		else:
			tax_string += line[8] + '\t' + 'unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;' + '\n'
			group_string += line[8] + '\t' + line[8].split(';')[0] + '\n'
	f.close()
	
	taxonomy_file = outputdir + sample_name + '_phylotype_STEP7.tax'
	check_it_and_remove_it(taxonomy_file)
	group_file = outputdir + sample_name + '_phylotype_STEP7.group'
	check_it_and_remove_it(group_file)
	f = open(taxonomy_file, 'w')
	f.write(tax_string)
	f.close()
	f = open(group_file, 'w')
	f.write(group_string)
	f.close()
	# step 2: create abundance file
	
	flag, stderr = execute_functions(mothur_summary_tax, processors, outputdir, 'multi', 'mothur', mothur_exec_path, taxonomy_file, group_file)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
	else:
		abundance_container = []
		extension_list = ['.tax.summary']
		flag = scandirs(outputdir, abundance_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(abundance_container[0], phylotype_file)
	return True


def phylotype_supervised(alignment_report, centroids_label, phylotype_file, mothur_exec_path, sample_name, processors, outputdir):
	# step 1: create taxonomy file
	tax_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax = line[1].split(';')
		species_name = tax[-2].split('_Strain')[0] + ';'
		tax_dict[line[0]] = ';'.join(tax[0:-2]) + ';' + species_name
	f.close()
	f = open(alignment_report, 'rU')
	tax_string = ''
	group_string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'H':
			tax_string += line[8] + '\t' + tax_dict[line[9]] + '\n'
			group_string += line[8] + '\t' + line[8].split(';')[0] + '\n'
		else:
			tax_string += line[8] + '\t' + 'unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;' + '\n'
			group_string += line[8] + '\t' + line[8].split(';')[0] + '\n'
	f.close()
	
	taxonomy_file = outputdir + sample_name + '_phylotype_STEP7.tax'
	check_it_and_remove_it(taxonomy_file)
	group_file = outputdir + sample_name + '_phylotype_STEP7.group'
	check_it_and_remove_it(group_file)
	f = open(taxonomy_file, 'w')
	f.write(tax_string)
	f.close()
	f = open(group_file, 'w')
	f.write(group_string)
	f.close()
	# step 2: create abundance file
	
	flag, stderr = execute_functions(mothur_summary_tax, processors, outputdir, 'multi', 'mothur', mothur_exec_path, taxonomy_file, group_file)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
	else:
		abundance_container = []
		extension_list = ['.tax.summary']
		flag = scandirs(outputdir, abundance_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(abundance_container[0], phylotype_file)
	return True


def create_sas_table(alignment_report, centroids_label, sas_file, sample_name, file_type, total_error_rate, processors, outputdir, mode=None):
	tax_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax = line[3].split(';')
		species_name = tax[-2]
		tax_dict[line[1]] = line[1] + ';;' + line[2] + ';;' + line[3] + ';;' + species_name
	f.close()
	# ###############################
	error_rate_flag = True
	error_dict = {}
	f = open(total_error_rate, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		error_dict[line[0]] = '\t'.join(line[1:])
	f.close()
	# ###########################################################
	sas_string = ''
	sas_string += 'Library_id' + '\t' + 'Read_id' + '\t' + 'Query_length' + '\t' + 'Reference_length' + '\t' + 'Compressed_alignment' + '\t' + 'Expected_error_rate' + '\t' + 'Average_error_probability' + '\t' + 'Average_Phred_score' + '\t'
	sas_string += 'Identity' + '\t' + 'Ref_id' + '\t' + 'Ref_acc' + '\t' + 'BEI' + '\t'
	sas_string += 'Kingdom_name' + '\t' + 'Phylum_name' + '\t' + 'Class_name' + '\t' + 'Order_name' + '\t'
	sas_string += 'Family_name' + '\t' + 'Genus_name' + '\t' + 'Species_name' + '\n'
	f = open(alignment_report, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tag = line[0]
		if tag == 'H':
			library_id = line[8].split(';')[0]
			temp_id = line[8].split(';')[1]
			temp_id = temp_id.split(';')[0]
			if 'M_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'M__' + count.zfill(9)
			elif 'S_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'S__' + count.zfill(9)
			elif 'R_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R2_' + count.zfill(9)
			elif 'F_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R1_' + count.zfill(9)
			hit_id = line[9]
			full_annotation = tax_dict[hit_id]
			target_id = full_annotation.split(';;')[0] + full_annotation.split(';;')[3]
			target_acc = hit_id
			target_length = full_annotation.split(';;')[1]
			seq_length = line[2]
			identity = line[3]
			alignment = line[7]
			if 'Bei' in line[9]:
				bei = '1'
			else:
				bei = '0'
			tax = full_annotation.split(';;')[2].split(';')
			sas_string += library_id + '\t' + read_id + '\t' + seq_length + '\t' + target_length + '\t' + alignment + '\t'
			if error_rate_flag is False:
				sas_string += '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t'
			elif error_rate_flag is True:
				sas_string += error_dict[line[8]] + '\t'
			sas_string += identity + '\t' + target_id + '\t' + target_acc + '\t' + bei + '\t'
			sas_string += tax[0] + '\t' + tax[1] + '\t' + tax[2] + '\t' + tax[3] + '\t' + tax[4] + '\t' + tax[5] + '\t' + tax[6] + '\n'
		elif tag == 'N':
			library_id = line[8].split(';')[0]
			temp_id = line[8].split(';')[1]
			temp_id = temp_id.split(';')[0]
			if 'M_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'M__' + count.zfill(9)
			elif 'R_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R2_' + count.zfill(9)
			elif 'F_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R1_' + count.zfill(9)
			hit_id = line[9]
			seq_length = line[2]
			sas_string += library_id + '\t' + read_id + '\t' + seq_length + '\t' + '.' + '\t' + '' + '\t'
			if error_rate_flag is False:
				sas_string += '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t'
			elif error_rate_flag is True:
				sas_string += error_dict[line[8]] + '\t'
			sas_string += '.' + '\t' + '.' + '\t' + '' + '\t' + '.' + '\t'
			sas_string += '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n'
	tax_dict = ''
	f.close()
	f = open(sas_file, 'w')
	f.write(sas_string)
	f.close()
	return True


def create_sas_table_supervised(alignment_report, centroids_label, sas_file, sample_name, file_type, total_error_rate, processors, outputdir, mode=None):
	tax_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax = line[1].split(';')
		species_name = tax[-2].split('_Strain')[0] + ';'
		acc = line[0].split(';')[0]
		acc_length = line[0].split(';')[1]
		tax_dict[line[0]] = acc + ';;' + acc_length + ';;' + ';'.join(tax[0:-2]) + ';' + species_name + ';;' + species_name
	f.close()
	# ###############################
	error_rate_flag = True
	error_dict = {}
	f = open(total_error_rate, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		error_dict[line[0]] = '\t'.join(line[1:])
	f.close()
	# ###########################################################
	sas_string = ''
	sas_string += 'Library_id' + '\t' + 'Read_id' + '\t' + 'Query_length' + '\t' + 'Reference_length' + '\t' + 'Compressed_alignment' + '\t' + 'Expected_error_rate' + '\t' + 'Average_error_probability' + '\t' + 'Average_Phred_score' + '\t'
	sas_string += 'Identity' + '\t' + 'Ref_id' + '\t' + 'Ref_acc' + '\t' + 'BEI' + '\t'
	sas_string += 'Kingdom_name' + '\t' + 'Phylum_name' + '\t' + 'Class_name' + '\t' + 'Order_name' + '\t'
	sas_string += 'Family_name' + '\t' + 'Genus_name' + '\t' + 'Species_name' + '\n'
	f = open(alignment_report, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tag = line[0]
		if tag == 'H':
			library_id = line[8].split(';')[0]
			temp_id = line[8].split(';')[1]
			temp_id = temp_id.split(';')[0]
			if 'M_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'M__' + count.zfill(9)
			elif 'S_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'S__' + count.zfill(9)
			elif 'R_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R2_' + count.zfill(9)
			elif 'F_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R1_' + count.zfill(9)
			hit_id = line[9]
			full_annotation = tax_dict[hit_id]
			target_id = full_annotation.split(';;')[0] + full_annotation.split(';;')[3]
			target_acc = hit_id.split(';')[0]
			target_length = full_annotation.split(';;')[1]
			seq_length = line[2]
			identity = line[3]
			alignment = line[7]
			if 'Bei' in line[9]:
				bei = '1'
			else:
				bei = '0'
			tax = full_annotation.split(';;')[2].split(';')
			sas_string += library_id + '\t' + read_id + '\t' + seq_length + '\t' + target_length + '\t' + alignment + '\t'
			if error_rate_flag is False:
				sas_string += '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t'
			elif error_rate_flag is True:
				sas_string += error_dict[line[8]] + '\t'
			sas_string += identity + '\t' + target_id + '\t' + target_acc + '\t' + bei + '\t'
			sas_string += tax[0] + '\t' + tax[1] + '\t' + tax[2] + '\t' + tax[3] + '\t' + tax[4] + '\t' + tax[5] + '\t' + tax[6] + '\n'
		elif tag == 'N':
			library_id = line[8].split(';')[0]
			temp_id = line[8].split(';')[1]
			temp_id = temp_id.split(';')[0]
			if 'M_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'M__' + count.zfill(9)
			elif 'R_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R2_' + count.zfill(9)
			elif 'F_' in temp_id:
				count = temp_id.split('_')[1]
				read_id = library_id + ';' + 'R1_' + count.zfill(9)
			hit_id = line[9]
			seq_length = line[2]
			sas_string += library_id + '\t' + read_id + '\t' + seq_length + '\t' + '.' + '\t' + '' + '\t'
			if error_rate_flag is False:
				sas_string += '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t'
			elif error_rate_flag is True:
				sas_string += error_dict[line[8]] + '\t'
			sas_string += '.' + '\t' + '.' + '\t' + '' + '\t' + '.' + '\t'
			sas_string += '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n'
	tax_dict = ''
	f.close()
	f = open(sas_file, 'w')
	f.write(sas_string)
	f.close()
	return True


def mothur_shared_correct(centroids_label, mothurout, corrected_mothurout, mode=None):
	label_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax = line[3].split(';')
		species_name = tax[-2]
		label_dict[line[1]] = line[1] + ';' + species_name
	f.close()
	f = open(mothurout, 'rU')
	string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'label':
			for j in line:
				if j.lower() in ['numotus', 'label', 'group']:
					string += j + '\t'
				else:
					if j not in label_dict:
						string += j
						string += '\t'
					else:
						tax = label_dict[j]
						#tax = j.split(';')
						string += tax
						string += '\t'
					
			string = string[:-1]
			string += '\n'
		else:
			string += '0.03' + '\t' + '\t'.join(line[1:]) + '\n'
	f.close()
	f = open(corrected_mothurout, 'w')
	f.write(string)
	f.close()
	return True


def mothur_shared_correct_supervised(centroids_label, mothurout, corrected_mothurout, mode=None):
	label_dict = {}
	f = open(centroids_label, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		tax = line[1].split(';')
		species_name = tax[-2]
		acc = line[0].split(';')[0]
		label_dict[acc] = line[0].split(';')[0] + ';' + species_name.split('_Strain')[0]
	f.close()
	f = open(mothurout, 'rU')
	string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'label':
			for j in line:
				if j.lower() in ['numotus', 'label', 'group']:
					string += j + '\t'
				else:
					if j not in label_dict:
						string += j
						string += '\t'
					else:
						tax = label_dict[j]
						#tax = j.split(';')
						string += tax
						string += '\t'
					
			string = string[:-1]
			string += '\n'
		else:
			string += '0.03' + '\t' + '\t'.join(line[1:]) + '\n'
	f.close()
	f = open(corrected_mothurout, 'w')
	f.write(string)
	f.close()
	return True


def error_rate_calculation(fasta_file, phred_error_rate, mothur_exec_path, processors, outputdir):
	flag, stderr = execute_functions(mothur_summary_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file)
	if flag is False:
		print "Execution of mothur_dist_seqs failed!!!"
	qual_container = []
	extension_list = ['.summary']
	flag = scandirs(outputdir, qual_container, extension_list)
	if flag is False:
		print "This extension is not availble: ", list_to_string(extension_list, ',')
		report("This extension is not availble: " + list_to_string(extension_list, ','))
		sys.exit(2)
	error_rate_string = ''
	for each_summary in qual_container:
		if 'STEP2' not in each_summary:
			continue
		else:
			exact_summary = each_summary

	f = open(exact_summary, 'rU')
	EE = 0.0
	avgP = 0.0
	avgQ = 0.0
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'seqname':
			continue
		if ';ee=' not in line[0]:
			error_rate_string += line[0] + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\n'
		else:
			ee_value = line[0].split(';ee=')[1]
			ee_value = ee_value[0:-1]
			if ee_value == '0.0':
				error_rate_string += line[0] + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\n'
				continue
			EE, avgP, avgQ = expected_error_rate_calc(ee_value, line[3])
			error_rate_string += line[0] + '\t' + EE + '\t' + avgP + '\t' + avgQ + '\n'
	f.close()
	f = open(phred_error_rate, 'w')
	f.write(error_rate_string)
	f.close()
	return True


def global_trimming_fastq(fastq_file, fasta_file, fastq_trimmed_file, fasta_trimmed_file, maxlength, name, usearch_exec_path, mothur_exec_path, processors, outputdir):
	if maxlength == 'auto':
		top_length = calculate_fasta_median_length(fasta_file, name, mothur_exec_path, processors, outputdir)
	elif maxlength is not False and maxlength != 'auto':
		top_length = maxlength
	elif maxlength is False:
		copy_file(fasta_file, fasta_trimmed_file)
		copy_file(fastq_file, fastq_trimmed_file)
		return True
	# ##########################
	flag, stderr = execute_functions(usearch_fastq_truncate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, fastq_trimmed_file, fasta_trimmed_file, top_length)
	if flag is False:
		print "Execution of usearch_fastx_filter failed!!!"
	return True


def global_trimming_fasta(fasta_file, fasta_trimmed_file, median_length, maxlength, name, usearch_exec_path, mothur_exec_path, processors, outputdir):
	if maxlength == 'auto':
		top_length = median_length
	elif maxlength is not False and maxlength != 'auto':
		top_length = maxlength
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "GLOBAL TRIMMING LENGTH:", top_length
		report("GLOBAL TRIMMING LENGTH:" + top_length)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	elif maxlength is False:
		copy_file(fasta_file, fasta_trimmed_file)
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "GLOBAL TRIMMING LENGTH: SKIPPED!!!"
		report("GLOBAL TRIMMING LENGTH: SKIPPED!!!")
		report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		return True
	# ##########################
	flag, stderr = execute_functions(usearch_fasta_truncate, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fasta_file, fasta_trimmed_file, top_length)
	if flag is False:
		print "Execution of usearch_fasta_filter failed!!!"
	return True


def calculate_fasta_median_length(fasta_file, absname, mothur_exec_path, processors, outputdir):
	summary_file = outputdir + absname + '_length_distributions_STEP5.txt'
	flag, stderr = execute_functions(mothur_summary_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file)
	if flag is False:
		print "Execution of mothur_summary_seqs failed!!!"
	mothur_output_container = []
	extension_list = ['.summary']
	flag = scandirs(outputdir, mothur_output_container, extension_list)
	if flag is False:
		print "This extension is not availble: ", extension_list
		sys.exit(2)
	else:
		for each_file in mothur_output_container:
			if 'current_files.summary' in each_file:
				continue
			os.rename(each_file, summary_file)
	the_list = []
	f = open(summary_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'seqname':
			continue
		the_list.append(int(line[3]))
	f.close()
	freq_dict = {}
	freq_dict = collections.Counter(the_list)
	freq_list = []
	freq_list = freq_dict.most_common(10)
	top_length, top_freq = grab_most_frequent_length(freq_list)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	freq_string = 'LENGTH FREQUENCY REPORT:\n'
	for i in range(len(freq_list)):
		c = i + 1
		freq_string += 'RANK_' + str(c) + ': ' + str(freq_list[i][0]) + '= ' + str(freq_list[i][1]) + '\n'
	print freq_string
	report(freq_string)
	print "MEDIAN LENGTH:", top_length, " With the frequency of ", top_freq
	report("MEDIAN LENGTH:" + top_length + " With the frequency of " + top_freq)
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	return top_length


def fasta_filtering_process(fasta_file, fasta_filtered_file, total_discarded_fasta, trunc_len, min_len, maxambig, maxhomop, usearch_exec_path, mothur_exec_path, absname, processors, outputdir):
	list_to_intersect = []
	list_to_discard = []
	#path, absname, ext = split_file_name(fasta_file)
	raw_read_count = fasta_read_count(fasta_file)
	report("1: MAXAMBIG FILTERING")
	print "1: MAXAMBIG FILTERING"
	if maxambig is not False:
		fasta_maxambig_filtered = outputdir + absname + '_maxambig_filtered_STEP3.fasta'
		fasta_maxambig_discarded = outputdir + absname + '_maxambig_discarded_STEP3.fasta'
		fasta_maxambig_discarded_temp = outputdir + absname + '_maxambig_discarded_STEP3.txt'
		fasta_maxambig_filtered_header = outputdir + absname + '_maxambig_filtered_header_STEP3.txt'
		fasta_maxambig_discarded_header = outputdir + absname + '_maxambig_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_maxambig = maxambig
		cluster_maxhomop = False
		flag, stderr = execute_functions(mothur_screen_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file, cluster_maxambig, cluster_maxhomop, cluster_min_len)
		if flag is False:
			print "Execution of mothur_screen_seqs failed!!!"
		mothur_output_container = []
		extension_list = ['.good.fasta']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", list_to_string(extension_list, ',')
			report("This extension is not availble: " + list_to_string(extension_list, ','))
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxambig_filtered)
		mothur_output_container = []
		extension_list = ['.bad.accnos']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxambig_discarded_temp)
		maxambig_read_count = fasta_read_count(fasta_maxambig_filtered)
		flag = check_and_verify(raw_read_count, maxambig_read_count, 5)
		if flag is False:
			print "maxambig trimming filter is too stringent, will disable its effect"
			report("maxambig trimming filter is too stringent, will disable its effect")
			report("MAXAMBIG FILTERING: " + FAILED_MARK)
			print "MAXAMBIG FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_maxambig_filtered, fasta_maxambig_filtered_header)
			header_list = []
			flag = file_to_list(fasta_maxambig_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_maxambig_filtered_header)
			
			grab_nth_column(fasta_maxambig_discarded_temp, 0, fasta_maxambig_discarded_header, '\t')
			check_it_and_remove_it(fasta_maxambig_discarded_temp)
			discard_list = []
			flag = file_to_list(fasta_maxambig_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			
			flag, stderr = execute_functions(mothur_get_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file, fasta_maxambig_discarded_header)
			if flag is False:
				print "Execution of mothur_screen_seqs failed!!!"
			mothur_output_container = []
			extension_list = ['.pick.fasta']
			flag = scandirs(outputdir, mothur_output_container, extension_list)
			if flag is False:
				print "This extension is not availble: ", extension_list
				sys.exit(2)
			else:
				#it has the file name in it
				os.rename(mothur_output_container[0], fasta_maxambig_discarded)
				check_it_and_remove_it(fasta_maxambig_discarded_header)

			print "NUMBER OF READS AFTER MAXAMBIG FILTERING: ", maxambig_read_count
			report("NUMBER OF READS AFTER MAXAMBIG FILTERING: " + maxambig_read_count)
			report("HOMOPOLYMER FILTERING: " + CHECK_MARK)
			print "HOMOPOLYMER FILTERING: ", CHECK_MARK
	else:
		report("HOMOPOLYMER FILTERING: " + FAILED_MARK)
		print "HOMOPOLYMER FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	report("2: HOMOPOLYMER FILTERING")
	print "2: HOMOPOLYMER FILTERING"
	if maxhomop is not False:
		fasta_maxhomop_filtered = outputdir + absname + '_maxhomop_filtered_STEP3.fasta'
		fasta_maxhomop_discarded = outputdir + absname + '_maxhomop_discarded_STEP3.fasta'
		fasta_maxhomop_discarded_temp = outputdir + absname + '_maxhomop_discarded_STEP3.txt'
		fasta_maxhomop_filtered_header = outputdir + absname + '_maxhomop_filtered_header_STEP3.txt'
		fasta_maxhomop_discarded_header = outputdir + absname + '_maxhomop_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_maxambig = False
		cluster_maxhomop = maxhomop
		flag, stderr = execute_functions(mothur_screen_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file, cluster_maxambig, cluster_maxhomop, cluster_min_len)
		if flag is False:
			print "Execution of mothur_screen_seqs failed!!!"
		mothur_output_container = []
		extension_list = ['.good.fasta']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", list_to_string(extension_list, ',')
			report("This extension is not availble: " + list_to_string(extension_list, ','))
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxhomop_filtered)
		mothur_output_container = []
		extension_list = ['.bad.accnos']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxhomop_discarded_temp)
		homopolymer_read_count = fasta_read_count(fasta_maxhomop_filtered)
		flag = check_and_verify(raw_read_count, homopolymer_read_count, 5)
		if flag is False:
			print "homopolymer trimming filter is too stringent, will disable its effect"
			report("homopolymer trimming filter is too stringent, will disable its effect")
			report("HOMOPOLYMER FILTERING: " + FAILED_MARK)
			print "HOMOPOLYMER FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_maxhomop_filtered, fasta_maxhomop_filtered_header)
			header_list = []
			flag = file_to_list(fasta_maxhomop_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_maxhomop_filtered_header)
			
			grab_nth_column(fasta_maxhomop_discarded_temp, 0, fasta_maxhomop_discarded_header, '\t')
			check_it_and_remove_it(fasta_maxhomop_discarded_temp)
			discard_list = []
			flag = file_to_list(fasta_maxhomop_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			
			flag, stderr = execute_functions(mothur_get_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file, fasta_maxhomop_discarded_header)
			if flag is False:
				print "Execution of mothur_screen_seqs failed!!!"
			mothur_output_container = []
			extension_list = ['.pick.fasta']
			flag = scandirs(outputdir, mothur_output_container, extension_list)
			if flag is False:
				print "This extension is not availble: ", extension_list
				sys.exit(2)
			else:
				#it has the file name in it
				os.rename(mothur_output_container[0], fasta_maxhomop_discarded)
				check_it_and_remove_it(fasta_maxhomop_discarded_header)

			print "NUMBER OF READS AFTER HOMOPOLYMER FILTERING: ", homopolymer_read_count
			report("NUMBER OF READS AFTER HOMOPOLYMER FILTERING: " + homopolymer_read_count)
			report("HOMOPOLYMER FILTERING: " + CHECK_MARK)
			print "HOMOPOLYMER FILTERING: ", CHECK_MARK
	else:
		report("HOMOPOLYMER FILTERING: " + FAILED_MARK)
		print "HOMOPOLYMER FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	if len(list_to_intersect) < 1:
		print "Filtering mode is disabled"
		copy_file(fasta_file, fasta_filtered_file)
	else:
		print "Find the intersection between filtering modes"
		report("Find the intersection between filtering modes")
		intersected_list = []
		intersected_list = list(set.intersection(*list_to_intersect))
		intersected_header = outputdir + 'filtered_intersected.txt'
		flag = list_to_file(intersected_list, intersected_header)
		flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, intersected_header, fasta_file, fasta_filtered_file)
		if flag is False:
			print "Execution of usearch_extract_fasta failed!!!"

		uniq_discarded_list = list(set(discard_list))
		uniq_discarded_list_file = outputdir + absname + '_total_discarded_while_filtering_header.txt'
		list_to_file(uniq_discarded_list, uniq_discarded_list_file)
		uniq_discarded_fasta = outputdir + absname + '_total_discarded_STEP3.fasta'
		flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, uniq_discarded_list_file, fasta_file, uniq_discarded_fasta)
		if flag is False:
			print "Execution of usearch_extract_fasta failed!!!"
		
		check_it_and_remove_it(uniq_discarded_list_file)
		flag = merge_files([uniq_discarded_fasta], total_discarded_fasta)
		if flag is not True:
			print "something went wrong on writing files during merging"
			error("something went wrong on writing files during merging")
			print "ABORTING!!!"
			sys.exit(2)

	read_count = fasta_read_count(fasta_filtered_file)
	print "The number of reads after filtering process is: ", read_count
	report("The number of reads after filtering process is: " + read_count)
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	
	return True


def fastq_filtering_process(fastq_file, fasta_file, fasta_filtered_file, total_discarded_fasta, trunc_len, min_len, maxambig, maxhomop, minqscore, ngs_qscore, expected_error_rate, usearch_exec_path, mothur_exec_path, absname, processors, outputdir):
	list_to_intersect = []
	list_to_discard = []
	#path, absname, ext = split_file_name(fastq_file)
	raw_read_count = fasta_read_count(fasta_file)
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("1: MAXAMBIG FILTERING")
	print "1: MAXAMBIG FILTERING"
	if maxambig is not False:
		fasta_maxambig_filtered = outputdir + absname + '_maxambig_filtered_STEP3.fasta'
		fasta_maxambig_discarded = outputdir + absname + '_maxambig_discarded_STEP3.fasta'
		fasta_maxambig_filtered_header = outputdir + absname + '_maxambig_filtered_header_STEP3.txt'
		fasta_maxambig_discarded_header = outputdir + absname + '_maxambig_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_trunc_len = False
		cluster_maxambig = maxambig
		cluster_expected_error_rate = False
		cluster_minqscore = False
		flag, stderr = execute_functions(usearch_fastq_filter, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, fasta_maxambig_filtered, fasta_maxambig_discarded, cluster_expected_error_rate, cluster_maxambig, cluster_minqscore, cluster_trunc_len, cluster_min_len, ngs_qscore)
		if flag is False:
			print "Execution of usearch_fastq_filter failed!!!"
		maxambig_read_count = fasta_read_count(fasta_maxambig_filtered)
		flag = check_and_verify(raw_read_count, maxambig_read_count, 5)
		if flag is False:
			print "MAXAMBIG PARAMETER IS TOO STRINGENT, SKIPPED!!!"
			report("MAXAMBIG PARAMETER IS TOO STRINGENT, SKIPPED!!!")
			report("MAXAMBIG FILTERING: " + FAILED_MARK)
			print "MAXAMBIG FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_maxambig_filtered, fasta_maxambig_filtered_header)
			header_list = []
			flag = file_to_list(fasta_maxambig_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_maxambig_filtered_header)
			
			discard_list = []
			flag = grab_fasta_header(fasta_maxambig_discarded, fasta_maxambig_discarded_header)
			flag = file_to_list(fasta_maxambig_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			check_it_and_remove_it(fasta_maxambig_discarded_header)

			print "NUMBER OF READS AFTER MAXAMBIG FILTERING: ", maxambig_read_count
			report("NUMBER OF READS AFTER MAXAMBIG FILTERING: " + maxambig_read_count)
			report("MAXAMBIG FILTERING: " + CHECK_MARK)
			print "MAXAMBIG FILTERING: ", CHECK_MARK
	else:
		report("MAXAMBIG FILTERING: " + FAILED_MARK)
		print "MAXAMBIG FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	report("2: EXPECTED ERROR RATE FILTERING")
	print "2: EXPECTED ERROR RATE FILTERING"
	if expected_error_rate is not False:
		fasta_expected_error_rate_filtered = outputdir + absname + '_expected_error_rate_filtered_STEP3.fasta'
		fasta_expected_error_rate_discarded = outputdir + absname + '_expected_error_rate_discarded_STEP3.fasta'
		fasta_expected_error_rate_filtered_header = outputdir + absname + '_expected_error_rate_filtered_header_STEP3.txt'
		fasta_expected_error_rate_discarded_header = outputdir + absname + '_expected_error_rate_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_trunc_len = False
		cluster_maxambig = False
		cluster_expected_error_rate = expected_error_rate
		cluster_minqscore = False
		flag, stderr = execute_functions(usearch_fastq_filter, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, fasta_expected_error_rate_filtered, fasta_expected_error_rate_discarded, cluster_expected_error_rate, cluster_maxambig, cluster_minqscore, cluster_trunc_len, cluster_min_len, ngs_qscore)
		if flag is False:
			print "Execution of usearch_fastq_filter failed!!!"
		
		error_rate_read_count = fasta_read_count(fasta_expected_error_rate_filtered)

		flag = check_and_verify(raw_read_count, error_rate_read_count, 5)
		if flag is False:
			print "expected error rate filter is too stringent, will disable its effect"
			report("expected error rate filter is too stringent, will disable its effect")
			report("EXPECTED ERROR RATE FILTERING: " + FAILED_MARK)
			print "EXPECTED ERROR RATE FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_expected_error_rate_filtered, fasta_expected_error_rate_filtered_header)
			header_list = []
			flag = file_to_list(fasta_expected_error_rate_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_expected_error_rate_filtered_header)
			
			discard_list = []
			flag = grab_fasta_header(fasta_expected_error_rate_discarded, fasta_expected_error_rate_discarded_header)
			flag = file_to_list(fasta_expected_error_rate_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			check_it_and_remove_it(fasta_expected_error_rate_discarded_header)

			print "NUMBER OF READS AFTER EXPECTED ERROR RATE FILTERING: ", error_rate_read_count
			report("NUMBER OF READS AFTER EXPECTED ERROR RATE FILTERING: " + error_rate_read_count)
			report("EXPECTED ERROR RATE FILTERING: " + CHECK_MARK)
			print "EXPECTED ERROR RATE FILTERING: ", CHECK_MARK
	else:
		report("EXPECTED ERROR RATE FILTERING: " + FAILED_MARK)
		print "EXPECTED ERROR RATE FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	report("3: MINQ SCORE FILTERING")
	print "3: MINQ SCORE FILTERING"
	if minqscore is not False:
		fasta_minqscore_filtered = outputdir + absname + '_minqscore_filtered_STEP3.fasta'
		fasta_minqscore_discarded = outputdir + absname + '_minqscore_discarded_STEP3.fasta'
		fasta_minqscore_filtered_header = outputdir + absname + '_minqscore_filtered_header_STEP3.txt'
		fasta_minqscore_discarded_header = outputdir + absname + '_minqscore_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_trunc_len = False
		cluster_maxambig = False
		cluster_expected_error_rate = False
		cluster_minqscore = minqscore
		flag, stderr = execute_functions(usearch_fastq_filter, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, fasta_minqscore_filtered, fasta_minqscore_discarded, cluster_expected_error_rate, cluster_maxambig, cluster_minqscore, cluster_trunc_len, cluster_min_len, ngs_qscore)
		if flag is False:
			print "Execution of usearch_fastq_filter failed!!!"
		
		minimum_quality_read_count = fasta_read_count(fasta_minqscore_filtered)
		flag = check_and_verify(raw_read_count, minimum_quality_read_count, 5)
		if flag is False:
			print "minqscore filter is too stringent, will disable its effect"
			report("minqscore filter is too stringent, will disable its effect")
			report("MINQ SCORE FILTERING: " + FAILED_MARK)
			print "MINQ SCORE FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_minqscore_filtered, fasta_minqscore_filtered_header)
			header_list = []
			flag = file_to_list(fasta_minqscore_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_minqscore_filtered_header)

			discard_list = []
			flag = grab_fasta_header(fasta_minqscore_discarded, fasta_minqscore_discarded_header)
			flag = file_to_list(fasta_minqscore_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			check_it_and_remove_it(fasta_minqscore_discarded_header)

			print "NUMBER OF READS AFTER MINIMUM QUALITY SCORE FILTERING: ", minimum_quality_read_count
			report("NUMBER OF READS AFTER MINIMUM QUALITY SCORE FILTERING: " + minimum_quality_read_count)
			report("MINQ SCORE FILTERING: " + CHECK_MARK)
			print "MINQ SCORE FILTERING: ", CHECK_MARK
	else:
		report("MINQ SCORE FILTERING: " + FAILED_MARK)
		print "MINQ SCORE FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	report("4: HOMOPOLYMER FILTERING")
	print "4: HOMOPOLYMER FILTERING"
	if maxhomop is not False:
		fasta_maxhomop_filtered = outputdir + absname + '_maxhomop_filtered_STEP3.fasta'
		fasta_maxhomop_discarded = outputdir + absname + '_maxhomop_discarded_STEP3.fasta'
		fasta_maxhomop_discarded_temp = outputdir + absname + '_maxhomop_discarded_STEP3.txt'
		fasta_maxhomop_filtered_header = outputdir + absname + '_maxhomop_filtered_header_STEP3.txt'
		fasta_maxhomop_discarded_header = outputdir + absname + '_maxhomop_discarded_header_STEP3.txt'
		cluster_min_len = False
		cluster_maxambig = False
		cluster_maxhomop = maxhomop
		flag, stderr = execute_functions(mothur_screen_seqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, fasta_file, cluster_maxambig, cluster_maxhomop, cluster_min_len)
		if flag is False:
			print "Execution of mothur_screen_seqs failed!!!"
		mothur_output_container = []
		extension_list = ['.good.fasta']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", list_to_string(extension_list, ',')
			report("This extension is not availble: " + list_to_string(extension_list, ','))
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxhomop_filtered)
		mothur_output_container = []
		extension_list = ['.bad.accnos']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], fasta_maxhomop_discarded_temp)
		homopolymer_read_count = fasta_read_count(fasta_maxhomop_filtered)
		flag = check_and_verify(raw_read_count, homopolymer_read_count, 5)
		if flag is False:
			print "homopolymer trimming filter is too stringent, will disable its effect"
			report("homopolymer trimming filter is too stringent, will disable its effect")
			report("HOMOPOLYMER FILTERING: " + FAILED_MARK)
			print "HOMOPOLYMER FILTERING: ", FAILED_MARK
		else:
			flag = grab_fasta_header(fasta_maxhomop_filtered, fasta_maxhomop_filtered_header)
			header_list = []
			flag = file_to_list(fasta_maxhomop_filtered_header, header_list)
			list_to_intersect.append(set(header_list))
			check_it_and_remove_it(fasta_maxhomop_filtered_header)
			
			grab_nth_column(fasta_maxhomop_discarded_temp, 0, fasta_maxhomop_discarded_header, '\t')
			check_it_and_remove_it(fasta_maxhomop_discarded_temp)
			discard_list = []
			flag = file_to_list(fasta_maxhomop_discarded_header, discard_list)
			list_to_discard.extend(discard_list)
			
			flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fasta_maxhomop_discarded_header, fastq_file, fasta_maxhomop_discarded)
			if flag is False:
				print "Execution of usearch_extract_fasta failed!!!"
			check_it_and_remove_it(fasta_maxhomop_discarded_header)

			print "NUMBER OF READS AFTER HOMOPOLYMER FILTERING: ", homopolymer_read_count
			report("NUMBER OF READS AFTER HOMOPOLYMER FILTERING: " + homopolymer_read_count)
			report("HOMOPOLYMER FILTERING: " + CHECK_MARK)
			print "HOMOPOLYMER FILTERING: ", CHECK_MARK
	else:
		report("HOMOPOLYMER FILTERING: " + FAILED_MARK)
		print "HOMOPOLYMER FILTERING: ", FAILED_MARK
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	# ###############################################################################
	if len(list_to_intersect) < 1:
		print "Filtering mode is disabled"
		report("Filtering mode is disabled")
		flag = execute_functions(fastq_to_fasta_converter, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, fasta_filtered_file, ngs_qscore)
		if flag is False:
			print "Execution of relable_fastq failed!!!"
	else:
		print "COLLECT AND CURATE ALL FINDINGS"
		report("COLLECT AND CURATE ALL FINDINGS")
		
		intersected_list = []
		intersected_list = list(set.intersection(*list_to_intersect))
		intersected_header = outputdir + 'filtered_intersected_STEP3.txt'
		flag = list_to_file(intersected_list, intersected_header)
		flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, intersected_header, fastq_file, fasta_filtered_file)
		if flag is False:
			print "Execution of usearch_extract_fasta failed!!!"

		uniq_discarded_list = list(set(discard_list))
		uniq_discarded_list_file = outputdir + absname + '_total_discarded_while_filtering_header.txt'
		list_to_file(uniq_discarded_list, uniq_discarded_list_file)
		uniq_discarded_fasta = outputdir + absname + '_total_discarded_STEP3.fasta'
		flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, uniq_discarded_list_file, fastq_file, uniq_discarded_fasta)
		if flag is False:
			print "Execution of usearch_extract_fasta failed!!!"

		check_it_and_remove_it(uniq_discarded_list_file)

		flag = merge_files([uniq_discarded_fasta], total_discarded_fasta)
		if flag is not True:
			print "something went wrong on writing files during merging"
			error("something went wrong on writing files during merging")
			print "ABORTING!!!"
			sys.exit(2)

	read_count = fasta_read_count(fasta_filtered_file)
	print "NUMBER OF READS AFTER DATA FILTERING: ", read_count
	report("NUMBER OF READS AFTER DATA FILTERING: " + read_count)
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	return True


def chimera_removal_process(fasta_file, chimera_reference, chimeric_free_fasta, chimeric_fasta, chimeric_report, usearch_exec_path, processors, outputdir):
	if isFileEmpty(fasta_file) is False:
		print "The ", fasta_file, " is empty!!!"
		error("The " + fasta_file + " is empty!!!")
		sys.exit(2)
	read_count = fasta_read_count(fasta_file)
	split_value = int(read_count) / 100000
	print "split value is :", split_value
	if split_value == 0:
		split_value = 1
	print "Splitting fasta file."
	report("Splitting fasta file")
	flag, stderr = execute_functions(usearch_fasta_split, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fasta_file, split_value)
	if flag is False:
		print "Execution of usearch_split_fasta failed!!!"
	print "Fasta file splitted."
	report("Fasta file splitted.")
	print "Removing chimeric reads.\n"
	report("Removing chimeric reads.")
	chimera_free_list = []
	chimeric_list = []
	chimeric_report_list = []
	for i in range(1, split_value + 2):
		split_fasta_file = outputdir + 'split' + str(i) + '.fastx'
		if isFileExist(split_fasta_file) is False:
			print "Split File is missing"
			break
		dechim_fasta = outputdir + 'split' + str(i) + '.good.fasta'
		check_it_and_remove_it(dechim_fasta)
	
		chim_fasta = outputdir + 'split' + str(i) + '.chimeric.fasta'
		check_it_and_remove_it(chim_fasta)
	
		chim_report = outputdir + 'split' + str(i) + '.chimera_report.txt'
		check_it_and_remove_it(chim_report)

		# remove chimera
		flag, stderr = execute_functions(usearch_chimera_removal, processors, outputdir, 'multi', 'usearch', usearch_exec_path, split_fasta_file, chimera_reference, chim_fasta, dechim_fasta, chim_report)
		if flag is False:
			print "Execution of usearch_uparse_cluster failed!!!"
		chimera_free_list.append(dechim_fasta)
		chimeric_list.append(chim_fasta)
		chimeric_report_list.append(chim_report)
		check_it_and_remove_it(split_fasta_file)
		#check_it_and_remove_it(chim_report)
	print "Chimeric reads removed."
	report("Chimeric reads removed.")
	flag = merge_files(chimera_free_list, chimeric_free_fasta)
	if flag is False:
		print "Error on merging files!!!"
	#removing remenant files
	for each_file in chimera_free_list:
		check_it_and_remove_it(each_file)
	flag = merge_files(chimeric_list, chimeric_fasta)
	if flag is False:
		print "Error on merging files!!!"
	for each_file in chimeric_list:
		check_it_and_remove_it(each_file)
	flag = merge_files(chimeric_report_list, chimeric_report)
	if flag is False:
		print "Error on merging files!!!"
	for each_file in chimeric_report_list:
		check_it_and_remove_it(each_file)
	return True


def data_prep(input_file_directory, filetype, filesystem, files_path, usearch_exec_path, mothur_exec_path, processors, outputdir):
	data_perp_container = []
	stderr = ''
	# ###################################################################################################################################### #
	if filetype == "fastq" and filesystem == "paired":
		flag, stderr = execute_functions(mothur_make_file, processors, outputdir, 'multi', 'mothur', mothur_exec_path, input_file_directory)
		if flag is False:
			print "Execution of mothur_make_file failed!!!"
		mothur_output_container = []
		extension_list = ['.paired.file']
		flag = scandirs(outputdir, mothur_output_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", list_to_string(extension_list, ',')
			error("This extension is not availble: " + list_to_string(extension_list, ','))
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			#it has the file name in it
			os.rename(mothur_output_container[0], files_path)
			temp_list = []
			file_to_list(files_path, temp_list)
			data_perp_container = temp_list
	# ###################################################################################################################################### #
	elif filetype == "fastq" and filesystem == "single":
		flag = scandirs(input_file_directory, data_perp_container, fastq_extensions)
		if flag is False:
			print "This extension is not availble: ", list_to_string(fastq_extensions, ',')
			error("This extension is not availble: " + list_to_string(fastq_extensions, ','))
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			flag = list_to_file(data_perp_container, files_path)
	# ###################################################################################################################################### #
	elif filetype == "fastq" and filesystem == "mixed":
		flag = scandirs(input_file_directory, data_perp_container, fastq_extensions)
		if flag is False:
			print "This extension is not availble: ", list_to_string(fastq_extensions, ',')
			error("This extension is not availble: " + list_to_string(fastq_extensions, ','))
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			flag = list_to_file(data_perp_container, files_path)
	# ###################################################################################################################################### #
	elif filetype == "fasta" and filesystem == "paired":
		print "Still in progress fork!!!!"
		print "your file type is fasta and sequencing mode is paired: Can not work with paired fasta files"
		error("your file type is fasta and sequencing mode is paired: Can not work with paired fasta files")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	# ###################################################################################################################################### #
	elif filetype == "fasta" and filesystem == 'single':
		flag = scandirs(input_file_directory, data_perp_container, fasta_extensions)
		if flag is False:
			print "This extension is not availble: ", fasta_extensions
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			flag = list_to_file(data_perp_container, files_path)
	# ###################################################################################################################################### #
	elif filetype == "fasta" and filesystem == 'mixed':
		flag = scandirs(input_file_directory, data_perp_container, fasta_extensions)
		if flag is False:
			print "This extension is not availble: ", fasta_extensions
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			flag = list_to_file(data_perp_container, files_path)
	# ###################################################################################################################################### #

	if len(data_perp_container) < 1:
		return False
	else:
		print len(data_perp_container), filesystem, filetype.upper(), "files found in input directory."
		report(str(len(data_perp_container)) + " " + filesystem + " " + filetype.upper() + " files found in input directory.")
	return True


def explode_fasta(fasta_file, separator, usearch_exec_path, processors, outputdir):
	header_file = outputdir + 'header.txt'
	flag = grab_fasta_header(fasta_file, header_file)
	fasta_handle = open(header_file, 'rU')
	group_dict = {}
	for line in fasta_handle:
		line = line.rstrip()
		if separator is False:
			fasta_header = line
		elif separator not in line:
			fasta_header = line
		else:
			fasta_header = line.split(separator)[0]
		fasta_header = slugify(fasta_header)
		if fasta_header not in group_dict:
			group_dict[fasta_header] = ''
			group_dict[fasta_header] = line + '\n'
		else:
			group_dict[fasta_header] += line + '\n'
	fasta_handle.close()
	for each_group in group_dict:
		group_string = ''
		group_string = group_dict[each_group]
		extract_label = outputdir + each_group + '_label.txt'
		f = open(extract_label, 'w')
		f.write(group_string)
		f.close()
		extracted_fasta = outputdir + each_group + '.fasta'
		flag, stderr = execute_functions(usearch_extract_fasta, processors, outputdir, 'multi', 'usearch', usearch_exec_path, extract_label, fasta_file, extracted_fasta)
		if flag is False:
			print "Execution of usearch_extract_fasta failed!!!"
	return True


def explode_fastq(fastq_file, separator, exploded_fastq_list, usearch_exec_path, processors, outputdir):
	# first we extract the fastq headers
	header_file = outputdir + 'headers.txt'
	flag, stderr = execute_functions(usearch_fastq_header_extract, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, header_file)
	if flag is False:
		print "Execution of usearch_fastq_header_extract failed!!!"
	fasta_handle = open(header_file, 'rU')
	sample_name_group_list = []
	for line in fasta_handle:
		line = line[1:]
		sample_name = line.split(separator)[0]
		if sample_name not in sample_name_group_list:
			sample_name_group_list.append(sample_name)

	for sample_name in sample_name_group_list:
		exploded_fastq = outputdir + sample_name + '.fastq'
		sample_target = sample_name + separator
		flag, stderr = execute_functions(usearch_extract_fastq_partial, processors, outputdir, 'multi', 'usearch', usearch_exec_path, sample_target, fastq_file, exploded_fastq)
		if flag is False:
			print "Execution of usearch_extract_fastq failed!!!"
		exploded_fastq_list.append(exploded_fastq)
	return True


def explode_fastq2(fastq_file, separator, usearch_exec_path, processors, outputdir):
	# first we extract the fastq headers
	print "EXPLODING MIXED FASTQ FILE."
	header_file = outputdir + 'header.txt'
	flag, stderr = execute_functions(usearch_fastq_header_extract, processors, outputdir, 'multi', 'usearch', usearch_exec_path, fastq_file, header_file)
	if flag is False:
		print "Execution of usearch_fastq_header_extract failed!!!"
	fasta_handle = open(header_file, 'rU')
	group_dict = {}
	for line in fasta_handle:
		if line[0] == '>':
			line = line.rstrip()
			line = line[1:]
			if separator is False:
				fasta_header = line
			elif separator not in line:
				fasta_header = line
			else:
				fasta_header = line.split(separator)[0]
			fasta_header = slugify(fasta_header)
			if fasta_header not in group_dict:
				group_dict[fasta_header] = ''
				group_dict[fasta_header] = line + '\n'
			else:
				group_dict[fasta_header] += line + '\n'
	fasta_handle.close()
	for each_group in group_dict:
		group_string = ''
		group_string = group_dict[each_group]
		extract_label = outputdir + each_group + '_label.txt'
		f = open(extract_label, 'w')
		f.write(group_string)
		f.close()
		extracted_fastq = outputdir + each_group + '.fastq'
		flag, stderr = execute_functions(usearch_extract_fastq, processors, outputdir, 'multi', 'usearch', usearch_exec_path, extract_label, fastq_file, extracted_fastq)
		if flag is False:
			print "Execution of usearch_extract_fastq failed!!!"
	return True


def execute(command, ** kwargs):
	assert isinstance(command, list), "Expected 'command' parameter to be a list containing the process/arguments to execute. Got %s of type %s instead" % (command, type(command))
	assert len(command) > 0, "Received empty list of parameters"
	retval = {
			"exitCode": -1,
			"stderr": u"",
			"stdout": u"",
			"execTime": datetime.timedelta(0),
			"command": None,
			"pid": None
		}
	retval["command"] = command
	log.info("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(command)))
	# print("::singleProcessExecuter > At %s, executing \"%s\"" % (datetime.datetime.now(), " ".join(parameter)))
	cwd = kwargs.get("cwd", os.getcwd())
	sheel = kwargs.get("shell", True)
	startDatetime = datetime.datetime.now()
	myPopy = subprocess.Popen(command, cwd=cwd, preexec_fn=os.seteuid(os.getuid()), shell=sheel, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	retval["pid"] = myPopy.pid
	log.debug("::singleProcessExecuter > Command \"%s\" got pid %s" % (" ".join(command), myPopy.pid))
	try:
		retval["stdout"], retval["stderr"] = myPopy.communicate()
		myPopy.wait()
	except OSError, osErr:
		log.debug("::singleProcessExecuter > Got %s %s in myPopy.communicate() when trying get output of command %s. It is probably a bug (more info: http://bugs.python.org/issue1731717)" % (osErr, type(osErr), command[0]))
	except Exception, e:
		log.warn("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s" % (type(e), e, " ".join(command)))
		log.debug("::singleProcessExecuter > Got %s %s when trying to get stdout/stderr outputs of %s. Showing traceback:\n%s" % (type(e), e, " ".join(command), traceback.format_exc()))
		raise
	retval["exitCode"] = myPopy.returncode
	retval["execTime"] = datetime.datetime.now() - startDatetime
	return retval


def execute_functions(function_name, processors, outputdir, thread_type, func_mode, *args):
	list_of_pid = []
	list_of_stderr = []
	list_of_stdout = []
	if thread_type == 'fork':
		threads = int(processors) + 1
		processors = '1'
	elif thread_type == 'multi':
		threads = 2
	for thread in range(1, threads):
		pid_file = 'run.pid' + str(thread)
		stderr_file = 'stderr' + str(thread)
		stdout_file = 'stdout' + str(thread)
		run_pid = outputdir + pid_file
		stderr = outputdir + stderr_file
		stdout = outputdir + stdout_file
		flag = function_name(processors, outputdir, stderr, stdout, run_pid, *args)
		if flag is False:
			sys.exit(2)
		list_of_pid.append(pid_file)
		list_of_stderr.append(stderr_file)
		list_of_stdout.append(stdout_file)
	flag, stderr = process_monitor(list_of_pid, list_of_stderr, list_of_stdout, outputdir, threads, func_mode)
	if flag is False:
		print "Process monitor failed."
	return (True, stderr)


def process_monitor(pid_list, stderr_list, stdout_list, outputdir, threads, mode):
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	finished_flag = False
	flag_list = {}
	for each_pid in pid_list:
		flag_list[each_pid] = False
	toolbar_width = threads
	sys.stdout.write("[%s]" % (" " * toolbar_width))
	sys.stdout.flush()
	sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['
	while finished_flag is False:
		for pid_file, stderr_file, stdout_file in itertools.izip(pid_list, stderr_list, stdout_list):
			f = open(outputdir + pid_file)
			the_pid = int(f.read().rstrip())
			if pid_exists(the_pid) is False:
				sys.stdout.write("OK")
				sys.stdout.flush()
				flag_list[pid_file] = True
				flag, stderr = validate_execution(stderr_file, stdout_file, outputdir, mode)
				if flag is False:
					sys.stdout.write(":(")
					report("[:()]")
					print "Error in result of this thread: ", str(the_pid)
					error("Error in result of this thread: " + str(the_pid))
					print "All generated threads killed."
					error("All generated threads killed.")
					kill_pid_list(pid_list, outputdir)
					#print stderr
					print "ABORTING!!!"
					report("ABORTING!!!")
					sys.exit(2)
			if False in flag_list.values():
				finished_flag = False
			else:
				finished_flag = True
		time.sleep(1)
	sys.stdout.write("\n")
	report("[OK]")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	return (True, stderr)


def pid_exists(pid):
	if pid < 0:
		return False
	if pid == 0:
		# According to "man 2 kill" PID 0 refers to every process
		# in the process group of the calling process.
		# On certain systems 0 is a valid PID but we have no way
		# to know that in a portable fashion.
		raise ValueError('invalid PID 0')
	try:
		os.kill(pid, 0)
	except OSError as err:
		if err.errno == errno.ESRCH:
			# ESRCH == No such process
			return False
		elif err.errno == errno.EPERM:
			# EPERM clearly means there's a process to deny access to
			return True
		else:
			# According to "man 2 kill" possible error values are
			# (EINVAL, EPERM, ESRCH)
			raise
	else:
		return True


def validate_execution(stderr, stdout, outputdir, mode):
	if mode == 'usearch':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		#print data
		if '---Fatal error---' in data or 'Invalid command line' in data:
			return (False, data)
		else:
			return(True, data)
	elif mode == 'mothur':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		if '[ERROR]:' in data:
			return (False, data)
		else:
			return(True, data)
	elif mode == 'mafft':
		f = open(outputdir + stdout, 'rU')
		data = f.read()
		if 'No such file or directory' in data or 'cannot open file' in data:
			return (False, data)
		else:
			return(True, data)


def kill_pid_list(pid_list, outputdir):
	for pid in pid_list:
		f = open(outputdir + pid)
		the_pid = int(f.read().rstrip())
		try:
			os.kill(the_pid, signal.SIGTERM)
		except OSError:
			pass
		print the_pid, "Killed!"
	return True


def set_qscore(ngs):
	if ngs in ['illumina']:
		ngs_qscore = '-fastq_ascii 33 -fastq_qmin 0 -fastq_qmax 41'
		SetIllumina()
	elif ngs in ['sanger']:
		ngs_qscore = '-fastq_ascii 64 -fastq_qmin 0 -fastq_qmax 41'
		SetSanger()
	elif ngs in ['pacbio']:
		ngs_qscore = '-fastq_ascii 33 -fastq_qmin 0 -fastq_qmax 42'
		SetPacBio()
	elif ngs in ['iontorrent']:
		ngs_qscore = '-fastq_ascii 33 -fastq_qmin 0 -fastq_qmax 45'
		SetIontorrent()

	return ngs_qscore


def SetSanger():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 64
	Q_Min = 0
	Q_Max = 41


def SetIllumina():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 33
	Q_Min = 0
	Q_Max = 41


def SetPacBio():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 33
	Q_Min = 0
	Q_Max = 42


def SetIontorrent():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 33
	Q_Min = 0
	Q_Max = 45

# ################################# UTITLITY FUNCTIONS ####################### #


def grab_nth_column(text_file, column_number, new_name, delimiter=None):
	string = ''
	if delimiter is None:
		delimiter = '\t'
	f = open(text_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split(delimiter)
		string += line[column_number] + '\n'
	f.close()
	o = open(new_name, 'w')
	o.write(string)
	o.close()
	return True


def relable_fasta(fasta_file, tag, suffix, new_name):
	
	if tag is None:
		copy_file(fasta_file, new_name)
	else:
		f = open(fasta_file, 'rU')
		string = ''
		c = 1
		for i in f:
			if i[0] == '>':
				string += '>' + tag + str(c) + suffix + '\n'
				c += 1
			else:
				string += i
		hdl = open(new_name, 'w')
		hdl.write(string)
		hdl.close()
		f.close()
		string = ''
	return True


def relable_single_fasta_reads(fasta_file, sample_name, new_name):
	f = open(fasta_file, 'rU')
	string = ''
	c = 1
	for i in f:
		if i[0] == '>':
			string += '>' + sample_name + ';S_' + str(c) + ';ee=0.0;' + '\n'
			c += 1
		else:
			string += i
	hdl = open(new_name, 'w')
	hdl.write(string)
	hdl.close()
	f.close()
	string = ''
	return True


def relable_mixed_fasta_reads(fasta_file, separator, new_name):
	f = open(fasta_file, 'rU')
	string = ''
	c = 1
	for i in f:
		if i[0] == '>':
			i = i.rstrip()
			i = i[1:]
			if separator is False:
				string += '>' + slugify(i) + ';S_' + str(c) + ';ee=0.0;' + '\n'
			else:
				i = i.split(separator)[0]
				string += '>' + slugify(i) + ';S_' + str(c) + ';ee=0.0;' + '\n'
			c += 1
		else:
			string += i
	hdl = open(new_name, 'w')
	hdl.write(string)
	hdl.close()
	f.close()
	string = ''
	return True


def expected_error_rate_calc(ee_value, length):
	ee_value = float(ee_value)
	length = float(length)
	avgP = ee_value / length
	avgQ = -10 * math.log10(avgP)
	return (str(round(ee_value, 2)), str(round(avgP, 4)), str(round(avgQ, 4)))


def grab_most_frequent_length(frequency_list):
	most_frequent_set = frequency_list[0]
	top_length = most_frequent_set[0]
	top_freq = most_frequent_set[1]
	top_limit = int(top_freq / 3)
	for i in range(len(frequency_list)):
		if top_length == frequency_list[i][0]:
			continue
		elif frequency_list[i][0] > top_length:
			continue
		elif frequency_list[i][0] < top_length and top_limit < frequency_list[i][1]:
			top_length = frequency_list[i][0]
			top_freq = frequency_list[i][1]
			top_limit = int(top_freq / 3)
	return (str(top_length), str(top_freq))


def copy_file(source, destination):
	shutil.copyfile(source, destination)
	return True


def check_and_verify(last_counter, next_counter, limit):
	# check two file and if size of next is less than 3% of last return False
	next_counter = int(next_counter)
	last_counter = int(last_counter)
	percentage = (next_counter * 100) / last_counter

	if percentage < limit:
		return False
	else:
		return True


def isFileEmpty(fname):
	# check if the path exist and have access
	if isFileExist(fname) is False:
		return False
	if os.path.getsize(fname) > 0:
		return True
	else:
		return False


def fasta_read_count(fasta_file):
	count = 0
	f = open(fasta_file)
	for i in f:
		if i[0] == '>':
			count += 1
	return str(count)


def fastq_read_count(fastq_file):
	print "Counting Reads!!!"
	counter = 0
	f = open(fastq_file, 'rU')
	for i in f:
		counter += 1
	f.close()
	size = counter / 4
	return str(size)


def merge_files(list_of_files, new_name):
	# open and merge list of files
	check_it_and_remove_it(new_name)
	merged_handle = open(new_name, 'a')
	for each_file in list_of_files:
		if isFileEmpty(each_file) is False:
			continue
		f = open(each_file, 'rU')
		merged_handle.write(f.read())
		f.close()
	merged_handle.close()
	return True


def split_file_name(file):
	path = os.path.dirname(file) + '/'
	name = os.path.basename(file)
	if '.' in name:
		ext = '.' + '.'.join(name.split('.')[1:])
		absname = name.split('.')[0]
	else:
		ext = 'no_extension'
		absname = name
	return (path, absname, ext)


def grab_fasta_header(fasta_file, header_file):
	string = ''
	f = open(fasta_file, 'rU')
	for i in f:
		if i[0] == '>':
			i = i.rstrip()
			string += i[1:] + '\n'
	f.close()
	write_string_down(string, header_file)
	string = ''
	return True


def write_string_down(new_string, new_name):
	f = open(new_name, 'w')
	f.write(new_string)
	f.close()
	return new_name


def slugify(filename):
	if '-' in filename:
		filename = filename.replace('-', '_')
	if '.' in filename:
		filename = filename.replace('.', '_')
	if ' ' in filename:
		filename = filename.replace(' ', '_')
	if ':' in filename:
		filename = filename.replace(' ', '_')
	return filename


def create_folder(path):
	# create folder in specified path
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(query_list)


def list_to_file(query_list, file_name, mode=None):
	# convert list to file by \n delimiter
	if mode is None:
		mode = 'w'
	file_hdl = open(file_name, mode)
	file_hdl.write(list_to_string(query_list, "\n"))
	file_hdl.close()
	return True


def file_to_list(query_file, list_name):
	# convert file to list
	file_hdl = open(query_file, 'rU')
	for i in file_hdl:
		i = i.rstrip()
		list_name.append(i)
	file_hdl.close()
	return True


def scandirs(path, container, ext_list, mode=None):
	# scan a spath and grab all files by specified extension
	for root, dirs, names in os.walk(path):
		for currentFile in names:
			path, absname, ext = split_file_name(currentFile)
			if mode is None:
				if ext in ext_list:
					container.append(os.path.join(root, currentFile))
			elif mode == 'partial':
				for each_ext in ext_list:
					if ext in each_ext:
						container.append(os.path.join(root, currentFile))
			elif mode == 'ex_partial':
				for each_ext in ext_list:
					if each_ext in ext:
						container.append(os.path.join(root, currentFile))
	if len(container) < 1:
		return False
	return True


def isfloat(x):
	try:
		float(x)
	except ValueError:
		return False
	else:
		return True


def isint(x):
	try:
		int(x)
	except ValueError:
		return False
	else:
		return True


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def error(report_string):
	f = open(report_file, "a")
	f.write('###############################  ERROR   ###################################\n')
	f.write(report_string)
	f.write("\n")
	f.write('############################################################################\n')
	f.close()


def report(report_string):
	f = open(report_file, "a")
	f.write(report_string.encode('utf8'))
	f.write("\n")
	f.close()


def check_it_and_remove_it(filename, noreport=False):
	try:
		os.remove(filename)
		if noreport is False:
			pass
	except OSError:
		pass


def isPathExist(path):
	# check if the path exist and have access
	if os.path.exists(path) and os.access(path, os.R_OK):
		return True
	else:
		return False
# ################################# END OF SCRIPT ############################ #


if __name__ == "__main__": main(sys.argv[1:])
