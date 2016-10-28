#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
16S ALPHA DIVERSITY ANALYSIS
"""
# ################################### IMPORTS ################################## #
import sys
import argparse
import time
import os
import multiprocessing
import datetime
import logging as log
import subprocess
import traceback
import time
import itertools
import errno
import shutil
import csv
import random
import decimal
import signal
import colorsys


# ################################### OBJECTS ################################## #
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
		print command
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


class entity:
	def __init__(self, info_hash):
		for key, value in info_hash.items():
			setattr(self, key, value)
	
	def response(self, key):
		return getattr(self, key)

	def get_taxon(self):
		return self.response('taxon')
	
	def get_taxlevel(self):
		return self.response('taxlevel')
	
	def get_total(self):
		return self.response('total')

	def get_rankid(self):
		return self.response('rankID')
		
	def get_values(self):
		return self.__dict__.values()
	
	def get_keys(self):
		return self.__dict__.keys()
	
	def get_abundance_values(self):
		abundance_values_list = []
		all_keys = self.get_keys()
		ignore_list = ['taxon', 'taxlevel', 'rankID', 'daughterlevels', 'total']
		for each_key in all_keys:
			if each_key not in ignore_list:
				abundance_values_list.append(self.response(each_key))
		return abundance_values_list

	def get_abundance_keys(self):
		abundance_keys_list = []
		all_keys = self.get_keys()
		ignore_list = ['taxon', 'taxlevel', 'rankID', 'daughterlevels', 'total']
		for each_key in all_keys:
			if each_key not in ignore_list:
				abundance_keys_list.append(each_key)
		return abundance_keys_list


# ################################### GLOBALS ################################## #
CHECK_MARK = "OK"
FAILED_MARK = ":("
outputdir = "/16S_data_analyser_output/"
report_file = outputdir + "16S_data_analyser_report.txt"
jslib_path = "/javascript/"
name = "microbiome_trial"
# ###################################   MAIN   ################################# #


def main(argv):
	report_string = ''

	# ############################# BEURACRATICS PROCEDURES
	report_string += "######################################################################################################################################\n"
	print "######################################################################################################################################"
	report_string += "16S SIMPLE PLOT EXECUTION HAS INITIATED" + '\n'
	print "16S SIMPLE PLOT EXECUTION HAS INITIATED"
	report_string += "Initiation time: " + time.strftime("%Y-%m-%d %H:%M") + '\n'
	print "Initiation time: ", time.strftime("%Y-%m-%d %H:%M")
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
	
	report_string += "16S_ALPHA DIVERSITY EXECUTION HAS INITIATED"
	report_string += time.strftime("%Y-%m-%d %H:%M")
	print "16S_ALPHA DIVERSITY EXECUTION HAS INITIATED"
	print time.strftime("%Y-%m-%d %H:%M")
	if isPathExist(outputdir) is False:
		print "[--outputdir]: Output file directory has Access/Exist issue!!!"
		report_string += "[--outputdir]: Output file directory has Access/Exist issue!!!" + '\n'
		print "ABORTING!!!"
		sys.exit(2)
	else:
		print "[--outputdir]: Output file directory is: ", outputdir
		report_string += "[--outputdir]: Output file directory is: " + outputdir
	check_it_and_remove_it(report_file, True)
	report(report_string)
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VERSION OF JAVASCRIPT LIBRARY IN JSLIB DIRECTORY"
	report("VERIFYING THE SANITY/VERSION OF JAVASCRIPT LIBRARY IN JSLIB DIRECTORY")
	print "###################################################################\n"
	report("###################################################################\n")
	
	jslib_list = ['jquery-2.2.3.min.js', 'jquery.tablesorter.js', 'plotly-latest.min.js', 'bootstrap.min.css', 'bootstrap.min.js', 'tether.min.js', 'fontawesome.min.js', 'fontawesome.min.css']
	#First Constant
	
	flag = test_javascript_library(jslib_path, jslib_list)
	if flag is False:
		report("JAVASCRIPT LIBRARY: " + FAILED_MARK)
		print "JAVASCRIPT LIBRARY: ", FAILED_MARK
		print "ABORTING!!!"
		sys.exit(2)
	
	js_path, alpha_path, data_path = visual_directory(outputdir, name, jslib_path, jslib_list)
	print "Done!!"
	report("Done!!")

# ################################### GENERAL FUNCTIONS ######################## #


def visual_directory(outputdir, name, jslib_directory, jslib_list):
	#design_outputdir, absname, ext = split_file_name(design_file)
	visual_path = outputdir + name + "_ANALYSIS/"
	flag = create_folder(visual_path)
	if flag is False:
		print "Can not create folder at specified path: ", visual_path
		print "ABORTING!!!"
		sys.exit(2)
	js_path = outputdir + name + "_ANALYSIS/" + name + "_JS/"
	flag = create_folder(js_path)
	if flag is False:
		print "Can not create folder at specified path: ", js_path
		print "ABORTING!!!"
		sys.exit(2)
	for js in jslib_list:
		shutil.copy2(jslib_directory + js, js_path + js)
	data_path = outputdir + name + "_ANALYSIS/" + name + "_DATA/"
	flag = create_folder(data_path)
	if flag is False:
		print "Can not create folder at specified path: ", data_path
		print "ABORTING!!!"
		sys.exit(2)
	return (js_path, visual_path, data_path)


def test_javascript_library(jslib_directory, jslib_list):
	for i in jslib_list:
		if isFileExist(jslib_directory + i):
			report("THE " + i + " LIBRARY IS FOUND")
			print "THE ", i, " LIBRARY IS FOUND"
		else:
			error("THIS JAVASCRIPT FILE IS MISSING:" + i)
			print "THIS JAVASCRIPT FILE IS MISSING:", i
			return False
	return True


# ################################### UTILITIES  ############################### #


def create_folder(path):
	# create folder in specified path
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False


def check_it_and_remove_it(filename, noreport=False):
	try:
		os.remove(filename)
		if noreport is False:
			pass
	except OSError:
		pass


def error(report_string):
	f = open(report_file, "a")
	f.write('###############################  ERROR   ###################################\n')
	f.write(report_string)
	f.write("\n")
	f.write('############################################################################\n')
	f.close()


def isPathExist(path):
	# check if the path exist and have access
	if os.path.exists(path) and os.access(path, os.R_OK):
		return True
	else:
		return False


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def report(report_string):
	f = open(report_file, "a")
	f.write(report_string)
	f.write("\n")
	f.close()
# ################################### END SCRIPT ############################### #
if __name__ == "__main__": main(sys.argv[1:])