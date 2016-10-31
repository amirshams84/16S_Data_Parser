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
execdir = "/exec/"
jslib = "/javascript/"
outputdir = "/16S_data_analyser_output/"
phylotype = "/test_data/ZAC_phylotype.txt"
design = "/test_data/ZAC_design.txt"
processors = "1"
name = "ZAC"
taxlevel = "genus"
remove_sample_file = "false"
keep_sample_file = "false"
normalize = "false"
# ###################################   MAIN   ################################# #


def main(argv):
	report_string = ''
	# ############################# PARSE INPUT ARGUMENTS
	parser = argparse.ArgumentParser()

	main_file = parser.add_argument_group('Main file parameters')
	main_file.add_argument("--phylotype", help="Phylotype file", action='store')
	main_file.add_argument("--taxlevel", help="taxonomy level on of [species, genus, family", action='store')
	main_file.add_argument("--design", help="design file", action='store')
	main_file.add_argument("--jslib", help="javascript library", action='store')
	main_file.add_argument("--outputdir", help="output directory", action='store')
	main_file.add_argument("--execdir", help="executables_directory", action='store')

	general = parser.add_argument_group('general parameters')
	general.add_argument("--processors", help="number of processors assigned", action='store')
	general.add_argument("--name", help="name of output files", action='store')

	beta = parser.add_argument_group('beta parameters')
	beta.add_argument("--remove_sample_file", help="list of samples to remove", action='store')
	beta.add_argument("--keep_sample_file", help="list of samples to keep", action='store')
	beta.add_argument("--normalize", help="Normalizing using totalgroup method", action='store')
	
	args = parser.parse_args()
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
	
	# ########### OUTPUT DIRECTORY
	args.outputdir = outputdir
	args.jslib = jslib
	args.execdir = execdir
	args.outputdir = outputdir
	args.phylotype = phylotype
	args.design = design
	args.processors = processors
	args.name = name
	args.remove_sample_file = remove_sample_file
	args.keep_sample_file = keep_sample_file
	args.taxlevel = taxlevel
	args.normalize = normalize
	# ########################

	if isPathExist(args.outputdir) is False:
		print "OUTPUTDIR IS NOT THERE"

	global report_file
	report_file = args.outputdir + "16S_simple_plot_report.txt"
	check_it_and_remove_it(report_file, True)
	report(report_string)

	# ########### SHARED FILE
	if isFileExist(args.phylotype) is False:
		error("[--phylotype]: phylotype file has Access/Exist issue")
		print "[--phylotype]: phylotypefile has Access/Exist issue"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		report("[--phylotype]: phylotype file is: " + args.phylotype)
		print "[--phylotype]: phylotype file is: ", args.phylotype

	# ########### DESIGN FILE
	if isFileExist(args.design) is False:
		print "[--design]: Design file has Access/Exist issue"
		error("[--design]: Design file has Access/Exist issue")
		print "ABORTING!!!"
		sys.exit(2)
	else:
		report("[--design]: Design file is: " + args.design)
		print "[--design]: Design file is: ", args.design

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

	# ########### TAX LEVEL CHECKING
	if args.taxlevel and args.taxlevel.strip():
		report("[--taxlevel]: Specified taxlevel for analysis is: " + args.taxlevel.lower())
		print "[--taxlevel]: Specified taxlevel for output files is: ", args.taxlevel.lower()
	else:
		error("[--taxlevel]: taxlevel can not be empty!!!")
		print "[--taxlevel]: taxlevel can not be empty!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	if args.taxlevel.lower() not in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
		error("[--taxlevel]: taxlevel value should be on of these ['phylum', 'class', 'order', 'family', 'genus', 'species']")
		print "[--taxlevel]: taxlevel value should be on of these ['phylum', 'class', 'order', 'family', 'genus', 'species']"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		args.taxlevel = args.taxlevel.lower()
	
	# ########### REMOVE SAMPLE FILE CHECKING
	if args.remove_sample_file.lower() == 'false':
		report("[--remove_sample_file]: Remove sample file is disbaled")
		print "[--remove_sample_file]: Remove sample file is disbaled"
		args.remove_sample_file = False
	elif isFileExist(args.remove_sample_file) is False:
		print "[--remove_sample_file]: Remove sample file has Access/Exist issue"
		error("[--remove_sample_file]: Remove sample file has Access/Exist issue")
		print "ABORTING!!!"
		sys.exit(2)
	else:
		report("Remove sample file(--remove_sample_file) is: " + args.remove_sample_file)
		print "Remove sample file(--remove_sample_file) is: ", args.remove_sample_file

	# ########### KEEP SAMPLE FILE CHECKING
	if args.keep_sample_file.lower() == 'false':
		report("[--keep_sample_file]: keep sample file is disbaled")
		print "[--keep_sample_file]: keep sample file is disbaled"
		args.keep_sample_file = False
	elif isFileExist(args.keep_sample_file) is False:
		print "[--keep_sample_file]: keep sample file has Access/Exist issue"
		error("[--keep_sample_file]: keep sample file has Access/Exist issue")
		print "ABORTING!!!"
		sys.exit(2)
	else:
		report("keep sample file(--keep_sample_file) is: " + args.keep_sample_file)
		print "keep sample file(--keep_sample_file) is: ", args.keep_sample_file

	# ########################
	if args.normalize.lower() == 'false':
		print "[--normalize]: Normalizing is disabled"
		args.normalize = False
	elif args.normalize.lower() == 'true':
		report("[--normalize]: Normalizing is enabled")
		print "[--normalize]: Normalizing is enabled"
		args.normalize = True
	else:
		report("[--normalize]: Normalizing should be true or false")
		print "[--normalize]: Normalizing should be true or false"
		print "Aborting!!!"
		sys.exit(2)
	
	# ########### EXECECUTIVE DIRECTORY CHECKING
	if args.execdir[-1] != '/':
		args.execdir += '/'
	if isPathExist(args.execdir) is False:
		error("[--execdir]: executables directory has Access/Exist issue!!!")
		print "[--execdir]: executables directory has Access/Exist issue!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
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
	# ########### JAVSCRIPT DIRECTORY CHECKING
	if args.jslib[-1] != '/':
		args.jslib += '/'
	if isPathExist(args.jslib) is False:
		error("[--jslib]: executables directory has Access/Exist issue!!!")
		print "[--jslib]: executables directory has Access/Exist issue!!!"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)

	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VERSION OF JAVASCRIPT LIBRARY IN JSLIB DIRECTORY"
	report("VERIFYING THE SANITY/VERSION OF JAVASCRIPT LIBRARY IN JSLIB DIRECTORY")
	print "###################################################################\n"
	report("###################################################################\n")

	jslib_list = ['jquery-2.2.3.min.js', 'jquery.tablesorter.js', 'plotly-latest.min.js', 'bootstrap.min.css', 'bootstrap.min.js', 'tether.min.js', 'fontawesome.min.js', 'font-awesome.min.css']
	flag = test_javascript_library(args.jslib, jslib_list)
	if flag is False:
		report("JAVASCRIPT LIBRARY: " + FAILED_MARK)
		print "JAVASCRIPT LIBRARY: ", FAILED_MARK
		print "ABORTING!!!"
		sys.exit(2)
	js_path, alpha_path, data_path = visual_directory(args.outputdir, args.name, args.jslib, jslib_list)

	entity_list = []
	entity_list = parse_abundance_file(args.phylotype)
	dict_of_entities = {}
	dict_of_entities = construct_object(entity, entity_list)
	shared_file_generated = create_shared_table(dict_of_entities, args.phylotype, args.taxlevel, args.name, args.outputdir)
	args.shared = shared_file_generated
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step1: Fix shared file label and future modification
	#for now we change the label from the distances to the specified args.name
	fixed_shared = args.outputdir + args.name + '_FIXED_shared_file_STEP1.txt'
	fixed_design = args.outputdir + args.name + '_FIXED_design_file_STEP1.txt'
	flag = fix_shared_design_file(args.shared, args.design, args.name, fixed_shared, fixed_design)
	if flag is True:
		print "Shared file label is fixed."
		print "# ##########################"
		print "Shared file replaced with Step1 Shared file(FIXED_shared_file_STEP1.txt)"
		args.shared = fixed_shared
		print "Design file replaced with Step1 Design file(FIXED_design_file_STEP1.txt)"
		args.design = fixed_design
		print "# ##########################"
	else:
		print "Some thing is wrong with fix_shared_design_file"
		sys.exit(2)
	print "# ####################################################################################"
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step2: Create shared file only contains design file samples.
	args.shared = fixed_shared
	
	filter_control = args.outputdir + args.name + '_UPDATED_control_file_STEP2.txt'
	design_to_control(args.design, filter_control)
	updated_shared_file = args.outputdir + args.name + '_UPDATED_shared_file_STEP2.txt'
	updated_design_file = args.outputdir + args.name + '_UPDATED_design_file_STEP2.txt'
	flag = update_shared_design_file(args.shared, args.design, filter_control, updated_shared_file, updated_design_file, args.name, mothur_exec_path, args.processors, args.outputdir)
	if flag is True:
		print "Shared file initial update is successfull."
		print "# ##########################"
		print "Shared file replaced with Step2 Shared file(UPDATED_shared_file_STEP2.txt)"
		args.shared = updated_shared_file
		print "Design file replaced with Step2 Design file(UPDATED_shared_file_STEP2.txt)"
		args.design = updated_design_file
		print "# ##########################"
	else:
		print "Something is wrong with update_shared_file."
		sys.exit(2)
	print "# ####################################################################################"

	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step3: Keep samples in shared file/Design file based on keep_sample_file.
	if args.keep_sample_file is not False and isFileExist(args.keep_sample_file) is True:
		samples_kept_shared_file = args.outputdir + args.name + '_SAMPLES_KEPT_shared_file_STEP3.txt'
		samples_kept_design_file = args.outputdir + args.name + '_SAMPLES_KEPT_design_file_STEP3.txt'
		print "keep_sample_file is set: so I am going to keep samples only listed."
		flag = keep_sample_shared_design(args.shared, args.design, args.keep_sample_file, samples_kept_shared_file, samples_kept_design_file, args.name, mothur_exec_path, args.processors, args.outputdir)
		if flag is True:
			print "keep_sample_shared_design is successfull."
			print "# ##########################"
			print "Shared file replaced with Step3 Shared file(SAMPLES_KEPT_shared_file_STEP3.txt)"
			args.shared = samples_kept_shared_file
			print "Design file replaced with Step3 Design file(SAMPLES_KEPT_design_file_STEP3.txt)"
			args.design = samples_kept_design_file
			print "# ##########################"
		else:
			print "Something is wrong with keep_sample_shared_design."
			sys.exit(2)
	else:
		print "Since you did not specifiy keep_sample_file, the STEP3 will be skipped."
	print "# ####################################################################################"

	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step4: Remove samples in shared file/Design file based on remove_sample_file.
	if args.remove_sample_file is not False and isFileExist(args.remove_sample_file) is True:
		samples_removed_shared_file = args.outputdir + args.name + '_SAMPLES_REMOVED_shared_file_STEP4.txt'
		samples_removed_design_file = args.outputdir + args.name + '_SAMPLES_REMOVED_design_file_STEP4.txt'
		print "remove_sample_file is set: so I am going to remove samples listed."
		flag = remove_sample_shared_design(args.shared, args.design, args.remove_sample_file, samples_removed_shared_file, samples_removed_design_file, args.name, mothur_exec_path, args.processors, args.outputdir)
		if flag is True:
			print "remove_sample_shared_design is successfull."
			print "# ##########################"
			print "Shared file replaced with Step4 Shared file(SAMPLES_REMOVED_shared_file_STEP4.txt)"
			args.shared = samples_removed_shared_file
			print "Design file replaced with Step4 Design file(SAMPLES_REMOVED_shared_file_STEP4.txt)"
			args.design = samples_removed_design_file
			print "# ##########################"
		else:
			print "Something is wrong with remove_sample_shared_design."
			sys.exit(2)
	else:
		print "Since you did not specifiy remove_sample_file, the STEP4 will be skipped."
	print "# ####################################################################################"
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step7: Normalize shared file.
	if args.normalize is True:
		normalized_shared_file = args.outputdir + args.name + '_NORMALIZED_shared_file_STEP7.txt'
		normalized_design_file = args.outputdir + args.name + '_NORMALIZED_design_file_STEP7.txt'
		flag = normalize_shared_design_maker(args.shared, args.design, args.name, normalized_shared_file, normalized_design_file, mothur_exec_path, args.processors, args.outputdir)
		if flag is True:
			print "normalize_shared_maker is successfull."
			print "# ##########################"
			print "Shared file replaced with Step7 Shared file(NORMALIZED_shared_file_STEP7.txt)"
			args.shared = normalized_shared_file
			print "Design file replaced with Step6 Design file(NORMALIZED_design_file_STEP7.txt)"
			args.design = normalized_design_file
			print "# ##########################"
		else:
			print "Something is wrong with normalize_shared_design_maker."
			sys.exit(2)
	else:
		print "Since you did not specifiy Normalized the STEP7 will be skipped."
	print "# ####################################################################################"

	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step7: RAREFACTION CURVE PLOT.
	rarefaction_design_dict = {}
	plotly_design_data_dict = {}
	rarefaction_file = args.outputdir + args.name + '_RAREFACTION_file_STEP7.txt'
	rarefaction_design_dict, plotly_design_data_dict = mothur_rarefaction_curve(args.shared, args.design, rarefaction_file, mothur_exec_path, args.processors, data_path, args.outputdir)
	rarefaction_javascript_plotly(rarefaction_design_dict, plotly_design_data_dict, args.name, args.design, js_path)
	rarefaction_file_name = args.name + '_RAREFACTION_file_STEP7.txt'
	#rarefaction_plotter(alpha_path, rarefaction_file_name, args.design, args.name)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step8: Sample_ABUNDANCE_FILE.
	abundance_columns_dict = {}
	abundance_otu_dict = {}
	simple_abundance_file = args.outputdir + args.name + '_SAMPLE_ABUNDANCE_file_STEP8.txt'
	abundance_columns_dict, abundance_otu_dict = abundance_table_plotly(args.shared, args.design, simple_abundance_file, data_path, args.outputdir)
	abundance_javascript_plotly(abundance_columns_dict, abundance_otu_dict, args.name, args.design, js_path)
	sample_abundance_file_name = args.name + '_SAMPLE_ABUNDANCE_file_STEP8.txt'
	#sample_abundance_plotter(alpha_path, sample_abundance_file_name, args.design, args.name)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step9: Bacterial_ABUNDANCE_FILE.
	bacterial_abundance_columns_dict = {}
	bacterial_abundance_category_dict = {}
	bacterial_abundance_file = args.outputdir + args.name + '_BACTERIAL_ABUNDANCE_file_STEP9.txt'
	bacterial_abundance_columns_dict, bacterial_abundance_category_dict, bact_list = bacterial_abundance_table_plotly(args.shared, args.design, bacterial_abundance_file, data_path, args.outputdir)
	bacterial_abundance_javascript_plotly(bacterial_abundance_columns_dict, bacterial_abundance_category_dict, bact_list, args.name, args.design, js_path)
	bacterial_abundance_file_name = args.name + '_BACTERIAL_ABUNDANCE_file_STEP9.txt'
	#bacterial_abundance_plotter(alpha_path, bacterial_abundance_file_name, args.design, args.name)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step10: ALPHA DIVERSITY SUMMARY FILE.
	alpha_diversity_summary_file = args.outputdir + args.name + '_ALPHA_DIVERSITY_SUMMARY_file_STEP10.txt'
	flag = summary_analysis(args.shared, alpha_diversity_summary_file, mothur_exec_path, args.processors, data_path, args.outputdir)
	summary_table_header, summary_table_body = summary_table(alpha_diversity_summary_file, args.design)
	alpha_diversity_summary_file_name = args.name + '_ALPHA_DIVERSITY_SUMMARY_file_STEP10.txt'
	#alpha_summary_plotter(alpha_path, alpha_diversity_summary_file_name, summary_table_header, summary_table_body, args.name, args.design)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step11: Biomarker Discovery file.
	# first we run the permutation of shared file and design files
	
	permuted_shared_dict, permuted_design_dict = shared_design_all_permutation(args.shared, args.design, mothur_exec_path, args.processors, args.outputdir, data_path)
	permuted_file_dict, permuted_table_dict = biomarker_discovery(permuted_shared_dict, permuted_design_dict, args.name, mothur_exec_path, args.processors, args.outputdir, data_path)
	biomarker_discovery_string = biomarker_discovery_html_maker(permuted_file_dict, permuted_table_dict)
	#biomarker_discovery_plotter(alpha_path, biomarker_discovery_string, args.name, args.design)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step12: Binomial PCA.
	# first we run the binomial permutation
	permuted_shared_list, permuted_design_list = shared_design_binomial_permutation(args.shared, args.design, mothur_exec_path, args.processors, args.outputdir, data_path)
	pca_html_string_dict = {}
	pca_html_string_dict = pca_3d_analyser(permuted_shared_list, permuted_design_list, args.name, mothur_exec_path, args.processors, args.outputdir, data_path, js_path)
	#pca_3d_plotter(alpha_path, pca_html_string_dict, args.name, args.design)
	#heatmap_analysis_plotly(args.shared, args.design, mothur_exec_path, args.processors, args.outputdir)
	#sys.exit(2)
	all_plotter(alpha_path, rarefaction_file_name, sample_abundance_file_name, bacterial_abundance_file_name, summary_table_header, summary_table_body, alpha_diversity_summary_file_name, biomarker_discovery_string, pca_html_string_dict, args.name, args.design)
	#sys.exit(2)

	"""
	#PCOA CALCULATIONS
	calc_list1 = ['braycurtis: Bray-Curtis index', 'manhattan: Manhattan distance(Taxicab geometry)', 'spearman: Spearman rank correlation coefficient', 'thetayc: the Yue & Clayton measure of dissimilarity']
	for calc in calc_list1:

		distance_matrix = distance_matrix_analyser(args.shared, calc, mothur_exec_path, args.name, args.processors, args.outputdir)
		pcoa_javastring, pcoa_sample_string, xmax, xmin, ymax, ymin, zmax, zmin = P3D_analyser_plotly(distance_matrix, args.design, calc, mothur_exec_path, args.processors, args.outputdir)
		pca_javascript_plotly(pcoa_javastring, pcoa_sample_string, args.design, xmax, xmin, ymax, ymin, zmax, zmin, calc, args.name, js_path)
		
	calc_list2 = ['jclass: Jaccard index', 'hamming: Hamming distance', 'anderberg: Anderberg index', 'whittaker: Whittaker index']
	for calc in calc_list2:

		distance_matrix = distance_matrix_analyser(args.shared, calc, mothur_exec_path, args.name, args.processors, args.outputdir)
		pcoa_javastring, pcoa_sample_string, xmax, xmin, ymax, ymin, zmax, zmin = P3D_analyser_plotly(distance_matrix, args.design, calc, mothur_exec_path, args.processors, args.outputdir)
		pca_javascript_plotly(pcoa_javastring, pcoa_sample_string, args.design, xmax, xmin, ymax, ymin, zmax, zmin, calc, args.name, js_path)
	"""
	copytree(args.outputdir + '/ZAC_ANALYSIS/', '/usr/share/nginx/html/')
	#html_maker(alpha_path, args.design, args.name, summary_table_header, summary_table_body, biomarker_discovery_string)
	# Beautiful Shared Table
	#shared_table_header, shared_table_body = shared_table_maker(args.shared, args.design)
	#shutil.copytree(args.outputdir + '/ZAC_ANALYSIS/', '/usr/share/nginx/html/')
	
	remove_mothur_log(os.getcwd())


# ################################### MAIN FUNCTIONS  ########################## #

def parse_abundance_file(abundance_file):
	delimiter = sniff(abundance_file)
	file_handle = open(abundance_file, 'rU')
	data_list = []
	header_list = []
	for line in file_handle:
		line = line.rstrip()
		line_list = []
		line_list = line.split(delimiter)
		if line_list[0].lower() == 'taxlevel':
			header_list = line_list
			#print header_list
		else:
			data_hash = {}
			for key, value in enumerate(line_list):
				if header_list[key] == 'taxlevel':
					data_hash[header_list[key]] = phylogeny_convert(value)
				else:
					data_hash[header_list[key]] = value
					
			#print data_hash['taxlevel']
			data_list.append(data_hash)
	return data_list


def phylogeny_convert(string):
	if string == '0':
		return 'root'
	elif string == '1':
		return 'kingdom'
	elif string == '2':
		return 'phylum'
	elif string == '3':
		return 'class'
	elif string == '4':
		return 'order'
	elif string == '5':
		return 'family'
	elif string == '6':
		return 'genus'
	elif string == '7':
		return 'species'


def construct_object(entity, entity_list):
	dict_of_entities = {}
	for each_entity in entity_list:
		#print each_entity
		a_class = entity(each_entity)
		dynamic_variable = ''
		dynamic_variable = a_class.get_taxon() + '_' + a_class.get_taxlevel()
		dict_of_entities[dynamic_variable] = a_class
	return dict_of_entities


def create_shared_table(dict_of_entities, abundance_file, taxlevel, name, outputdir):
	all_entities = dict_of_entities.keys()
	group_list = ['Groups']
	group_list.extend(dict_of_entities[dict_of_entities.keys()[0]].get_abundance_keys())
	label_list = ['label', taxlevel]
	num_temp = []
	shared_list = []
	for each_entity in all_entities:
		if dict_of_entities[each_entity].get_taxlevel() == taxlevel.lower():
			temp_list = []
			temp_list.append(dict_of_entities[each_entity].get_taxon())
			num_temp.append(dict_of_entities[each_entity].get_taxon())
			temp_list.extend(dict_of_entities[each_entity].get_abundance_values())
			shared_list.append(temp_list)
	num_otus_list = ['numOtus', str(len(num_temp))]
	transposed_shared_list = transpose(shared_list)
	shared_string = ''
	line_counter = 0
	for group, entity in itertools.izip(group_list, transposed_shared_list):
		if line_counter == 0:
			shared_string += label_list[0] + '\t' + group + '\t' + num_otus_list[0] + '\t' + list_to_string(list(entity)) + '\n'
			
		else:
			shared_string += label_list[1] + '\t' + group + '\t' + num_otus_list[1] + '\t' + list_to_string(list(entity)) + '\n'
			
		line_counter += 1

	shared_name = outputdir + name + '_' + taxlevel + '.shared'
	write_string_down(shared_string, shared_name)
	return shared_name


def fix_shared_design_file(shared_file, design_file, absname, new_shared_name, new_design_file):
	f = open(shared_file, 'rU')
	new_shared_string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] == 'label':
			for otu in line:
				if otu in ['Group', 'label', 'numOtus', 'group']:
					new_shared_string += otu + '\t'
				else:
					#new_shared_string += 'OTU' + str(otu_counter).zfill(zfill_value) + ';' + otu + '\t'
					new_shared_string += otu + '\t'
			new_shared_string = new_shared_string[:-1]
			new_shared_string += '\n'
		else:
			line = line[1:]
			new_shared_string += absname + '\t' + '\t'.join(line) + '\n'
	f.close()
	write_string_down(new_shared_string, new_shared_name)
	new_shared_string = ''
	copy_file(design_file, new_design_file)
	return True


def update_shared_design_file(shared_file, design_file, control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# we update shared file based on this control file
	f = open(control_file, 'rU')
	control_string = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_string += i + '-'
	control_string = control_string[:-1]
	f.close()

	flag, stderr = execute_functions(mothur_get_groups, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, control_string)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file)
		scanned_container = []
		extension_list = ['.' + name + '.pick' + extension]
		flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], new_shared_file)
	copy_file(design_file, new_design_file)
	return True


def keep_sample_shared_design(shared_file, design_file, control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# First we are reading sample listed in control file in string and list
	f = open(control_file, 'rU')
	control_string = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_string += i + '-'
	control_string = control_string[:-1]
	f.close()
	# Second we are filtering design file based on what is listed on control list
	new_design_string = 'Sample\tDesign\n'
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] in control_data:
			new_design_string += i + '\n'
	f.close()
	o = open(new_design_file, 'w')
	o.write(new_design_string)
	o.close()
	new_design_string = ''
	control_data = []
	#Third we are going to keep samples that only available in list at shared file
	flag, stderr = execute_functions(mothur_get_groups, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, control_string)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file)
		scanned_container = []
		extension_list = ['.' + name + '.pick' + extension]
		flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], new_shared_file)
	return True


def remove_sample_shared_design(shared_file, design_file, control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# First we are reading sample listed in control file in string and list
	f = open(control_file, 'rU')
	control_string = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_string += i + '-'
	control_string = control_string[:-1]
	f.close()
	# Second we are filtering design file based on what is listed on control list
	new_design_string = ''
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] not in control_data:
			new_design_string += i + '\n'
	f.close()
	o = open(new_design_file, 'w')
	o.write(new_design_string)
	o.close()
	new_design_string = ''
	control_data = []
	#Third we are going to keep samples that only available in list at shared file
	flag, stderr = execute_functions(mothur_remove_groups, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, control_string)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file)
		scanned_container = []
		extension_list = ['.' + name + '.pick' + extension]
		flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], new_shared_file)
	return True


def keep_otu_shared_design(shared_file, design_file, otu_control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# First we are reading sample listed in control file in string and list
	f = open(otu_control_file, 'rU')
	otu_control_list = []
	for i in f:
		i = i.rstrip()
		otu_control_list.append(i)
	f.close()
	# #####################################################
	keep_otu_list = list(set(otu_control_list))
	print "Keeping Following OTUs:", keep_otu_list
	header = []
	columns = {}
	
	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	sim_head = ''
	for head in header:
		if ';' in head:
			sim_head_list = head.split(';')
			sim_head = ';'.join(sim_head_list[1:])
			#print sim_head
		if sim_head.lower() in ['group', 'numotus', 'label']:
			continue
		elif sim_head not in keep_otu_list:
			continue
		elif sim_head in keep_otu_list:
			OTU_dict[head] = columns[head]
		else:
			print "shared file is not sane."
	#print OTU_dict
	#sys.exit(2)
	
	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	copy_file(design_file, new_design_file)
	return True


def remove_otu_shared_design(shared_file, design_file, otu_control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# First we are reading sample listed in control file in string and list
	f = open(otu_control_file, 'rU')
	otu_control_list = []
	for i in f:
		i = i.rstrip()
		otu_control_list.append(i)
	f.close()
	# #####################################################
	remove_otu_list = list(set(otu_control_list))
	print "Removing Following OTUs:", remove_otu_list
	header = []
	columns = {}
	
	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	sim_head = ''
	for head in header:
		if ';' in head:
			sim_head_list = head.split(';')
			sim_head = ';'.join(sim_head_list[1:])
			#print sim_head
		elif head.lower() in ['group', 'numotus', 'label']:
			continue
		if sim_head in remove_otu_list:
			continue
		elif sim_head not in remove_otu_list:
			OTU_dict[head] = columns[head]
		else:
			print "shared file is not sane."
	#print OTU_dict
	#sys.exit(2)
	
	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	copy_file(design_file, new_design_file)
	return True


def filter_shared_design_file(shared_file, design_file, control_file, name, new_shared_file, new_design_file, mothur_exec_path, processors, outputdir):
	# ########################################################################
	path, absname, ext = split_file_name(shared_file)
	f = open(control_file, 'rU')
	control_string = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_string += i + '-'
	f.close()
	# ########################################################################
	f = open(design_file, 'rU')
	new_design_string = ''
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] in control_data:
			continue
		new_design_string += line[0] + '\t' + line[1] + '\n'
	f.close()
	controlled_design_file = outputdir + absname + '_controlled_design.txt'
	write_string_down(new_design_string, controlled_design_file)
	os.rename(controlled_design_file, new_design_file)
	# ########################################################################
	removing_groups = control_string[:-1]
	flag, stderr = execute_functions(mothur_remove_groups, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, removing_groups)
	if flag is False:
		print "Execution of mothur_rarefaction_single failed!!!"
	else:
		scanned_container = []
		extension_list = ['.' + name + '.pick' + ext]
		flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	removing_groups = ''
	os.rename(scanned_container[0], new_shared_file)
	return True


def remove_OTUs_from_shared(shared_file, otu_control, new_shared_file):
	removing_otu_list = []
	f = open(otu_control, 'rU')
	for i in f:
		i = i.rstrip()
		removing_otu_list.append(i)
	removing_otu_list = list(set(removing_otu_list))
	#print "Removing Following OTUs:", removing_otu_list
	header = []
	columns = {}
	
	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		if head in removing_otu_list:
			continue
		else:
			OTU_dict[head] = columns[head]

	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	return True


def keep_OTUs_in_shared(shared_file, otu_control_list, new_shared_file):
	
	removing_otu_list = list(set(otu_control_list))
	#print "Removing Following OTUs:", removing_otu_list
	header = []
	columns = {}
	
	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		if head not in removing_otu_list:
			continue
		else:
			OTU_dict[head] = columns[head]

	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	return True


def mothur_shared_parser(mothur_shared_file):
	sharedfile = open(mothur_shared_file, 'rU')
	header = []
	header = sharedfile.readline().rstrip().split('\t')
	#print header
	sharedfile.close()
	sharedfile = open(mothur_shared_file, 'rU')
	reader = csv.DictReader(sharedfile, header, delimiter="\t")
	columns = {}
	for row in reader:
		for key, value in row.items():
			if key == value:
				continue
			elif key in columns:
				columns[key].append(value)
			else:
				columns[key] = [value]
	return (header, columns)


def normalize_shared_design_maker(shared_file, design_file, name, new_shared_file, new_design_file, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	method = 'totalgroup'
	flag, stderr = execute_functions(mothur_normalize, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, method)
	if flag is False:
		print "Execution of mothur_normalized failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file)
		scanned_container = []
		extension_list = ['.' + name + '.norm.shared']
		flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], new_shared_file)
	# we filter
	f = open(new_shared_file, 'rU')
	group_list = []
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[1].lower() == 'group':
			continue
		else:
			group_list.append(line[1])
	f.close()
	new_design_string = 'Sample\tDesign\n'
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] in group_list:
			new_design_string += i + '\n'
	f.close()
	write_string_down(new_design_string, new_design_file)
	new_design_string = ''
	return True


def subsample_shared_design_maker(shared_file, design_file, name, new_shared_file, new_design_file, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_subsample, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file)
	if flag is False:
		print "Execution of mothur_subsample failed!!!"
	else:
		scanned_container = []
		extension_list = ['.' + name + '.subsample' + ext]
		flag = scandirs(outputdir, scanned_container, extension_list)
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], new_shared_file)
	# we filter
	f = open(new_shared_file, 'rU')
	group_list = []
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[1].lower() == 'group':
			continue
		else:
			group_list.append(line[1])
	f.close()
	new_design_string = 'Sample\tDesign\n'
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] in group_list:
			new_design_string += i + '\n'
	f.close()
	write_string_down(new_design_string, new_design_file)
	new_design_string = ''
	return True


def mothur_rarefaction_curve(shared_file, design_file, rarefaction_file, mothur_exec_path, processors, data_path, outputdir):
	#min_sh, max_sh, ave_sh, median_sh = shared_max_min(shared_file)
	ave_sh = 0
	flag, stderr = execute_functions(mothur_rarefaction_single, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, ave_sh)
	if flag is False:
		print "Execution of mothur_rarefaction_single failed!!!"
		sys.exit(2)
	else:
		scanned_container = []
		extension_list = ['.groups.r_chao']
		flag = scandirs(outputdir, scanned_container, extension_list)
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	flag = remove_extension_files(outputdir, '.rabund')
	os.rename(scanned_container[0], rarefaction_file)
	copy_file(rarefaction_file, data_path)
	header = []
	columns = {}
	header, columns = mothur_result_parser(rarefaction_file)
	
	design_dictionary = {}
	design_dictionary = design_dict_maker(design_file)
	rarefaction_series_dict = {}
	rarefaction_xaxis_string = ''
	
	for head in header:
		if 'lci' in head or 'hci' in head:
			continue
		elif 'numsampled' in head:
			rarefaction_xaxis_string = list_to_string(columns[head], ',')

	for head in header:
		if 'lci' in head or 'hci' in head or 'numsampled' in head:
			continue
		else:
			rarefaction_series_dict[head] = "		var " + slugify(head) + "= {\n"
			rarefaction_series_dict[head] += "			x:" + "[" + rarefaction_xaxis_string + "],\n"
			rarefaction_series_dict[head] += "			y:" + "[" + list_to_string(list_variance(columns[head], 0.01), ',') + "],\n"
			rarefaction_series_dict[head] += "			mode: 'scatter',\n"
			rarefaction_series_dict[head] += "			mode: 'lines',\n"
			rarefaction_series_dict[head] += "			hoverinfo: 'text',\n"
			rarefaction_series_dict[head] += "			text:'" + head.split('-')[1] + "',\n"
			rarefaction_series_dict[head] += "			name: '" + head.split('-')[1] + "'\n};\n"
	rarefaction_design_dict = {}
	plotly_design_data_dict = {}
	for each_sample in rarefaction_series_dict.keys():
		each_sample_design = each_sample.split('-')[1]

		if design_dictionary[each_sample_design] not in rarefaction_design_dict:
			rarefaction_design_dict[design_dictionary[each_sample_design]] = rarefaction_series_dict[each_sample]
			plotly_design_data_dict[design_dictionary[each_sample_design]] = slugify(each_sample) + ','
		else:
			rarefaction_design_dict[design_dictionary[each_sample_design]] += rarefaction_series_dict[each_sample]
			plotly_design_data_dict[design_dictionary[each_sample_design]] += slugify(each_sample) + ','
	
	return (rarefaction_design_dict, plotly_design_data_dict)


def remove_extension_files(outputdir, extension):
	extension_list = []
	extension_list.append(extension)
	scanned_container = []
	flag = scandirs(outputdir, scanned_container, extension_list, 'ex_partial')
	print "Scanning to find", extension, "Files.."
	if flag is False:
		print "Failed :("
		print "This extension is not available: ", extension_list
		sys.exit(2)
	else:
		print "NICE :) Found them"
		counter = 1
		for file in scanned_container:
			#print "Removing File#", str(counter), ":", file
			counter += 1
			check_it_and_remove_it(file)
	return True


def copy_file(source_file, destination):
	shutil.copy2(source_file, destination)
	return True


def abundance_table_plotly(shared_file, design_file, simple_abundance_file, data_path, outputdir):
	# ###############################################################################
	# first we aggregated the similar OTU content
	#aggregated_shared = outputdir + 'sample_abundance_aggregated_shared_file.txt'
	#aggregate_otu_shared_file(shared_file, aggregated_shared, simple_abundance_file)
	#copy_file(simple_abundance_file, data_path)
	# ###############################################################################
	# Second we parse it
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)
	#shared_dict = parse_shared_file(aggregated_shared)
	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	abundance_otu_dict = {}
	abundance_columns_dict = {}

	sample_name_list = shared_dict.keys()
	otu_name_list = shared_dict[sample_name_list[0]].keys()
	otu_name_color_dict = symbol_color_dict(otu_name_list, 'color_only')
	otu_counter = 0
	for each_design in reverse_design_dict:
		sample_list = reverse_design_dict[each_design]
		otu_value_list = []
		
		xaxis = list_to_string(sample_list, "','")
		if each_design not in abundance_columns_dict:
			abundance_columns_dict[each_design] = ''
		for each_otu in otu_name_list:
			otu_value_list = []
			for each_sample in sample_list:
				otu_value_list.append(float(shared_dict[each_sample][each_otu]))
			if all(value == 0 for value in otu_value_list) is True:
				continue
			otu_json_string = ''
			otu_temp = each_otu.split(';')
			otu_hover_name = ''
			
			if len(otu_temp) > 2:
				otu = otu_temp[0] + ';' + otu_temp[-3] + ';' + otu_temp[-2]
				otu_hover_name = otu_temp[-2]
			
			elif len(otu_temp) > 1:
				otu = otu_temp[0] + ';' + otu_temp[1]
				otu_hover_name = otu_temp[1]
			else:
				otu = otu_temp[0]
				otu_hover_name = otu
			otu_counter += 1
			if each_design not in abundance_otu_dict:
				abundance_otu_dict[each_design] = 'OTU' + str(otu_counter).zfill(10) + ','
			else:
				abundance_otu_dict[each_design] += 'OTU' + str(otu_counter).zfill(10) + ','
			
			otu_json_string += "\tvar " + 'OTU' + str(otu_counter).zfill(10) + " = {\n"
			otu_json_string += "\tx:" + "['" + xaxis + "'],\n"
			otu_json_string += "\t\ty:" + "["
			for each_sample in sample_list:
				otu_json_string += shared_dict[each_sample][each_otu] + ','
			otu_json_string += "],\n"
			otu_json_string += "\t\ttype: 'bar',\n"

			otu_json_string += "\t\thoverinfo: 'all',\n"
			otu_json_string += "\t\tautobinx: false,\n"
			if each_otu in otu_name_color_dict:
				otu_json_string += "\t\tmarker: {color: '" + otu_name_color_dict[each_otu] + "'},\n"
			else:
				otu_json_string += "\t\tmarker: {color: '" + random_color() + "'},\n"
			otu_json_string += "\t\txbins: {size: 0.2},\n"
			otu_json_string += "\t\torientation: 'v',\n"
			
			#otu_json_string += "\t\ttext:'" + slugify(otu) + "',\n"
			otu_json_string += "\t\tname: '" + slugify(otu_hover_name) + "'\n};"

			abundance_columns_dict[each_design] += otu_json_string
	
	return (abundance_columns_dict, abundance_otu_dict)


def aggregate_otu_shared_file(shared_file, new_shared_file, simple_file):
	header = []
	columns = {}
	header, columns = mothur_shared_parser(shared_file)
	bact_dict = {}
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		bact_split = head.split(';')
		if len(bact_split) > 1:
			bact_name = bact_split[-2]
		else:
			bact_name = bact_split[0]
		if bact_name not in bact_dict:
			bact_dict[bact_name] = columns[head]
		else:
			bact_dict[bact_name] = [sum(x) for x in zip(list_string_to_float(bact_dict[bact_name]), list_string_to_float(columns[head]))]

	numotus = str(len(bact_dict.keys()))
	new_shared_string = ''
	simple_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + bact_dict.keys())
	simple_string = list_to_string(['Group'] + bact_dict.keys())
	new_shared_string += '\n'
	simple_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		simple_string += columns['Group'][i]
		simple_string += '\t'
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in bact_dict.keys():
			new_shared_string += str(bact_dict[each_otu][i])
			new_shared_string += '\t'
			simple_string += str(bact_dict[each_otu][i])
			simple_string += '\t'
		new_shared_string += '\n'
		simple_string += '\n'

	write_string_down(new_shared_string, new_shared_file)
	write_string_down(simple_string, simple_file)
	return True


def bacterial_abundance_table_plotly(shared_file, design_file, bacterial_abundance_file, data_path, outputdir):
	# ###############################################################################
	# first we aggregated the similar OTU content
	#aggregated_shared = outputdir + 'bacterial_abundance_aggregated_shared_file.txt'
	#aggregate_otu_shared_file(shared_file, aggregated_shared, bacterial_abundance_file)
	#copy_file(bacterial_abundance_file, data_path)
	# ###############################################################################
	# Second we parse it
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)
	#shared_dict = parse_shared_file(aggregated_shared)
	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	sample_name_list = shared_dict.keys()
	bacterial_list = shared_dict[sample_name_list[0]].keys()
	#bacterial_list = sort_otus(shared_file)
	bact_list = []
	for each_bact in bacterial_list:
		
		bact = each_bact.split(';')
		if len(bact) > 1:
			bact_list.append(bact[-1])
		else:
			bact_list.append(bact[0])
	
	#print bacterial_abundance_category_dict
	xaxis = list_to_string(bact_list, "','")
	otu_name_color_dict = symbol_color_dict(bact_list, 'color_only')
	#print bacterial_abundance_category_dict['Adenoma']
	bacterial_abundance_sample_dict = {}
	bacterial_abundance_columns_dict = {}
	for each_design in reverse_design_dict:
		sample_list = reverse_design_dict[each_design]
		if each_design not in bacterial_abundance_columns_dict:
			bacterial_abundance_columns_dict[each_design] = ''
		for each_sample in sample_list:
			if each_design not in bacterial_abundance_sample_dict:
				bacterial_abundance_sample_dict[each_design] = each_sample + ','
			else:
				bacterial_abundance_sample_dict[each_design] += each_sample + ','
			otu_value_list = []
			for each_bact in bacterial_list:
				otu_value_list.append(float(shared_dict[each_sample][each_bact]))
			if all(value == 0 for value in otu_value_list) is True:
				continue
			#print each_sample
			otu_json_string = ''
			#otu_json_string += '{\n' + 'name:"' + each_sample + '",\n'
			otu_json_string += "\tvar " + slugify(each_sample) + "= {\n"
			otu_json_string += "\tx:" + "['" + xaxis + "'],\n"
			otu_json_string += "\t\ty:" + "["
			otu_json_string += list_to_string(shared_dict[each_sample].values(), ',')
			otu_json_string += "],\n"
			otu_json_string += "\t\ttype: 'bar',\n"
			otu_json_string += "\t\thoverinfo: 'all',\n"
			otu_json_string += "\t\tautobinx: false,\n"
			otu_json_string += "\t\tmarker: {color: '"
			for each_bact in bact_list:
				if each_bact in otu_name_color_dict:
					otu_json_string += otu_name_color_dict[each_bact] + ', '
				else:
					otu_json_string += random_color() + ', '
			otu_json_string += "'},\n"
			otu_json_string += "\t\txbins: {size: 0.2},\n"
			otu_json_string += "\t\torientation: 'v',\n"
			otu_json_string += "\t\tname: '" + slugify(each_sample) + "'\n};"
			
			bacterial_abundance_columns_dict[each_design] += otu_json_string
	
	return(bacterial_abundance_columns_dict, bacterial_abundance_sample_dict, bact_list)


def summary_analysis(shared_file, summary_file_name, mothur_exec_path, processors, data_path, outputdir):
	flag, stderr = execute_functions(mothur_summary_single, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file)
	if flag is False:
		print "Execution of summary_single failed!!!"
	else:
		scanned_container = []
		extension_list = ['.groups.summary']
		flag = scandirs(outputdir, scanned_container, extension_list)
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], summary_file_name)
	flag = remove_extension_files(outputdir, '.rabund')
	copy_file(summary_file_name, data_path)
	return True


def biomarker_discovery(permuted_shared_dict, permuted_design_dict, name, mothur_exec_path, processors, outputdir, data_path):
	permuted_keys = []
	permuted_keys = permuted_shared_dict.keys()
	permuted_file_dict = {}
	permuted_table_dict = {}
	permuted_kruskal_wallis_dict = {}
	permuted_indicator_dict = {}
	permuted_lefse_dict = {}
	permuted_metastats_dict = {}
	for permuted_key in permuted_keys:
		print "calculating ", permuted_key
		permuted_file_dict[permuted_key] = {}
		permuted_table_dict[permuted_key] = {}
		permuted_kruskal_wallis_dict[permuted_key] = {}
		
		permuted_indicator_dict[permuted_key] = {}
		
		permuted_lefse_dict[permuted_key] = {}
		
		permuted_metastats_dict[permuted_key] = {}
		
		for each_shared, each_design in itertools.izip(permuted_shared_dict[permuted_key], permuted_design_dict[permuted_key]):
			temp_name = each_shared.split('Shared_file_')[1]
			temp_name = temp_name.split('.txt')[0]
			significant_OTUs_dict = {}
			significant_OTUs_list = []
			#kruskal_name = 'kruskal_wallis_' + temp_name
			significant_OTUs_dict['kruskal_wallis'], kw_list = kw_analyser(each_shared, each_design, name, mothur_exec_path, processors, outputdir)
			significant_OTUs_list.extend(kw_list)
			
			#////////////////////////////////////
			#lefse_name = 'lefse_' + temp_name
			significant_OTUs_dict['lefse'], lefse_list = lefse_analyser(each_shared, each_design, name, mothur_exec_path, processors, outputdir)
			significant_OTUs_list.extend(lefse_list)

			#////////////////////////////////////
			#indicator_name = 'indicator_' + temp_name
			significant_OTUs_dict['indicator'], indicator_list = indicator_analyser(each_shared, each_design, name, mothur_exec_path, processors, outputdir)
			significant_OTUs_list.extend(indicator_list)

			#////////////////////////////////////
			#metastats_name = 'metastats_' + temp_name
			significant_OTUs_dict['metastats'], metastats_list = metastats_analyser(each_shared, each_design, name, mothur_exec_path, processors, outputdir)
			significant_OTUs_list.extend(metastats_list)
			
			permuted_table_dict[permuted_key][temp_name] = {}
			signature_header_table, signature_body_table = signature_table(significant_OTUs_dict, each_design, name)
			permuted_table_dict[permuted_key][temp_name] = (signature_header_table, signature_body_table)
			permuted_file_dict[permuted_key][temp_name] = {}
			if len(significant_OTUs_list) > 0:
				ranked_shared = data_path + name + '_ranked_shared_file.txt'
				flag = keep_OTUs_in_shared(each_shared, significant_OTUs_list, ranked_shared)
				if flag is False:
					print "Warning"
				new_orange = data_path + 'Orange_file_filtered_' + temp_name + '.txt'
				orange_shared = orange_juice(ranked_shared, each_design, name, mothur_exec_path, processors, outputdir)
				os.rename(orange_shared, new_orange)

				new_orange = data_path + 'Orange_file_raw_' + temp_name + '.txt'
				orange_shared = orange_juice(each_shared, each_design, name, mothur_exec_path, processors, outputdir)
				os.rename(orange_shared, new_orange)
				
				new_shared = data_path + 'Abundance_file_filtered_' + temp_name + '.txt'
				os.rename(ranked_shared, new_shared)
				
				new_shared = data_path + 'Abundance_file_raw_' + temp_name + '.txt'
				os.rename(each_shared, new_shared)
				permuted_file_dict[permuted_key][temp_name] = (new_shared, new_orange)
			else:
				print "EMPTY LIST DETECTED"

	return(permuted_file_dict, permuted_table_dict)


def pca_3d_analyser(permuted_shared_list, permuted_design_list, name, mothur_exec_path, processors, outputdir, data_path, js_path):
	calc_list = ['braycurtis: Bray-Curtis index', 'manhattan: Manhattan distance(Taxicab geometry)']
	pca_html_string_dict = {}
	for each_shared, each_design in itertools.izip(permuted_shared_list, permuted_design_list):
		PCA_name = each_shared.split('binomial_')[1]
		PCA_name = PCA_name.split('.txt')[0]
		distance_matrix = outputdir + name + '_' + PCA_name + '_braycurtis_distance_matrix.txt'
		flag = distance_matrix_analyser(each_shared, calc_list[0], distance_matrix, mothur_exec_path, name, processors, outputdir)
		if flag is False:
			print "Something is wrong with distance_matrix_analyser."
		pcoa_javastring, pcoa_sample_string, xmax, xmin, ymax, ymin, zmax, zmin = P3D_analyser_plotly(distance_matrix, each_design, calc_list[0], mothur_exec_path, processors, outputdir)
		pca_javascript_plotly(pcoa_javastring, pcoa_sample_string, PCA_name, each_design, xmax, xmin, ymax, ymin, zmax, zmin, calc_list[0], name, js_path)
		pca_html_maker(pca_html_string_dict, PCA_name, calc_list[0], name)
		print "##################################################################################################################################"
	return pca_html_string_dict


def distance_matrix_analyser(shared_file, calc, matrix_name, mothur_exec_path, absname, processors, outputdir):
	calc = calc.split(':')[0]
	flag, stderr = execute_functions(mothur_distance_shared, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, calc)
	if flag is False:
		print "Execution of mothur_distance_matrix failed!!!"
	else:
		scanned_container = []
		extension_list = ['.' + calc + '.' + absname + '.square.dist']
		flag = scandirs(outputdir, scanned_container, extension_list)
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container[0], matrix_name)
	return True


def P3D_analyser_plotly(dist_file, design_file, calc, mothur_exec_path, processors, outputdir):
	flag, stderr = execute_functions(mothur_PCOA, processors, outputdir, 'multi', 'mothur', mothur_exec_path, dist_file)
	if flag is False:
		print "Execution of mothur_rarefaction_single failed!!!"
	else:
		scanned_container = []
		extension_list = ['.pcoa.axes']
		flag = scandirs(outputdir, scanned_container, extension_list)
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container:
				print "File#", str(counter), ":", file
				counter += 1
	dict_of_design = {}
	pcoa_dict = {}
	dict_of_design = design_dict_maker(design_file)
	pcoafile = open(scanned_container[0], 'rU')
	
	xmax = ymax = zmax = 0
	xmin = ymin = zmin = 0
	for line in pcoafile:
		line = line.rstrip()
		line_split = line.split('\t')
		if line_split[0].lower() == 'group':
			continue
		elif line_split[0] in dict_of_design.keys():
			if len(line_split) < 4:
				zero_list = ['0.0', '0.0', '0.0']
				line_split.extend(zero_list)
			#pcoa_dict[line_split[0]] = (dict_of_design[line_split[0]], '[' + line_split[1] + ',' + line_split[2] + ',' + line_split[3] + ']')
			coord_list = []
			coord_list.extend([line_split[1], line_split[2], line_split[3]])
			pcoa_dict[line_split[0]] = (dict_of_design[line_split[0]], coord_list)
			if isfloat(line_split[1]) and float(line_split[1]) > xmax: xmax = float(line_split[1])
			if isfloat(line_split[2]) and float(line_split[2]) > ymax: ymax = float(line_split[2])
			if isfloat(line_split[3]) and float(line_split[3]) > zmax: zmax = float(line_split[3])
			if isfloat(line_split[1]) and float(line_split[1]) < xmin: xmin = float(line_split[1])
			if isfloat(line_split[2]) and float(line_split[2]) < ymin: ymin = float(line_split[2])
			if isfloat(line_split[3]) and float(line_split[3]) < zmin: zmin = float(line_split[3])
	pcoafile.close()
	dict_of_symbol_color = symbol_color_dict(list(set(dict_of_design.values())))
	#print pcoa_dict
	#sys.exit(2)
	pcoa_sample_string = ''
	pcoa_javastring = ''
	#pcoa_javastring += "\n\tvar " + calc + "= {\n"
	for key, value in pcoa_dict.iteritems():
		pcoa_sample_string += key + ','
		pcoa_javastring += "\n\tvar " + key + "= {\n"
		pcoa_javastring += "\t\tx:" + "['" + value[1][0] + "'],\n"
		pcoa_javastring += "\t\ty:" + "['" + value[1][1] + "'],\n"
		pcoa_javastring += "\t\tz:" + "['" + value[1][2] + "'],\n"
		pcoa_javastring += "\t\ttype: 'scatter3d',\n"
		pcoa_javastring += "\t\thoverinfo: 'text+name',\n"
		pcoa_javastring += "\t\thovermode: 'closest',\n"
		pcoa_javastring += "\t\tmode: 'markers',\n"
		pcoa_javastring += "\t\t'legendgroup': '" + value[0] + "',\n"
		pcoa_javastring += "\t\tmarker: {color: '" + dict_of_symbol_color[value[0]][1] + "', size: 5, symbol: '" + dict_of_symbol_color[value[0]][0] + "', opacity: 0.9 },\n"
		pcoa_javastring += "\t\tname: '" + value[0] + "',\n"
		pcoa_javastring += "\t\ttext: '" + slugify(key) + "'\n"
		pcoa_javastring += "\t};"
		"""
		pcoa_javastring += '\n\t\t\t{\n'
		pcoa_javastring += '\t\t\t\tname: "' + key + '[' + value[0] + ']",\n'
		pcoa_javastring += '\t\t\t\tdata: [' + value[1] + '],\n'
		pcoa_javastring += '\t\t\t\tcolor: '" + dict_of_symbol_color[value[0]][1] + '",\n'
		pcoa_javastring += '\t\t\t\tmarker:{symbol:"' + dict_of_symbol_color[value[0]][0] + '"}'
		pcoa_javastring += '\t\t\t},\n'
		"""
	check_it_and_remove_it(scanned_container[0])
	flag = remove_extension_files(outputdir, '.pcoa.loadings')
	return (pcoa_javastring, pcoa_sample_string, xmax, xmin, ymax, ymin, zmax, zmin)


def shared_design_all_permutation(shared_file, design_file, mothur_exec_path, processors, outputdir, data_path):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	permuted_shared_dict = {}
	permuted_design_dict = {}
	rev_keys = design_dictionary_rev.keys()
	for step in range(2, len(rev_keys) + 1):
		permuted_shared_dict['TIER_' + str(step)] = []
		permuted_design_dict['TIER_' + str(step)] = []
		#print list(itertools.combinations(rev_keys, step))
		design_list = list(itertools.combinations(rev_keys, step))
		for design in design_list:
			design_sample_list = []
			#first we extract desired samples
			new_shared_file_name = data_path + 'Shared_file_tier' + str(step) + '_'
			new_design_file_name = data_path + 'Design_file_tier' + str(step) + '_'
			for each_class in design:
				new_shared_file_name += each_class + '_vs_'
				new_design_file_name += each_class + '_vs_'
				sample_list = design_dictionary_rev[each_class]
				design_sample_list.extend(sample_list)
			new_shared_file_name = new_shared_file_name[:-4] + '.txt'
			new_design_file_name = new_design_file_name[:-4] + '.txt'
			design_sample_string = ''
			design_sample_string = list_to_string(design_sample_list, '\n')
			#print all_sample_string
			#sys.exit(2)
			design_sample_string_file = outputdir + 'design_sample_string_file.txt'
			write_string_down(design_sample_string, design_sample_string_file)
			#sys.exit(2)
			control_shared, control_design = mothur_update(shared_file, design_file, design_sample_string_file, mothur_exec_path, processors, outputdir)
			os.rename(control_shared, new_shared_file_name)
			os.rename(control_design, new_design_file_name)
			check_it_and_remove_it(design_sample_string_file)
			permuted_shared_dict['TIER_' + str(step)].append(new_shared_file_name)
			permuted_design_dict['TIER_' + str(step)].append(new_design_file_name)
	
	return (permuted_shared_dict, permuted_design_dict)


def shared_design_binomial_permutation(shared_file, design_file, mothur_exec_path, processors, outputdir, data_path):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	permuted_shared_list = []
	permuted_design_list = []
	rev_keys = design_dictionary_rev.keys()
	design_list = list(itertools.combinations(rev_keys, 2))
	for design in design_list:
		design_sample_list = []
		#first we extract desired samples
		new_shared_file_name = data_path + 'Shared_file_binomial' + '_'
		new_design_file_name = data_path + 'Design_file_binomial' + '_'
		for each_class in design:
			new_shared_file_name += each_class + '_vs_'
			new_design_file_name += each_class + '_vs_'
			sample_list = design_dictionary_rev[each_class]
			design_sample_list.extend(sample_list)
		new_shared_file_name = new_shared_file_name[:-4] + '.txt'
		new_design_file_name = new_design_file_name[:-4] + '.txt'
		design_sample_string = ''
		design_sample_string = list_to_string(design_sample_list, '\n')
		#print all_sample_string
		#sys.exit(2)
		design_sample_string_file = outputdir + 'design_sample_string_file.txt'
		write_string_down(design_sample_string, design_sample_string_file)
		#sys.exit(2)
		control_shared, control_design = mothur_update(shared_file, design_file, design_sample_string_file, mothur_exec_path, processors, outputdir)
		os.rename(control_shared, new_shared_file_name)
		os.rename(control_design, new_design_file_name)
		check_it_and_remove_it(design_sample_string_file)
		permuted_shared_list.append(new_shared_file_name)
		permuted_design_list.append(new_design_file_name)
	
	return (permuted_shared_list, permuted_design_list)


def shared_table_maker(shared_file, design_file):
	dict_of_shared, header_dict = shared_dict(shared_file)
	dict_of_design = design_dict_maker(design_file)
	rev_design_dict = {}
	for k, v in dict_of_design.iteritems():
		rev_design_dict.setdefault(v, []).append(k)
	#print rev_design_dict.items()
	#('WT', ['A9482-L1', 'A9477-L1', 'A9478-L1', 'A9486-L1', 'A9484-L1'])
	#print dict_of_shared.items()
	#('A9486-L1', [('label', 'family'), ('Group', 'A9486-L1'), ('numOtus', '4'), ('*o_RF32', '2'), ('f_Deferribacteraceae', '5336'), ('f_Desulfovibrionaceae', '966'), ('f_F16', '3096')])
	#print header_dict.items()
	#[(0, 'label'), (1, 'Group'), (2, 'numOtus'), (3, '*o_RF32'), (4, 'f_Deferribacteraceae'), (5, 'f_Desulfovibrionaceae'), (6, 'f_F16')]
	
	string_list = ''

	table_header = ''
	table_header += '\t\t<thead>\n'
	table_header += '\t\t\t<tr>\n'
	table_header += '\t\t\t\t<th scope="rowgroup">State</th>\n'
	for key, value in header_dict.iteritems():
		if value in ['numOtus', 'label']:
			continue
		else:
			table_header += '\t\t\t\t<th scope="col">' + value + '</th>\n'

	header_keys = header_dict.keys()
	table_header += '\t\t\t</tr>\n'
	table_header += '\t\t</thead>\n'
	
	table_body = ''
	for key, value_list in rev_design_dict.iteritems():
		if key in ['class', 'sample', 'treatment']:
			continue
		else:
			table_body += '\t\t<tbody align="center">\n'
			table_body += '\t\t\t\t<th rowspan="' + str(len(value_list) + 1) + '" scope="rowgroup" class="metadata">' + key + '</th>\n'
			for each_value in value_list:
				list_of_values_of_each_value = dict_of_shared[each_value]
				#[('label', 'family'), ('Group', 'A9486-L1'), ('numOtus', '4'), ('*o_RF32', '2'), ('f_Deferribacteraceae', '5336'), ('f_Desulfovibrionaceae', '966'), ('f_F16', '3096')])
				table_body += '\t\t\t<tr>\n'
				
				for each_key in header_keys:
					if each_key == 2 or each_key == 0:
						continue
					else:
						
						table_body += '\t\t\t\t<td>' + list_of_values_of_each_value[each_key][1] + '</td>\n'

				table_body += '\t\t\t</tr>\n'
	table_body += '\t\t</tbody>\n'
	return (table_header, table_body)


def orange_juice(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	flag, stderr = execute_functions(mothur_get_relabund, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file)
	if flag is False:
		print "Execution of mothur_get_relabund failed!!!"
	else:
		relabund_container = []
		extension_list = ['.relabund']
		flag = scandirs(outputdir, relabund_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	f = open(relabund_container[0], 'rU')
	design_hash = design_dict_maker(design_file, 'orange')
	shared_string = ''
	for line in f:
		line = line.rstrip()
		line_list = line.split('\t')
		if line_list[1] == 'Group':
			#shared_string += line_list[1] + '\t' + '\t'.join(line_list[3:]) + '\t' + 'class' + '\n'
			shared_string += line_list[1] + '\t'
			for entities in line_list[3:]:
				#shared_string += ';'.join(entities.split(';')[1:])
				shared_string += entities
				shared_string += '\t'

			shared_string += 'class' + '\n'
		else:
			shared_string += line_list[1] + '\t' + '\t'.join(line_list[3:]) + '\t' + design_hash[line_list[1]] + '\n'
	orange = outputdir + name + '_orange_file.tsv'
	o = open(orange, 'w')
	o.write(shared_string)
	o.close()
	check_it_and_remove_it(relabund_container[0])
	return orange


def design_to_control(design_file, control_file):
	string = ''
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		string += line[0] + '\n'
	write_string_down(string, control_file)
	return True


def share_to_control(shared_file, control_file):
	string = ''
	f = open(shared_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[1].lower() == 'group':
			continue
		string += line[1] + '\n'
	write_string_down(string, control_file)
	return True


def kw_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_kruskal_wallis, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, design_file)
	if flag is False:
		print "Execution of mothur_kruskal_wallis failed!!!"
	else:
		kw_container = []
		extension_list = ['.' + name + '.kruskall_wallis']
		flag = scandirs(outputdir, kw_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	kw_dict = {}
	kw_list = []
	header = []
	columns = {}
	header, columns = mothur_result_parser(kw_container[0])
	this_list = list_string_to_float(columns[header[2]])
	if this_list is False:
		empty_dict = {}
		empty_list = []
		print "EMPTY LIST DETECTED"
		return (empty_dict, empty_list)
	maxl, minl, avel, medl = list_stat(this_list)
	flag, pw = significant_value(this_list, minl)
	if flag is False:
		print "Set significant limit to ", minl
		sig_value = minl
	else:
		print "Set significant limit to ", pw
		sig_value = pw
	f = open(kw_container[0], 'rU')
	
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if isfloat(line[2]) is False:
			continue
		if float(line[2]) < sig_value:
			kw_dict[line[0]] = (line[1], line[2])
			kw_list.append(line[0])
	#check_it_and_remove_it(kw_container[0])
	return (kw_dict, kw_list)


def indicator_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_indicator, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, design_file)
	if flag is False:
		print "Execution of mothur_kruskal_wallis failed!!!"
	else:
		indicator_container = []
		extension_list = ['.indicator.summary']
		flag = scandirs(outputdir, indicator_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	indicator_dict = {}
	indicator_list = []
	header = []
	columns = {}
	header, columns = mothur_result_parser(indicator_container[0])
	this_list = list_string_to_float(columns[header[3]])
	if this_list is False:
		empty_dict = {}
		empty_list = []
		return (empty_dict, empty_list)
	maxl, minl, avel, medl = list_stat(this_list)
	flag, pw = significant_value(this_list, minl)
	if flag is False:
		print "Set significant limit to ", minl
		sig_value = minl
	else:
		print "Set significant limit to ", pw
		sig_value = pw
	f = open(indicator_container[0], 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if '<0.001000' == line[3]:
			indicator_dict[line[0]] = (line[2], '0.001')
			indicator_list.append(line[0])
		if isfloat(line[3]) is False:
			continue
		if float(line[3]) < sig_value:
			indicator_dict[line[0]] = (line[2], line[3])
			indicator_list.append(line[0])
	#check_it_and_remove_it(indicator_container[0])
	return (indicator_dict, indicator_list)


def metastats_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_metastats, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, design_file)
	if flag is False:
		print "Execution of mothur_kruskal_wallis failed!!!"
	else:
		metastats_container = []
		extension_list = ['.metastats.logfile']
		flag = scandirs(outputdir, metastats_container, extension_list, 'ex_partial')
		if flag is False:
			print "This extension is not availble: ", extension_list
			
		else:
			for each_log in metastats_container:
				check_it_and_remove_it(each_log)

		metastats_container = []
		extension_list = ['.metastats']
		flag = scandirs(outputdir, metastats_container, extension_list, 'ex_partial')
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	metastats_dict = {}
	metastats_list = []
	for each_metastats_file in metastats_container:
		path, absname, ext = split_file_name(each_metastats_file)
		compare_name = ext.split('.')[2]
		new_metastats = outputdir + absname + '_new_metastats' + ext
		flag = remove_line_from_file(each_metastats_file, 6, new_metastats)
		header = []
		columns = {}
		header, columns = mothur_result_parser(new_metastats)
		this_list = list_string_to_float(columns[header[7]])
		if this_list is False:
			empty_dict = {}
			empty_list = []
			return (empty_dict, empty_list)
		new_list = []
		new_list = remove_value_from_list(this_list, '0.0')
		
		maxl, minl, avel, medl = list_stat(new_list)
		flag, pw = significant_value(new_list, minl)
		if flag is False:
			print "Set significant limit to ", minl
			sig_value = minl
		else:
			print "Set significant limit to ", pw
			sig_value = pw
		f = open(new_metastats, 'rU')
		for i in f:
			i = i.rstrip()
			line = i.split('\t')
			if len(line) < 7:
				continue
			elif line[7] in ['-', ' ', '', None]:
				continue
			elif isfloat(line[7]) is False:
				continue
			elif float(line[7]) < sig_value and float(line[2]) > 0.0 and float(line[5]) > 0.0:
				metastats_dict[line[0]] = (compare_name, line[7])
				metastats_list.append(line[0])
		check_it_and_remove_it(new_metastats)
		check_it_and_remove_it(each_metastats_file)
	return (metastats_dict, metastats_list)


def lefse_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_lefse, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, design_file)
	if flag is False:
		print "Execution of mothur_lefse failed!!!"
	else:
		lefse_container = []
		extension_list = ['.' + name + '.lefse_summary']
		flag = scandirs(outputdir, lefse_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			empty_dict = {}
			empty_list = []
			return (empty_dict, empty_list)

	lefse_dict = {}
	lefse_list = []
	header = []
	columns = {}
	header, columns = mothur_result_parser(lefse_container[0])
	this_list = list_string_to_float(columns[header[4]])
	if this_list is False:
		empty_dict = {}
		empty_list = []
		check_it_and_remove_it(lefse_container[0])
		return (empty_dict, empty_list)
	maxl, minl, avel, medl = list_stat(this_list)
	flag, pw = significant_value(this_list, minl)
	if flag is False:
		print "Set significant limit to ", minl
		sig_value = minl
	else:
		print "Set significant limit to ", pw
		sig_value = pw
	f = open(lefse_container[0], 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if len(line) < 5:
			continue
		elif line[4] in ['-', ' ', '', None]:
			continue
		elif isfloat(line[4]) is False:
			continue
		if float(line[4]) < sig_value:
			lefse_dict[line[0]] = (line[3], line[4])
			lefse_list.append(line[0])
	f.close()
	check_it_and_remove_it(lefse_container[0])
	return (lefse_dict, lefse_list)


def mothur_update(shared_file, design_file, control_file, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	f = open(control_file, 'rU')
	control_string = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_string += i + '-'
	f.close()
	#print control_data
	#print control_string
	f = open(design_file, 'rU')
	new_design_string = 'Sample\tDesign\n'
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[0] not in control_data:
			continue
		new_design_string += line[0] + '\t' + line[1] + '\n'
	f.close()
	controlled_design_file = outputdir + absname + '_controlled_design.txt'
	write_string_down(new_design_string, controlled_design_file)
	keep_groups = control_string[:-1]
	flag, stderr = execute_functions(mothur_get_groups, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, keep_groups)
	if flag is False:
		print "Execution of mothur_rarefaction_single failed!!!"
	else:
		control_container = []
		extension_list = ['.pick' + ext]
		flag = scandirs(outputdir, control_container, extension_list, 'ex_partial')
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	controlled_shared_file = outputdir + absname + '_updated' + ext
	os.rename(control_container[0], controlled_shared_file)
	return (controlled_shared_file, controlled_design_file)


def P3D_analyser(dist_file, design_file, calc, mothur_exec_path, processors, outputdir):
	flag, stderr = execute_functions(mothur_PCOA, processors, outputdir, 'multi', 'mothur', mothur_exec_path, dist_file)
	if flag is False:
		print "Execution of mothur_rarefaction_single failed!!!"
	else:
		pca_container = []
		extension_list = ['.pcoa.axes']
		flag = scandirs(outputdir, pca_container, extension_list)
		if flag is False:
			print "This extension is not availble: ", extension_list
			sys.exit(2)
	dict_of_design = {}
	pcoa_dict = {}
	dict_of_design = design_dict_maker(design_file)
	pcoafile = open(pca_container[0], 'rU')
	pcoa_javastring = ''
	xmax = ymax = zmax = 0
	xmin = ymin = zmin = 0
	for line in pcoafile:
		line = line.rstrip()
		line_split = line.split('\t')
		if line_split[0].lower() == 'group':
			continue
		elif line_split[0] in dict_of_design.keys():
			if len(line_split) < 4:
				zero_list = ['0.0', '0.0', '0.0']
				line_split.extend(zero_list)
			pcoa_dict[line_split[0]] = (dict_of_design[line_split[0]], '[' + line_split[1] + ',' + line_split[2] + ',' + line_split[3] + ']')
			if isfloat(line_split[1]) and float(line_split[1]) > xmax: xmax = float(line_split[1])
			if isfloat(line_split[2]) and float(line_split[2]) > ymax: ymax = float(line_split[2])
			if isfloat(line_split[3]) and float(line_split[3]) > zmax: zmax = float(line_split[3])
			if isfloat(line_split[1]) and float(line_split[1]) < xmin: xmin = float(line_split[1])
			if isfloat(line_split[2]) and float(line_split[2]) < ymin: ymin = float(line_split[2])
			if isfloat(line_split[3]) and float(line_split[3]) < zmin: zmin = float(line_split[3])
	
	dict_of_symbol_color = symbol_color_dict(list(set(dict_of_design.values())))
	for key, value in pcoa_dict.iteritems():
		pcoa_javastring += '\n\t\t\t{\n'
		pcoa_javastring += '\t\t\t\tname: "' + key + '[' + value[0] + ']",\n'
		pcoa_javastring += '\t\t\t\tdata: [' + value[1] + '],\n'
		pcoa_javastring += '\t\t\t\tcolor: "' + dict_of_symbol_color[value[0]][1] + '",\n'
		pcoa_javastring += '\t\t\t\tmarker:{symbol:"' + dict_of_symbol_color[value[0]][0] + '"}'
		pcoa_javastring += '\t\t\t},\n'
	
	return (pcoa_javastring, xmax, xmin, ymax, ymin, zmax, zmin)


def summary_table(shared_file, design_file):
	dict_of_shared, header_dict = shared_dict(shared_file)
	rev_design_dict = design_dict_maker(design_file, 'reverse')
	
	table_header = ''
	table_header += '\t\t\t\t\t\t<thead>\n'
	table_header += '\t\t\t\t\t\t\t<tr>\n'
	table_header += '\t\t\t\t\t\t\t\t<th class="text-center">State</th>\n'
	for key, value in header_dict.iteritems():
		if value in ['numOtus', 'label'] or 'lci' in value or 'hci' in value:
			continue
		elif value == 'nseqs':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Number of Sequences</th>\n'
		elif value == 'simpson':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Simpson(Diversity)</th>\n'
		elif value == 'shannon':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Shannon(Diversity)</th>\n'
		elif value == 'invsimpson':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Inverse Simpson(Diversity)</th>\n'
		elif value == 'npshannon':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">non-Parametric Shannon(Diversity)</th>\n'
		elif value == 'sobs':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Observed richness(Richness)</th>\n'
		elif value == 'jackknife':
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">Jackknife(Richness)</th>\n'
		else:
			table_header += '\t\t\t\t\t\t\t\t\t<th class="text-center">' + value + '</th>\n'
	header_keys = header_dict.keys()
	table_header += '\t\t\t\t\t\t\t</tr>\n'
	table_header += '\t\t\t\t\t\t</thead>\n'
	table_body = ''
	table_body += '\t\t\t\t\t<tbody align="center">\n'

	for key, value_list in rev_design_dict.iteritems():
		if key in ['class', 'sample', 'treatment']:
			continue
		else:
			#table_body += '\t\t\t\t<th rowspan="' + str(len(value_list)+1) + '" scope="rowgroup" class="metadata">' + key + '</th>\n'
			for each_value in value_list:

				list_of_values_of_each_value = dict_of_shared[each_value]
				"""[('label', 'zac'), ('group', 'Ademona1'), ('nseqs', '74568.000000'), ('simpson', '0.062251'), ('simpson_lci', '0.061383'), ('simpson_hci', '0.063118'), ('invsimpson', '16.064047'), ('invsimpson_lci', '15.843278'), ('invsimpson_hci', '16.291056'), ('shannon', '3.725330'), ('shannon_lci', '3.712430'), ('shannon_hci', '3.738230'), ('npshannon', '3.736071'), ('sobs', '534.000000'), ('jackknife', '628.000000'), ('jackknife_lci', '601.125834'), ('jackknife_hci', '654.874166')]"""

				table_body += '\t\t\t\t\t\t\t\t\t\t<tr name="div_alpha_' + key + '">\n'
				table_body += '\t\t\t\t\t\t\t\t\t\t\t<td>' + key + '</td>\n'
				for each_key in header_keys:
					if each_key == -1 or 'lci' in header_dict[each_key] or 'hci' in header_dict[each_key] or 'label' in header_dict[each_key]:
						continue
					else:
						table_body += '\t\t\t\t\t\t\t\t\t\t\t<td>' + round_float(list_of_values_of_each_value[each_key][1]) + '</td>\n'
				table_body += '\t\t\t\t\t\t\t\t\t\t</tr>\n'
	table_body += '\t\t\t\t\t\t\t\t\t</tbody>\n'
	
	return (table_header, table_body)


def abundance_table(shared_file, design_file, top, bottom, outputdir):
	min_sh, max_sh, ave_sh, median_sh = shared_max_min(shared_file)
	new_shared = outputdir + 'abundance_new_shared.txt'
	if bottom is not False or top is not False:
		filter_shared(shared_file, top, bottom, new_shared)
		shared_file = new_shared
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)
	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	abundance_category_dict = {}
	abundance_columns_dict = {}

	sample_name_list = shared_dict.keys()
	otu_name_list = shared_dict[sample_name_list[0]].keys()
	#print sample_name_list
	#print otu_name_list
	#colors_list = ['#7cb5ec', '#434348', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#5645b3', '#a23456', '#babbab', '#fffdd1', '#5645b3', '#800000', '#a23456', '#ba3ba3', '#00bb33', '#1b2d45', '#05377b', '#e5e5e5', '#021631', '#002448', '#042c62', '#032149', '#949494', '#123456', '#696969', '#444444', '#c0ffee', '#098765', '#567890', '#654321', '#123456', '#a23456', '#bababa', '#babbab', '#ba3ba3', '#baebae', '#b3b3b3', '#dbb234', '#b3a271', '#0b03b3', '#0303bb', '#00bb33', '#3300bb', '#0033bb', '#bb0033', '#cc0033', '#d3ad33']
	for each_design in reverse_design_dict:
		if each_design not in abundance_category_dict:
			abundance_category_dict[each_design] = list_to_string(reverse_design_dict[each_design], "','")
	#print abundance_category_dict
	#print abundance_category_dict
	for each_design in reverse_design_dict:
		sample_list = reverse_design_dict[each_design]
		if each_design not in abundance_columns_dict:
			abundance_columns_dict[each_design] = ''
		for each_otu in otu_name_list:
			otu_json_string = ''
			otu_temp = each_otu.split(';')
			if len(otu_temp) > 1:

				otu = otu_temp[0] + ';' + otu_temp[-3] + ';' + otu_temp[-2]
			else:
				otu = otu_temp[0]
			otu_json_string += '{\n' + 'name:"' + otu + '",\n'
			otu_json_string += 'data: ['
			for each_sample in sample_list:
				otu_json_string += shared_dict[each_sample][each_otu] + ','
			otu_json_string += ']\n},'
			abundance_columns_dict[each_design] += otu_json_string
	#print abundance_columns_dict['Cancer']
	#sys.exit(2)
	return(abundance_columns_dict, abundance_category_dict)


def bacterial_abundance_table(shared_file, design_file, top, bottom, outputdir):
	new_shared = outputdir + 'collapsed.txt'
	collapse_bacterial_otu(shared_file, new_shared)
	shared_file = new_shared
	min_sh, max_sh, ave_sh, median_sh = shared_max_min(shared_file)
	new_shared = outputdir + 'bacterial_new_shared.txt'
	if bottom is not False or top is not False:
		filter_shared(shared_file, top, bottom, new_shared)
		shared_file = new_shared
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	sample_name_list = shared_dict.keys()
	
	bacterial_list = shared_dict[sample_name_list[0]].keys()
	#bacterial_list = sort_otus(shared_file)
	bact_list = []
	for each_bact in bacterial_list:
		
		bact = each_bact.split(';')
		if len(bact) > 1:
			bact_list.append(bact[-2])
		else:
			bact_list.append(bact[0])
	
	bacterial_abundance_category_dict = {}
	#colors_list = ['#7cb5ec', '#434348', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#5645b3', '#a23456', '#babbab', '#fffdd1', '#5645b3', '#800000', '#a23456', '#ba3ba3', '#00bb33', '#1b2d45', '#05377b', '#e5e5e5', '#021631', '#002448', '#042c62', '#032149', '#949494', '#123456', '#696969', '#444444', '#c0ffee', '#098765', '#567890', '#654321', '#123456', '#a23456', '#bababa', '#babbab', '#ba3ba3', '#baebae', '#b3b3b3', '#dbb234', '#b3a271', '#0b03b3', '#0303bb', '#00bb33', '#3300bb', '#0033bb', '#bb0033', '#cc0033', '#d3ad33']
	for each_design in reverse_design_dict:
		if each_design not in bacterial_abundance_category_dict:
			bacterial_abundance_category_dict[each_design] = list_to_string(bact_list, "','")
	#print bacterial_abundance_category_dict
	bacterial_abundance_columns_dict = {}
	for each_design in reverse_design_dict:
		sample_list = reverse_design_dict[each_design]
		if each_design not in bacterial_abundance_columns_dict:
			bacterial_abundance_columns_dict[each_design] = ''
		for each_sample in sample_list:
			#print each_sample
			otu_json_string = ''
			otu_json_string += '{\n' + 'name:"' + each_sample + '",\n'
			otu_json_string += 'data: ['
			#print list_to_string(shared_dict[each_sample].values(), ',')
			#otu_json_string += list_to_string(list_float_to_string(sorted(list_string_to_float(shared_dict[each_sample].values()), reverse=True)), ',')
			otu_json_string += list_to_string(shared_dict[each_sample].values(), ',')
			otu_json_string += ']\n},'
			bacterial_abundance_columns_dict[each_design] += otu_json_string

	#print bacterial_abundance_columns_dict
	
	return(bacterial_abundance_columns_dict, bacterial_abundance_category_dict)


def heatmap_analysis_plotly(shared_file, design_file, mothur_exec_path, processors, outputdir):
	new_shared_file = outputdir + 'temp.shared'
	collapse_bacterial_otu(shared_file, new_shared_file)
	reduced_shared_dict, otu_list = shared_to_design_reduction(new_shared_file, design_file, mothur_exec_path, processors, outputdir)
	
	heatmap_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########     SCATTERPLOT PLOTLY 				      #########
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:
		heatmap_javascript_string += 'var data_heatmap_' + each_design + ' = [' + list_to_string(reduced_shared_dict[each_design], ',') + '];\n'
	
	print heatmap_javascript_string
	sys.exit(2)
	'''
	pca_javascript_string += pca_string
	pca_javascript_string += "\tvar data_pca_" + PCA_name + '_' + calc_name + " = [" + pca_sample_string + "];\n"
	
	pca_javascript_string += """
	var layout_pca_""" + PCA_name + '_' + calc_name + """ = {
		title: 'PCA - """ + calc + """ ',
		titlefont: {
			family: 'Avenir',
			size: 20
		},
		//showlegend: """ + legend_status + """,
		showlegend: false,
		traceorder:'normal',
		autosize: false,
		height: 500,
		width: 1000,
		dragmode: false,
		
		margin: {
			l: 25,
			r: 25,
			b: 50,
			t: 50,
			pad: 1
		},
		scene: {
			xaxis: {
				title: 'PC1',
				titlefont: {
					family: 'Avenir',
					size: 15 },
				range:[""" + xmin + """, """ + xmax + """],
				autorange: false,
				showbackground: true,
				nticks:10,
				autotick:true,
				fixedrange: false,
				gridcolor: "white",
				linecolor: "white",
				zerolinecolor: "white",
				gridwidth:4,
				backgroundcolor: "rgba(0,0,0,0.02)",
				showspikes: false,
			},
			yaxis: {
				title: 'PC2',
				titlefont: {
					family: 'Avenir',
					size: 15
				},
				range:[""" + ymin + """, """ + ymax + """],
				showbackground: true,
				fixedrange: true,
				nticks:5,
				autorange: false,
				gridcolor: "white",
				linecolor: "white",
				gridwidth:4,
				zerolinecolor: "white",
				backgroundcolor: "rgba(0,0,0,0.03)",
				showspikes: false,
				type:"linear"
			},
			zaxis: {
				title: 'PC3',
				titlefont: {
					family: 'Avenir',
					size: 15},
				range:[""" + zmin + """, """ + zmax + """],
				showbackground: true,
				fixedrange: true,
				gridcolor: "white",
				linecolor: "white",
				zerolinecolor: "white",
				gridwidth:4,
				backgroundcolor: "rgba(0,0,0,0.04)",
				nticks:5,
				autorange: false,
				showspikes: false,
			},
			//cameraposition:[[0.4, 0.5, -0.3, -0.4], [0.1, 0, -0.1], 4],
			aspectratio:{x:1,y:1,z:1},
			aspectmode:'manual'
		}

	};
	"""
	pca_javascript_string += """\n
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	
	f = open(js_path + absname + '_' + PCA_name + '_' + calc_name + '_pca_plotly.js', 'w')
	f.write(pca_javascript_string)
	f.close()
	'''
	return True


def shared_to_design_reduction(shared_file, design_file, mothur_exec_path, processors, outputdir):
	reduced_shared_dict = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	print design_dictionary_rev
	for each_design in design_dictionary_rev.keys():
		sample_list = design_dictionary_rev[each_design]
		#print sample_list
		selected_shared_list = []
		f = open(shared_file, 'rU')
		for i in f:
			i = i.rstrip()
			line = i.split('\t')
			if line[1].lower() == 'group':
				otu_list = line[3:]
			elif line[1] in sample_list:
				abundance_list = line[3:]

				abundance_float_list = [float(x) for x in abundance_list]
				#print abundance_float_list
				if len(selected_shared_list) == 0:
					selected_shared_list = [0.0] * len(abundance_float_list)
				#print selected_shared_list
				temp = [x + y for x, y in zip(selected_shared_list, abundance_float_list)]
				selected_shared_list = temp
				#print selected_shared_list
		#print selected_shared_list
		f.close()
		reduced_shared_dict[each_design] = [str(x) for x in selected_shared_list]
	return (reduced_shared_dict, otu_list)


# ################################### MOTHUR  ################################## #


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


def mothur_rarefaction_single(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, freq_value):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'rarefaction.single'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'iters=100, calc=chao')
	#parameter_list.append(',' + space + 'freq=' + str(freq_value))
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_summary_single(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'summary.single'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'size=100000, iters=10000, calc=nseqs-simpson-invsimpson-shannon-npshannon-sobs-jack')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_distance_shared(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, calc):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'dist.shared'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'calc=' + calc)
	parameter_list.append(',' + space + 'iters=10000')
	parameter_list.append(',' + space + 'output=square')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_PCOA(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, dist_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'pcoa'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('phylip=' + dist_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_normalize(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, method):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'normalize.shared'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'method=' + method)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_remove_groups(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, removing_groups):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'remove.groups'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'groups=' + removing_groups)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_groups(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, removing_groups):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.groups'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'groups=' + removing_groups)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_kruskal_wallis(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, design_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'kruskal.wallis'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'design=' + design_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_indicator(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, design_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'indicator'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'design=' + design_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_lefse(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, design_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'lefse'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'design=' + design_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_metastats(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, design_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'metastats'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'design=' + design_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_relabund(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.relabund'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'scale=totalgroup')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_subsample(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'sub.sample'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


# ################################### UTILITY ################################## #
def copytree(src, dst, symlinks=False, ignore=None):
	for item in os.listdir(src):
		s = os.path.join(src, item)
		d = os.path.join(dst, item)
		if os.path.isdir(s):
			shutil.copytree(s, d, symlinks, ignore)
		else:
			shutil.copy2(s, d)
	return True
def min_abundance(shared_file, new_shared_file, minimum_abundance):
	header = []
	columns = {}
	abundance_list = []
	OTU_dict = {}
	header, columns = mothur_result_parser(shared_file)
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		else:
			abundance_list = list_string_to_float(columns[head])
			if sum(abundance_list) < minimum_abundance:
				print "LOW ABUNDANCE OTU", head
				continue
			else:
				OTU_dict[head] = columns[head]

	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	return True


def permutation_of_list(permute_list, size):
	return list(itertools.permutations(permute_list, int(size)))


def transpose_shared_table(shared_file):
	shared_dict = parse_shared_file(shared_file)
	shared_header_list = ['OTU']
	shared_header_list.extend(shared_dict.keys())
	sample_name_list = shared_dict.keys()
	bacterial_dict = {}
	for each_sample in sample_name_list:
		for bacts, values in shared_dict[each_sample].iteritems():
			if bacts not in bacterial_dict:
				bacterial_dict[bacts] = [values]
			else:
				bacterial_dict[bacts].append(values)
	shared_body_list = []
	for key, values in bacterial_dict.iteritems():
		shared_body_list.append(key)
		shared_body_list.extend(values)
	print shared_header_list
	print shared_body_list


def transpose(one_list):
	trans_data = list(itertools.izip_longest(*one_list, fillvalue=''))
	return trans_data


def remove_mothur_log(outputdir):
	# scan directory for files with logfile extension and remove it
	template = ['.logfile', '.rabund']
	list_of_files = []
	flag = scandirs(outputdir, list_of_files, template, 'ex_partial')
	if flag is True:
		for i in list_of_files:
			check_it_and_remove_it(i)
	else:
		pass


def sort_otus(shared_file):
	header = []
	columns = {}
	header, columns = mothur_result_parser(shared_file)
	sorted_shared_dict = {}

	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		sum_value = sum(list_string_to_float(columns[head]))
		sorted_shared_dict[head] = sum_value
	#then we sort them
	sorted_list = sorted(sorted_shared_dict, key=sorted_shared_dict.get, reverse=True)

	return sorted_list


def collapse_otus(otu_list, collapsed_otu_list):
	uniq_otu = []
	for otu in otu_list:
		split_otu = otu.split(';')
		if len(split_otu) > 3:
			if split_otu[-2] in uniq_otu:
				continue
			else:
				uniq_otu.append(split_otu[-2])
				collapsed_otu_list.append(otu)
		else:
			collapsed_otu_list.append(otu)
		
	return True


def collapse_bacterial_otu(shared_file, new_shared_file):
	header = []
	columns = {}
	
	header, columns = mothur_result_parser(shared_file)
	
	bact_dict = {}
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		bact_split = head.split(';')
		if len(bact_split) > 1:
			bact_name = bact_split[-2]
		else:
			bact_name = bact_split[0]
		if bact_name not in bact_dict:
			bact_dict[bact_name] = [head]
		else:
			bact_dict[bact_name] = [sum(x) for x in zip(list_string_to_float(bact_dict[bact_name]), list_string_to_float(columns[head]))]
	numotus = str(len(bact_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + bact_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in bact_dict.keys():
			new_shared_string += str(bact_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'

	write_string_down(new_shared_string, new_shared_file)
	return True


def remove_value_from_list(the_list, the_value):
	new_list = []
	for i in the_list:
		if i == the_value:
			continue
		new_list.append(i)
	return new_list


def remove_line_from_file(the_file, line_number, new_file):
	string = ''
	f = open(the_file, 'rU')
	for idx, i in enumerate(f):
		if idx < line_number:
			continue
		string += i
	write_string_down(string, new_file)
	return True


def significant_value(list_float, minimum_value):
	sig_list = []
	pw = 0.05
	if pw > minimum_value:
		return (True, pw)
	while pw > minimum_value:
		for i in list_float:
			if i < pw:
				sig_list.append(i)
		if len(sig_list) > 0:
			break
		pw += 0.01
	if len(sig_list) > 0:
		return (True, pw)
	else:
		return (False, 0.0)


def list_stat(this_list):
	maxl = max(this_list)
	minl = min(this_list)
	avel = sum(this_list) / len(this_list)
	medl = median(this_list)
	return (maxl, minl, avel, medl)


def mothur_shared_ranking(shared, ranking_list, new_shared_file):
	if len(ranking_list) < 1:
		return shared
	unique_ranking_list = list(set(ranking_list))
	
	if len(unique_ranking_list) < 1:
		print "your ranking list is empty"
		return False
	label_list = ['label', 'Group', 'numOtus'] + unique_ranking_list
	shared_file = open(shared, 'rU')
	reader = csv.DictReader(shared_file, delimiter='\t')
	new_shared = ''
	new_shared = list_to_string(label_list)
	new_shared += '\n'
	for row in reader:
		for key in label_list:
			if key == 'numOtus':
				new_shared += str(len(unique_ranking_list)) + '\t'
			elif key in row:
				new_shared += row[key] + '\t'
		new_shared += '\n'
	write_string_down(new_shared, new_shared_file)
	return True


def mothur_remove_otu(shared_file, otu_control, new_shared_file):
	otu_list = []
	f = open(otu_control, 'rU')
	for i in f:
		i = i.rstrip()
		if i != '':
			otu_list.append(i + ';')
	otu_list = list(set(otu_list))
	print "Removing Following OTUs:", otu_list
	header = []
	columns = {}
	
	header, columns = mothur_result_parser(shared_file)
	bact_dict = {}
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		flag = True
		for otu in otu_list:
			if otu in head:
				flag = False
		if flag:
			bact_dict[head] = columns[head]

	numotus = str(len(bact_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Group', 'numOtus', ] + bact_dict.keys())
	new_shared_string += '\n'
	
	for i in range(len(columns['label'])):
		
		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Group'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in bact_dict.keys():
			new_shared_string += str(bact_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'

	write_string_down(new_shared_string, new_shared_file)
	return True


def write_string_down(new_string, new_name):
	f = open(new_name, 'w')
	f.write(new_string)
	f.close()
	return new_name


def check_it_and_remove_it(filename, noreport=False):
	try:
		os.remove(filename)
		if noreport is False:
			pass
	except OSError:
		pass


def filter_shared(shared_file, top_value, bottom_value, new_shared):
	header = []
	columns = {}
	header, columns = mothur_result_parser(shared_file)
	ok_list = []
	bottom_list = []
	top_list = []
	sum_value = 0.0
	
	#numOtus = len(list_string_to_float(columns['numOtus']))
	#print numOtus
	for head in header:
		if head.lower() in ['group', 'numotus', 'label']:
			continue
		sum_value = sum(list_string_to_float(columns[head]))

		if bottom_value is not False and sum_value >= bottom_value:
			bottom_list.append(head)
		if top_value is not False and sum_value <= top_value:
			top_list.append(head)
			#print "less_than_shared_ave"
			#print sum_value
			#print head
			#print columns[head]
	ok_list = intersect(top_list, bottom_list)
	#print len(ok_list)
	#print len(header)
	mothur_shared_ranking(shared_file, ok_list, new_shared)
	return True


def intersect(a, b):
	if len(a) == 0:
		return list(set(b))
	elif len(b) == 0:
		return list(set(a))
	return list(set(a) & set(b))


def list_string_to_float(list_string):
	new_list = []
	for i in list_string:
		if i in ['-', ' ', '', None]:
			continue
		if isfloat(i) is False:
			continue
		else:
			new_list.append(float(i))
	if len(new_list) < 1:
		print "list does not have any digits value"
		return False
	return new_list


def list_float_to_string(float_string):
	list_string = [str(i) for i in float_string]
	return list_string


def round_float(value):
	if isfloat(value) is False:
	
		return value
	else:
		return str(round(decimal.Decimal(value), 3))


def shared_dict(shared_file):
	header_dict = {}
	dict_of_shared = {}
	delimiter = sniff(shared_file)
	f = open(shared_file, 'rU')
	for line in f:
		line = line.rstrip()
		line_list = line.split(delimiter)
		#print line_list
		if line_list[0] == 'label':
			for each_cell in line_list:
				header_dict[line_list.index(each_cell)] = each_cell
		else:
			dict_of_shared[line_list[1]] = []
			for each_cell in line_list:
				
				dict_of_shared[line_list[1]].append((header_dict[line_list.index(each_cell)], each_cell))
	
	return (dict_of_shared, header_dict)


def sniff(template_file):

	template = open(template_file, 'rU')
	try:
		dialect = csv.Sniffer().sniff(template.readline(), [',', '\t'])
	except csv.Error:
		print "file with bad delimiter, please use comma or tab as delimiter "
		template.close()
		sys.exit(2)
	template.close()
	return dialect.delimiter


def parse_shared_file(shared_file):
	header = []
	columns = {}
	header, columns = mothur_result_parser(shared_file)
	shared_dict = {}
	for idx, each_sample in enumerate(columns['Group']):
		if columns['Group'][idx] not in shared_dict:
			shared_dict[each_sample] = {}
	for head in header:
		if head in ['label', 'numOtus', 'Group']:
			continue
		for idx, each_value in enumerate(columns[head]):
			sample_name = columns['Group'][idx]
			shared_dict[sample_name][head] = each_value

	return shared_dict


def random_color():
	golden_ratio_conjugate = 0.618033988749895
	randvalue = random.randint(1, 256)
	randvalue += golden_ratio_conjugate
	randvalue %= 1
	#hsv_to_rgb(randvalue, 0.5, 0.95)
	color = "#%06x" % random.randint(1, 0xFFFFFE)
	#print color
	return color


def color_generate(color_count):
	max_value = 16581375  #255**3
	interval = int(max_value / color_count)
	colors_list = [hex(I)[2:].zfill(6) for I in range(20, max_value, interval)]
	return colors_list


def hsv_color_generate(color_count):
	HSV_tuples = [(x * 1.0 / color_count, 0.5, 0.5) for x in xrange(color_count)]
	hex_out = []
	for rgb in HSV_tuples:
		rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
		hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))

	return hex_out


def symbol_color_dict(target_list, mode=None):
	print len(target_list)
	dict_of_symbol_color = {}
	#symbol_list = ["circle", "square", "triangle", "diamond", "triangle-down", "cross", "circle", "square", "diamond", "triangle", "triangle-down", "circle", "square", "diamond", "triangle", "triangle-down", "circle", "square", "diamond", "triangle", "triangle-down", "circle", "square", "diamond", "triangle", "triangle-down", "circle", "square", "diamond", "triangle", "triangle-down", "circle", "square", "diamond", "triangle", "triangle-down"]
	symbol_list = ["circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x", "circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x", "circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x", "circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x", "circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x", "circle", "circle-open", "square", "square-open", "diamond", "diamond-open", "cross", "x"]
	#symbol_list = ["circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "triangle-left", "triangle-right", "triangle-ne", "triangle-se", "triangle-sw"]
	#color_list = ["#000000", "#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE", "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#FFEEE8", "#774D00", "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D", "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400", "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F", "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231", "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C", "#DEFF74", "#00FFC6", "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78", "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167", "#009BFF", "#E85EBE"]
	#color_list = ["#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80","#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100","#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F","#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09","#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66","#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C","#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81","#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00","#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700","#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329","#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]
	#colors_list = [ "#e0eeee", "#c1cdcd", "#838b8b", "#f5f5dc", "#ffe4c4", "#eed5b7", "#cdb79e", "#8b7d6b", "#000000", "#ffebcd", "#0000ff", "#0000ee", "#00008b", "#8a2be2", "#a52a2a", "#ff4040", "#ee3b3b", "#cd3333", "#8b2323", "#deb887", "#ffd39b", "#eec591", "#cdaa7d", "#8b7355", "#5f9ea0", "#98f5ff", "#8ee5ee", "#7ac5cd", "#53868b", "#7fff00", "#76ee00", "#66cd00", "#458b00", "#d2691e", "#ff7f24", "#ee7621", "#cd661d", "#ff7f50", "#ff7256", "#ee6a50", "#cd5b45", "#8b3e2f", "#6495ed", "#fff8dc", "#eee8cd", "#cdc8b1", "#8b8878", "#00ffff", "#00eeee", "#00cdcd", "#008b8b", "#b8860b", "#ffb90f", "#eead0e", "#cd950c", "#8b6508", "#006400", "#bdb76b", "#556b2f", "#caff70", "#bcee68", "#a2cd5a", "#6e8b3d", "#ff8c00", "#ff7f00", "#ee7600", "#cd6600", "#8b4500", "#9932cc", "#bf3eff", "#b23aee", "#9a32cd", "#68228b", "#e9967a", "#8fbc8f", "#c1ffc1", "#b4eeb4", "#9bcd9b", "#698b69", "#483d8b", "#2f4f4f", "#97ffff", "#8deeee", "#79cdcd", "#528b8b", "#00ced1", "#9400d3", "#ff1493", "#ee1289", "#cd1076", "#8b0a50", "#00bfff", "#00b2ee", "#009acd", "#00688b", "#696969", "#1e90ff", "#1c86ee", "#1874cd", "#104e8b", "#b22222", "#ff3030", "#ee2c2c", "#cd2626", "#8b1a1a", "#fffaf0", "#228b22", "#dcdcdc", "#f8f8ff", "#ffd700", "#eec900", "#cdad00", "#8b7500", "#daa520", "#ffc125", "#eeb422", "#cd9b1d", "#8b6914", "#bebebe", "#030303", "#1a1a1a", "#1c1c1c", "#1f1f1f", "#212121", "#242424", "#262626", "#292929", "#2b2b2b", "#2e2e2e", "#303030", "#050505", "#333333", "#363636", "#383838", "#3b3b3b", "#3d3d3d", "#404040", "#424242", "#454545", "#474747", "#4a4a4a", "#080808", "#4d4d4d", "#4f4f4f", "#525252", "#545454", "#575757", "#595959", "#5c5c5c", "#5e5e5e", "#616161", "#636363", "#0a0a0a", "#666666", "#696969", "#6b6b6b", "#6e6e6e", "#707070", "#737373", "#757575", "#787878", "#7a7a7a", "#7d7d7d", "#0d0d0d", "#7f7f7f", "#828282", "#858585", "#878787", "#8a8a8a", "#8c8c8c", "#8f8f8f", "#919191", "#949494", "#969696", "#0f0f0f", "#999999", "#9c9c9c", "#9e9e9e", "#a1a1a1", "#a3a3a3", "#a6a6a6", "#a8a8a8", "#ababab", "#adadad", "#b0b0b0", "#121212", "#b3b3b3", "#b5b5b5", "#b8b8b8", "#bababa", "#bdbdbd", "#bfbfbf", "#c2c2c2", "#c4c4c4", "#c7c7c7", "#c9c9c9", "#141414", "#cccccc", "#cfcfcf", "#d1d1d1", "#d4d4d4", "#d6d6d6", "#d9d9d9", "#dbdbdb", "#dedede", "#e0e0e0", "#e3e3e3", "#171717", "#e5e5e5", "#e8e8e8", "#ebebeb", "#ededed", "#f0f0f0", "#f2f2f2", "#f7f7f7", "#fafafa", "#fcfcfc", "#00ff00", "#00ee00", "#00cd00", "#008b00", "#adff2f", "#f0fff0", "#e0eee0", "#c1cdc1", "#838b83", "#ff69b4", "#ff6eb4", "#ee6aa7", "#cd6090", "#8b3a62", "#cd5c5c", "#ff6a6a", "#ee6363", "#cd5555", "#8b3a3a", "#fffff0", "#eeeee0", "#cdcdc1", "#8b8b83", "#f0e68c", "#fff68f", "#eee685", "#cdc673", "#8b864e", "#e6e6fa", "#fff0f5", "#eee0e5", "#cdc1c5", "#8b8386", "#7cfc00", "#fffacd", "#eee9bf", "#cdc9a5", "#8b8970", "#eedd82", "#add8e6", "#bfefff", "#b2dfee", "#9ac0cd", "#68838b", "#f08080", "#e0ffff", "#d1eeee", "#b4cdcd", "#7a8b8b", "#ffec8b", "#eedc82", "#cdbe70", "#8b814c", "#fafad2", "#d3d3d3", "#ffb6c1", "#ffaeb9", "#eea2ad", "#cd8c95", "#8b5f65", "#ffa07a", "#ee9572", "#cd8162", "#8b5742", "#20b2aa", "#87cefa", "#b0e2ff", "#a4d3ee", "#8db6cd", "#607b8b", "#8470ff", "#778899", "#b0c4de", "#cae1ff", "#bcd2ee", "#a2b5cd", "#6e7b8b", "#ffffe0", "#eeeed1", "#cdcdb4", "#8b8b7a", "#32cd32", "#faf0e6", "#ff00ff", "#ee00ee", "#cd00cd", "#8b008b", "#b03060", "#ff34b3", "#ee30a7", "#cd2990", "#8b1c62", "#66cdaa", "#66cdaa", "#0000cd", "#ba55d3", "#e066ff", "#d15fee", "#b452cd", "#7a378b", "#9370db", "#ab82ff", "#9f79ee", "#8968cd", "#5d478b", "#3cb371", "#7b68ee", "#00fa9a", "#48d1cc", "#c71585", "#191970", "#f5fffa", "#ffe4e1", "#eed5d2", "#cdb7b5", "#8b7d7b", "#ffe4b5", "#ffdead", "#eecfa1", "#cdb38b", "#8b795e", "#000080", "#fdf5e6", "#6b8e23", "#c0ff3e", "#b3ee3a", "#698b22", "#ffa500", "#ee9a00", "#cd8500", "#8b5a00", "#ff4500", "#ee4000", "#cd3700", "#8b2500", "#da70d6", "#ff83fa", "#ee7ae9", "#cd69c9", "#8b4789", "#db7093", "#eee8aa", "#98fb98", "#9aff9a", "#90ee90", "#7ccd7c", "#548b54", "#afeeee", "#bbffff", "#aeeeee", "#96cdcd", "#668b8b", "#db7093", "#ff82ab", "#ee799f", "#cd6889", "#8b475d", "#ffefd5", "#ffdab9", "#eecbad", "#cdaf95", "#8b7765", "#ffc0cb", "#ffb5c5", "#eea9b8", "#cd919e", "#8b636c", "#dda0dd", "#ffbbff", "#eeaeee", "#cd96cd", "#8b668b", "#b0e0e6", "#a020f0", "#663399", "#9b30ff", "#912cee", "#7d26cd", "#551a8b", "#ff0000", "#ee0000", "#cd0000", "#8b0000", "#bc8f8f", "#ffc1c1", "#eeb4b4", "#cd9b9b", "#8b6969", "#4169e1", "#4876ff", "#436eee", "#3a5fcd", "#27408b", "#8b4513", "#fa8072", "#ff8c69", "#ee8262", "#cd7054", "#8b4c39", "#f4a460", "#54ff9f", "#4eee94", "#43cd80", "#2e8b57", "#fff5ee", "#eee5de", "#cdc5bf", "#8b8682", "#a0522d", "#ff8247", "#ee7942", "#cd6839", "#8b4726", "#87ceeb", "#87ceff", "#7ec0ee", "#6ca6cd", "#4a708b", "#6a5acd", "#836fff", "#7a67ee", "#6959cd", "#473c8b", "#708090", "#c6e2ff", "#b9d3ee", "#9fb6cd", "#6c7b8b", "#fffafa", "#eee9e9", "#cdc9c9", "#8b8989", "#00ff7f", "#00ee76", "#00cd66", "#008b45", "#4682b4", "#63b8ff", "#5cacee", "#4f94cd", "#36648b", "#d2b48c", "#ffa54f", "#ee9a49", "#cd853f", "#8b5a2b", "#d8bfd8", "#ffe1ff", "#eed2ee", "#cdb5cd", "#8b7b8b", "#ff6347", "#ee5c42", "#cd4f39", "#8b3626", "#40e0d0", "#00f5ff", "#00e5ee", "#00c5cd", "#00868b", "#ee82ee", "#d02090", "#ff3e96", "#ee3a8c", "#cd3278", "#8b2252", "#f5deb3", "#ffe7ba", "#eed8ae", "#cdba96", "#8b7e66", "#ffffff", "#f5f5f5", "#ffff00", "#eeee00", "#cdcd00", "#8b8b00", "#9acd32", "#B0171F", "#DC143C", "#FFB6C1", "#FFAEB9", "#EEA2AD", "#CD8C95", "#8B5F65", "#FFC0CB", "#FFB5C5", "#EEA9B8", "#CD919E", "#8B636C", "#DB7093", "#FF82AB", "#EE799F", "#CD6889", "#8B475D", "#FFF0F5", "#EEE0E5", "#CDC1C5", "#8B8386", "#FF3E96", "#EE3A8C", "#CD3278", "#8B2252", "#FF69B4", "#FF6EB4", "#EE6AA7", "#CD6090", "#8B3A62", "#872657", "#FF1493", "#EE1289", "#CD1076", "#8B0A50", "#FF34B3", "#EE30A7", "#CD2990", "#8B1C62", "#C71585", "#D02090", "#DA70D6", "#FF83FA", "#EE7AE9", "#CD69C9", "#8B4789", "#D8BFD8", "#FFE1FF", "#EED2EE", "#CDB5CD", "#8B7B8B", "#FFBBFF", "#EEAEEE", "#CD96CD", "#8B668B", "#DDA0DD", "#EE82EE", "#FF00FF", "#EE00EE", "#CD00CD", "#8B008B", "#800080", "#BA55D3", "#E066FF", "#D15FEE", "#B452CD", "#7A378B", "#9400D3", "#9932CC", "#BF3EFF", "#B23AEE", "#9A32CD", "#68228B", "#4B0082", "#8A2BE2", "#9B30FF", "#912CEE", "#7D26CD", "#551A8B", "#9370DB", "#AB82FF", "#9F79EE", "#8968CD", "#5D478B", "#483D8B", "#8470FF", "#7B68EE", "#6A5ACD", "#836FFF", "#7A67EE", "#6959CD", "#473C8B", "#F8F8FF", "#E6E6FA", "#0000FF", "#0000EE", "#0000CD", "#00008B", "#000080", "#191970", "#3D59AB", "#4169E1", "#4876FF", "#436EEE", "#3A5FCD", "#27408B", "#6495ED", "#B0C4DE", "#CAE1FF", "#BCD2EE", "#A2B5CD", "#6E7B8B", "#778899", "#708090", "#C6E2FF", "#B9D3EE", "#9FB6CD", "#6C7B8B", "#1E90FF", "#1C86EE", "#1874CD", "#104E8B", "#F0F8FF", "#4682B4", "#63B8FF", "#5CACEE", "#4F94CD", "#36648B", "#87CEFA", "#B0E2FF", "#A4D3EE", "#8DB6CD", "#607B8B", "#87CEFF", "#7EC0EE", "#6CA6CD", "#4A708B", "#87CEEB", "#00BFFF", "#00B2EE", "#009ACD", "#00688B", "#33A1C9", "#ADD8E6", "#BFEFFF", "#B2DFEE", "#9AC0CD", "#68838B", "#B0E0E6", "#98F5FF", "#8EE5EE", "#7AC5CD", "#53868B", "#00F5FF", "#00E5EE", "#00C5CD", "#00868B", "#5F9EA0", "#00CED1", "#F0FFFF", "#E0EEEE", "#C1CDCD", "#838B8B", "#E0FFFF", "#D1EEEE", "#B4CDCD", "#7A8B8B", "#BBFFFF", "#AEEEEE", "#96CDCD", "#668B8B", "#2F4F4F", "#97FFFF", "#8DEEEE", "#79CDCD", "#528B8B", "#00FFFF", "#00EEEE", "#00CDCD", "#008B8B", "#008080", "#48D1CC", "#20B2AA", "#03A89E", "#40E0D0", "#808A87", "#00C78C", "#7FFFD4", "#76EEC6", "#66CDAA", "#458B74", "#00FA9A", "#F5FFFA", "#00FF7F", "#00EE76", "#00CD66", "#008B45", "#3CB371", "#54FF9F", "#4EEE94", "#43CD80", "#2E8B57", "#00C957", "#BDFCC9", "#3D9140", "#F0FFF0", "#E0EEE0", "#C1CDC1", "#838B83", "#8FBC8F", "#C1FFC1", "#B4EEB4", "#9BCD9B", "#698B69", "#98FB98", "#9AFF9A", "#90EE90", "#7CCD7C", "#548B54", "#32CD32", "#228B22", "#00FF00", "#00EE00", "#00CD00", "#008B00", "#008000", "#006400", "#308014", "#7CFC00", "#7FFF00", "#76EE00", "#66CD00", "#458B00", "#ADFF2F", "#CAFF70", "#BCEE68", "#A2CD5A", "#6E8B3D", "#556B2F", "#6B8E23", "#C0FF3E", "#B3EE3A", "#9ACD32", "#698B22", "#FFFFF0", "#EEEEE0", "#CDCDC1", "#8B8B83", "#F5F5DC", "#FFFFE0", "#EEEED1", "#CDCDB4", "#8B8B7A", "#FAFAD2", "#FFFF00", "#EEEE00", "#CDCD00", "#8B8B00", "#808069", "#808000", "#BDB76B", "#FFF68F", "#EEE685", "#CDC673", "#8B864E", "#F0E68C", "#EEE8AA", "#FFFACD", "#EEE9BF", "#CDC9A5", "#8B8970", "#FFEC8B", "#EEDC82", "#CDBE70", "#8B814C", "#E3CF57", "#FFD700", "#EEC900", "#CDAD00", "#8B7500", "#FFF8DC", "#EEE8CD", "#CDC8B1", "#8B8878", "#DAA520", "#FFC125", "#EEB422", "#CD9B1D", "#8B6914", "#B8860B", "#FFB90F", "#EEAD0E", "#CD950C", "#8B6508", "#FFA500", "#EE9A00", "#CD8500", "#8B5A00", "#FFFAF0", "#FDF5E6", "#F5DEB3", "#FFE7BA", "#EED8AE", "#CDBA96", "#8B7E66", "#FFE4B5", "#FFEFD5", "#FFEBCD", "#FFDEAD", "#EECFA1", "#CDB38B", "#8B795E", "#FCE6C9", "#D2B48C", "#9C661F", "#FF9912", "#FAEBD7", "#FFEFDB", "#EEDFCC", "#CDC0B0", "#8B8378", "#DEB887", "#FFD39B", "#EEC591", "#CDAA7D", "#8B7355", "#FFE4C4", "#EED5B7", "#CDB79E", "#8B7D6B", "#E3A869", "#ED9121", "#FF8C00", "#FF7F00", "#EE7600", "#CD6600", "#8B4500", "#FF8000", "#FFA54F", "#EE9A49", "#CD853F", "#8B5A2B", "#FAF0E6", "#FFDAB9", "#EECBAD", "#CDAF95", "#8B7765", "#FFF5EE", "#EEE5DE", "#CDC5BF", "#8B8682", "#F4A460", "#C76114", "#D2691E", "#FF7F24", "#EE7621", "#CD661D", "#8B4513", "#292421", "#FF7D40", "#FF6103", "#8A360F", "#A0522D", "#FF8247", "#EE7942", "#CD6839", "#8B4726", "#FFA07A", "#EE9572", "#CD8162", "#8B5742", "#FF7F50", "#FF4500", "#EE4000", "#CD3700", "#8B2500", "#5E2612", "#E9967A", "#FF8C69", "#EE8262", "#CD7054", "#8B4C39", "#FF7256", "#EE6A50", "#CD5B45", "#8B3E2F", "#8A3324", "#FF6347", "#EE5C42", "#CD4F39", "#8B3626", "#FA8072", "#FFE4E1", "#EED5D2", "#CDB7B5", "#8B7D7B", "#FFFAFA", "#EEE9E9", "#CDC9C9", "#8B8989", "#BC8F8F", "#FFC1C1", "#EEB4B4", "#CD9B9B", "#8B6969", "#F08080", "#CD5C5C", "#FF6A6A", "#EE6363", "#8B3A3A", "#CD5555", "#A52A2A", "#FF4040", "#EE3B3B", "#CD3333", "#8B2323", "#B22222", "#FF3030", "#EE2C2C", "#CD2626", "#8B1A1A", "#FF0000", "#EE0000", "#CD0000", "#8B0000", "#800000", "#8E388E", "#7171C6", "#7D9EC0", "#388E8E", "#71C671", "#8E8E38", "#C5C1AA", "#C67171", "#555555", "#1E1E1E", "#282828", "#515151", "#5B5B5B", "#848484", "#8E8E8E", '#7cb5ec', '#434348', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#5645b3', '#a23456', '#babbab', '#fffdd1', '#5645b3', '#800000', '#a23456', '#ba3ba3', '#00bb33', '#1b2d45', '#05377b', '#e5e5e5', '#021631', '#002448', '#042c62', '#032149', '#949494', '#123456', '#696969', '#444444', '#c0ffee', '#098765', '#567890', '#654321', '#123456', '#a23456', '#bababa', '#babbab', '#ba3ba3', '#baebae', '#b3b3b3', '#dbb234', '#b3a271', '#0b03b3', '#0303bb', '#00bb33', '#3300bb', '#0033bb', '#bb0033', '#cc0033', '#d3ad33']
	#color_list = ["#52b2b5", "#842541", "#d93c6a", "#776153", "#fbccae", "#371404", "#702a09", "#3f0523", "#7c0b45", "#494632", "#948e64", "#0b2077", "#1238cf", "#0f4033", "#1e8369", "#641207", "#677110", "#6e603b", "#d1b771", "#d6e822", "#b21f0e", "#00ffff", "#f6cece", "#191200", "#cbdad5", "#a00093", "#a00090", "#a00082", "#a00076", "#a00073", "#582233", "#67313f", "#793c43", "#000181", "#000183", "#000186", "#ff00fb", "#cd6ae8", "#ff00ff", "#ff00fb", "#ff00ff", "#cd6ae8", "#B88183", "#922329", "#5A0007", "#D7BFC2", "#D86A78", "#FF8A9A", "#3B000A", "#E20027", "#943A4D", "#5B4E51", "#B05B6F", "#FEB2C6", "#D83D66", "#895563", "#FF1A59", "#FFDBE5", "#CC0744", "#CB7E98", "#997D87", "#6A3A4C", "#FF2F80", "#6B002C", "#A74571", "#C6005A", "#FF5DA7", "#300018", "#B894A6", "#FF90C9", "#7C6571", "#A30059", "#DA007C", "#5B113C", "#402334", "#D157A0", "#DDB6D0", "#885578", "#962B75", "#A97399", "#D20096", "#E773CE", "#AA5199", "#E704C4", "#6B3A64", "#FFA0F2", "#6F0062", "#B903AA", "#C895C5", "#FF34FF", "#320033", "#DBD5DD", "#EEC3FF", "#BC23FF", "#671190", "#201625", "#F5E1FF", "#BC65E9", "#D790FF", "#72418F", "#4A3B53", "#9556BD", "#B4A8BD", "#7900D7", "#A079BF", "#958A9F", "#837393", "#64547B", "#3A2465", "#353339", "#BCB1E5", "#9F94F0", "#9695C5", "#0000A6", "#000035", "#636375", "#00005F", "#97979E", "#7A7BFF", "#3C3E6E", "#6367A9", "#494B5A", "#3B5DFF", "#C8D0F6", "#6D80BA", "#8FB0FF", "#0045D2", "#7A87A1", "#324E72", "#00489C", "#0060CD", "#789EC9", "#012C58", "#99ADC0", "#001325", "#DDEFFF", "#59738A", "#0086ED", "#75797C", "#BDC9D2", "#3E89BE", "#8CD0FF", "#0AA3F7", "#6B94AA", "#29607C", "#404E55", "#006FA6", "#013349", "#0AA6D8", "#658188", "#5EBCD1", "#456D75", "#0089A3", "#B5F4FF", "#02525F", "#1CE6FF", "#001C1E", "#203B3C", "#A3C8C9", "#00A6AA", "#00C6C8", "#006A66", "#518A87", "#E4FFFC", "#66E1D3", "#004D43", "#809693", "#15A08A", "#00846F", "#00C2A0", "#00FECF", "#78AFA1", "#02684E", "#C2FFED", "#47675D", "#00D891", "#004B28", "#8ADBB4", "#0CBD66", "#549E79", "#1A3A2A", "#6C8F7D", "#008941", "#63FFAC", "#1BE177", "#006C31", "#B5D6C3", "#3D4F44", "#4B8160", "#66796D", "#71BB8C", "#04F757", "#001E09", "#D2DCD5", "#00B433", "#9FB2A4", "#003109", "#A3F3AB", "#456648", "#51A058", "#83A485", "#7ED379", "#D1F7CE", "#A1C299", "#061203", "#1E6E00", "#5EFF03", "#55813B", "#3B9700", "#4FC601", "#1B4400", "#C2FF99", "#788D66", "#868E7E", "#83AB58", "#374527", "#98D058", "#C6DC99", "#A4E804", "#76912F", "#8BB400", "#34362D", "#4C6001", "#DFFB71", "#6A714A", "#222800", "#6B7900", "#3A3F00", "#BEC459", "#FEFFE6", "#A3A489", "#9FA064", "#FFFF00", "#61615A", "#FFFFFE", "#9B9700", "#CFCDAC", "#797868", "#575329", "#FFF69F", "#8D8546", "#F4D749", "#7E6405", "#1D1702", "#CCAA35", "#CCB87C", "#453C23", "#513A01", "#FFB500", "#A77500", "#D68E01", "#B79762", "#7A4900", "#372101", "#886F4C", "#A45B02", "#E7AB63", "#FAD09F", "#C0B9B2", "#938A81", "#A38469", "#D16100", "#A76F42", "#5B4534", "#5B3213", "#CA834E", "#FF913F", "#953F00", "#D0AC94", "#7D5A44", "#BE4700", "#FDE8DC", "#772600", "#A05837", "#EA8B66", "#391406", "#FF6832", "#C86240", "#29201D", "#B77B68", "#806C66", "#FFAA92", "#89412E", "#E83000", "#A88C85", "#F7C9BF", "#643127", "#E98176", "#7B4F4B", "#1E0200", "#9C6966", "#BF5650", "#BA0900", "#FF4A46", "#F4ABAA", "#000000", "#452C2C", "#C8A1A1"]
	#color_list = ["#0048BA", "#4C2F27", "#B0BF1A", "#7CB9E8", "#C9FFE5", "#B284BE", "#5D8AA8", "#00308F", "#72A0C1", "#AF002A", "#F0F8FF", "#84DE02", "#E32636", "#C46210", "#EFDECD", "#E52B50", "#9F2B68", "#F19CBB", "#AB274F", "#D3212D", "#3B7A57", "#FFBF00", "#FF7E00", "#FF033E", "#9966CC", "#A4C639", "#F2F3F4", "#CD9575", "#665D1E", "#915C83", "#841B2D", "#FAEBD7", "#008000", "#8DB600", "#FBCEB1", "#00FFFF", "#7FFFD4", "#D0FF14", "#4B5320", "#3B444B", "#8F9779", "#E9D66B", "#B2BEB5", "#87A96B", "#FF9966", "#A52A2A", "#FDEE00", "#6E7F80", "#568203", "#C39953", "#007FFF", "#F0FFFF", "#F0FFFF", "#DBE9F4", "#2E5894", "#89CFF0", "#A1CAF1", "#F4C2C2", "#FEFEFA", "#FF91AF", "#21ABCD", "#FAE7B5", "#FFE135", "#006A4E", "#E0218A", "#7C0A02", "#1DACD6", "#848482", "#98777B", "#BCD4E6", "#9F8170", "#FA6E79", "#F5F5DC", "#9C2542", "#E88E5A", "#FFE4C4", "#3D2B1F", "#967117", "#CAE00D", "#BFFF00", "#FE6F5E", "#BF4F51", "#000000", "#3D0C02", "#54626F", "#253529", "#3B3C36", "#BFAFB2", "#FFEBCD", "#A57164", "#318CE7", "#ACE5EE", "#FAF0BE", "#0000FF", "#1F75FE", "#0093AF", "#0087BD", "#0018A8", "#333399", "#0247FE", "#A2A2D0", "#00B9FB", "#5DADEC", "#ACE5EE", "#126180", "#5072A7", "#6699CC", "#0D98BA", "#553592", "#8A2BE2", "#4F86F7", "#1C1CF0", "#DE5D83", "#79443B", "#0095B6", "#E3DAC9", "#DDE26A", "#CC0000", "#006A4E", "#873260", "#0070FF", "#B5A642", "#CB4154", "#1DACD6", "#66FF00", "#BF94E4", "#D891EF", "#C32148", "#1974D2", "#FF007F", "#08E8DE", "#D19FE8", "#FFAA1D", "#3399FF", "#F4BBFF", "#FF55A3", "#FB607F", "#004225", "#CD7F32", "#737000", "#964B00", "#A52A2A", "#AF6E4D", "#cc9966", "#6B4423", "#1B4D3E", "#FFC1CC", "#E7FEFF", "#7BB661", "#F0DC82", "#480607", "#800020", "#DEB887", "#A17A74", "#CC5500", "#E97451", "#8A3324", "#BD33A4", "#702963", "#536872", "#5F9EA0", "#91A3B0", "#006B3C", "#ED872D", "#E30022", "#FFF600", "#A67B5B", "#4B3621", "#1E4D2B", "#A3C1AD", "#C19A6B", "#EFBBCC", "#78866B", "#FFFF99", "#FFEF00", "#FF0800", "#E4717A", "#00BFFF", "#592720", "#C41E3A", "#00CC99", "#960018", "#D70040", "#EB4C42", "#FF0038", "#FFA6C9", "#B31B1B", "#56A0D3", "#ED9121", "#00563F", "#062A78", "#703642", "#C95A49", "#92A1CF", "#ACE1AF", "#007BA7", "#2F847C", "#B2FFFF", "#4997D0", "#DE3163", "#EC3B83", "#007BA7", "#2A52BE", "#6D9BC3", "#007AA5", "#E03C31", "#A0785A", "#F7E7CE", "#F1DDCF", "#36454F", "#232B2B", "#E68FAC", "#DFFF00", "#7FFF00", "#DE3163", "#FFB7C5", "#954535", "#DE6FA1", "#A8516E", "#AA381E", "#856088", "#4AFF00", "#7B3F00", "#D2691E", "#FFA700", "#98817B", "#E34234", "#CD607E", "#D2691E", "#E4D00A", "#9FA91F", "#7F1734", "#FBCCE7", "#0047AB", "#D2691E", "#965A3E", "#6F4E37", "#C4D8E2", "#F88379", "#002E63", "#8C92AC", "#B87333", "#DA8A67", "#AD6F69", "#CB6D51", "#996666", "#FF3800", "#FF7F50", "#F88379", "#FF4040", "#FD7C6E", "#893F45", "#FBEC5D", "#B31B1B", "#6495ED", "#FFF8DC", "#2E2D88", "#FFF8E7", "#FFBCD9", "#81613C", "#FFFDD0", "#DC143C", "#BE0032", "#990000", "#F5F5F5", "#00FFFF", "#00B7EB", "#4E82B4", "#28589C", "#188BC2", "#4682BF", "#58427C", "#FFD300", "#F56FA1", "#FFFF31", "#F0E130", "#00008B", "#666699", "#654321", "#88654E", "#5D3954", "#A40000", "#08457E", "#986960", "#CD5B45", "#008B8B", "#536878", "#B8860B", "#A9A9A9", "#013220", "#006400", "#1F262A", "#00416A", "#6E6EF9", "#1A2421", "#BDB76B", "#483C32", "#734F96", "#534B4F", "#543D37", "#8B008B", "#A9A9A9", "#003366", "#4A5D23", "#556B2F", "#FF8C00", "#9932CC", "#779ECB", "#03C03C", "#966FD6", "#C23B22", "#E75480", "#003399", "#4F3A3C", "#301934", "#872657", "#8B0000", "#E9967A", "#560319", "#8FBC8F", "#3C1414", "#8CBED6", "#483D8B", "#2F4F4F", "#177245", "#918151", "#FFA812", "#483C32", "#CC4E5C", "#00CED1", "#D1BEA8", "#9400D3", "#9B870C", "#00703C", "#555555", "#D70A53", "#40826D", "#A9203E", "#EF3038", "#E9692C", "#DA3287", "#FAD6A5", "#B94E48", "#704241", "#C154C1", "#056608", "#0E7C61", "#004B49", "#333366", "#F5C71A", "#9955BB", "#CC00CC", "#820000", "#D473D4", "#355E3B", "#FFCBA4", "#FF1493", "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80","#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100","#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F","#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09","#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66","#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C","#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81","#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00","#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700","#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]
	color_list = ["#DB3340", "#E8B71A", "#1FDA9A", "#28ABE3", "#665D1E", "#915C83", "#841B2D", "#008000", "#FF3424", "#3284FF", "#FFBB00", "#26B14C", "#8DB600", "#5E412F", "#84DE02", "#E32636", "#3396fb", "#9F2B68", "#F19CBB", "#AB274F", "#D3212D", "#3B7A57", "#FFBF00", "#FF7E00", "#FF033E", "#9966CC", "#A4C639", "#0f4033", "#B2BEB5", "#87A96B", "#FF9966", "#A52A2A", "#FDEE00", "#6E7F80", "#568203", "#C39953", "#007FFF", "#DBE9F4", "#2E5894", "#89CFF0", "#A1CAF1", "#0048BA", "#B0BF1A", "#7CB9E8", "#B284BE", "#00308F", "#72A0C1", "#AF002A", "#CD9575", "#665D1E", "#915C83", "#841B2D", "#FAEBD7", "#008000", "#8DB600", "#FBCEB1", "#00FFFF", "#7FFFD4", "#D0FF14", "#4B5320", "#3B444B", "#8F9779", "#E9D66B", "#B2BEB5", "#87A96B", "#FF9966", "#A52A2A", "#FDEE00", "#6E7F80", "#568203", "#C39953", "#007FFF", "#F0FFFF", "#F0FFFF", "#DBE9F4", "#2E5894", "#89CFF0", "#A1CAF1", "#F4C2C2", "#FEFEFA", "#FF91AF", "#21ABCD", "#FAE7B5", "#FFE135", "#006A4E", "#E0218A", "#7C0A02", "#1DACD6", "#848482", "#98777B", "#BCD4E6", "#9F8170", "#FA6E79", "#F5F5DC", "#9C2542", "#E88E5A", "#FFE4C4", "#3D2B1F", "#967117", "#CAE00D", "#BFFF00", "#FE6F5E", "#BF4F51", "#000000", "#3D0C02", "#54626F", "#253529", "#3B3C36", "#BFAFB2", "#FFEBCD", "#A57164", "#318CE7", "#ACE5EE", "#FAF0BE", "#0000FF", "#1F75FE", "#0093AF", "#0087BD", "#0018A8", "#333399", "#0247FE", "#A2A2D0", "#00B9FB", "#5DADEC", "#ACE5EE", "#126180", "#5072A7", "#6699CC", "#0D98BA", "#553592", "#8A2BE2", "#4F86F7", "#1C1CF0", "#DE5D83", "#79443B", "#0095B6", "#E3DAC9", "#DDE26A", "#CC0000", "#006A4E", "#873260", "#0070FF", "#B5A642", "#CB4154", "#1DACD6", "#66FF00", "#BF94E4", "#D891EF", "#C32148", "#1974D2", "#FF007F", "#08E8DE", "#D19FE8", "#FFAA1D", "#3399FF", "#F4BBFF", "#FF55A3", "#FB607F", "#004225", "#CD7F32", "#737000", "#964B00", "#A52A2A", "#AF6E4D", "#cc9966", "#6B4423", "#1B4D3E", "#FFC1CC", "#E7FEFF", "#7BB661", "#F0DC82", "#480607", "#800020", "#DEB887", "#A17A74", "#CC5500", "#E97451", "#8A3324", "#BD33A4", "#702963", "#536872", "#5F9EA0", "#91A3B0", "#006B3C", "#ED872D", "#E30022", "#FFF600", "#A67B5B", "#4B3621", "#1E4D2B", "#A3C1AD", "#C19A6B", "#EFBBCC", "#78866B", "#FFFF99", "#FFEF00", "#FF0800", "#E4717A", "#00BFFF", "#592720", "#C41E3A", "#00CC99", "#960018", "#D70040", "#EB4C42", "#FF0038", "#FFA6C9", "#B31B1B", "#56A0D3", "#ED9121", "#00563F", "#062A78", "#703642", "#C95A49", "#92A1CF", "#ACE1AF", "#007BA7", "#2F847C", "#B2FFFF", "#4997D0", "#DE3163", "#EC3B83", "#007BA7", "#2A52BE", "#6D9BC3", "#007AA5", "#E03C31", "#A0785A", "#F7E7CE", "#F1DDCF", "#36454F", "#232B2B", "#E68FAC", "#DFFF00", "#7FFF00", "#DE3163", "#FFB7C5", "#954535", "#DE6FA1", "#A8516E", "#AA381E", "#856088", "#4AFF00", "#7B3F00", "#D2691E", "#FFA700", "#98817B", "#E34234", "#CD607E", "#D2691E", "#E4D00A", "#9FA91F", "#7F1734", "#FBCCE7", "#0047AB", "#D2691E", "#965A3E", "#6F4E37", "#C4D8E2", "#F88379", "#002E63", "#8C92AC", "#B87333", "#DA8A67", "#AD6F69", "#CB6D51", "#996666", "#FF3800", "#FF7F50", "#F88379", "#FF4040", "#FD7C6E", "#893F45", "#FBEC5D", "#B31B1B", "#6495ED", "#FFF8DC", "#2E2D88", "#FFF8E7", "#FFBCD9", "#81613C", "#FFFDD0", "#DC143C", "#BE0032", "#990000", "#F5F5F5", "#00FFFF", "#00B7EB", "#4E82B4", "#28589C", "#188BC2", "#4682BF", "#58427C", "#FFD300", "#F56FA1", "#FFFF31", "#F0E130", "#00008B", "#666699", "#654321", "#88654E", "#5D3954", "#A40000", "#08457E", "#986960", "#CD5B45", "#008B8B", "#536878", "#B8860B", "#A9A9A9", "#013220", "#006400", "#1F262A", "#00416A", "#6E6EF9", "#1A2421", "#BDB76B", "#483C32", "#734F96", "#534B4F", "#543D37", "#8B008B", "#A9A9A9", "#003366", "#4A5D23", "#556B2F", "#FF8C00", "#9932CC", "#779ECB", "#03C03C", "#966FD6", "#C23B22", "#E75480", "#003399", "#4F3A3C", "#301934", "#872657", "#8B0000", "#E9967A", "#560319", "#8FBC8F", "#3C1414", "#8CBED6", "#483D8B", "#2F4F4F", "#177245", "#918151", "#FFA812", "#483C32", "#CC4E5C", "#00CED1", "#D1BEA8", "#9400D3", "#9B870C", "#00703C", "#555555", "#D70A53", "#40826D", "#A9203E", "#EF3038", "#E9692C", "#DA3287", "#FAD6A5", "#B94E48", "#704241", "#C154C1", "#056608", "#0E7C61", "#004B49", "#333366", "#F5C71A", "#9955BB", "#CC00CC", "#820000", "#D473D4", "#355E3B", "#FFCBA4", "#FF1493", "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]

	#color_list = color_generate(color_count)
	#color_list = hsv_color_generate(color_count)
	i = 0
	for target, color in itertools.izip(target_list, color_list):
		if mode is None:
			dict_of_symbol_color[target] = (symbol_list[i], color)
			i += 1
		elif mode == 'color_only':
			dict_of_symbol_color[target] = color

	"""
	for target in target_list:
		if mode is None:
			dict_of_symbol_color[target] = (symbol_list[i], random_color())
		elif mode == 'color_only':
			dict_of_symbol_color[target] = random_color()
		i += 1
	"""
	return dict_of_symbol_color


def list_variance(list, tolerance):
	filtered_list = []
	temp = 0
	#tolerance=0.05
	for i in list:
		if i == 'NA':
			continue
		if isfloat(i) and float(i) > float(temp) + tolerance:
			filtered_list.append(i)
		temp = i
	return filtered_list


def design_dict_maker(design_file, mode=None):
	design = open(design_file, 'rU')
	dict_of_design = {}
	data = design.readlines()
	if mode == 'orange':
		for line in data:
			line = line.rstrip()
			line_split = line.split('\t')
			dict_of_design[line_split[0]] = line_split[1]
		return dict_of_design
	elif mode is None:
		for line in data:
			if data.index(line) == 0:
				continue
			line = line.rstrip()
			line_split = line.split('\t')
			dict_of_design[line_split[0]] = line_split[1]
		return dict_of_design
	elif mode == 'reverse':
		for line in data:
			if data.index(line) == 0:
				continue
			line = line.rstrip()
			line_split = line.split('\t')
			if line_split[1] not in dict_of_design:
				dict_of_design[line_split[1]] = [line_split[0]]
			elif line_split[1] in dict_of_design:
				dict_of_design[line_split[1]].append(line_split[0])
		return dict_of_design


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(query_list)


def mothur_result_parser(mothur_result_file):
	rarefile = open(mothur_result_file, 'rU')
	header = []
	header = rarefile.readline().rstrip().split('\t')
	#print header
	rarefile.close()
	rarefile = open(mothur_result_file, 'rU')
	reader = csv.DictReader(rarefile, header, delimiter="\t")
	columns = {}
	for row in reader:
		for key, value in row.items():
			if key == value:
				continue
			elif key in columns:
				columns[key].append(value)
			else:
				columns[key] = [value]
	return (header, columns)


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


def split_file_name(file):
	path = os.path.dirname(file) + '/'
	name = os.path.basename(file)
	if '.' in name:
		ext = '.' + '.'.join(name.split('.')[1:])
		absname = name.split('.')[0]
	else:
		ext = '.no_extension'
		absname = name
	return (path, absname, ext)


def shared_max_min(shared_file):
	list_number = []
	f = open(shared_file, 'rU')
	for line in f:
		line = line.rstrip()
		line_list = line.split('\t')
		for idx, i in enumerate(line_list):
			if idx in [0, 1, 2]:
				continue
			if isfloat(i):
				i = float(i)
				list_number.append(i)
	min_sh = min(list_number)
	max_sh = max(list_number)
	ave_sh = sum(list_number) / 1000
	median_sh = median(list_number)

	return (min_sh, max_sh, ave_sh, median_sh)


def median(lst):
	sortedLst = sorted(lst)
	lstLen = len(lst)
	index = (lstLen - 1) // 2

	if (lstLen % 2):
		return sortedLst[index]
	else:
		return (sortedLst[index] + sortedLst[index + 1]) / 2.0


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


def process_monitor(pid_list, stderr_list, stdout_list, outputdir, threads, mode):
	print "Execution started."
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
				sys.stdout.write("#")
				sys.stdout.flush()
				flag_list[pid_file] = True
				flag, stderr = validate_execution(stderr_file, stdout_file, outputdir,  mode)
				if flag is False:
					print "Error in result of this thread: ", the_pid
					print "All generated threads will be killed"
					kill_pid_list(pid_list, outputdir)
					#print stderr
					print "ABORTING!!!"
					sys.exit(2)
			if False in flag_list.values():
				finished_flag = False
			else:
				finished_flag = True
		time.sleep(1)
	sys.stdout.write("\n")
	#remove_list_of_file(stderr_list, outputdir)
	#remove_list_of_file(stdout_list, outputdir)
	#remove_list_of_file(pid_list, outputdir)
	print "Execution completed."

	return (True, stderr)


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


def report(report_string):
	f = open(report_file, "a")
	f.write(report_string)
	f.write("\n")
	f.close()


def error(report_string):
	f = open(report_file, "a")
	f.write('###############################  ERROR   ###################################\n')
	f.write(report_string)
	f.write("\n")
	f.write('############################################################################\n')
	f.close()


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


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def isPathExist(path):
	# check if the path exist and have access
	if os.path.exists(path) and os.access(path, os.R_OK):
		return True
	else:
		return False


def isFileEmpty(fname):
	# check if the path exist and have access
	if isFileExist(fname) is False:
		return False
	if os.path.getsize(fname) > 0:
		return True
	else:
		return False


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


def create_folder(path):
	# create folder in specified path
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False


# ################################### JAVASCRIPT ############################### #
def pca_3d_plotter(alpha_path, pca_html_string_dict, absname, design_file):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	check_box_string = ''
	for each_design in design_dictionary_rev:
		############################################################################################################################################################################
		check_box_string += """
				<div class="checkbox">
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>
				</div>
		"""
	pca_html_string = ''
	pca_script_string = ''
	pca_plotly_script_string = ''
	for pca_key in pca_html_string_dict.keys():
		pca_script_string += pca_html_string_dict[pca_key][1] + '\n'
		pca_html_string += pca_html_string_dict[pca_key][0] + '\n'
		pca_plotly_script_string += pca_html_string_dict[pca_key][2] + '\n'
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
				""" + pca_script_string + """
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">

				

						<!--PCA 3D PLOT-->
						<div id="PCA_3D_PLOTTER" class="panel panel-default">
							<div id="PCA_3D_PLOTTER_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_PCA_3D_PLOTTER_BODY">PCA_3D_PLOTTER</a>
							</h4>
							</div><!-- END of PCA_3D_PLOTTER_HEADING -->
							<div id="COLLAPSE_PCA_3D_PLOTTER_BODY" class="panel-collapse collapse">
								<div id="PCA_3D_PLOTTER_BODY" class="panel-body">
									<div id="PCA_3D_PLOTTER_ACCORDION" class="panel-group">
										""" + pca_html_string + """
									</div><!-- END of PCA_3D_PLOTTER_ACCORDION -->
								</div><!-- END of PCA_3D_PLOTTER_BODY -->
							</div><!-- END of COLLAPSE_PCA_3D_PLOTTER_BODY -->
						</div><!-- END of PCA_3D_PLOTTER_GROUP -->
						<!--END OF PCA_3D_PLOTTER-->
						


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					console.log(selectedbox)
					console.log(unselectedbox)

					var all_divs = document.getElementsByTagName("div");
					var biomarker_discovery_divs = []

					for(var div_counter=0; div_counter< all_divs.length; div_counter++){
						if (all_divs[div_counter].id.indexOf("biomarker_discovery") != -1 ){
							biomarker_discovery_divs.push(all_divs[div_counter])
						}
					}
					
					
					var bio_div_to_hide_list = [];
					var bio_div_to_keep_list = [];
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
						for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
							if (biomarker_discovery_divs[div_counter].id.indexOf(selectedbox[i]) != -1){
								var div_to_hide = document.getElementById(biomarker_discovery_divs[div_counter].id);
								bio_div_to_hide_list.push(div_to_hide);
								
							}
						}

					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
					for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
						var div_to_keep = document.getElementById(biomarker_discovery_divs[div_counter].id)
						if(bio_div_to_hide_list.indexOf(div_to_keep) == -1){

							bio_div_to_keep_list.push(div_to_keep);
						}
					}
					console.log(bio_div_to_hide_list);
					for (var j=0; j < bio_div_to_hide_list.length; j++){
						if (bio_div_to_hide_list[j] != null){
							if (bio_div_to_hide_list[j].style.display !== 'none') {
								bio_div_to_hide_list[j].style.display = 'none';
							};
						}
					}
					console.log(bio_div_to_keep_list);
					for (var j=0; j < bio_div_to_keep_list.length; j++){
						if (bio_div_to_keep_list[j] != null){
							if (bio_div_to_keep_list[j].style.display == 'none') {
								bio_div_to_keep_list[j].style.display = 'block';
							};
						}
					}
					
				}
			</script>
			<script type="text/javascript">
			""" + pca_plotly_script_string + """
			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def all_plotter(alpha_path, rarefaction_file, sample_abundance_file, bacterial_abundance_file, summary_table_header, summary_table_body, alpha_diversity_summary_file_name, biomarker_discovery_string, pca_html_string_dict, absname, design_file):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	
	rarefaction_html = ''
	rarefaction_script = ''
	
	sample_abundance_html = ''
	sample_abundance_script = ''

	bacterial_abundance_html = ''
	bacterial_abundance_script = ''

	check_box_string = ''
	for each_design in design_dictionary_rev:
		rarefaction_html += """<div id='div_rarefaction_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		rarefaction_script += """\tPlotly.newPlot('div_rarefaction_""" + each_design + """', data_rarefaction_""" + each_design + """, layout_rarefaction_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""

		sample_abundance_html += """<div id='div_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		sample_abundance_script += """\tPlotly.newPlot('div_abundance_""" + each_design + """', data_abundance_""" + each_design + """, layout_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""

		bacterial_abundance_html += """<div id='div_bacterial_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		bacterial_abundance_script += """\tPlotly.newPlot('div_bacterial_abundance_""" + each_design + """', data_bacterial_abundance_""" + each_design + """, layout_bacterial_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""

		check_box_string += """
				<div class="checkbox">\n
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>\n
				</div>\n
		"""
	pca_html_string = ''
	pca_script_string = ''
	pca_plotly_script_string = ''
	for pca_key in pca_html_string_dict.keys():
		pca_script_string += pca_html_string_dict[pca_key][1] + '\n'
		pca_html_string += pca_html_string_dict[pca_key][0] + '\n'
		pca_plotly_script_string += pca_html_string_dict[pca_key][2] + '\n'
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
				""" + pca_script_string + """
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
						<!--RAREFACTION-->
						<div id="RAREFACTION" class="panel panel-default">
							<div id="RAREFACTION_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_RAREFACTION_BODY">RAREFACTION CURVE PLOT</a>
								<a data-toggle="collapse" href="#INFO_RAREFACTION"><span class="fa fa-download fa-1x pull-right"></span></a>
								<div id="INFO_RAREFACTION" class="panel-collapse collapse" role="tabpanel">
									<div class="form-inline" role="form">
										<br>
										<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + rarefaction_file + """' download='""" + rarefaction_file + """'>Rarefaction_plot_file <span class="fa fa-download fa-1x right"></span></a>
									</div>
								</div>
							</h4>
							</div><!-- END of RAREFACTION_HEADING -->
							<div id="COLLAPSE_RAREFACTION_BODY" class="panel-collapse collapse">
								<div id="RAREFACTION_BODY" class="panel-body">\n
									""" + rarefaction_html + """
								</div><!-- END of RAREFACTION_BODY -->
							</div><!-- END of COLLAPSE_RAREFACTION_BODY -->
						</div><!-- END of RAREFACTION_GROUP -->
						<!--END OF RAREFACTION-->

						<!--SAMPLE ABUNDANCE-->
						<div id="SAMPLE_ABUNDANCE" class="panel panel-default">
							<div id="SAMPLE_ABUNDANCE_HEADING" class="panel-heading">
								<h4 class="panel-title">
									<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_SAMPLE_ABUNDANCE_BODY">SAMPLE_ABUNDANCE BAR PLOT</a>
									<a data-toggle="collapse" href="#INFO_SAMPLE_ABUNDANCE"><span class="fa fa-download fa-1x pull-right"></span></a>
										<div id="INFO_SAMPLE_ABUNDANCE" class="panel-collapse collapse" role="tabpanel">
											<div class="form-inline" role="form">
												<br>
											<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + sample_abundance_file + """' download='""" + sample_abundance_file + """'>Sample_abundance_file <span class="fa fa-download fa-1x right"></span></a>
										</div>
								</div>
								</h4>
							</div><!-- END of SAMPLE_ABUNDANCE_HEADING -->
							<div id="COLLAPSE_SAMPLE_ABUNDANCE_BODY" class="panel-collapse collapse">
								<div id="SAMPLE_ABUNDANCE_BODY" class="panel-body">\n
									""" + sample_abundance_html + """
								</div><!-- END of SAMPLE_ABUNDANCE_BODY -->
							</div><!-- END of COLLAPSE_SAMPLE_ABUNDANCE_BODY -->
						</div><!-- END of SAMPLE_ABUNDANCE_GROUP -->
						<!--END OF SAMPLE ABUNDANCE-->


						<!--BACTERIAL ABUNDANCE-->
						<div id="BACTERIAL_ABUNDANCE" class="panel panel-default">
							<div id="BACTERIAL_ABUNDANCE_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BACTERIAL_ABUNDANCE_BODY">BACTERIAL_ABUNDANCE BAR PLOT</a>
								<a data-toggle="collapse" href="#INFO_BACTERIAL_ABUNDANCE"><span class="fa fa-download fa-1x pull-right"></span></a>
									<div id="INFO_BACTERIAL_ABUNDANCE" class="panel-collapse collapse" role="tabpanel">
										<div class="form-inline" role="form">
											<br>
										<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + bacterial_abundance_file + """' download='""" + bacterial_abundance_file + """'>Bacterial_abundance_file <span class="fa fa-download fa-1x right"></span></a>
									</div>
								</div>
							</h4>
							</div><!-- END of BACTERIAL_ABUNDANCE_HEADING -->
							<div id="COLLAPSE_BACTERIAL_ABUNDANCE_BODY" class="panel-collapse collapse">
								<div id="BACTERIAL_ABUNDANCE_BODY" class="panel-body">\n
									""" + bacterial_abundance_html + """
								</div><!-- END of BACTERIAL_ABUNDANCE_BODY -->
							</div><!-- END of COLLAPSE_BACTERIAL_ABUNDANCE_BODY -->
						</div><!-- END of BACTERIAL_ABUNDANCE_GROUP -->
						<!--END OF BACTERIAL ABUNDANCE-->


						<!--ALPHA DIVERTSITY SUMMARY-->
						<div id="ALPHA_DIVERTSITY_SUMMARY" class="panel panel-default">
							<div id="ALPHA_DIVERTSITY_SUMMARY_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY">ALPHA DIVERSITY SUMMARY</a>
								<a data-toggle="collapse" href="#INFO_ALPHA_DIVERTSITY_SUMMARY"><span class="fa fa-download fa-1x pull-right"></span></a>
									<div id="INFO_ALPHA_DIVERTSITY_SUMMARY" class="panel-collapse collapse" role="tabpanel">
										<div class="form-inline" role="form">
											<br>
										<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + alpha_diversity_summary_file_name + """' download='""" + alpha_diversity_summary_file_name + """'>Alpha_diversity_summary_file_name <span class="fa fa-download fa-1x right"></span></a>
									</div>
								</div>
							</h4>
							</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_HEADING -->
							<div id="COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-collapse collapse">
								<div id="ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-body">
									<!-- Table -->
									<div class="table-responsive">
										<table id="alpha_diversity_summary" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%">
											""" + summary_table_header + """
											""" + summary_table_body + """
										</table>
									</div>
									<!-- END of TABLE -->
									
								</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_BODY -->
							</div><!-- END of COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY -->
						</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_GROUP -->
						<!--END OF ALPHA DIVERTSITY SUMMARY-->


						<!--PCA 3D PLOT-->
						<div id="PCA_3D_PLOTTER" class="panel panel-default">
							<div id="PCA_3D_PLOTTER_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_PCA_3D_PLOTTER_BODY">PCA_3D_PLOTTER</a>
							</h4>
							</div><!-- END of PCA_3D_PLOTTER_HEADING -->
							<div id="COLLAPSE_PCA_3D_PLOTTER_BODY" class="panel-collapse collapse">
								<div id="PCA_3D_PLOTTER_BODY" class="panel-body">
									<div id="PCA_3D_PLOTTER_ACCORDION" class="panel-group">
										""" + pca_html_string + """
									</div><!-- END of PCA_3D_PLOTTER_ACCORDION -->
								</div><!-- END of PCA_3D_PLOTTER_BODY -->
							</div><!-- END of COLLAPSE_PCA_3D_PLOTTER_BODY -->
						</div><!-- END of PCA_3D_PLOTTER_GROUP -->
						<!--END OF PCA_3D_PLOTTER-->



						<!--BIOMARKER_DISCOVERY-->
						<div id="BIOMARKER_DISCOVERY_GROUP" class="panel panel-default">
							<div id="BIOMARKER_DISCOVERY_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BIOMARKER_DISCOVERY_BODY">BIOMARKER_DISCOVERY</a>
							</h4>
							</div><!-- END of BIOMARKER_DISCOVERY_HEADING -->
							<div id="COLLAPSE_BIOMARKER_DISCOVERY_BODY" class="panel-collapse collapse">
								<div id="BIOMARKER_DISCOVERY_BODY" class="panel-body">
									<div id="BIOMARKER_DISCOVERY_ACCORDION" class="panel-group">
												""" + biomarker_discovery_string + """
									</div><!-- END of BIOMARKER_DISCOVERY_ACCORDION -->
								</div><!-- END of BIOMARKER_DISCOVERY_BODY -->
							</div><!-- END of COLLAPSE_BIOMARKER_DISCOVERY_BODY -->
						</div><!-- END of BIOMARKER_DISCOVERY_GROUP -->
						<!--END OF BIOMARKER_DISCOVERY-->

				


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					console.log(selectedbox)
					console.log(unselectedbox)

					var all_divs = document.getElementsByTagName("div");
					var biomarker_discovery_divs = []

					for(var div_counter=0; div_counter< all_divs.length; div_counter++){
						if (all_divs[div_counter].id.indexOf("biomarker_discovery") != -1 ){
							biomarker_discovery_divs.push(all_divs[div_counter])
						}
					}
					
					
					var bio_div_to_hide_list = [];
					var bio_div_to_keep_list = [];
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
						for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
							if (biomarker_discovery_divs[div_counter].id.indexOf(selectedbox[i]) != -1){
								var div_to_hide = document.getElementById(biomarker_discovery_divs[div_counter].id);
								bio_div_to_hide_list.push(div_to_hide);
								
							}
						}

					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
					for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
						var div_to_keep = document.getElementById(biomarker_discovery_divs[div_counter].id)
						if(bio_div_to_hide_list.indexOf(div_to_keep) == -1){

							bio_div_to_keep_list.push(div_to_keep);
						}
					}
					console.log(bio_div_to_hide_list);
					for (var j=0; j < bio_div_to_hide_list.length; j++){
						if (bio_div_to_hide_list[j] != null){
							if (bio_div_to_hide_list[j].style.display !== 'none') {
								bio_div_to_hide_list[j].style.display = 'none';
							};
						}
					}
					console.log(bio_div_to_keep_list);
					for (var j=0; j < bio_div_to_keep_list.length; j++){
						if (bio_div_to_keep_list[j] != null){
							if (bio_div_to_keep_list[j].style.display == 'none') {
								bio_div_to_keep_list[j].style.display = 'block';
							};
						}
					}
					
				}
			</script>
			<script type="text/javascript">
				""" + rarefaction_script + """
				""" + bacterial_abundance_script + """
				""" + sample_abundance_script + """
				""" + pca_plotly_script_string + """

			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def biomarker_discovery_plotter(alpha_path, biomarker_discovery_string, absname, design_file):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	check_box_string = ''
	for each_design in design_dictionary_rev:
		############################################################################################################################################################################
		check_box_string += """
				<div class="checkbox">
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>
				</div>
		"""
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
				

						<!--BIOMARKER_DISCOVERY-->
						<div id="BIOMARKER_DISCOVERY_GROUP" class="panel panel-default">
							<div id="BIOMARKER_DISCOVERY_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BIOMARKER_DISCOVERY_BODY">BIOMARKER_DISCOVERY</a>
							</h4>
							</div><!-- END of BIOMARKER_DISCOVERY_HEADING -->
							<div id="COLLAPSE_BIOMARKER_DISCOVERY_BODY" class="panel-collapse collapse">
								<div id="BIOMARKER_DISCOVERY_BODY" class="panel-body">
									<div id="BIOMARKER_DISCOVERY_ACCORDION" class="panel-group">
												""" + biomarker_discovery_string + """
									</div><!-- END of BIOMARKER_DISCOVERY_ACCORDION -->
								</div><!-- END of BIOMARKER_DISCOVERY_BODY -->
							</div><!-- END of COLLAPSE_BIOMARKER_DISCOVERY_BODY -->
						</div><!-- END of BIOMARKER_DISCOVERY_GROUP -->
						<!--END OF BIOMARKER_DISCOVERY-->
						


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					console.log(selectedbox)
					console.log(unselectedbox)

					var all_divs = document.getElementsByTagName("div");
					var biomarker_discovery_divs = []

					for(var div_counter=0; div_counter< all_divs.length; div_counter++){
						if (all_divs[div_counter].id.indexOf("biomarker_discovery") != -1 ){
							biomarker_discovery_divs.push(all_divs[div_counter])
						}
					}
					
					
					var bio_div_to_hide_list = [];
					var bio_div_to_keep_list = [];
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
						for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
							if (biomarker_discovery_divs[div_counter].id.indexOf(selectedbox[i]) != -1){
								var div_to_hide = document.getElementById(biomarker_discovery_divs[div_counter].id);
								bio_div_to_hide_list.push(div_to_hide);
								
							}
						}

					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
					for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
						var div_to_keep = document.getElementById(biomarker_discovery_divs[div_counter].id)
						if(bio_div_to_hide_list.indexOf(div_to_keep) == -1){

							bio_div_to_keep_list.push(div_to_keep);
						}
					}
					console.log(bio_div_to_hide_list);
					for (var j=0; j < bio_div_to_hide_list.length; j++){
						if (bio_div_to_hide_list[j] != null){
							if (bio_div_to_hide_list[j].style.display !== 'none') {
								bio_div_to_hide_list[j].style.display = 'none';
							};
						}
					}
					console.log(bio_div_to_keep_list);
					for (var j=0; j < bio_div_to_keep_list.length; j++){
						if (bio_div_to_keep_list[j] != null){
							if (bio_div_to_keep_list[j].style.display == 'none') {
								bio_div_to_keep_list[j].style.display = 'block';
							};
						}
					}
					
				}
			</script>
			<script type="text/javascript">

			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def alpha_summary_plotter(alpha_path, alpha_diversity_summary_file_name, summary_table_header, summary_table_body, absname, design_file):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	check_box_string = ''
	for each_design in design_dictionary_rev:
		############################################################################################################################################################################
		check_box_string += """
				<div class="checkbox">
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>
				</div>
		"""
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
				

						<!--ALPHA DIVERTSITY SUMMARY-->
						<div id="ALPHA_DIVERTSITY_SUMMARY" class="panel panel-default">
							<div id="ALPHA_DIVERTSITY_SUMMARY_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY">ALPHA DIVERSITY SUMMARY</a>
								<a data-toggle="collapse" href="#INFO_ALPHA_DIVERTSITY_SUMMARY"><span class="fa fa-download fa-1x pull-right"></span></a>
									<div id="INFO_ALPHA_DIVERTSITY_SUMMARY" class="panel-collapse collapse" role="tabpanel">
										<div class="form-inline" role="form">
											<br>
										<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + alpha_diversity_summary_file_name + """' download='""" + alpha_diversity_summary_file_name + """'>Alpha_diversity_summary_file_name <span class="fa fa-download fa-1x right"></span></a>
									</div>
								</div>
							</h4>
							</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_HEADING -->
							<div id="COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-collapse collapse">
								<div id="ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-body">
									<!-- Table -->
									<div class="table-responsive">
										<table id="alpha_diversity_summary" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%">
											""" + summary_table_header + """
											""" + summary_table_body + """
										</table>
									</div>
									<!-- END of TABLE -->
									
								</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_BODY -->
							</div><!-- END of COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY -->
						</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_GROUP -->
						<!--END OF ALPHA DIVERTSITY SUMMARY-->
						


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];

					var all_divs = document.getElementsByTagName("div");
					var biomarker_discovery_divs = []

					for(var div_counter=0; div_counter< all_divs.length; div_counter++){
						if (all_divs[div_counter].id.indexOf("biomarker_discovery") != -1 ){
							biomarker_discovery_divs.push(all_divs[div_counter])
						}
					}

					
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						console.log(div_to_hide_list);
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
						for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
							if (biomarker_discovery_divs[div_counter].id.indexOf(selectedbox[i]) != -1 ){
								var div_to_hide = document.getElementById(biomarker_discovery_divs[div_counter].id);
								console.log(div_to_hide)
								
								if (div_to_hide.style.display !== 'none') {
									div_to_hide.style.display = 'none';
								
								}
							}
						}
					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						console.log(div_to_keep_list);
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
						for(var div_counter = 0; div_counter < biomarker_discovery_divs.length; div_counter++){
							if (biomarker_discovery_divs[div_counter].id.indexOf(unselectedbox[i]) != -1 ){
								var div_to_hide = document.getElementById(biomarker_discovery_divs[div_counter].id);
								console.log(div_to_hide)
								if (biomarker_discovery_divs[div_counter] != null){
									if (biomarker_discovery_divs[div_counter].style.display == 'none') {
										biomarker_discovery_divs[div_counter].style.display = 'block';
									};
								}
							}
						}
					}
				}
			</script>
			<script type="text/javascript">

			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def bacterial_abundance_plotter(alpha_path, bacterial_abundance_file, design_file, absname):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	bacterial_abundance_html = ''
	bacterial_abundance_script = ''
	check_box_string = ''
	for each_design in design_dictionary_rev:
		############################################################################################################################################################################
		bacterial_abundance_html += """<div id='div_bacterial_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		bacterial_abundance_script += """\tPlotly.newPlot('div_bacterial_abundance_""" + each_design + """', data_bacterial_abundance_""" + each_design + """, layout_bacterial_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""
		check_box_string += """
				<div class="checkbox">
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>
				</div>
		"""
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
				

						<!--BACTERIAL ABUNDANCE-->
						<div id="BACTERIAL_ABUNDANCE" class="panel panel-default">
							<div id="BACTERIAL_ABUNDANCE_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BACTERIAL_ABUNDANCE_BODY">BACTERIAL_ABUNDANCE BAR PLOT</a>
								<a data-toggle="collapse" href="#INFO_BACTERIAL_ABUNDANCE"><span class="fa fa-download fa-1x pull-right"></span></a>
									<div id="INFO_BACTERIAL_ABUNDANCE" class="panel-collapse collapse" role="tabpanel">
										<div class="form-inline" role="form">
											<br>
										<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + bacterial_abundance_file + """' download='""" + bacterial_abundance_file + """'>Bacterial_abundance_file <span class="fa fa-download fa-1x right"></span></a>
									</div>
								</div>
							</h4>
							</div><!-- END of BACTERIAL_ABUNDANCE_HEADING -->
							<div id="COLLAPSE_BACTERIAL_ABUNDANCE_BODY" class="panel-collapse collapse">
								<div id="BACTERIAL_ABUNDANCE_BODY" class="panel-body">\n
									""" + bacterial_abundance_html + """
								</div><!-- END of BACTERIAL_ABUNDANCE_BODY -->
							</div><!-- END of COLLAPSE_BACTERIAL_ABUNDANCE_BODY -->
						</div><!-- END of BACTERIAL_ABUNDANCE_GROUP -->
						<!--END OF BACTERIAL ABUNDANCE-->
						


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						console.log(div_to_hide_list);
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						console.log(div_to_keep_list);
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
				}
			</script>
			<script type="text/javascript">
				""" + bacterial_abundance_script + """
			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def sample_abundance_plotter(alpha_path, sample_abundance_file, design_file, absname):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	sample_abundance_html = ''
	sample_abundance_script = ''
	check_box_string = ''
	for each_design in design_dictionary_rev:
		############################################################################################################################################################################
		sample_abundance_html += """<div id='div_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		sample_abundance_script += """\tPlotly.newPlot('div_abundance_""" + each_design + """', data_abundance_""" + each_design + """, layout_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""
		check_box_string += """
				<div class="checkbox">
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>
				</div>
		"""
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
				

						<!--SAMPLE ABUNDANCE-->
						<div id="SAMPLE_ABUNDANCE" class="panel panel-default">
							<div id="SAMPLE_ABUNDANCE_HEADING" class="panel-heading">
								<h4 class="panel-title">
									<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_SAMPLE_ABUNDANCE_BODY">SAMPLE_ABUNDANCE BAR PLOT</a>
									<a data-toggle="collapse" href="#INFO_SAMPLE_ABUNDANCE"><span class="fa fa-download fa-1x pull-right"></span></a>
										<div id="INFO_SAMPLE_ABUNDANCE" class="panel-collapse collapse" role="tabpanel">
											<div class="form-inline" role="form">
												<br>
											<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + sample_abundance_file + """' download='""" + sample_abundance_file + """'>Sample_abundance_file <span class="fa fa-download fa-1x right"></span></a>
										</div>
								</div>
								</h4>
							</div><!-- END of SAMPLE_ABUNDANCE_HEADING -->
							<div id="COLLAPSE_SAMPLE_ABUNDANCE_BODY" class="panel-collapse collapse">
								<div id="SAMPLE_ABUNDANCE_BODY" class="panel-body">\n
									""" + sample_abundance_html + """
								</div><!-- END of SAMPLE_ABUNDANCE_BODY -->
							</div><!-- END of COLLAPSE_SAMPLE_ABUNDANCE_BODY -->
						</div><!-- END of SAMPLE_ABUNDANCE_GROUP -->
					<!--END OF SAMPLE ABUNDANCE-->
						


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						console.log(div_to_hide_list);
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						console.log(div_to_keep_list);
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
				}
			</script>
			<script type="text/javascript">
				""" + sample_abundance_script + """
			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def rarefaction_plotter(alpha_path, rarefaction_file, design_file, absname):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	rarefaction_html = ''
	rarefaction_css = ''
	rarefaction_script = ''
	check_box_string = ''
	for each_design in design_dictionary_rev:
		rarefaction_html += """<div id='div_rarefaction_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		rarefaction_script += """\tPlotly.newPlot('div_rarefaction_""" + each_design + """', data_rarefaction_""" + each_design + """, layout_rarefaction_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n"""
		rarefaction_css += """
		#div_rarefaction_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		check_box_string += """
				<div class="checkbox">\n
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>\n
				</div>\n
		"""
	javastring = ''
	javastring += """
		<!DOCTYPE html>
		<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
		<head>
			<meta charset="utf-8">
			<meta name="viewport" content="width=device-width, initial-scale=1">
			<!--PLOTLY LIBRARAY MAIN LIBRARY-->
			<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
			<!--############################-->
			<!--BOOTSTRAP LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
				<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
			<!--############################-->
			<!--DATATABLE LIBRARAY-->
				<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
				<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
				<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
			<!--############################-->
			<!--CUSTOMIZED SCRIPTS-->
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
				<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
			<!--############################-->
			<!--OPEN SCRIPTS-->
				<script type="text/javascript">
					$(document).ready(function() {
						$('#alpha_diversity_summary').DataTable( {
							"order": [[ 1, "desc" ]],
							"paging": false,
							"info": false,
							"ordering": true,
							"searching": false
						} );
						
					});
				</script>
			<!--############################-->
			<!--STYLES-->
				<style type="text/css">
					.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
					a:link {color: #669;}      /* unvisited link */
					a:visited {color: #669;}   /* visited link */
					a:hover {color: #669;}     /* mouse over link */
					a:active {color: #669;}    /* selected link */
				</style>
			<!--############################-->
		</head>
		<body>
			<div id="CONTAINER" class="container">
				<h2>16S data analysis results</h2>
				<div class="row">
					<dl class="dl-horizontal">
					<dt>Request Name</dt>
					<!--dd>kraemerls-20160216-H1</dd-->
					<dt>Project Title</dt>
					<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
					<dt>Project Description</dt>
					<!--dd>
					<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
					<dt>Experiment</dt>
					<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
					<dt>PI</dt>
					<!--dd>Helene Rosenberg</dd-->
					</dl>
				</div>
				<div class="row">
					<form role="form" method="post">
						<fieldset class="form-group">
							""" + check_box_string + """
							<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
						</fieldset>
					</form>
				</div>
				<div class="row">
					<div id="MAIN_PANEL_GROUP" class="panel-group">
				<!--RAREFACTION-->
				<div id="RAREFACTION" class="panel panel-default">
					<div id="RAREFACTION_HEADING" class="panel-heading">
					<h4 class="panel-title">
						<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_RAREFACTION_BODY">RAREFACTION CURVE PLOT</a>
						<a data-toggle="collapse" href="#INFO_RAREFACTION"><span class="fa fa-download fa-1x pull-right"></span></a>
						<div id="INFO_RAREFACTION" class="panel-collapse collapse" role="tabpanel">
							<div class="form-inline" role="form">
								<br>
								<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + rarefaction_file + """' download='""" + rarefaction_file + """'>Rarefaction_plot_file <span class="fa fa-download fa-1x right"></span></a>
							</div>
						</div>
					</h4>
					</div><!-- END of RAREFACTION_HEADING -->
					<div id="COLLAPSE_RAREFACTION_BODY" class="panel-collapse collapse">
						<div id="RAREFACTION_BODY" class="panel-body">\n
							""" + rarefaction_html + """
						</div><!-- END of RAREFACTION_BODY -->
					</div><!-- END of COLLAPSE_RAREFACTION_BODY -->
				</div><!-- END of RAREFACTION_GROUP -->
				<!--END OF RAREFACTION-->

				


					</div><!--END OF MAIN_PANEL_GROUP-->
				</div>
			</div><!--END OF MAIN CONTAINER-->
		</body>
			<script type="text/javascript">
				function getSelectedChbox(frm) {
					var selectedbox = [];
					var unselectedbox = [];
					var inpfields = frm.getElementsByTagName('input');
					var nr_inpfields = inpfields.length;
					for(var i = 0; i < nr_inpfields; i++) {
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
					if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
					}
					return [selectedbox, unselectedbox];
				}
				document.getElementById('submit').onclick = function(){
					var box = getSelectedChbox(this.form);
					var selectedbox = box[0];
					
					var unselectedbox = box[1];
					
					for(var i=0; i < selectedbox.length; i++){

						var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
						if (div_to_hide != null){
							if (div_to_hide.style.display !== 'none') {
							div_to_hide.style.display = 'none';
							};
						}
						var div_to_hide_list = [];
						div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
						console.log(div_to_hide_list);
						for (var j=0; j < div_to_hide_list.length; j++){
							if (div_to_hide_list[j] != null){
								if (div_to_hide_list[j].style.display !== 'none') {
									div_to_hide_list[j].style.display = 'none';
								};
							}
						}
					}
					for(var i=0; i < unselectedbox.length; i++){
						var div_to_keep = document.getElementById("div_abundance_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						if (div_to_keep != null){
						var div_to_keep = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
						if (div_to_keep != null){
							if (div_to_keep.style.display == 'none') {
							div_to_keep.style.display = 'block';
							};
						}
						var div_to_keep_list = [];
						div_to_keep_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
						console.log(div_to_keep_list);
						for (var j=0; j < div_to_keep_list.length; j++){
							if (div_to_keep_list[j] != null){
								if (div_to_keep_list[j].style.display == 'none') {
									div_to_keep_list[j].style.display = 'table-row';
								};
							}
						}
					}
				}
			</script>
			<script type="text/javascript">
				""" + rarefaction_script + """
			</script>
		</html>
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def html_maker(alpha_path, design_file, absname, summary_table_header, summary_table_body, biomarker_discovery_string):
	#Default mode bar ['sendDataToCloud', 'zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian' ]
	# Default mode bar 3d['zoom3d', 'pan3d', 'orbitRotation', 'tableRotation', 'resetCameraDefault3d', 'resetCameraLastSave3d', 'hoverClosest3d', 'resetViews' ]
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	rarefaction_html = ''
	rarefaction_css = ''
	rarefaction_script = ''
	abundance_html = ''
	abundance_css = ''
	abundance_script = ''
	bacterial_abundance_html = ''
	bacterial_abundance_css = ''
	bacterial_abundance_script = ''
	for each_design in design_dictionary_rev:
		rarefaction_html += """<div id='div_rarefaction_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		rarefaction_script += """\tPlotly.newPlot('div_rarefaction_""" + each_design + """', data_rarefaction_""" + each_design + """, layout_rarefaction_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
	  });\n"""
		rarefaction_css += """
		#div_rarefaction_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		############################################################################################################################################################################
		abundance_html += """<div id='div_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		bacterial_abundance_html += """<div id='div_bacterial_abundance_""" + each_design + """'></div>\n\t\t\t\t\t\t"""
		abundance_script += """\tPlotly.newPlot('div_abundance_""" + each_design + """', data_abundance_""" + each_design + """, layout_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
	  });\n"""

		bacterial_abundance_script += """\tPlotly.newPlot('div_bacterial_abundance_""" + each_design + """', data_bacterial_abundance_""" + each_design + """, layout_bacterial_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'pan2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
	  });\n"""
		abundance_css += """
		#div_abundance_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		bacterial_abundance_css += """
		#div_bacterial_abundance_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
	
	javastring = ''
	javastring += """
	<!DOCTYPE html>
	<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
	<head>
		<meta charset="utf-8">
		<meta name="viewport" content="width=device-width, initial-scale=1">
		<!--PLOTLY LIBRARAY MAIN LIBRARY-->
		<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
		<script type="text/javascript" src="./""" + absname + """_js/plotly-latest.min.js"></script>
		<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
		<!--############################-->
		<!--BOOTSTRAP LIBRARAY-->
			<link rel="stylesheet" type="text/css" href="./""" + absname + """_js/bootstrap.css">
			<script type="text/javascript" src="./""" + absname + """_js/bootstrap.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/jquery.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/tether.min.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/fontawesome.min.js"></script>
		<!--############################-->
		<!--DATATABLE LIBRARAY-->
			<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"/>
			<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
			<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
			<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
		<!--############################-->
		<!--CUSTOMIZED SCRIPTS-->
			<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_plotly.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance_plotly.js"></script>
			<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance_plotly.js"></script>
		<!--############################-->
		<!--OPEN SCRIPTS-->
			<script type="text/javascript">
				$(document).ready(function() {
					$('#alpha_diversity_summary').DataTable( {
						"order": [[ 1, "desc" ]],
						"paging": false,
						"info": false,
						"ordering": true,
						"searching": false
					} );
					
				});
			</script>
		<!--############################-->
		<!--STYLES-->
			<style type="text/css">
				.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
					background-color: yellow;
				}
				a:link {color: #669;}      /* unvisited link */
				a:visited {color: #669;}   /* visited link */
				a:hover {color: #669;}     /* mouse over link */
				a:active {color: #669;}    /* selected link */
			</style>
		<!--############################-->
	</head>
	<body>
		<div id="CONTAINER" class="container">
			<h2>16S data analysis results</h2>
			<div class="row">
				<dl class="dl-horizontal">
				<dt>Request Name</dt>
				<!--dd>kraemerls-20160216-H1</dd-->
				<dt>Project Title</dt>
				<!--dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd-->
				<dt>Project Description</dt>
				<!--dd>
				<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
				<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naïve mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
				<dt>Experiment</dt>
				<!--dd>Gut Microbiome 16S rRNA sequencing</dd-->
				<dt>PI</dt>
				<!--dd>Helene Rosenberg</dd-->
				</dl>
			</div>
			<div class="row">
				<form role="form" method="post">
					<fieldset class="form-group">
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="Control" >Control</label>
						</div>
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="C_naïve_NOD2_TLR2">C_naïve_NOD2_TLR2</label>
						</div>
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="D_NOD2_TLR2_L_plantarum_challenge">D_NOD2_TLR2_L_plantarum_challenge</label>
						</div>
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="A_naïve_C57BL_6">A_naïve_C57BL_6</label>
						</div>
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="B_C57BL_6_L_plantarum_challenge">B_C57BL_6_L_plantarum_challenge</label>
						</div>
						<div class="checkbox">
							<label><input type="checkbox" name="chb[]" value="MOCK">MOCK</label>
						</div>
						<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
					</fieldset>
				</form>
			</div>
			<div class="row">
				<div id="MAIN_PANEL_GROUP" class="panel-group">
				<!--SAMPLE ABUNDANCE-->
					<div id="SAMPLE_ABUNDANCE" class="panel panel-default">
						<div id="SAMPLE_ABUNDANCE_HEADING" class="panel-heading">
							<h4 class="panel-title">
								<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_SAMPLE_ABUNDANCE_BODY">SAMPLE_ABUNDANCE BAR PLOT</a>
							</h4>
						</div><!-- END of SAMPLE_ABUNDANCE_HEADING -->
						<div id="COLLAPSE_SAMPLE_ABUNDANCE_BODY" class="panel-collapse collapse">
							<div id="SAMPLE_ABUNDANCE_BODY" class="panel-body">\n
								""" + abundance_html + """
							</div><!-- END of SAMPLE_ABUNDANCE_BODY -->
						</div><!-- END of COLLAPSE_SAMPLE_ABUNDANCE_BODY -->
					</div><!-- END of SAMPLE_ABUNDANCE_GROUP -->
				<!--END OF SAMPLE ABUNDANCE-->


				<!--BACTERIAL ABUNDANCE-->
				<div id="BACTERIAL_ABUNDANCE" class="panel panel-default">
					<div id="BACTERIAL_ABUNDANCE_HEADING" class="panel-heading">
					<h4 class="panel-title">
						<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BACTERIAL_ABUNDANCE_BODY">BACTERIAL_ABUNDANCE BAR PLOT</a>
						<a data-toggle="collapse" href="#INFO_BACTERIAL_ABUNDANCE"><span class="fa fa-download fa-1x pull-right"></span></a>
							<div id="INFO_BACTERIAL_ABUNDANCE" class="panel-collapse collapse" role="tabpanel">
								<div class="form-inline" role="form">
									<br>
								<a type="button" class="btn btn-default btn-sm" href='./""" + absname + """_DATA/""" + bacterial_abundance_file + """' download='""" + bacterial_abundance_file + """'>Sample_abundance_file <span class="fa fa-download fa-1x right"></span></a>
							</div>
						</div>
					</h4>
					</div><!-- END of BACTERIAL_ABUNDANCE_HEADING -->
					<div id="COLLAPSE_BACTERIAL_ABUNDANCE_BODY" class="panel-collapse collapse">
						<div id="BACTERIAL_ABUNDANCE_BODY" class="panel-body">\n
							""" + bacterial_abundance_html + """
						</div><!-- END of BACTERIAL_ABUNDANCE_BODY -->
					</div><!-- END of COLLAPSE_BACTERIAL_ABUNDANCE_BODY -->
				</div><!-- END of BACTERIAL_ABUNDANCE_GROUP -->
				<!--END OF BACTERIAL ABUNDANCE-->


			<!--RAREFACTION-->
			<div id="RAREFACTION" class="panel panel-default">
				<div id="RAREFACTION_HEADING" class="panel-heading">
				<h4 class="panel-title">
					<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_RAREFACTION_BODY">RAREFACTION CURVE PLOT</a>
				</h4>
				</div><!-- END of RAREFACTION_HEADING -->
				<div id="COLLAPSE_RAREFACTION_BODY" class="panel-collapse collapse">
					<div id="RAREFACTION_BODY" class="panel-body">\n
						""" + rarefaction_html + """
					</div><!-- END of RAREFACTION_BODY -->
				</div><!-- END of COLLAPSE_RAREFACTION_BODY -->
			</div><!-- END of RAREFACTION_GROUP -->
			<!--END OF RAREFACTION-->

			<!--ALPHA DIVERTSITY SUMMARY-->
			<div id="ALPHA_DIVERTSITY_SUMMARY" class="panel panel-default">
				<div id="ALPHA_DIVERTSITY_SUMMARY_HEADING" class="panel-heading">
				<h4 class="panel-title">
					<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY">ALPHA DIVERSITY SUMMARY</a>
				</h4>
				</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_HEADING -->
				<div id="COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-collapse collapse">
					<div id="ALPHA_DIVERTSITY_SUMMARY_BODY" class="panel-body">
						<!-- Table -->
						<div class="table-responsive">
							<table id="alpha_diversity_summary" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%">
								""" + summary_table_header + """
								""" + summary_table_body + """
							</table>
						</div>
						<!-- END of TABLE -->
						
					</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_BODY -->
				</div><!-- END of COLLAPSE_ALPHA_DIVERTSITY_SUMMARY_BODY -->
			</div><!-- END of ALPHA_DIVERTSITY_SUMMARY_GROUP -->
			<!--END OF ALPHA DIVERTSITY SUMMARY-->

			<!--BIOMARKER_DISCOVERY-->
			<div id="BIOMARKER_DISCOVERY_GROUP" class="panel panel-default">
				<div id="BIOMARKER_DISCOVERY_HEADING" class="panel-heading">
				<h4 class="panel-title">
					<a data-toggle="collapse" data-parent="#MAIN_PANEL_GROUP" href="#COLLAPSE_BIOMARKER_DISCOVERY_BODY">BIOMARKER_DISCOVERY</a>
				</h4>
				</div><!-- END of BIOMARKER_DISCOVERY_HEADING -->
				<div id="COLLAPSE_BIOMARKER_DISCOVERY_BODY" class="panel-collapse collapse">
					<div id="BIOMARKER_DISCOVERY_BODY" class="panel-body">
						<div id="BIOMARKER_DISCOVERY_ACCORDION" class="panel-group">
									""" + biomarker_discovery_string + """
						</div><!-- END of BIOMARKER_DISCOVERY_ACCORDION -->
					</div><!-- END of BIOMARKER_DISCOVERY_BODY -->
				</div><!-- END of COLLAPSE_BIOMARKER_DISCOVERY_BODY -->
			</div><!-- END of BIOMARKER_DISCOVERY_GROUP -->
			<!--END OF BIOMARKER_DISCOVERY-->


				</div><!--END OF MAIN_PANEL_GROUP-->
			</div>
		</div><!--END OF MAIN CONTAINER-->
	</body>
		<script type="text/javascript">
			function getSelectedChbox(frm) {
				var selectedbox = [];
				var unselectedbox = [];
				var inpfields = frm.getElementsByTagName('input');
				var nr_inpfields = inpfields.length;
				for(var i = 0; i < nr_inpfields; i++) {
				if(inpfields[i].type == 'checkbox' && inpfields[i].checked == true) selectedbox.push(inpfields[i].value);
				if(inpfields[i].type == 'checkbox' && inpfields[i].checked == false) unselectedbox.push(inpfields[i].value);
				}
				return [selectedbox, unselectedbox];
			}
			document.getElementById('submit').onclick = function(){
				var box = getSelectedChbox(this.form);
				var selectedbox = box[0];
				
				var unselectedbox = box[1];
				
				for(var i=0; i < selectedbox.length; i++){
					var div_to_hide = document.getElementById("div_abundance_".concat(selectedbox[i]));
					if (div_to_hide.style.display !== 'none') {
					div_to_hide.style.display = 'none';
					};
					var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(selectedbox[i]));
					if (div_to_hide.style.display !== 'none') {
					div_to_hide.style.display = 'none';
					};
					var div_to_hide = document.getElementById("div_rarefaction_".concat(selectedbox[i]));
					if (div_to_hide.style.display !== 'none') {
					div_to_hide.style.display = 'none';
					};
					var div_to_hide_list = [];
					div_to_hide_list = document.getElementsByName("div_alpha_".concat(selectedbox[i]));
					console.log(div_to_hide_list);
					for (var j=0; j < div_to_hide_list.length; j++){
						if (div_to_hide_list[j].style.display !== 'none') {
							div_to_hide_list[j].style.display = 'none';
						};
					}
				}
				for(var i=0; i < unselectedbox.length; i++){
					var div_to_hide = document.getElementById("div_abundance_".concat(unselectedbox[i]));
					console.log(div_to_hide);
					if (div_to_hide.style.display == 'none') {
					div_to_hide.style.display = 'block';
					};
					var div_to_hide = document.getElementById("div_bacterial_abundance_".concat(unselectedbox[i]));
					if (div_to_hide.style.display == 'none') {
					div_to_hide.style.display = 'block';
					};
					var div_to_hide = document.getElementById("div_rarefaction_".concat(unselectedbox[i]));
					if (div_to_hide.style.display == 'none') {
					div_to_hide.style.display = 'block';
					};
					var div_to_hide_list = [];
					div_to_hide_list = document.getElementsByName("div_alpha_".concat(unselectedbox[i]));
					console.log(div_to_hide_list);
					for (var j=0; j < div_to_hide_list.length; j++){
						if (div_to_hide_list[j].style.display == 'none') {
							div_to_hide_list[j].style.display = 'table-row';
						};
					}
				}
			}
		</script>
		<script type="text/javascript">
			""" + abundance_script + """
			""" + bacterial_abundance_script + """
			""" + rarefaction_script + """
		</script>
	</html>
		"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def alpha_html(alpha_path, absname, design_file, pca_calc_list_1, pca_calc_list_2, summary_table_header, summary_table_body, signature_header_table, signature_body_table):
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')

	javastring = ''
	rarefaction_css = ''
	rarefaction_html = '<td>\n<div id="wrapper">'
	abundance_css = ''
	abundance_html = '<td>\n<div id="wrapper">'
	bacterial_abundance_css = ''
	bacterial_abundance_html = '<td>\n<div id="wrapper">'
	pca_css = ''
	pca_html_1 = ''
	pca_html_2 = ''
	plotly_script = '<script>'

	for each_design in design_dictionary_rev:
		plotly_script += "Plotly.newPlot('myDiv', rarefaction_" + each_design + ");"
		rarefaction_css += """
		#rarefaction_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		rarefaction_html += """
		<div id='rarefaction_""" + each_design + """'>
		</div>
		"""
		abundance_css += """
		#abundance_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		abundance_html += """

		<div id='abundance_""" + each_design + """'>
		</div>
		"""
		bacterial_abundance_css += """
		#bacterial_abundance_""" + each_design + """{
		min-width: 1000px;
		max-width: 1000px;
		
		border: 0.01px outset grey;
		margin: 1px auto;
		float:left;

		}
		"""
		bacterial_abundance_html += """

		<div id='bacterial_abundance_""" + each_design + """'>
		</div>
		"""

	rarefaction_html += '</div>\n</td>\n'
	abundance_html += '</div>\n</td>\n'
	bacterial_abundance_html += '</div>\n</td>\n'
	plotly_script = '</script>'
	#collectors_css = ''
	#collectors_html = ''
	pca_jscript_1 = ''
	for calc in pca_calc_list_1:
		calc_name = calc.split(':')[0]
		pca_jscript_1 += """<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_""" + calc_name + """_pca.js"></script>\n"""
		pca_css += """
		#PCOA_"""+calc_name+""" {
			min-width: 1000px;
			max-width: 1000px;
			border: 0.01px outset black;
			margin: 1px auto;
			float:left;
			}\n"""
		pca_html_1 += """<td>
		
								<div id='PCOA_""" + calc_name + """'></div>
								<button id="2D_PCOA_""" + calc_name + """"class="autocompare">switch to 2D</button>
								
		</td>\n"""
	pca_jscript_2 = ''
	for calc in pca_calc_list_2:
		calc_name = calc.split(':')[0]
		pca_jscript_2 += """<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_""" + calc_name + """_pca.js"></script>\n"""
		pca_css += """
		#PCOA_"""+calc_name+""" {
			min-width: 1000px;
			max-width: 1000px;
			border: 0.01px outset black;
			margin: 1px auto;
			float:left;
			}\n"""
		pca_html_2 += """<td>
		
								<div id='PCOA_""" + calc_name + """'></div>
								<button id="2D_PCOA_""" + calc_name + """"class="autocompare">switch to 2D</button>
								
		</td>\n"""
	
	javastring += """
		<!DOCTYPE html>
			<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
			<meta charset="utf-8"/>
				<head>
					<script type="text/javascript" src="./""" + absname + """_js/jquery-2.2.3.min.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/highcharts.src.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/highcharts-3d.src.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/exporting.src.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/export-csv.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/jquery.tablesorter.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_rarefaction_ploty.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_abundance.js"></script>
					<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_bacterial_abundance.js"></script>
					""" + pca_jscript_1 + """
					""" + pca_jscript_2 + """

					<script type="text/javascript">
						$(document).ready(function() { 
						// call the tablesorter plugin 
						$("table").tablesorter(); 
					});
					""" + plotly_script + """

					</script>
					<style type="text/css">
					#wrapper{
						max-width: 5000px;
						max-height: 5000px;
						min-width: 1000px;
						border: 1px solid black;
						overflow-x: auto;
						overflow-y: auto;
					}
					#wrapper_pca{
						max-width: 15000px;
						max-height: 5000px;
						min-width: 1000px;
						border: 1px solid black;
						overflow-x: auto;
						overflow-y: auto;
					}
					
					#shared_table
					{
						font-family: "Avenir";
						font-size: 12px;
						background: #fff;
						margin-top: 50px;
						width: 1800px;
						height: 400px;
						border: 1px solid #6678b1;

						border-collapse: collapse;
						text-align: center;
						
					}
					#shared_table th
					{
						font-size: 14px;
						font-weight: normal;
						color: #039;
						padding: 10px 8px;
						border-bottom: 2px solid #6678b1;
						text-align: center;
					}
					#shared_table td
					{
						border-bottom: 1px solid #ccc;
						color: #669;
						padding: 6px 8px;
						text-align: center;
					}
					#shared_table tbody tr:hover td
					{
						color: #009;
						background-color: yellow;
					}
					#shared_table .odd
					{
						background: #e8edff; 
					}
					#shared_table .metadata
					{
						border-bottom: 2px solid #ccc;
						color: #669;
						padding: 6px 8px;
						text-align: center; 
					}
					""" + bacterial_abundance_css + """
					""" + abundance_css + """
					""" + rarefaction_css + """
					""" + pca_css + """
					
					
					
					</style>
				</head>
				<body>
					<table align="center">
						<tr>
						""" + bacterial_abundance_html + """
						""" + abundance_html + """
						""" + rarefaction_html + """
						

						</tr>
					</table>
					<div id="wrapper_pca">
					<table align="center">
						<tr>
						
						""" + pca_html_1 + """
						
						</tr>
					</table>
					</div>
					<div id="wrapper_pca">
					<table align="center">
						<tr>
						
						""" + pca_html_2 + """
						
						</tr>
					</table>
					</div>
					
					<table id="shared_table" summary="shared_table">
					
					""" + signature_header_table + """
					""" + signature_body_table + """
					
					
					</table>
					

					<table id="shared_table" summary="shared_table" class="tablesorter">
					
					""" + summary_table_header + """
					""" + summary_table_body + """
					
					
					</table>
					

					
				</body>
			</html>
		
	"""
	f = open(alpha_path + absname + '_ANALYSIS.html', 'w')
	f.write(javastring)
	f.close()


def rarefaction_javascript_plotly(rarefaction_design_dict, plotly_design_data_dict, name, design_file, js_path):
	rarefaction_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########     RAREFACTION USING PLOTLY      #######################
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:
		rarefaction_javascript_string += rarefaction_design_dict[each_design]
		rarefaction_javascript_string += "\tvar data_rarefaction_" + each_design + " = [" + plotly_design_data_dict[each_design] + "];\n"
		rarefaction_javascript_string += """
		var layout_rarefaction_""" + each_design + """ = {
			title: "Rarefaction curve plot of """ + each_design + """ ",
			titlefont: {
						family: 'Avenir',
						size: 20 },
			showlegend: true,
			legend: {
				font: {
					family: 'Avenir',
							},
				borderwidth: 0,
				traceorder: 'normal'

			},
			autosize: true,
			width: 1000,
			height: 500,
			hovermode:'closest',

			xaxis: {
				title: 'Number of sequences sampled',
				titlefont: {
					family: 'Avenir',
					size: 15 },
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			yaxis: {
				title: 'Number of detected species',
				titlefont: {
					family: 'Avenir',
					size: 15 },
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			

			};
		"""
		rarefaction_javascript_string += """\n
		//####################################################################
		//####################################################################
		//####################################################################
		"""
	
	f = open(js_path + name + '_' + 'rarefaction_plotly.js', 'w')
	f.write(rarefaction_javascript_string)
	f.close()


def abundance_javascript_plotly(abundance_design_dict, abundance_otu_dict, name, design_file, js_path):
	abundance_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########     Abundance table USING PLOTLY      ###################
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	minimumwidth = '1000'

	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:
		if len(design_dictionary_rev[each_design]) <= 30:
			minimumwidth = '1000'
		else:
			minimumwidth = str(len(design_dictionary_rev[each_design]) * 30)

		abundance_javascript_string += abundance_design_dict[each_design]
		abundance_javascript_string += "\tvar data_abundance_" + each_design + " = [" + abundance_otu_dict[each_design] + "];\n"
		abundance_javascript_string += """
		var layout_abundance_""" + each_design + """ = {
			title: "Sample abundance of """ + each_design + """ ",
			titlefont: {
						family: 'Avenir',
						size: 20 },
			showlegend: false,
			autosize: false,
			width: """ + minimumwidth + """,
			height: 500,
			hovermode:'closest',
			bargap: 0.25,
			bargroupgap: 0.3,
			barmode: 'stack',
			
			xaxis: {
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			yaxis: {
				title: 'Total abundance value',
				titlefont: {
					family: 'Avenir',
					size: 15},
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			margin: {
				l: 50,
				r: 50,
				b: 200,
				t: 100,
				pad: 0

			},
			

			};
		"""
		abundance_javascript_string += """\n
		//####################################################################
		//####################################################################
		//####################################################################
		"""
	
	f = open(js_path + name + '_abundance_plotly.js', 'w')
	f.write(abundance_javascript_string)
	f.close()


def bacterial_abundance_javascript_plotly(abundance_design_dict, abundance_otu_dict, bact_list, name, design_file, js_path):
	abundance_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########     BACTERIAL Abundance table USING PLOTLY      #########
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	minimumwidth = '1000'

	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:
		if len(bact_list) <= 30:
			minimumwidth = '1000'
		else:
			minimumwidth = str(len(bact_list) * 30)

		abundance_javascript_string += abundance_design_dict[each_design]
		abundance_javascript_string += "\tvar data_bacterial_abundance_" + each_design + " = [" + abundance_otu_dict[each_design] + "];\n"
		abundance_javascript_string += """
		var layout_bacterial_abundance_""" + each_design + """ = {
			title: "Bacterial abundance of """ + each_design + """ ",
			titlefont: {
						family: 'Avenir',
						size: 20 },
			showlegend: false,
			autosize: false,
			width: """ + minimumwidth + """,
			height: 500,
			hovermode:'closest',
			bargap: 0.25,
			bargroupgap: 0.3,
			barmode: 'stack',
			
			xaxis: {
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			yaxis: {
				title: 'Total abundance value',
				titlefont: {
					family: 'Avenir',
					size: 15},
				autorange: true,
				tickfont: {
					family: 'Avenir',
					},
			},
			margin: {
				l: 50,
				r: 50,
				b: 200,
				t: 100,
				pad: 0

			},
			

			};
		"""
		abundance_javascript_string += """\n
		//####################################################################
		//####################################################################
		//####################################################################
		"""
	
	f = open(js_path + name + '_bacterial_abundance_plotly.js', 'w')
	f.write(abundance_javascript_string)
	f.close()


def abundance_javascript(abundance_columns_dict, abundance_category_dict, design_file, name, js_path):
	minimumwidth = 0
	sample_dict = {}
	sample_dict = design_dict_maker(design_file)
	sample_count = sample_dict.keys()
	legend_status = 'false'
	"""
	if len(sample_count) > 5:
		legend_status = 'false'
	else:
		legend_status = 'true'
	"""
	abundance_javascript_string = ''
	abundance_javascript_string += """
	//###################################################################
	//###################################################################
	//###################################################################
	//########################## Abundance      #########################
	//###################################################################
	//###################################################################
	//###################################################################
	//###################################################################

	"""
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:
		if len(abundance_category_dict) < 10:
			minimumwidth = 1000
		else:
			minimumwidth = len(abundance_category_dict) * 5

		abundance_javascript_string += """
	$(function () {
					$('#abundance_""" + each_design + """').highcharts({
						chart: {
							width: """ + str(minimumwidth) + """,
							type: 'column',
							style: {
							fontFamily: '"Avenir"', // default font
						}
						},
						title: {
							text: 'Bacterial abundance of """ + each_design + """',
							align: 'left'
						},
						
						credits: {enabled: false},
						legend: {
							enabled: """ + legend_status + """,
							layout: "horizontal"
						},
						xAxis: {
							categories:['""" + abundance_category_dict[each_design] + """']
							
						},
						yAxis: {
							visible : true,
							allowDecimals: false,
							title: {
								text: 'Total abundance value'
							}
						},
						tooltip: {
									
							hideDelay: 10,
							shared: false
						},
						plotOptions: {
							column: {
								stacking: 'normal'
							},
							series: {
								stickyTracking: false,
								pointWidth: 20
							}
						},
				

						series: [""" + abundance_columns_dict[each_design] + """]
					});
				});
	"""
	f = open(js_path + name + '_' + 'abundance.js', 'w')
	f.write(abundance_javascript_string)
	f.close()


def bacterial_abundance_javascript(abundance_columns_dict, abundance_category_dict, design_file, name, js_path):
	
	sample_dict = {}
	sample_dict = design_dict_maker(design_file)
	sample_count = sample_dict.keys()
	legend_status = 'false'
	"""
	if len(sample_count) > 5:
		legend_status = 'false'
	else:
		legend_status = 'true'
	"""
	abundance_javascript_string = ''
	abundance_javascript_string += """
	//###################################################################
	//###################################################################
	//###################################################################
	//########################## Abundance      #########################
	//###################################################################
	//###################################################################
	//###################################################################
	//###################################################################

	"""
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	for each_design in design_dictionary_rev:

		abundance_javascript_string += """
	$(function () {
					$('#bacterial_abundance_""" + each_design + """').highcharts({
						chart: {
							width: """ + str(len(abundance_category_dict[each_design])*5) + """,
							type: 'column',
							style: {
							fontFamily: '"Avenir"', // default font
						}
						},
						title: {
							text: 'Bacterial abundance of """ + each_design + """',
							align: 'left'
						},
						
						credits: {enabled: false},
						legend: {
							enabled: """ + legend_status + """,
							layout: "horizontal"
						},
						xAxis: {
							crosshair: true,
							categories:['""" + abundance_category_dict[each_design] + """']
			
						},
						yAxis: {
							visible : true,
							allowDecimals: false,
							title: {
								text: 'Total abundance value'
							}
						},
						tooltip: {
									
							hideDelay: 10,
							shared: false
						},
						plotOptions: {
							column: {
								stacking: 'normal'
							},
							series: {
								stickyTracking: false,
								pointWidth: 20
							}
						},
				

						series: [""" + abundance_columns_dict[each_design] + """]
					});
				});
	"""
	f = open(js_path + name + '_' + 'bacterial_abundance.js', 'w')
	f.write(abundance_javascript_string)
	f.close()


def pca_javascript(pca_string, design_file, xmax, xmin, ymax, ymin, zmax, zmin, calc, absname, js_path):
	sample_dict = {}
	sample_dict = design_dict_maker(design_file)
	sample_count = sample_dict.keys()
	if len(sample_count) > 5:
		legend_status = 'false'
	else:
		legend_status = 'true'
	calc_name = calc.split(':')[0]
	xmax += (xmax/2)
	ymax += (ymax/2)
	zmax += (zmax/2)
	xmin += (xmin/2)
	ymin += (ymin/2)
	zmin += (zmin/2)
	pca_javascript_string = ''
	pca_javascript_string += '''

	//####################################################################
	//####################################################################
	//####################################################################
	//###############################     PCA      #######################
	//####################################################################
	//####################################################################
	//####################################################################
				$(function () {

				// Give the points a 3D feel by adding a radial gradient
				Highcharts.getOptions().colors = $.map(Highcharts.getOptions().colors, function (color) {
				return {
						radialGradient: {
							cx: 0.4,
							cy: 0.3,
							r: 0.5
						},
						stops: [
							[0, color],
							[1, Highcharts.Color(color).brighten(-0.2).get('rgb')]
						]
					};
				});

				// Set up the chart
				var defaultOptions={
					chart: {
						type: "scatter",
						reflow: true,
						//margin: 100,
						//width: 1000,
						//height: 400,
						style: {
							fontFamily: '"Avenir"', // default font
							},
						options3d: {
							enabled: true,
							alpha: 10,
							beta: 30,
							depth: 250,
							viewDistance: 5,

							frame: {
								bottom: { size: 1, color: 'rgba(0,0,0,0.02)' },
								back: { size: 1, color: 'rgba(0,0,0,0.04)' },
								side: { size: 1, color: 'rgba(0,0,0,0.06)' }
							}
						}
					},
					title: {
						text: ' Principle component analysis of ''' + absname + ''' '
					},
					subtitle: {
						text: 'Method: ''' + calc + ''' '
					},
					plotOptions: {
						scatter: {
							width: 10,
							height: 10,
							depth: 10
						}
					},
					yAxis: {
						min: ''' + str(ymin) + ''',
						max: ''' + str(ymax) + ''',
						gridLineWidth: 1,
						title: {
									enabled: true,
									text: "Y axis"
								},
					},
					xAxis: {
						min: ''' + str(xmin) + ''',
						max: ''' + str(xmax) + ''',
						gridLineWidth: 1,
						title: {
									enabled: true,
									text: "X axis"
								},
					},
					zAxis: {
						min: ''' + str(zmin) + ''',
						max: ''' + str(zmax) + ''',
						gridLineWidth: 1,
						title: {
									enabled: true,
									text: "Z axis"
								},
					},
					legend: {
						enabled: ''' + legend_status + ''',
						layout: "horizontal"
					},
					credits: {
						enabled: false
					},
					tooltip: {
									
						headerFormat: '<b><span style="color:{series.color}">\u25CF</span><b>{series.name}</b><br>',
						//pointFormat: '<b>{point.x}</b>: {point.y}',
						valueDecimals: 5,
						shared: false
						
										
					},
					series: [''' + pca_string + ''']
				}
				$("#PCOA_''' + calc_name + '''").highcharts(defaultOptions);
				var chart = $("#PCOA_''' + calc_name + '''").highcharts();


				$("#PCOA_''' + calc_name + '''").bind("mousedown.hc touchstart.hc", function (e) {
					e = chart.pointer.normalize(e);
			
					var posX = e.pageX,
						posY = e.pageY,
						alpha = chart.options.chart.options3d.alpha,
						beta = chart.options.chart.options3d.beta,
						newAlpha,
						newBeta,
						sensitivity = 5; // lower is more sensitive
			
					$(document).bind({
						"mousemove.hc touchdrag.hc": function (e) {
							// Run beta
							newBeta = beta + (posX - e.pageX) / sensitivity;
							newBeta = Math.min(100, Math.max(-100, newBeta));
							chart.options.chart.options3d.beta = newBeta;
			
							// Run alpha
							newAlpha = alpha + (e.pageY - posY) / sensitivity;
							newAlpha = Math.min(100, Math.max(-100, newAlpha));
							chart.options.chart.options3d.alpha = newAlpha;
			
							chart.redraw(false);
						},
						"mouseup touchend": function () {
							$(document).unbind(".hc");
						}
					});
				});
				//end of mouse event rotation
				//###############switch  Dimension
				var dimension=false;
				
				$("#2D_PCOA_''' + calc_name + '''").click(function () {
					if(!dimension){
						chart.destroy();
						defaultOptions.chart.zoomType="x";
						defaultOptions.chart.panning = true;
						defaultOptions.chart.panKey = "shift";
						defaultOptions.legend.layout="horizontal";
						//defaultOptions.subtitle.text='tax level:family';

						//defaultOptions.legend.floating=true;
						defaultOptions.chart.options3d.enabled = dimension;
						chart = new Highcharts.Chart(defaultOptions);
				$("#2D_PCOA_''' + calc_name + '''").html("switch to 3D");//change the title of button
					}else if(dimension){
						chart.destroy();
						defaultOptions.chart.zoomType="";
						defaultOptions.chart.options3d.enabled = dimension;
						//defaultOptions.subtitle.text="PCA 3D plot";
						chart = new Highcharts.Chart(defaultOptions);

						
						$("#2D_PCOA_''' + calc_name + '''").html("switch to 2D");
						

						}
							 dimension = !dimension
					
					
					
					});
				});
		'''
	f = open(js_path + absname + '_' + calc_name + '_pca.js', 'w')
	f.write(pca_javascript_string)
	f.close()


def ncbi_otu_linker(otu):
	ncbi_link = "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=10000&dbkey=all&keywords="
	search_term = ''
	otu_list = otu.split(';')
	otu_term = otu_list[-2]
	if 's_' in otu_term:
		search_term += otu_list[-3].split('g_')[1]
		search_term += '+' + otu_term.split('s_')[1]
	else:
		search_term += otu_term.split('_')[1]
	
	ncbi_link += search_term
	return ncbi_link


def signature_table(signature_dict, design_file, absname):
	dict_of_design = design_dict_maker(design_file, 'reverse')
	group_list = dict_of_design.keys()
	group_string = ''
	for i in group_list:
		group_string += i + '//'
	group_string = slugify(group_string)
	tab_name = group_string.replace('//', ']-[')
	tab_name = tab_name[:-2]
	tab_name = '[' + tab_name
	data_parent = 'TIER_'
	Tier_name = 'tier' + str(tab_name.count(']-[') + 1)
	data_parent += str(tab_name.count(']-[') + 1) + '_BODY_ACCORDION'
	group_string = group_string.replace('//', '_vs_')
	group_string = group_string[:-4]
	group_string = Tier_name + '_' + group_string
	table_header = """
			<!--""" + group_string + """-->
				<div id='biomarker_discovery_""" + group_string + """' class="panel panel-default">
					<div id='""" + group_string + """_HEADING' class="panel-heading">
						<h4 class="panel-title">
							<a data-toggle="collapse" data-parent="#""" + data_parent + """" href="#COLLAPSE_""" + group_string + """_BODY">""" + tab_name + """</a>
							<a data-toggle="collapse" href="#INFO_""" + group_string + """"><span class="fa fa-download fa-1x pull-right"></span></a>
								<div id="INFO_""" + group_string + """" class="panel-collapse collapse" role="tabpanel">
								<div class="form-inline" role="form">
								<br>
								<a type="button" class="btn btn-default btn-sm" href="./""" + absname + """_DATA/Abundance_file_filtered_""" + group_string + """.txt" download="Abundance_file_filtered_""" + group_string + """.txt">Abundance table filtered <span class="fa fa-download fa-1x right"></span></a>
								<a type="button" class="btn btn-default btn-sm" href="./""" + absname + """_DATA/Abundance_file_raw_""" + group_string + """.txt" download="Abundance_file_raw_""" + group_string + """.txt">Abundance table raw <span class="fa fa-download fa-1x right"></span></a>
								<a type="button" class="btn btn-default btn-sm" href="./""" + absname + """_DATA/Orange_file_filtered_""" + group_string + """.txt" download="Orange_file_filtered_""" + group_string + """.txt">Orange table filtered <span class="fa fa-download fa-1x right"></span></a>
								<a type="button" class="btn btn-default btn-sm" href="./""" + absname + """_DATA/Orange_file_raw_""" + group_string + """.txt" download="Orange_file_raw_""" + group_string + """.txt">Orange table raw <span class="fa fa-download fa-1x right"></span></a>
								</div>
								</div>
						</h4>
					</div><!-- END of """ + group_string + """_HEADING -->
					<div id="COLLAPSE_""" + group_string + """_BODY" class="panel-collapse collapse">
						<div id='""" + group_string + """_BODY' class="panel-body">\n
				"""

	table_header += """
									<!-- Table -->
										<div class="table-responsive">
	"""
	table_header += '						<table id="BIOMARKER_DISCOVERY_TABLE" class="table table-bordered table-hover small" cellspacing="0" width="100%">\n'
	table_header += '							<thead>\n'
	table_header += '								<tr>\n'
	table_header += '									<th><p class="text-center">Analysis Name</p></th>\n'
	#table_header += '\t\t\t\t<th>Group name</th>\n'
	table_header += '							<th><p class="text-center">Significant OTUs</p></th>\n'
	table_header += '							<th><p class="text-center">Test value</p></th>\n'
	table_header += '							<th><p class="text-center">P value</p></th>\n'
	table_header += '								</tr>\n'
	table_header += '							</thead>\n'
	#//////////////////////////////////////////////////////

	table_body = ''
	table_body += '								<tbody>\n'
	
	for test, bact in signature_dict.iteritems():
		if len(bact) == 0:
			table_body += '								<tr>\n'
			table_body += '									<th><p class="text-center">' + test + '</p></th>\n'
			#table_body += '\t\t\t\t<td>' + group_string + '</td>\n'
			table_body += '									<td><p class="text-center">No Significant OTUs</p></td>\n'
			table_body += '									<td><p class="text-center">N/A</p></td>\n'
			table_body += '									<td><p class="text-center">N/A</p></td>\n'
			table_body += '								</tr>\n'
		else:
			bact_list = sorted(signature_dict[test], key=signature_dict[test].get, reverse=True)
			collapsed_bact_list = []
			collapse_otus(bact_list, collapsed_bact_list)
			#table_body += '\t\t<tr>\n'
			table_body += '								<th rowspan="' + str(len(collapsed_bact_list) + 1) + '" scope="rowgroup" class="metadata"><p class="text-center">' + test + '</p></th>\n'
			for each_otu in collapsed_bact_list:
				table_body += '								<tr>\n'
				#table_body += '\t\t\t\t<td>' + group_string + '</td>\n'
				table_body += '									<td><p class="text-left"><a href="' + each_otu + '" style="text-decoration:none;" target="_blank">' + each_otu + '</p></td>\n'
				#table_body += '									<td><p class="text-left"><a href="' + ncbi_otu_linker(each_otu) + '" style="text-decoration:none;" target="_blank">' + each_otu + '</p></td>\n'
				table_body += '									<td><p class="text-center">' + test + ':' + signature_dict[test][each_otu][0] + '</p></td>\n'
				table_body += '									<td><p class="text-center">' + signature_dict[test][each_otu][1] + '</p></td>\n'
				table_body += '								</tr>\n'
	table_body += '								</tbody>\n'
	table_body += '							</table>\n'
	table_body += '						</div><!-- END of responsive table -->\n'
	#table_body = ''
	
	table_footer = """
			</div><!-- END of """ + group_string + """_BODY -->
		</div><!-- END of COLLAPSE_""" + group_string + """_BODY -->
	</div><!-- END of """ + group_string + """_GROUP -->
	<!--END OF """ + group_string + """-->
	"""
	table_body += table_footer
	return(table_header, table_body)


def biomarker_discovery_html_maker(permuted_file_dict, permuted_table_dict):
	
	biomarker_discovery_string = ''
	tier_group = permuted_table_dict.keys()
	print tier_group
	for each_tier in tier_group:
		biomarker_discovery_string += """
			<div id='""" + each_tier + """' class="panel panel-primary">
				<div id='""" + each_tier + """_HEADING' class="panel-heading">
					<h4 class="panel-title">
						<a data-toggle="collapse" data-parent="#BIOMARKER_DISCOVERY_ACCORDION" href="#COLLAPSE_""" + each_tier + """_BODY" style="color: rgb(200,200,200)">TIER """ + each_tier.split('_')[1] + """</a>
					</h4>
				</div><!-- END of """ + each_tier + """_HEADING -->
				<div id="COLLAPSE_""" + each_tier + """_BODY" class="panel-collapse collapse">
					<div id='""" + each_tier + """_BODY' class="panel-body">
						<div id='""" + each_tier + """_BODY_ACCORDION' class="panel-group">
		"""
		print permuted_table_dict[each_tier].keys()
		for each_set_of_analysis in permuted_table_dict[each_tier].keys():
			biomarker_discovery_string += permuted_table_dict[each_tier][each_set_of_analysis][0]
			biomarker_discovery_string += '\n'
			biomarker_discovery_string += permuted_table_dict[each_tier][each_set_of_analysis][1]
		biomarker_discovery_string += """
					</div><!-- <!-- END of """ + each_tier + """_BODY_ACCORDION -->
				</div><!-- <!-- END of """ + each_tier + """_BODY -->
			</div><!-- <!-- END of COLLAPSE_""" + each_tier + """_BODY-->
		</div><!-- <!-- END of """ + each_tier + """-->
		"""

	return biomarker_discovery_string


def pca_javascript_plotly(pca_string, pca_sample_string, PCA_name, design_file, xmax, xmin, ymax, ymin, zmax, zmin, calc, absname, js_path):
	calc_name = calc.split(':')[0]
	xmax += (xmax / 2)
	ymax += (ymax / 2)
	zmax += (zmax / 2)
	xmin += (xmin / 2)
	ymin += (ymin / 2)
	zmin += (zmin / 2)
	xmax = str(xmax)
	ymax = str(ymax)
	zmax = str(zmax)
	xmin = str(xmin)
	ymin = str(ymin)
	zmin = str(zmin)
	sample_count = pca_sample_string.split(',')
	if len(sample_count) > 25:
		legend_status = 'false'
	else:
		legend_status = 'true'
	pca_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########     BACTERIAL Abundance table USING PLOTLY      #########
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	pca_javascript_string += pca_string
	pca_javascript_string += "\tvar data_pca_" + PCA_name + '_' + calc_name + " = [" + pca_sample_string + "];\n"
	
	pca_javascript_string += """
	var layout_pca_""" + PCA_name + '_' + calc_name + """ = {
		title: 'PCA - """ + calc + """ ',
		titlefont: {
			family: 'Avenir',
			size: 20
		},
		//showlegend: """ + legend_status + """,
		showlegend: false,
		traceorder:'normal',
		autosize: false,
		height: 500,
		width: 1000,
		dragmode: false,
		
		margin: {
			l: 25,
			r: 25,
			b: 50,
			t: 50,
			pad: 1
		},
		scene: {
			xaxis: {
				title: 'PC1',
				titlefont: {
					family: 'Avenir',
					size: 15 },
				range:[""" + xmin + """, """ + xmax + """],
				autorange: false,
				showbackground: true,
				nticks:10,
				autotick:true,
				fixedrange: false,
				gridcolor: "white",
				linecolor: "white",
				zerolinecolor: "white",
				gridwidth:4,
				backgroundcolor: "rgba(0,0,0,0.02)",
				showspikes: false,
			},
			yaxis: {
				title: 'PC2',
				titlefont: {
					family: 'Avenir',
					size: 15
				},
				range:[""" + ymin + """, """ + ymax + """],
				showbackground: true,
				fixedrange: true,
				nticks:5,
				autorange: false,
				gridcolor: "white",
				linecolor: "white",
				gridwidth:4,
				zerolinecolor: "white",
				backgroundcolor: "rgba(0,0,0,0.03)",
				showspikes: false,
				type:"linear"
			},
			zaxis: {
				title: 'PC3',
				titlefont: {
					family: 'Avenir',
					size: 15},
				range:[""" + zmin + """, """ + zmax + """],
				showbackground: true,
				fixedrange: true,
				gridcolor: "white",
				linecolor: "white",
				zerolinecolor: "white",
				gridwidth:4,
				backgroundcolor: "rgba(0,0,0,0.04)",
				nticks:5,
				autorange: false,
				showspikes: false,
			},
			//cameraposition:[[0.4, 0.5, -0.3, -0.4], [0.1, 0, -0.1], 4],
			aspectratio:{x:1,y:1,z:1},
			aspectmode:'manual'
		}

	};
	"""
	pca_javascript_string += """\n
	//####################################################################
	//####################################################################
	//####################################################################
	"""
	
	f = open(js_path + absname + '_' + PCA_name + '_' + calc_name + '_pca_plotly.js', 'w')
	f.write(pca_javascript_string)
	f.close()


def pca_html_maker(pca_string_dict, PCA_name, calc, absname):
	calc_name = calc.split(':')[0]
	tag_name = '['
	tag_name += PCA_name.replace('_vs_', ']-[')
	tag_name += ']'
	pca_html = ''
	pca_html += """
				<div id='biomarker_discovery_""" + PCA_name + """' class="panel panel-primary">
				<div id='""" + PCA_name + """_HEADING' class="panel-heading">
				<h4 class="panel-title">
				<a data-toggle="collapse" data-parent="#PCA_3D_PLOTTER_ACCORDION" href="#COLLAPSE_""" + PCA_name + """_BODY" style="color: rgb(200,200,200)">""" + tag_name + """</a>
				</h4>
				</div><!-- END of """ + PCA_name + """_HEADING -->
				<div id="COLLAPSE_""" + PCA_name + """_BODY" class="panel-collapse collapse">
					<div id='PCOA_""" + PCA_name + '_' + calc_name + """'></div>\n
					</div>
					</div>\n"""
	pca_script = """<script type="text/javascript" src="./""" + absname + """_JS/""" + absname + """_""" + PCA_name + '_' + calc_name + """_pca_plotly.js"></script>\n"""
	pca_plotly_script = """
	Plotly.newPlot('PCOA_""" + PCA_name + '_' + calc_name + """', data_pca_""" + PCA_name + '_' + calc_name + """, layout_pca_""" + PCA_name + '_' + calc_name + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['sendDataToCloud', 'tableRotation', 'resetCameraDefault3d', 'resetCameraLastSave3d','hoverClosest3d'],
		});
	"""
	pca_string_dict[PCA_name] = (pca_html, pca_script, pca_plotly_script)
	return True




# ################################### END SCRIPT ############################### #
if __name__ == "__main__": main(sys.argv[1:])