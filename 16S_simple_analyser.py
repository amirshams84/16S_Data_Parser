# 16S simple parser
# ################################### IMPORTS ################################## #
import sys
import argparse
import time
import os
import multiprocessing
import itertools
import platform
import csv
import random
import errno
import signal
import datetime
import logging as log
import subprocess
import traceback
import shutil
import collections
import math
import decimal



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
#outputdir = "/16S_simple_analyser_results/"
outputdir = "/"
phylotype = "/test_data/ZAC_phylotype.txt"
design = "/test_data/ZAC_design.txt"
processors = multiprocessing.cpu_count()
name = "ZAC"
taxlevel = "family"
remove_sample_file = "false"
keep_sample_file = "false"
normalize = "false"
# ###################################   MAIN   ################################# #


def main(argv):
	report_string = ''
	# ############################# PARSE INPUT ARGUMENTS
	parser = argparse.ArgumentParser()
	# ##################### MAIN FILE PARAMETERS
	main_file = parser.add_argument_group('Main file parameters')
	main_file.add_argument("--phylotype", help="Phylotype file", action='store')
	main_file.add_argument("--shared", help="Phylotype file", action='store')
	main_file.add_argument("--biom", help="Phylotype file", action='store')
	main_file.add_argument("--taxlevel", help="taxonomy level on of [species, genus, family]", action='store')
	main_file.add_argument("--design", help="design file", action='store')
	main_file.add_argument("--jslib", help="javascript library", action='store')
	main_file.add_argument("--outputdir", help="output directory", action='store')
	main_file.add_argument("--execdir", help="executables_directory", action='store')

	# #################### GENERAL PARAMETERS
	general = parser.add_argument_group('general parameters')
	general.add_argument("--processors", help="number of processors assigned", action='store')
	general.add_argument("--name", help="name of output files", action='store')

	# #################### BETA PARAMETERS
	beta = parser.add_argument_group('beta parameters')
	beta.add_argument("--remove_sample_file", help="list of samples to remove", action='store')
	beta.add_argument("--keep_sample_file", help="list of samples to keep", action='store')
	beta.add_argument("--normalize", help="Normalizing using totalgroup method", action='store')

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
	# ########### EXECECUTIVE DIRECTORY CHECKING
	if args.outputdir is None:
		args.outputdir = outputdir
		report_string += "Using Default output directory\n"
	global report_file
	report_file = args.outputdir + "16S_simple_analyser_report.txt"
	check_it_and_remove_it(report_file, True)
	report(report_string)
	args.jslib = jslib
	args.execdir = execdir
	if args.phylotype is None:
		args.phylotype = phylotype
		report("Using Default phylotype file")
		print "Using Default phylotype file"
	if args.design is None:
		args.design = design
		report("Using Default design file")
		print "Using Default design file"

	args.processors = str(processors)
	args.name = name
	args.remove_sample_file = remove_sample_file
	args.keep_sample_file = keep_sample_file
	args.taxlevel = taxlevel
	args.normalize = normalize
	
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
		#report("mothur execution file is: " + mothur_exec_path)
		print "mothur execution file is: ", mothur_exec_path
		#report("Testing mothur executables: ")
		print "Testing mothur executables: "
		flag, stderr = execute_functions(test_mothur, args.processors, args.outputdir, 'multi', 'mothur', mothur_exec_path)
		if flag is False:
			#report("[" + FAILED_MARK + "]")
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

	# ###############################
	# Filtering Shared file based on Design samples
	design_control_file = args.outputdir + args.name + '_design_control_file.txt'
	flag = design_to_control(args.design, design_control_file)
	if flag is False:
		print "ABORTING!!!"
		sys.exit(2)
	updated_shared_file = args.outputdir + args.name + '_UPDATED_shared_file_STEP2.txt'
	updated_design_file = args.outputdir + args.name + '_UPDATED_design_file_STEP2.txt'
	flag = update_shared_design_file(args.shared, args.design, design_control_file, updated_shared_file, updated_design_file, args.name, mothur_exec_path, args.processors, args.outputdir)
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
	# Step7: RAREFACTION CURVE PLOT.
	rarefaction_design_dict = {}
	plotly_design_data_dict = {}
	rarefaction_file = args.outputdir + args.name + '_RAREFACTION_file_STEP7.txt'
	rarefaction_design_dict, plotly_design_data_dict = mothur_rarefaction_curve(args.shared, args.design, rarefaction_file, mothur_exec_path, args.processors, data_path, args.outputdir)
	rarefaction_javascript_plotly(rarefaction_design_dict, plotly_design_data_dict, args.name, args.design, js_path)
	#rarefaction_file_name = args.name + '_RAREFACTION_file_STEP7.txt'
	#rarefaction_plotter(alpha_path, rarefaction_file_name, args.design, args.name)
	print "# ####################################################################################"
	
	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step8: Sample_ABUNDANCE_FILE.
	sample_abundance_columns_dict = {}
	sample_abundance_otu_dict = {}
	simple_abundance_file = args.outputdir + args.name + '_SAMPLE_ABUNDANCE_file_STEP8.txt'
	sample_abundance_columns_dict, sample_abundance_otu_dict = sample_abundance_table_plotly(args.shared, args.design, simple_abundance_file, data_path, args.outputdir)
	sample_abundance_javascript_plotly(sample_abundance_columns_dict, sample_abundance_otu_dict, args.name, args.design, js_path)
	#sample_abundance_file_name = args.name + '_SAMPLE_ABUNDANCE_file_STEP8.txt'
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
	#bacterial_abundance_file_name = args.name + '_BACTERIAL_ABUNDANCE_file_STEP9.txt'
	#bacterial_abundance_plotter(alpha_path, bacterial_abundance_file_name, args.design, args.name)
	print "# ####################################################################################"

	# ##############################################################################################################################################################
	# ##############################################################################################################################################################
	# Step10: ALPHA DIVERSITY SUMMARY FILE.
	alpha_diversity_summary_file = args.outputdir + args.name + '_ALPHA_DIVERSITY_SUMMARY_file_STEP10.txt'
	flag = summary_analysis(args.shared, alpha_diversity_summary_file, mothur_exec_path, args.processors, data_path, args.outputdir)
	alpha_diversity_summary_table_header, alpha_diversity_summary_table_body = summary_table(alpha_diversity_summary_file, args.design)
	#alpha_diversity_summary_file_name = args.name + '_ALPHA_DIVERSITY_SUMMARY_file_STEP10.txt'
	#alpha_summary_plotter(alpha_path, alpha_diversity_summary_file_name, summary_table_header, summary_table_body, args.name, args.design)
	print "# ####################################################################################"

	pca_html_string_dict = {}
	html_plotter(alpha_path, args.name, args.design, alpha_diversity_summary_table_header, alpha_diversity_summary_table_body, pca_html_string_dict)
	remove_mothur_log(os.getcwd())
	#make_archive(args.outputdir)
	print "16S SIMPLE ANALYSER EXECUTION COMPLETED AT ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("16S SIMPLE ANALYSER EXECUTION COMPLETED AT " + time.strftime("%Y-%m-%d %H:%M:%S"))
	#all_plotter(alpha_path, rarefaction_file_name, sample_abundance_file_name, bacterial_abundance_file_name, summary_table_header, summary_table_body, alpha_diversity_summary_file_name, biomarker_discovery_string, pca_html_string_dict, args.name, args.design)
# ################################### ALPHA DIVERSITY SUMMARY FUNCTIONS ################### #


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


# ################################### BACTERIAL ABUNDANCE FUNCTIONS ####################### #


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


# ################################### SAMPLE ABUNDANCE FUNCTIONS ########################## #


def sample_abundance_table_plotly(shared_file, design_file, simple_abundance_file, data_path, outputdir):
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


def sample_abundance_javascript_plotly(abundance_design_dict, abundance_otu_dict, name, design_file, js_path):
	sample_abundance_javascript_string = """
	//####################################################################
	//####################################################################
	//####################################################################
	//##########    SAMPLE Abundance table USING PLOTLY      ###################
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

		sample_abundance_javascript_string += abundance_design_dict[each_design]
		sample_abundance_javascript_string += "\tvar data_sample_abundance_" + each_design + " = [" + abundance_otu_dict[each_design] + "];\n"
		sample_abundance_javascript_string += """
		var layout_sample_abundance_""" + each_design + """ = {
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
		sample_abundance_javascript_string += """\n
		//####################################################################
		//####################################################################
		//####################################################################
		"""
	
	f = open(js_path + name + '_sample_abundance_plotly.js', 'w')
	f.write(sample_abundance_javascript_string)
	f.close()

# ################################### RAREFACTION FUNCTIONS ############################### #


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


# ################################### MOTHUR_FUNCTIONS #################################### #

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


# ################################### SPECIFIC_FUNCTIONS ################################## #


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


def random_color():
	golden_ratio_conjugate = 0.618033988749895
	randvalue = random.randint(1, 256)
	randvalue += golden_ratio_conjugate
	randvalue %= 1
	#hsv_to_rgb(randvalue, 0.5, 0.95)
	color = "#%06x" % random.randint(1, 0xFFFFFE)
	#print color
	return color


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


def construct_object(entity, entity_list):
	dict_of_entities = {}
	for each_entity in entity_list:
		#print each_entity
		a_class = entity(each_entity)
		dynamic_variable = ''
		dynamic_variable = a_class.get_taxon() + '_' + a_class.get_taxlevel()
		dict_of_entities[dynamic_variable] = a_class
	return dict_of_entities


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


def design_to_control(design_file, control_file):
	string = ''
	f = open(design_file, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		string += line[0] + '\n'
	write_string_down(string, control_file)
	return True

# ################################### EXECUTION_FUNCTIONS #################################### #


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


# ################################### UTILITIES_FUNCTIONS #################################### #

def make_archive(target_file_path):
	shutil.make_archive(target_file_path, 'zip', target_file_path)
	return True


def round_float(value):
	if isfloat(value) is False:
	
		return value
	else:
		return str(round(decimal.Decimal(value), 3))


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


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(query_list)


def transpose(one_list):
	trans_data = list(itertools.izip_longest(*one_list, fillvalue=''))
	return trans_data


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


def copy_file(source_file, destination):
	shutil.copy2(source_file, destination)
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


def write_string_down(new_string, new_name):
	f = open(new_name, 'w')
	f.write(new_string)
	f.close()
	return new_name


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def create_folder(path):
	# create folder in specified path
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False


def isPathExist(path):
	# check if the path exist and have access
	if os.path.exists(path) and os.access(path, os.R_OK):
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


def report(report_string):
	f = open(report_file, "a")
	f.write(report_string.encode('utf8'))
	f.write("\n")
	f.close()


# ################################### HTML & JAVASCRIPT ######################## #

def html_plotter(alpha_path, absname, design_file, alpha_diversity_summary_table_header, alpha_diversity_summary_table_body, pca_html_string_dict):
	# #################################################
	design_dictionary_rev = {}
	design_dictionary_rev = design_dict_maker(design_file, 'reverse')
	sorted_design_dictionary_rev = sorted(design_dictionary_rev, key=design_dictionary_rev.get, reverse=False)
	design_dictionary_rev = {}
	design_dictionary_rev = sorted_design_dictionary_rev
	# #################################################
	pca_html_string = ''
	pca_script_string = ''
	pca_plotly_script_string = ''
	for pca_key in pca_html_string_dict.keys():
		pca_script_string += pca_html_string_dict[pca_key][1] + '\n'
		pca_html_string += pca_html_string_dict[pca_key][0] + '\n'
		pca_plotly_script_string += pca_html_string_dict[pca_key][2] + '\n'
	# #################################################
	alpha_diversity_summary_html = ''
	alpha_diversity_summary_html += """
	<div id="alpha_diversity" class="container-fluid">
		<h2>Alpha diversity summary</h2>
		<div class="row">
		<dl class="dl-horizontal">
			<dt>Description</dt>
				<dd>The analysis of diversity index of each sample individually in an area or landscape of interest.</dd>
			<dt>Diversity index</dt>
			<dd>A diversity index is a quantitative measure that reflects how many different types (such as species) there are in a dataset, and simultaneously takes into account how evenly the basic entities (such as individuals) are distributed among those types. The value of a diversity index increases both when the number of types increases and when evenness increases. For a given number of types, the value of a diversity index is maximized when all types are equally abundant.</dd>
			
		</dl>
		</div>
		<div class="row">
			<div class="col-md-12">
				<div class="panel panel-default">
					<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
					<h3 class="panel-title">ALPHA Diversity Summary</h3>
					</div>
					<div class="panel-body">
						<div class="table-responsive">
											<table id="alpha_diversity_summary" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%">
											""" + alpha_diversity_summary_table_header + """
											""" + alpha_diversity_summary_table_body + """
											</table>
						</div>
					</div>
				</div>
			</div>
		</div>
	</div>
	"""
	# #################################################
	rarefaction_html = ''
	rarefaction_html = """
	<div id="rarefaction" class="container-fluid">
			<h2>Rarefaction curve plot</h2>
			<div class="row">
			<dl class="dl-horizontal">
				<dt>Description</dt>
					<dd>
		Rarefaction is a technique used in numerical ecology. Most commonly, rarefaction is used to determine whether all the species in an ecosystem have been observed, as discussed in this Wikipedia article.
		More generally, the goal of rarefaction is determine whether sufficient observations have been made to get a reasonable estimate of a quantity (call it R) that has been measured by sampling. The most commonly considered quantities are species richness (the number of different species in an environment or ecosystem) and alpha diversity (a measure of species diversity that may attempt to extrapolate to larger numbers of samples and/or take into account the abundance distribution -- there are many different definitions).

					</dd>
				
			</dl>
		</div>
		<div class="row">
	"""
	rarefaction_script = ''
	# #################################################
	sample_abundance_html = ''
	sample_abundance_html = """
	<div id="sample_abundance" class="container-fluid">
			<h2>Sample Abundance BarPlot</h2>
			<div class="row">
			<dl class="dl-horizontal">
				<dt>Description</dt>
					<dd>
		A bar chart based on the true abundance of each Sample's entities
					</dd>
				
			</dl>
		</div>
		<div class="row">

	"""
	sample_abundance_script = ''
	# #################################################
	bacterial_abundance_html = ''
	bacterial_abundance_html += """
	<div id="bacterial_abundance" class="container-fluid">
			<h2>Bacterial Abundance BarPlot</h2>
			<div class="row">
			<dl class="dl-horizontal">
				<dt>Description</dt>
					<dd>
		A bar chart based on the true abundance of each Sample's entities
					</dd>
				
			</dl>
		</div>
		<div class="row">

	"""
	bacterial_abundance_script = ''
	# #################################################
	check_box_string = ''
	for each_design in design_dictionary_rev:
		# ###########################################
		check_box_string += """
				<div class="checkbox">\n
					<label><input type="checkbox" name="chb[]" value='""" + each_design + """' >""" + each_design + """</label>\n
				</div>\n
		"""
		# ###########################################
		rarefaction_html += """
		<div class="col-md-6">
			<div class="panel panel-default">
				<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
					<h3 class="panel-title">""" + each_design + """</h3>
				</div>
				<div class="panel-body">
					<div id='div_rarefaction_""" + each_design + """'></div>
				</div>
			</div>
		</div>
		"""
		rarefaction_script += """\tPlotly.newPlot('div_rarefaction_""" + each_design + """', data_rarefaction_""" + each_design + """, layout_rarefaction_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n
		"""
		# ###########################################
		sample_abundance_html += """
		<div class="col-md-6">
			<div class="panel panel-default">
				<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
					<h3 class="panel-title">""" + each_design + """</h3>
				</div>
				<div class="panel-body">
					<div id='div_sample_abundance_""" + each_design + """'></div>
				</div>
			</div>
		</div>
		"""
		sample_abundance_script += """\tPlotly.newPlot('div_sample_abundance_""" + each_design + """', data_sample_abundance_""" + each_design + """, layout_sample_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n
		"""
		# ###########################################
		bacterial_abundance_html += """
		<div class="col-md-6">
			<div class="panel panel-default">
				<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
					<h3 class="panel-title">""" + each_design + """</h3>
				</div>
				<div class="panel-body">
					<div id='div_bacterial_abundance_""" + each_design + """'></div>
				</div>
			</div>
		</div>
		"""
		bacterial_abundance_script += """\tPlotly.newPlot('div_bacterial_abundance_""" + each_design + """', data_bacterial_abundance_""" + each_design + """, layout_bacterial_abundance_""" + each_design + """, {
		displaylogo: false,
		displayModeBar: true,
		modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],
		});\n
		"""
		# ###########################################
	rarefaction_html += """</div><!--End of Rarefaction class=row div-->"""
	rarefaction_html += """</div><!--End of Rarefaction class=container-fluid div-->"""
	sample_abundance_html += """</div><!--End of sample_abundance class=row div-->"""
	sample_abundance_html += """</div><!--End of sample_abundance class=container-fluid div-->"""
	bacterial_abundance_html += """</div><!--End of bacterial_abundance class=row div-->"""
	bacterial_abundance_html += """</div><!--End of bacterial_abundance class=container-fluid div-->"""
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
		<script type="text/javascript" src="./""" + absname + """_js/""" + absname + """_sample_abundance_plotly.js"></script>
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
		body {
		position: relative;
		}
		.table-hover tbody tr:hover td, .table-hover tbody tr:hover th {
						background-color: yellow;
					}
		a:link {color: #669;}		/* unvisited link */
		a:visited {color: #669;}	/* visited link */
		a:hover {color: #669;}		/* mouse over link */
		a:active {color: #669;}		/* selected link */
		#design {padding-top:50px;}
		#alpha_diversity {padding-top:50px;}
		#rarefaction {padding-top:50px;}
		#sample_abundance {padding-top:50px;}
		#bacterial_abundance {padding-top:50px;}
		#pcoa {padding-top:50px;}
		#biomarker_discovery {padding-top:50px;}
	</style>
	<!--############################-->
</head>
<body>
	<nav class="navbar navbar-inverse navbar-fixed-top" style="border:0px; background-color: #E6E6FF;">
		<div class="container-fluid">
			<div class="navbar-header">
					<button type="button" class="navbar-toggle" data-toggle="collapse" data-target="#myNavbar">
					<span class="icon-bar"></span>
					<span class="icon-bar"></span>
					<span class="icon-bar"></span>
					</button>
			</div>
			<div class="collapse navbar-collapse" id="myNavbar">
				<ul class="nav navbar-nav">
					<li><a href="#design"><span class="fa fa-cogs fa-2x"></span> Design</a></li>
					<li><a href="#alpha_diversity"><span class="fa fa-table fa-2x"></span> ALPHA Diversity summary</a></li>
					<li><a href="#rarefaction"><span class="fa fa-area-chart fa-2x"></span> Rarefaction</a></li>
					<li><a href="#sample_abundance"><span class="fa fa-bar-chart fa-2x"></span> Sample abundance</a></li>
					<li><a href="#bacterial_abundance"><span class="fa fa-bar-chart fa-2x"></span> Bacterial abundance</a></li>
					<li><a href="#pcoa"><span class="fa fa-cube fa-2x"></span> PCoA</a></li>
					<li><a  href="#biomarker_discovery"><span class="fa fa-th-list fa-2x"></span> Biomarker discovery</a></li>
				</ul>
			</div>
		</div>
	</nav>
	<div id="design" class="container-fluid">
	
		<h2>16S data analysis results</h2>
		<div class="row">
		<dl class="dl-horizontal">
			<dt>Request Name</dt>
				<dd>kraemerls-20160216-H1</dd>
			<dt>Project Title</dt>
				<dd>Microbiome composition of NOD2-/-TLR2-/- mice</dd>
			<dt>Project Description</dt>
				<dd>
				<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
				<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naiive mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd>
			<dt>Experiment</dt>
				<dd>Gut Microbiome 16S rRNA sequencing</dd>
			<dt>PI</dt>
				<dd>Helene Rosenberg</dd>
		</dl>
		</div>
		<div class="col-md-4 col-md-offset-4">
				<form role="form" method="post">
					<fieldset class="form-group">
						""" + check_box_string + """
						<button type="button" class="btn btn-sm btn-info" id="submit">FILTER</button>
					</fieldset>
				</form>
		</div>
	</div>
	""" + rarefaction_html + """
	""" + sample_abundance_html + """
	""" + bacterial_abundance_html + """
	""" + alpha_diversity_summary_html + """
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
	""" + sample_abundance_script + """
	""" + bacterial_abundance_script + """
	
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
					<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naiive mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
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


# ################################### END SCRIPT ############################### #
if __name__ == "__main__": main(sys.argv[1:])