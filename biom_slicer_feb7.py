# ################################### INFO ##################################### #
# 16S biom Slicer 1.0
# Author: Amir Shams
# Email: amir.shams84@gmail.com
# ################################### IMPORTED LIBRARY ######################### #


import sys
import os
import errno
import datetime
import signal
import logging as log
import time
import subprocess
import traceback
import itertools
import argparse
import multiprocessing
import platform
import pandas
import numpy
import csv
import shutil
import decimal
import zipfile
import random
import collections
import operator
import plotly
import plotly.graph_objs as PLOTLY_GO
import biom
# ################################### GLOBALS ################################## #


CHECK_MARK = "OK"
FAILED_MARK = ":("
DEFAULT_OUTPUTDIR = "/BIOM_SLICER_OUTPUTDIR/"
DEFAULT_TESTDIR = "/BIOM_SLICER_TESTDIR/"
DEFAULT_EXECDIR = "/BIOM_SLICER_EXECDIR/"
DEFAULT_PROCESSORS = str(multiprocessing.cpu_count())
DEFAULT_PREFIX = "BIOM_SLICER"
NO_BETA = True
CURRENT_PATH = "./"
TAXONOMY_EXIST = False
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
# ################################### EXECUTIONS ############################### #


def execute_functions(function_name, processors, outputdir, thread_type, func_mode, *args):
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
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
	else:
		print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	return (True, stderr)


def process_monitor(pid_list, stderr_list, stdout_list, outputdir, threads, mode):

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
# ################################### LOGGING & REPORTING ##################### #


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
# ################################### VISUALIZER ############################# #


def html_visualizer(final_string, html_string=None, javascript_string=None):

	if html_string is None and javascript_string is None:
		if final_string == '':

			print "Initial html script appending"
			final_string += """
				<!DOCTYPE html>
							<html lang='en' xml:lang='en' xmlns='http://www.w3.org/1999/xhtml'>
								<head>
									<meta charset="utf-8">
									<meta name="viewport" content="width=device-width, initial-scale=1">
									<!--JQUERY LIBRARAY MAIN LIBRARY-->
									<script type="text/javascript" src="https://code.jquery.com/jquery-1.12.3.js"></script>
									<!--PLOTLY LIBRARAY MAIN LIBRARY-->
									<script type="text/javascript" async src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_SVG"></script>
									<script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
									<!--BOOTSTRAP LIBRARAY MAIN LIBRARY-->
									<script type="text/javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
									<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
									<!--BOOTSTRAP DATATABLES LIBRARAY MAIN LIBRARY-->
									<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css"/>
									<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
									<script type="text/javascript" src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>
									<!--FONTAWESOME LIBRARAY MAIN LIBRARY-->
									<link rel="stylesheet" type="text/css" href="ttps://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css"/>
									<!--############################-->

									<!--OPEN SCRIPTS-->
										<script type="text/javascript">
											$(document).ready(function() {
												$('#statistics_table').DataTable( {
													"order": [[ 0, "asc" ]],
													"paging": false,
													"info": false,
													"ordering": true,
													"searching": false,
													"scrollY": "600px",
													"lengthMenu": [[25, 50, -1], [25, 50, "All"]]
												});
												$('#biomarker_table').DataTable( {
													"order": [[ 1, "asc" ]],
													"paging": false,
													"info": false,
													"ordering": true,
													"searching": false,
													"scrollY": "600px",
													"lengthMenu": [[25, 50, -1], [25, 50, "All"]]
												});

											});
									//###################################################
									//#######RELOAD ON RESIZE
									window.addEventListener('resize', function () {
									"use strict";
									window.location.reload();
									});

									//###################################################
											function linechart_filter(target_div){
											var myPlot = document.getElementById(target_div);
											var new_data = [];
											myPlot.on('plotly_click', function(targetdata){
											Plotly.deleteTraces(myPlot, targetdata.points[0].curveNumber)
												});

											};
									//###################################################
											function barchart_filter(target_div){
											var myPlot = document.getElementById(target_div);
											var new_data = [];
											myPlot.on('plotly_click', function(targetdata){

											var index_value = '';
											index_value = targetdata.points[0].x;
											var index_place = 0;
											index_place  = targetdata.points[0].fullData.x.indexOf(index_value);
											//console.log(index_value);
											//console.log(index_place);
											//console.log(targetdata);
											main_data_size = myPlot.data.length;
											new_data = myPlot.data
											for (var main_data_index = 0; main_data_index < main_data_size; main_data_index++){
												new_data[main_data_index].x.splice(index_place, 1);
												new_data[main_data_index].y.splice(index_place, 1);

											}
											//console.log(data_sample_abundance_Cancer[0].x)
											//Plotly.newPlot('myDiv', data, layout);
											//Plotly.deleteTraces(myPlot, data.points[0].x)
											//Plotly.newPlot('myDiv', data, layout);
											Plotly.redraw(target_div, new_data)
											//Plotly.relayout('myDiv',data)
										});

										};
										//##############################################
										function add_annotate(target_div){
											var myPlot = document.getElementById(target_div);
											var new_data = [];
											myPlot.on('plotly_click', function(target_data){
												var index_value = '';
												index_value = target_data.points[0].x;
												var index_place = 0;
												index_place  = target_data.points[0].fullData.x.indexOf(index_value);
											//console.log(index_value);
											//console.log(index_place);
											//console.log(targetdata);
												main_data_size = myPlot.data.length;
												new_data = myPlot.data;
												var annotation_string = 'Max:';
												var abundance_value_dict = {};
												for (var main_data_index = 0; main_data_index < main_data_size; main_data_index++){
													abundance_value_dict[myPlot.data[main_data_index].name] = myPlot.data[main_data_index].y[index_place];
												}
												sorted_abundance_value_array = sort_dictionary_by_value_descendingly(abundance_value_dict);
												//console.log(sorted_abundance_value_array);
												yPosition = 0;
												for (var i = 0; i < 1; i++){

													annotation_string += sorted_abundance_value_array[i] + '(' + abundance_value_dict[sorted_abundance_value_array[i]] + ')';

												}
												console.log(target_data.points[0].fullData.y);
												yPosition = target_data.points[0].fullData.y.reduce((a, b) => a + b, 0);
												console.log(yPosition);
												annotation = {
													text: annotation_string,
													x: index_value,
													y: yPosition,
													showarrow: true,
													font: {
														family: 'Avenir',
														size: 10,
														//color: '#ffffff'
													},
													align: 'center',
													arrowhead: 2,
													arrowsize: 1,
													arrowwidth: 2,
													arrowcolor: '#636363',
													ax: 20,
													ay: -30,
													bordercolor: '#c7c7c7',
													borderwidth: 2,
													borderpad: 4,
													//bgcolor: '#ff7f0e',
													opacity: 0.8
												}
												annotations = myPlot.layout.annotations || [];
												annotations.push(annotation);
												Plotly.relayout(target_div,{annotations: annotations});

											});
										};
										//##############################################
										function add_comment(target_div, target_comment){
											var myPlot = document.getElementById(target_div);
											var myComment = document.getElementById(target_comment).value;
											console.log(myComment);
											var new_data = [];
											myPlot.on('plotly_click', function(target_data){
												var index_value = '';
												index_value = target_data.points[0].x;
												var index_place = 0;
												index_place  = target_data.points[0].fullData.x.indexOf(index_value);
												annotation = {
													text: myComment,
													x: index_value,
													//y: yPosition,
													showarrow: true,
													font: {
														family: 'Avenir',
														size: 10,
														//color: '#ffffff'
													},
													align: 'center',
													arrowhead: 2,
													arrowsize: 1,
													arrowwidth: 2,
													arrowcolor: '#636363',
													ax: 20,
													ay: -30,
													bordercolor: '#c7c7c7',
													borderwidth: 2,
													borderpad: 4,
													//bgcolor: '#ff7f0e',
													opacity: 0.8
												}
												annotations = myPlot.layout.annotations || [];
												annotations.push(annotation);
												Plotly.relayout(target_div,{annotations: annotations});

											});
										};
										//##############################################

										function sort_dictionary_by_value_descendingly(myDict){
											var keys = [];
											for(var key in myDict){
												keys.push(key);
											}
											return keys.sort(function(a,b){return myDict[b] - myDict[a]});

										};
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
								<div class="container-fluid">
								<h2>biom Slicer results</h2>
								<div class="row">
									<dl class="dl-horizontal">
									<dt>Request Name</dt>
									<dd>20160216-H1</dd>
									<dt>Project Title</dt>
									<dd>biom composition</dd>
									<dt>Project Description</dt>
									<dd>
									<strong> Brief background: </strong><small>Having bred a strain of double gene-deleted mice that lack the pattern recognition receptors nucleotide-binding and oligomerization domain-containing 2 (NOD2) and toll-like receptor 2 (TLR2), these mice will likely be valuable in experiments demonstrating altered host response to probiotic bacteria and respiratory viruses in the airways.  Dysregulated capacity to sense microbial patterns by these two receptors may contribute to downstream immunomodulation, altered inflammation in the airways.</small><br>
									<strong>Goal/purpose or rationale: </strong><small>We would like to know if the lack of these two receptors leave mice with an altered composition of commensal microbes at baseline (untouched/naive mice), both in the gut and in the lung. Sequence comparison of microbial 16S rRNA genes within these tissues (knockout vs. WT) offers a powerful analytical tool.</small></dd-->
									<dt>Experiment</dt>
									<dd>Gut biom 16S rRNA sequencing</dd>
									<dt>PI</dt>
									<dd>...</dd>
									</dl>
								</div>
				"""
		else:
			final_string += '</body>\n</html>\n'
	elif html_string is not None and javascript_string is None:
		final_string += html_string
	elif html_string is not None and javascript_string is not None:
		final_string += html_string
		final_string += javascript_string
	return final_string


def plotly_html_maker(main_div_name, plot_type, figure_object_dict, mode_bar=None):
	temp_plotly_string = ''
	for each_figure in figure_object_dict:
		print each_figure
		plotly_div = plotly.offline.plot(figure_object_dict[each_figure], include_plotlyjs=False, output_type='div', show_link=False)
		temp_string = plotly_div.split('Plotly.newPlot(', 1)[1]
		temp2_string = temp_string.split(',', 1)[1]
		plotly_javascript_string = temp2_string.split('})', 1)[0]
		if mode_bar is None:
			plotly_javascript_string += """, displaylogo: false, displayModeBar: true, modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],}"""
		elif mode_bar is False:
			plotly_javascript_string += """, displaylogo: false, displayModeBar: false, modeBarButtonsToRemove: ['zoom2d', 'select2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'sendDataToCloud', 'orbitRotation', 'tableRotation', 'pan3d', 'zoom3d','resetCameraDefault3d','resetCameraLastSave3d','hoverClosest3d'],}"""
		temp_plotly_string += """Plotly.newPlot("PLOTLY_""" + main_div_name + '_' + plot_type + '_' + each_figure + """_DIV", """ + plotly_javascript_string + """);"""

	panel_count = len(figure_object_dict.keys())
	if panel_count == 1:
		bootstrap_column_string = 'col-xs-12 col-sm-12 col-md-12 col-lg-12'
	if panel_count == 2:
		bootstrap_column_string = 'col-xs-6 col-sm-6 col-md-6 col-lg-6'
	elif panel_count == 3:
		bootstrap_column_string = 'col-xs-4 col-sm-4 col-md-4 col-lg-4'
	elif panel_count == 4:
		bootstrap_column_string = 'col-xs-3 col-sm-3 col-md-3 col-lg-3'
	elif panel_count > 4:
		bootstrap_column_string = 'col-xs-2 col-sm-2 col-md-2 col-lg-2'

	plotly_div_string = """
		<div class="row">
	"""

	for each_figure in figure_object_dict:
		plotly_div_string += """
					<div id="PLOTLY_""" + main_div_name + '_' + plot_type + '_' + each_figure + """_DIV" class=" """ + bootstrap_column_string + """ "></div>

		"""

	plotly_div_string += """
		</div>
	"""

	plotly_html_string = """
	<h4>""" + main_div_name + ' ' + plot_type + """</h4>
	<div class="row">
		<dl class="dl-horizontal">
		<dt>Brief Description</dt>
		<dd>
		<strong> Brief background: </strong><small>Some description over here is useful</small><br>
		<strong>Goal/purpose or rationale: </strong><small>Alittle detail of the purpose over here</small></dd-->
		</dl>
	</div>
	<div id="PLOTLY_""" + main_div_name + '_' + plot_type + """_MAIN_DIV" class="container-fluid">
		<div class="row">

			<div class="col-xs-12 col-sm-12 col-md-12 col-lg-12">
				<div class="panel panel-default" >
					<div id="PLOTLY_""" + main_div_name + '_' + plot_type + """_HEADER_DIV" class="panel-heading" style="border:0px; background-color: #F4F9FD;">
						<h3 class="panel-title">""" + main_div_name + """</h3>
					</div>
					<div id="PLOTLY_""" + main_div_name + '_' + plot_type + """_BODY_DIV" class="panel-body" align="center">

							""" + plotly_div_string + """
					</div>
					<div id="PLOTLY_""" + main_div_name + '_' + plot_type + """_FOOTER_DIV" class="panel-footer">
						<form class="form-inline" role="form">
							<label class="form-check-inline" style="font-family:'Avenir';margin-left: 5px;">
								<input type="radio" name="options" id="ALPHA_DIVERSITY_INDEX_DIV_normal_mode" autocomplete="off" checked > Normal mode
							</label>
							<!--label class="form-check-inline" style="font-family:'Avenir';margin-left: 5px;">
								<input type="radio" name="options" id="ALPHA_DIVERSITY_INDEX_DIV_delete_sample_mode" autocomplete="off" onclick="barchart_filter('ALPHA_DIVERSITY_INDEX_DIV');"> Click to Delete mode
							</label-->
						</form>
					</div>
				</div>
			</div>
		</div>
	</div>
	"""
	plotly_script = """
				<script type="text/javascript">
					//####################################################################
					//####################################################################
					//####################################################################
					//#####   PLOTLY_""" + main_div_name + '_' + plot_type + """  ########
					//####################################################################
					//####################################################################
					//####################################################################
					//####################################################################
	"""
	plotly_script += temp_plotly_string

	plotly_script += """
					//####################################################################
					//####################################################################
					//####################################################################
				</script>
	"""
	return (plotly_html_string, plotly_script)


def plotly_Natural_Abundance_Barplot_NORMAL_old(shared_file, design_file):
	headers_list = []
	columns_dict = {}

	headers_list, columns_dict = mothur_shared_parser(shared_file)
	sample_name_list = columns_dict['Groups']

	OTU_name_list = []
	for each_head in headers_list:
		if each_head.lower() in ['label', 'groups', 'numotus']:
			continue
		else:
			OTU_name_list.append(each_head)

	OTU_color_dict = color_lover_generate(OTU_name_list)

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	# ######################## PLOTLY DATA CREATE
	Natural_Abundance_barplot_TRACE_list = []
	Natural_Abundance_barplot_FIGURE_objects_dict = {}
	Natural_Abundance_barplot_LAYOUT = {}

	for each_design in reverse_design_dict:
		sample_list = reverse_design_dict[each_design]
		#design_text_name = []
		for each_OTU in OTU_name_list:
			OTU_value_list = []
			for each_sample in sample_list:
				each_sample_index = columns_dict['Groups'].index(each_sample)
				OTU_value_list.append(columns_dict[each_OTU][each_sample_index])
			#design_text_name.append(each_design)
			Natural_Abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
				x=sample_list,
				y=OTU_value_list,
				name=each_OTU,
				text=[each_design] * len(sample_list),
				hoverinfo='name+y+text',
				marker=dict(
					color=OTU_color_dict[each_OTU],
				)
			)
			)

	Natural_Abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
		barmode='stack',
		height=600,
		autosize=True,
		showlegend=True,
		legend=dict(
			font=dict(family='Avenir', size=10),
			orientation="v",
			traceorder="normal"
		),

		title='Bacterial Taxonomic Composition',
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
		),
		font=dict(family='Avenir', size=10),
		yaxis=dict(
			showgrid=True,
			title='Natural Abundance Value',
			titlefont=dict(
				family='Avenir',
				size=15,
			)
		),
		margin=dict(
			l=100,
			r=100,
			b=100,
			t=100,
		)
	)

	#print Natural_Abundance_barplot_TRACE_objects_dict[each_design].keys()
	Natural_Abundance_barplot_FIGURE_objects_dict['natural_abundance'] = PLOTLY_GO.Figure(data=Natural_Abundance_barplot_TRACE_list, layout=Natural_Abundance_barplot_LAYOUT)

	# ###################################  PLOTLY PLOT CREATE ###############################

	plotly_html_string, plotly_script = plotly_html_maker('NATURAL_ABUNDANCE', 'STACKED_BARPLOT', Natural_Abundance_barplot_FIGURE_objects_dict, mode_bar=False)

	return (plotly_html_string, plotly_script)


def plotly_Rarefaction_Curve_Linechart_DROPDOWN_old(shared_file, design_file, mothur_exec_path, processors, outputdir):
	# ######################## Rarefaction curve CALCULATION
	frequency_value = '1000'
	flag, stderr = execute_functions(mothur_rarefaction_single, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, frequency_value)
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
	rarefaction_file = scanned_container[0]
	# ######################## PARSING
	headers_list = []
	columns_dict = {}
	raw_sample_name_list = []
	headers_list, columns_dict = mothur_shared_parser(rarefaction_file)
	number_of_sequences_sampled_list = columns_dict['numsampled']
	for each_head in headers_list:
		if 'lci' in each_head:
			continue
		elif 'hci' in each_head:
			continue
		elif 'numsampled' in each_head:
			continue
		else:
			raw_sample_name_list.append(each_head)
	rarefaction_sample_dict_values = {}
	for each_sample in raw_sample_name_list:
		clean_sample_name = each_sample.split('-')[1]
		rarefaction_sample_dict_values[clean_sample_name] = list_variance(columns_dict[each_sample], 0.01)  # removing NA element and variant factors

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	visibility_list_length = len(rarefaction_sample_dict_values.keys())
	visibility_flag_dict = {}
	visibility_flag_dict['ALL'] = [True] * visibility_list_length
	VISIBLE_flag_list = []
	VISIBLE_flag_list.append(
		dict(
			args=[
				dict(
					#title=title_dictionary[each_index],
					visible=visibility_flag_dict['ALL']
				)
			],
			label='ALL',
			method='restyle',

		),
	)
	for each_design in reverse_design_dict:
		visibility_flag_dict[each_design] = [False] * visibility_list_length
		sample_list = reverse_design_dict[each_design]
		for each_sample in rarefaction_sample_dict_values.keys():
			if each_sample in sample_list:
				each_sample_index = rarefaction_sample_dict_values.keys().index(each_sample)
				visibility_flag_dict[each_design][each_sample_index] = True
		VISIBLE_flag_list.append(
			dict(
				args=[
					dict(
						#title=title_dictionary[each_index],
						visible=visibility_flag_dict[each_design]
					)
				],
				label=each_design,
				method='restyle',

			),
		)

	# ######################## PLOTLY DATA CREATE
	RAREFACTION_curve_linechart_TRACE_list = []
	RAREFACTION_curve_linechart_LAYOUT_objects_dict = {}
	RAREFACTION_curve_linechart_FIGURE_objects_dict = {}
	for each_sample in rarefaction_sample_dict_values.keys():
		RAREFACTION_curve_linechart_TRACE_list.append(PLOTLY_GO.Scatter(
				x=number_of_sequences_sampled_list,
				y=rarefaction_sample_dict_values[each_sample],
				name=each_sample,
				hoverinfo='name',
				mode='lines',
			)
		)

	RAREFACTION_curve_linechart_LAYOUT_objects_dict = PLOTLY_GO.Layout(
		height=600,
		autosize=True,
		title='Rarefaction curve plot',
		showlegend=True,
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=12
		),
		xaxis=dict(
			showgrid=False,
			title='Number of sampling effort',
			titlefont=dict(
				family='Avenir',
				size=10,
			)

		),
		font=dict(family='Avenir', size=10),
		yaxis=dict(
			showgrid=True,
			title='Number of detected species',
			titlefont=dict(
				family='Avenir',
				size=10,
			)
		),
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0
			)
		]
		)

	)
	RAREFACTION_curve_linechart_FIGURE_objects_dict['rarefaction'] = PLOTLY_GO.Figure(data=RAREFACTION_curve_linechart_TRACE_list, layout=RAREFACTION_curve_linechart_LAYOUT_objects_dict)
	# ###################################  PLOTLY PLOT CREATE
	plotly_html_string, plotly_script = plotly_html_maker('RAREFACTION_CURVE', 'LINE_CHART', RAREFACTION_curve_linechart_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Alpha_Diversity_Index_Boxplot_DROPDOWN_old(shared_file, design_file, mothur_exec_path, processors, outputdir):
	# ######################## ALPHA DIVERSITY INDEX CALCULATION
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
	alpha_diversity_index_file = scanned_container[0]
	# ######################## PARSING
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_result_parser(alpha_diversity_index_file)
	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')
	# ######################## PLOTLY DATA CREATE
	ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict = {}
	ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict = {}
	ALPHA_Diversity_Index_boxplot_TRACE_list = []
	alpha_index_list = ['shannon', 'simpson', 'invsimpson', 'sobs', 'ace']
	visibility_list_length = len(alpha_index_list) * len(reverse_design_dict.keys())
	alpha_visibility_list = []
	# ####################### TO GENERATE SWITCHING STRING
	VISIBLE_flag_list = []
	visibility_flag_dict = {}
	alpha_visibility_list = [False] * visibility_list_length
	step_list = range(len(reverse_design_dict.keys()))
	visible_index = 0
	for i in xrange(0, visibility_list_length, len(reverse_design_dict.keys())):
		for each_step in step_list:
			alpha_visibility_list[i + each_step] = True
		visibility_flag_dict[visible_index] = alpha_visibility_list
		alpha_visibility_list = [False] * visibility_list_length
		visible_index += 1
	#print visibility_flag_dict
	# #######################
	title_dictionary = {}
	title_dictionary['simpson'] = 'Simpson diversity index'
	title_dictionary['invsimpson'] = 'Inverse-Simpson diversity index'
	title_dictionary['shannon'] = 'Shannon diversity index'
	title_dictionary['sobs'] = 'Sobs(observed number of species) richness index'
	title_dictionary['ace'] = 'ACE(Abundance-base Coverage Estimator) richness index'
	for each_index in alpha_index_list:
		for each_design in reverse_design_dict:
			sample_list = reverse_design_dict[each_design]
			index_value_list = []
			for each_sample in sample_list:
				each_sample_index = columns_dict['group'].index(each_sample)
				index_value_list.append(columns_dict[each_index][each_sample_index])

			ALPHA_Diversity_Index_boxplot_TRACE_list.append(PLOTLY_GO.Box(
				x=each_design,
				y=index_value_list,
				#text=each_design,
				hoverinfo='y',
				name=each_design,
				boxmean='sd',
				visible=visibility_flag_dict[alpha_index_list.index(each_index)][reverse_design_dict.keys().index(each_design) * alpha_index_list.index(each_index)],
				legendgroup=title_dictionary[each_index]
			)
			)
		VISIBLE_flag_list.append(
			dict(
				args=[
					dict(
						#title=title_dictionary[each_index],
						visible=visibility_flag_dict[alpha_index_list.index(each_index)]
					)
				],
				label=title_dictionary[each_index],
				method='restyle',

			),
		)

	ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict = PLOTLY_GO.Layout(
		boxmode='group',
		height=600,
		autosize=True,
		showlegend=False,
		title='ALPHA Diversity index',
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
		),
		font=dict(family='Avenir', size=12),
		yaxis=dict(
			showgrid=True,
			#title=each_index + ' Index Value',
			titlefont=dict(
				family='Avenir',
				size=16,
			)
		),
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0
			)
		])
	)
	#print VISIBLE_flag_list
	ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict['alpha_diversity_index'] = PLOTLY_GO.Figure(data=ALPHA_Diversity_Index_boxplot_TRACE_list, layout=ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict)
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('ALPHA_DIVERSITY_INDEX', 'BOX_PLOT', ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def brief_statistics_table_function(shared_file_PATH, design_MATRIX_reverse):

	shared_DICT, shared_LIST = any_file_to_dict_converter(shared_file_PATH)
	#design_dictionary = design_dict_maker(design_file_PATH)
	biom_matrix_header_list = ['Groups', 'Sample Name', 'Total number', 'Classified number', 'Unclassified number', 'Classified percentage', 'Unclassified percentage']
	no_unclassified_flag = True
	if 'unclassified' in shared_LIST:
		no_unclassified_flag = False
	else:
		no_unclassified_flag = True
	thead_string = '								<thead>\n								<tr>\n'
	for thead in biom_matrix_header_list:
		thead_string += '									<th class="text-center">' + thead + '</th>\n'
	thead_string += '								</tr>\n								</thead>\n'
	tbody_string = '								<tbody>\n'
	bacterial_list = []

	for each_head in shared_LIST:
		if each_head.lower() in ['label', 'numotus', 'groups']:
			continue
		else:
			bacterial_list.append(each_head)
	biom_matrix_length = len(shared_DICT['numOtus'])

	for each_row in range(biom_matrix_length):
		each_row_list = []
		for each_bacteria in bacterial_list:
			each_row_list.append(shared_DICT[each_bacteria][each_row])
		each_row_list_int = list(map(int, each_row_list))
		total_read_count = int(sum(each_row_list_int))
		#we test for being outliers

		tbody_string += '									<tr>\n'
		#tbody_string += '									<td class="text-center">' + str(design_dictionary[columns_dict['Groups'][each_row]]) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(design_MATRIX_reverse[shared_DICT['Groups']][0]) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(shared_DICT['Groups'][each_row]) + '</td>\n'

		tbody_string += '									<td class="text-center">' + str(total_read_count) + '</td>\n'
		if no_unclassified_flag is False:
			unclassified_read_count = int(shared_DICT['unclassified'][each_row])
		else:
			unclassified_read_count = 0
		classified_read_count = int(total_read_count - unclassified_read_count)
		tbody_string += '									<td class="text-center">' + str(classified_read_count) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(unclassified_read_count) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(round_float(percentage(classified_read_count, total_read_count))) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(round_float(percentage(unclassified_read_count, total_read_count))) + '</td>\n'
		tbody_string += '								</tr>\n'
	tbody_string += '								</tbody>\n'

	brief_statistics_table_html_string = """
								<div id="BRIEF_STATISTICS" class="container-fluid">
									<div class="row">
										<div class="col-xs-12 col-sm-12 col-md-12 col-lg-12">
											<div class="panel panel-default">
												<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
													<h3 class="panel-title">SEQUENCING BRIEF STATISTICS</h3>
												</div>
												<div class="panel-body" align="center">
													<!-- Table -->
													<div class="table-responsive">
													<table id="brief_statistics_table" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%" style="font-family:'Avenir';">
													""" + thead_string + """
													""" + tbody_string + """
													</table>

												</div>
												</div>
												<div class="panel-footer">
												</div>
											</div>
										</div>
									</div>
								</div>
	"""

	return brief_statistics_table_html_string


def plotly_Bacterial_relative_abundance_Barplot_old(shared_file, design_file, OTU_metadata_file_PATH):


	condensed_shared_dict = {}
	condensed_shared_dict = condense_shared_file_by_design(shared_file, design_file)
	#print condensed_shared_dict
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)

	design_name_list = condensed_shared_dict.keys()
	OTU_name_list = condensed_shared_dict[design_name_list[0]].keys()
	#print sample_name_list
	#print OTU_name_list
	#print shared_dict[sample_name_list[0]]
	real_OTU_name_LIST = []
	for each_OTU in OTU_name_list:
		real_OTU_name_LIST.append(OTU_MATRIX_reverse[each_OTU]['OTU_ID'])


	Bacterial_relative_abundance_barplot_TRACE_list = []
	Bacterial_relative_abundance_barplot_LAYOUT = {}
	Bacterial_relative_abundance_barplot_FIGURE_objects_dict = {}

	for each_design in design_name_list:
		Total_value_of_each_OTU_list = []
		for each_OTU in OTU_name_list:
			Total_value_of_each_OTU_list.append(sum(condensed_shared_dict[each_design][each_OTU]))
		Bacterial_relative_abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
			x=real_OTU_name_LIST,
			y=Total_value_of_each_OTU_list,
			name=each_design,
			hoverinfo='all',
		)
		)
	Bacterial_relative_abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
		barmode='group',
		height=600,
		autosize=True,
		title='Bacterial Relative abundance',
		showlegend=True,
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
			#color='#7f7f7f'
		),
		font=dict(family='Avenir', size=10),
		xaxis=dict(
			autorange=True,
		),
		yaxis=dict(
			title='Relative abundance',
			titlefont=dict(family='Avenir', size=12),
		),
		margin=dict(
			l=100,
			r=100,
			b=100,
			t=100,
		)
	)

	Bacterial_relative_abundance_barplot_FIGURE_objects_dict['bacterial_relative_abundance'] = PLOTLY_GO.Figure(data=Bacterial_relative_abundance_barplot_TRACE_list, layout=Bacterial_relative_abundance_barplot_LAYOUT)

	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('BACTERIAL_RELATIVE_ABUNDANCE', 'STACKED_BARPLOT', Bacterial_relative_abundance_barplot_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Relative_Abundance_Piechart_old(shared_file, design_file):
	# ############################### CONDENSE SHARED FILE
	condensed_shared_dict = {}
	condensed_shared_dict = condense_shared_file_by_design(shared_file, design_file)
	design_name_list = condensed_shared_dict.keys()
	OTU_name_list = condensed_shared_dict[design_name_list[0]].keys()
	# ###################################  PLOTLY DATA CREATE
	Relative_Abundance_Piechart_TRACE_objects_dict = {}
	#Relative_Abundance_Piechart_LAYOUT_objects_dict = {}
	Relative_Abundance_Piechart_FIGURE_objects_dict = {}
	Relative_Abundance_Piechart_ANNOTATIONS_objects_dict = {}

	Domain_step = 1.0 / len(design_name_list)
	Domain_tolerance = 0.01
	#Domain_step = 0.2
	Domain_min = 0
	Domain_max = Domain_step
	Domain_list = [Domain_min, Domain_max]
	Domain_mid = (Domain_max + Domain_min) / 2

	for each_design in design_name_list:
		Relative_Abundance_Piechart_TRACE_objects_dict[each_design] = {}
		#Relative_Abundance_Piechart_LAYOUT_objects_dict[each_design] = {}
		Relative_Abundance_Piechart_ANNOTATIONS_objects_dict[each_design] = {}
		Total_value_of_each_OTU_list = []
		for each_OTU in OTU_name_list:
			Total_value_of_each_OTU_list.append(sum(condensed_shared_dict[each_design][each_OTU]))
		Relative_Abundance_Piechart_TRACE_objects_dict[each_design] = PLOTLY_GO.Pie(
			labels=OTU_name_list,
			values=Total_value_of_each_OTU_list,
			name=each_design,
			type='pie',
			hoverinfo="label+percent+name",
			textinfo="all",
			hole=0.4,
			textposition="inside",
			domain=dict(
				x=Domain_list,
				y=[0.0, 1.0]
			)
		)
		Relative_Abundance_Piechart_ANNOTATIONS_objects_dict[each_design] = dict(
			showarrow=False,
			text=each_design,
			xanchor="center",
			xref="paper",
			yanchor="bottom",
			yref="paper",
			y=0.5,
			x=Domain_mid
		)

		Domain_tolerance = 0.01
		Domain_step = 1.0 / len(design_name_list)
		Domain_min = Domain_max + Domain_tolerance
		Domain_max = Domain_max + Domain_step + Domain_tolerance
		Domain_list = [Domain_min, Domain_max]
		Domain_mid = (Domain_min + Domain_max) / 2

	Relative_Abundance_Piechart_LAYOUT = PLOTLY_GO.Layout(
		title='Bacterial Taxonomic Composition',
		autosize=True,
		height=700,
		titlefont=dict(
			family='Avenir',
			size=16
		),
		showlegend=True,
		legend=dict(
			font=dict(family='Avenir', size=12),
			orientation="v",
			traceorder="normal"
		),
		font=dict(family='Avenir', size=15),
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
		annotations=PLOTLY_GO.Annotations(Relative_Abundance_Piechart_ANNOTATIONS_objects_dict.values())
	)
	Relative_Abundance_Piechart_FIGURE_objects_dict['pie'] = PLOTLY_GO.Figure(data=Relative_Abundance_Piechart_TRACE_objects_dict.values(), layout=Relative_Abundance_Piechart_LAYOUT)
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('Relative_Abundance', 'PIECHART', Relative_Abundance_Piechart_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Principal_Coordinate_Analysis(shared_file, design_file, mothur_exec_path, prefix, processors, outputdir):
	# ##################### DISTANCE MATRIX CALCULATION
	distance_matrix_file = outputdir + prefix + '_Braycurtis_distance_matrix.txt'
	calc = 'braycurtis'
	flag = distance_matrix_analyser(shared_file, calc, distance_matrix_file, mothur_exec_path, prefix, processors, outputdir)
	if flag is False:
		print "Something is wrong with distance_matrix_analyser."
	# ##################### PCOA CALCULATION
	pcoa_axis_file = outputdir + prefix + '_PCOA_axis_file.txt'
	flag, stderr = execute_functions(mothur_PCOA, processors, outputdir, 'multi', 'mothur', mothur_exec_path, distance_matrix_file)
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
	os.rename(scanned_container[0], pcoa_axis_file)
	# ##################### PARSING AXIS FILE
	AXIS_dict = {}
	AXIS_dict = parse_pcoa_axis_file(pcoa_axis_file)

	sample_name_list = AXIS_dict.keys()
	AXIS_list = AXIS_dict[sample_name_list[0]].keys()
	#print AXIS_dict
	condensed_AXIS_dict = {}
	variance_AXIS_dict = {}
	for each_AXIS in AXIS_list:
		AXIS_total_count_list = []
		for each_sample in sample_name_list:
			AXIS_total_count_list.append(float(AXIS_dict[each_sample][each_AXIS]))
		variance_AXIS_dict[each_AXIS] = numpy_variance(AXIS_total_count_list)
		condensed_AXIS_dict[each_AXIS] = AXIS_total_count_list
	#print variance_AXIS_dict
	#sorted_variance_AXIS_dict = sorted(variance_AXIS_dict.items(), key=operator.itemgetter(1), reverse=True)
	#sorted_variance_keys = sorted(variance_AXIS_dict.items(), key=operator.itemgetter(0), reverse=True)
	#sorted_variance_values = sorted(variance_AXIS_dict.items(), key=operator.itemgetter(1), reverse=True)
	varinace_XAXIS_list = []
	varinace_YAXIS_list = []
	for i in range(1, 11):
		varinace_XAXIS_list.append('PC' + str(i))
		varinace_YAXIS_list.append(variance_AXIS_dict['axis' + str(i)])
	#print varinace_XAXIS_list
	#print varinace_YAXIS_list

	# ##################### VARIANCE PLOT
	PCOA_VARIANCE_barplot_TRACE = {}
	PCOA_VARIANCE_barplot_LAYOUT = {}
	PCOA_VARIANCE_barplot_FIGURE_objects_dict = {}
	PCOA_VARIANCE_barplot_TRACE = PLOTLY_GO.Bar(
		x=varinace_XAXIS_list,
		y=varinace_YAXIS_list,
		name='Varinace'
	)
	PCOA_VARIANCE_SCATTER_TRACE = PLOTLY_GO.Scatter(
		x=varinace_XAXIS_list,
		y=varinace_YAXIS_list,
		name='Scree Plot'
	)
	PCOA_VARIANCE_barplot_LAYOUT = PLOTLY_GO.Layout(
		title='Explained variance by different principal coordinates',
		height=600,
		autosize=True,
		showlegend=True,
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
	)
	PCOA_VARIANCE_barplot_FIGURE_objects_dict['a_scree'] = PLOTLY_GO.Figure(data=[PCOA_VARIANCE_barplot_TRACE, PCOA_VARIANCE_SCATTER_TRACE], layout=PCOA_VARIANCE_barplot_LAYOUT)
	# ##################### 3D Scatter plot

	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_result_parser(pcoa_axis_file)

	design_dict = {}
	design_dict = design_dict_maker(design_file)

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	#design_color_dict = color_lover_generate(list(set(design_dict.values())))
	sample_name_list = columns_dict['group']

	PCOA_3D_SCATTERPLOT_TRACE_objects_dict = {}
	for each_design in reverse_design_dict:
		PCOA_3D_SCATTERPLOT_TRACE_objects_dict[each_design] = {}
		sample_list = reverse_design_dict[each_design]
		AXIS1_list = []
		AXIS2_list = []
		AXIS3_list = []
		for each_sample in sample_list:

			each_sample_index = columns_dict['group'].index(each_sample)
			AXIS1_list.append(columns_dict['axis1'][each_sample_index])
			AXIS2_list.append(columns_dict['axis2'][each_sample_index])
			AXIS3_list.append(columns_dict['axis3'][each_sample_index])
		PCOA_3D_SCATTERPLOT_TRACE_objects_dict[each_design] = PLOTLY_GO.Scatter3d(
			x=AXIS1_list,
			y=AXIS2_list,
			z=AXIS3_list,
			type='scatter3d',
			mode='markers',
			name=each_design,
			marker=dict(
				symbol='circle',
				size=8,
				opacity=0.9
			),
			text=sample_list,
			hoverinfo='name+text'
		)

	PCOA_3D_scatterplot_LAYOUT = PLOTLY_GO.Layout(
		title='PCOA PLOT',
		autosize=True,
		showlegend=True,
		hovermode='closest',
		height=700,
		titlefont=dict(
			family='Avenir',
			size=16
		),
		dragmode=True,
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
		legend=dict(
			font=dict(family='Avenir', size=15),
			orientation="v",
			traceorder="normal"
		),
		scene=dict(
			xaxis=dict(
				title='PC1',
				titlefont=dict(
					family='Avenir',
					size=16
				),
				autorange=True,
				showbackground=True,
				nticks=10,
				#autotick=True,
				fixedrange=False,
				gridcolor="white",
				linecolor="white",
				zerolinecolor="white",
				gridwidth=4,
				backgroundcolor="rgba(0,0,0,0.02)",
				showspikes=False,
			),
			yaxis=dict(
				title='PC2',
				titlefont=dict(
					family='Avenir',
					size=16
				),
				autorange=True,
				showbackground=True,
				nticks=10,
				#autotick=True,
				fixedrange=False,
				gridcolor="white",
				linecolor="white",
				zerolinecolor="white",
				gridwidth=4,
				backgroundcolor="rgba(0,0,0,0.02)",
				showspikes=False,

			),
			zaxis=dict(
				title='PC3',
				titlefont=dict(
					family='Avenir',
					size=16
				),
				autorange=True,
				showbackground=True,
				nticks=10,
				#autotick=True,
				fixedrange=False,
				gridcolor="white",
				linecolor="white",
				zerolinecolor="white",
				gridwidth=4,
				backgroundcolor="rgba(0,0,0,0.02)",
				showspikes=False,

			),
			aspectratio=dict(
				x=1,
				y=1,
				z=1
			),
			aspectmode='manual'
		)
	)
	PCOA_VARIANCE_barplot_FIGURE_objects_dict['b_3d'] = PLOTLY_GO.Figure(data=PCOA_3D_SCATTERPLOT_TRACE_objects_dict.values(), layout=PCOA_3D_scatterplot_LAYOUT)
	# ##################### HIERARCHICAL CLUSTERING CALCULATION



	# ##################### CORR AXES CALCULATION
	corr_axis_file = outputdir + prefix + '_CORR_axis_file.txt'
	flag, stderr = execute_functions(mothur_corr_axes, processors, outputdir, 'multi', 'mothur', mothur_exec_path, pcoa_axis_file, shared_file)
	if flag is False:
		print "Execution of mothur_corr_axes failed!!!"
	else:
		scanned_container = []
		extension_list = ['.spearman.corr.axes']
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
	os.rename(scanned_container[0], corr_axis_file)
	# ##################### PARSE CORR AXES RESULT
	"""
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_result_parser(corr_axis_file)
	AXIS1_list = columns_dict['axis1']
	AXIS2_list = columns_dict['axis2']
	OTU_list = columns_dict['OTU']
	length_list = map(float, columns_dict['length'])
	length_list = [int(x*20) for x in length_list]
	print length_list
	plotly_CORR_AXES_SCATTER_TRACE = PLOTLY_GO.Scatter(
		x=AXIS1_list,
		y=AXIS2_list,
		text=OTU_list,
		hoverinfo='text',
		mode='markers',
		marker=dict(
			symbol='circle',
			size=length_list,
			#opacity=0.9
		),

	)
	plotly_CORR_AXES_SCATTER_LAYOUT = PLOTLY_GO.Layout(
		title='blah bklah',
		height=700

	)
	#PCOA_VARIANCE_barplot_FIGURE_objects_dict['c_CORR'] = PLOTLY_GO.Figure(data=[plotly_CORR_AXES_SCATTER_TRACE], layout=plotly_CORR_AXES_SCATTER_LAYOUT)
	"""
	# ###################################  PLOTLY PLOT CREATE
	plotly_html_string, plotly_script = plotly_html_maker('PCOA_VARIANCE', 'SCREE_BARPLOT', PCOA_VARIANCE_barplot_FIGURE_objects_dict, mode_bar=False)

	return (plotly_html_string, plotly_script)


def plotly_Dominant_Entities_Distribution_Boxplot(shared_file, design_file):
	# ######################## PARSING DATA
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_shared_parser(shared_file)

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

	OTU_name_list = []
	for each_head in headers_list:
		if each_head.lower() in ['label', 'groups', 'numotus']:
			continue
		else:
			OTU_name_list.append(each_head)
	# ######################## PLOTLY DATA CREATE
	Dominant_Entities_Distribution_boxplot_TRACE_objects_dict = {}
	Dominant_Entities_Distribution_boxplot_LAYOUT_objects_dict = {}
	Dominant_Entities_Distribution_boxplot_FIGURE_objects_dict = {}
	for each_OTU in OTU_name_list:
		Dominant_Entities_Distribution_boxplot_TRACE_objects_dict[each_OTU] = {}
		Dominant_Entities_Distribution_boxplot_LAYOUT_objects_dict[each_OTU] = {}
		for each_design in reverse_design_dict:
			sample_list = reverse_design_dict[each_design]
			OTU_value_list = []
			for each_sample in sample_list:
				each_sample_index = columns_dict['Groups'].index(each_sample)
				OTU_value_list.append(columns_dict[each_OTU][each_sample_index])
			Dominant_Entities_Distribution_boxplot_TRACE_objects_dict[each_OTU][each_design] = PLOTLY_GO.Box(
				x=each_design,
				y=OTU_value_list,
				name=each_design,
				hoverinfo='y+name',
				boxmean=True,
				#orientation='h',
			)
		Dominant_Entities_Distribution_boxplot_LAYOUT_objects_dict[each_OTU] = PLOTLY_GO.Layout(
			boxmode='group',
			height=400,
			autosize=True,
			showlegend=False,
			title=each_OTU,
			hovermode='closest',
			titlefont=dict(
				family='Avenir',
				size=12
			),
			font=dict(family='Avenir', size=10),
			yaxis=dict(
				showgrid=True,
				#title=each_OTU + ' Index Value',
				titlefont=dict(
					family='Avenir',
					size=10,
				)
			),
			margin=dict(
				l=40,
				r=30,
				b=80,
				t=100,
			),
		)
		Dominant_Entities_Distribution_boxplot_FIGURE_objects_dict[each_OTU] = PLOTLY_GO.Figure(data=Dominant_Entities_Distribution_boxplot_TRACE_objects_dict[each_OTU].values(), layout=Dominant_Entities_Distribution_boxplot_LAYOUT_objects_dict[each_OTU])
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('Dominant_Entities_Distribution', 'BOX_PLOT', Dominant_Entities_Distribution_boxplot_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def statistics_table_function(shared_file_PATH, design_MATRIX_forward, design_MATRIX_reverse):
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	design_header_LIST = []
	design_header_LIST = design_MATRIX_forward.keys()
	if len(design_header_LIST) > 2 and 'DEFAULT_DESIGN' in design_header_LIST:
		design_header_LIST.pop(design_header_LIST.index('DEFAULT_DESIGN'))
	sample_id = design_header_LIST.pop(design_header_LIST.index('SAMPLE_ID'))
	design_header_LIST.insert(0, sample_id)
	#print design_header_LIST
	#design_header_LIST[0], design_header_LIST[design_header_LIST.index('SAMPLE_ID')] = design_header_LIST[design_header_LIST.index('SAMPLE_ID')], design_header_LIST[0]
	#print design_header_LIST
	statistics_table_header_list = []
	statistics_table_header_list.extend(design_header_LIST)
	statistics_table_header_list.extend(['Total number', 'Classified number', 'Unclassified number', 'Classified percentage', 'Unclassified percentage'])

	#statistics_table_header_list = design_header_LIST['Groups', 'Sample Name', 'Total number', 'Classified number', 'Unclassified number', 'Classified percentage', 'Unclassified percentage']
	no_unclassified_flag = True
	if 'unclassified' in shared_file_vertical_DICT:
		no_unclassified_flag = False
	else:
		no_unclassified_flag = True
	table_head_string = ''
	table_head_string += '<thead>\n'
	table_head_string += '\t<tr>\n'
	for each_head in statistics_table_header_list:
		table_head_string += '\t\t<th class="text-center">' + each_head + '</th>\n'
	table_head_string += '\t</tr>\n'
	table_head_string += '</thead>\n'

	shared_file_horizontal_DICT, shared_file_horizontal_LIST = any_file_to_dict_converter_horizontal(shared_file_PATH, 1)
	#print shared_file_horizontal_DICT

	table_body_string = ''
	table_body_string += '\t<tbody>\n'
	for each_sample in shared_file_horizontal_LIST:
		OTU_value_list = map(int, shared_file_horizontal_DICT[each_sample][3:])
		total_read_count = sum(OTU_value_list)
		if no_unclassified_flag is False:
			unclassified_read_count = int(shared_file_vertical_DICT['unclassified'][each_sample])
		else:
			unclassified_read_count = 0
		classified_read_count = int(total_read_count - unclassified_read_count)
		
		table_body_string += '\t\t<tr>\n'
		for each_header in design_header_LIST:
			table_body_string += '\t\t\t<td class="text-center">' + str(design_MATRIX_reverse[each_sample][each_header]) + '</td>\n'
		#table_body_string += '\t\t\t<td class="text-center">' + str(each_sample) + '</td>\n'
		table_body_string += '\t\t\t<td class="text-center">' + str(total_read_count) + '</td>\n'
		table_body_string += '\t\t\t<td class="text-center">' + str(classified_read_count) + '</td>\n'
		table_body_string += '\t\t\t<td class="text-center">' + str(unclassified_read_count) + '</td>\n'
		table_body_string += '\t\t\t<td class="text-center">' + str(round_float(percentage(classified_read_count, total_read_count))) + '</td>\n'
		table_body_string += '\t\t\t<td class="text-center">' + str(round_float(percentage(unclassified_read_count, total_read_count))) + '</td>\n'
		table_body_string += '\t\t</tr>\n'
	table_body_string += '\t</tbody>\n'

	statistics_table_html_string = """
		<div id="STATISTICS_TABLE" class="container-fluid">
			<div class="row">
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-12">
					<div class="panel panel-default">
						<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
							<h3 class="panel-title">STATISTICS TABLE</h3>
						</div>
						<div class="panel-body" align="center">
							<!-- Table -->
							<div class="table-responsive">
							<table id="statistics_table" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%" style="font-family:'Avenir';">
							""" + table_head_string + """
							""" + table_body_string + """
							</table>

						</div>
						</div>
						<div class="panel-footer">
						</div>
					</div>
				</div>
			</div>
		</div>
	"""
	#print statistics_table_html_string
	return statistics_table_html_string


def plotly_Natural_Abundance_Barplot(shared_file_PATH, design_file_PATH, OTU_metadata_file_PATH):
	#LOAD UP DICT
	# ++++++++++++++++++++++++++++++ DESIGN MATRIX REVERSE AND FORWARD
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	#Forward: {'phinchID': {'24': ['SAMPLE_ALIAS_025'], '25': ['SAMPLE_ALIAS_026'], '26': ['SAMPLE_ALIAS_027'],
	#Reverse: {'SAMPLE_ALIAS_024': {'phinchID': '23', 'temp': '20.8 C', 'collection_date': '2012-12-12T17:00:00+08:00'
	# ----------------------------- END OF DESIGN MATRIX REVERSE AND FORWARD

	OTU_percentile_DICT, OTU_threshold_DICT = shared_file_to_OTU_percentile_dict_converter(shared_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_name_list = []
	for each_head in shared_file_vertical_LIST:
		if each_head.lower() in ['label', 'groups', 'numotus', 'group']:
			continue
		else:
			OTU_name_list.append(each_head)
	OTU_color_dict = color_lover_generate(OTU_name_list)

	#shared_file_horizontal_DICT, shared_file_horizontal_LIST = any_file_to_dict_converter_horizontal(shared_file_PATH, 1)
	sample_name_LIST = shared_file_vertical_DICT['Groups']
	design_sample_name_LIST = []
	for each_sample in sample_name_LIST:
		design_sample_name_LIST.append(design_MATRIX_reverse[each_sample]['SAMPLE_ID'])

	#print OTU_MATRIX_reverse
	#print OTU_MATRIX_forward
	# ######################## VISIBILE
	visibility_list_length = len(OTU_name_list)
	visibility_flag_dict = {}
	visibility_flag_dict['All'] = [True] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'] = [False] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')] = [False] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')] = [False] * visibility_list_length
	visibility_flag_dict['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')] = [False] * visibility_list_length
	# ######################## PLOTLY DATA CREATE
	Natural_Abundance_barplot_TRACE_list = []
	Natural_Abundance_barplot_FIGURE_objects_dict = {}
	Natural_Abundance_barplot_LAYOUT = {}

	for each_OTU in OTU_name_list:
		Natural_Abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
				x=design_sample_name_LIST,
				y=shared_file_vertical_DICT[each_OTU],
				name=OTU_MATRIX_reverse[each_OTU]['OTU_ID'],
				hoverinfo='y+name',
				marker=dict(
					color=OTU_color_dict[each_OTU],
				)
			)
		)
		each_OTU_index = OTU_name_list.index(each_OTU)
		# #############################################
		if each_OTU in OTU_percentile_DICT['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')]:
			visibility_flag_dict['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')]:
			visibility_flag_dict[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')]:
			visibility_flag_dict[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance']:
			visibility_flag_dict[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'][each_OTU_index] = True
		# #############################################
	VISIBLE_flag_list = []
	for each_percentage_dict in visibility_flag_dict:
			VISIBLE_flag_list.append(
				dict(
					args=[
						dict(
							#title=title_dictionary[each_index],
							visible=visibility_flag_dict[each_percentage_dict]
						)
					],
					label=each_percentage_dict,
					method='restyle',
				),
			)

	Natural_Abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
		barmode='stack',
		height=600,
		autosize=True,
		showlegend=True,
		legend=dict(
			font=dict(family='Avenir', size=10),
			orientation="v",
			traceorder="normal"
		),

		title='Bacterial Taxonomic Composition',
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
		),
		font=dict(family='Avenir', size=10),
		yaxis=dict(
			showgrid=True,
			title='Natural Abundance Value',
			titlefont=dict(
				family='Avenir',
				size=15,
			)
		),
		margin=dict(
			l=100,
			r=100,
			b=200,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0,
				visible=True
			)
		]
		)

	)
	Natural_Abundance_barplot_FIGURE_objects_dict['natural_abundance'] = PLOTLY_GO.Figure(data=Natural_Abundance_barplot_TRACE_list, layout=Natural_Abundance_barplot_LAYOUT)
	
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('NATURAL_ABUNDANCE', 'STACKED_BARPLOT', Natural_Abundance_barplot_FIGURE_objects_dict, mode_bar=False)

	return (plotly_html_string, plotly_script)


def plotly_Rarefaction_Curve_Linechart_DROPDOWN(shared_file_PATH, design_MATRIX_forward, design_MATRIX_reverse, mothur_exec_path, processors, outputdir):

	#1 heuristic method to find the nseqs value
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_total_count_LIST = []
	for each_OTU in shared_file_vertical_LIST[3:]:
		OTU_total_count_LIST.append(sum(map(int, shared_file_vertical_DICT[each_OTU])))
	min_threshold = 50
	#print OTU_total_count_LIST
	rare_value = numpy_percentile(OTU_total_count_LIST, min_threshold)
	#sys.exit(2)
	#2 remove rare
	frequency_value = str(int(rare_value))

	# ######################## Rarefaction curve CALCULATION
	#frequency_value = '0.25'
	flag, stderr = execute_functions(mothur_rarefaction_single, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file_PATH, frequency_value)
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
	rarefaction_file = add_extension(scanned_container[0], '.txt')
	# ######################## PARSING
	rarefaction_file_vertical_DICT, rarefaction_file_vertical_LIST = any_file_to_dict_converter_vertical(rarefaction_file)
	number_of_sequences_sampled_list = rarefaction_file_vertical_DICT['numsampled']
	
	raw_sample_name_list = []
	for each_head in rarefaction_file_vertical_LIST:
		if 'lci' in each_head:
			continue
		elif 'hci' in each_head:
			continue
		elif 'numsampled' in each_head:
			continue
		else:
			raw_sample_name_list.append(each_head)

	rarefaction_sample_dict_values = {}
	for each_sample in raw_sample_name_list:
		clean_sample_name = each_sample.split('-')[1]
		rarefaction_sample_dict_values[clean_sample_name] = list_variance(rarefaction_file_vertical_DICT[each_sample], 0.01)  # removing NA element and variant factors

	visibility_list_length = len(rarefaction_sample_dict_values.keys())
	visibility_flag_dict = {}
	visibility_flag_dict['ALL'] = [True] * visibility_list_length
	VISIBLE_flag_list = []
	VISIBLE_flag_list.append(
		dict(
			args=[
				dict(
					#title=title_dictionary[each_index],
					visible=visibility_flag_dict['ALL']
				)
			],
			label='ALL',
			method='restyle',

		),
	)
	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)
	reverse_design_dict = design_MATRIX_forward[selected_design]
	
	for each_design in reverse_design_dict:
		visibility_flag_dict[each_design] = [False] * visibility_list_length
		sample_list = reverse_design_dict[each_design]
		for each_sample in rarefaction_sample_dict_values.keys():
			if each_sample in sample_list:
				each_sample_index = rarefaction_sample_dict_values.keys().index(each_sample)
				visibility_flag_dict[each_design][each_sample_index] = True
		VISIBLE_flag_list.append(
			dict(
				args=[
					dict(
						#title=title_dictionary[each_index],
						visible=visibility_flag_dict[each_design]
					)
				],
				label=each_design,
				method='restyle',

			),
		)
	# ######################## PLOTLY DATA CREATE
	RAREFACTION_curve_linechart_TRACE_list = []
	RAREFACTION_curve_linechart_LAYOUT_objects_dict = {}
	RAREFACTION_curve_linechart_FIGURE_objects_dict = {}
	for each_sample in rarefaction_sample_dict_values.keys():
		RAREFACTION_curve_linechart_TRACE_list.append(PLOTLY_GO.Scatter(
				x=number_of_sequences_sampled_list,
				y=rarefaction_sample_dict_values[each_sample],
				name=design_MATRIX_reverse[each_sample]['SAMPLE_ID'],
				hoverinfo='name',
				mode='lines',
			)
		)
	RAREFACTION_curve_linechart_LAYOUT_objects_dict = PLOTLY_GO.Layout(
		height=600,
		autosize=True,
		title='Rarefaction curve plot',
		showlegend=True,
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=12
		),
		xaxis=dict(
			showgrid=False,
			title='Number of sampling effort',
			titlefont=dict(
				family='Avenir',
				size=10,
			)

		),
		font=dict(family='Avenir', size=10),
		yaxis=dict(
			showgrid=True,
			title='Number of detected species',
			titlefont=dict(
				family='Avenir',
				size=10,
			)
		),
		margin=dict(
			l=100,
			r=100,
			b=100,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0,
				visible=True
			)
		]
		)
	)

	RAREFACTION_curve_linechart_FIGURE_objects_dict['rarefaction'] = PLOTLY_GO.Figure(data=RAREFACTION_curve_linechart_TRACE_list, layout=RAREFACTION_curve_linechart_LAYOUT_objects_dict)
	# ###################################  PLOTLY PLOT CREATE
	plotly_html_string, plotly_script = plotly_html_maker('RAREFACTION_CURVE', 'LINE_CHART', RAREFACTION_curve_linechart_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Alpha_Diversity_Index_Boxplot_DROPDOWN(shared_file_PATH, design_file_PATH, mothur_exec_path, processors, outputdir):
	#LOAD UP THE DICT
	# ++++++++++++++++++++++++++++++ DESIGN MATRIX REVERSE AND FORWARD
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	#Forward: {'phinchID': {'24': ['SAMPLE_ALIAS_025'], '25': ['SAMPLE_ALIAS_026'], '26': ['SAMPLE_ALIAS_027'],
	#Reverse: {'SAMPLE_ALIAS_024': {'phinchID': '23', 'temp': '20.8 C', 'collection_date': '2012-12-12T17:00:00+08:00'
	# ----------------------------- END OF DESIGN MATRIX REVERSE AND FORWARD

	# ######################## ALPHA DIVERSITY INDEX CALCULATION
	flag, stderr = execute_functions(mothur_summary_single, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file_PATH)
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
	flag = remove_extension_files(outputdir, '.rabund')
	alpha_diversity_index_file = add_extension(scanned_container[0], '.txt')
	# ######################## PARSING
	alpha_diversity_index_vertical_DICT, alpha_diversity_index_vertical_LIST = any_file_to_dict_converter_vertical(alpha_diversity_index_file)
	#Select a goo design

	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)

	reverse_design_dict = design_MATRIX_forward[selected_design]
	# ######################## PLOTLY DATA CREATE
	ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict = {}
	ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict = {}
	ALPHA_Diversity_Index_boxplot_TRACE_list = []
	alpha_index_list = ['shannon', 'simpson', 'invsimpson', 'sobs', 'ace']
	visibility_list_length = len(alpha_index_list) * len(reverse_design_dict.keys())
	alpha_visibility_list = []
	# ####################### TO GENERATE SWITCHING STRING
	VISIBLE_flag_list = []
	visibility_flag_dict = {}
	alpha_visibility_list = [False] * visibility_list_length
	step_list = range(len(reverse_design_dict.keys()))
	visible_index = 0
	for i in xrange(0, visibility_list_length, len(reverse_design_dict.keys())):
		for each_step in step_list:
			alpha_visibility_list[i + each_step] = True
		visibility_flag_dict[visible_index] = alpha_visibility_list
		alpha_visibility_list = [False] * visibility_list_length
		visible_index += 1
	# #######################
	title_dictionary = {}
	title_dictionary['simpson'] = 'Simpson diversity index'
	title_dictionary['invsimpson'] = 'Inverse-Simpson diversity index'
	title_dictionary['shannon'] = 'Shannon diversity index'
	title_dictionary['sobs'] = 'Sobs(observed number of species) richness index'
	title_dictionary['ace'] = 'ACE(Abundance-base Coverage Estimator) richness index'
	for each_index in alpha_index_list:
		for each_design in reverse_design_dict:
			sample_list = reverse_design_dict[each_design]
			index_value_list = []
			for each_sample in sample_list:
				each_sample_index = alpha_diversity_index_vertical_DICT['group'].index(each_sample)
				index_value_list.append(alpha_diversity_index_vertical_DICT[each_index][each_sample_index])

			ALPHA_Diversity_Index_boxplot_TRACE_list.append(PLOTLY_GO.Box(
				x=each_design,
				y=index_value_list,
				#text=each_design,
				hoverinfo='y',
				name=each_design,
				boxmean='sd',
				visible=visibility_flag_dict[alpha_index_list.index(each_index)][reverse_design_dict.keys().index(each_design) * alpha_index_list.index(each_index)],
				legendgroup=title_dictionary[each_index]
			)
			)
		VISIBLE_flag_list.append(
			dict(
				args=[
					dict(
						#title=title_dictionary[each_index],
						visible=visibility_flag_dict[alpha_index_list.index(each_index)]
					)
				],
				label=title_dictionary[each_index],
				method='restyle',

			),
		)

	ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict = PLOTLY_GO.Layout(
		boxmode='group',
		height=600,
		autosize=True,
		showlegend=False,
		title='ALPHA Diversity index',
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
		),
		font=dict(family='Avenir', size=12),
		yaxis=dict(
			showgrid=True,
			#title=each_index + ' Index Value',
			titlefont=dict(
				family='Avenir',
				size=16,
			)
		),
		margin=dict(
			l=40,
			r=30,
			b=80,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0
			)
		])
	)
	#print VISIBLE_flag_list
	ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict['alpha_diversity_index'] = PLOTLY_GO.Figure(data=ALPHA_Diversity_Index_boxplot_TRACE_list, layout=ALPHA_Diversity_Index_boxplot_LAYOUT_objects_dict)
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('ALPHA_DIVERSITY_INDEX', 'BOX_PLOT', ALPHA_Diversity_Index_boxplot_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def biomarker_discovery_html_table_function(significant_OTUs_DICT, OTU_metadata_file_PATH):
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	OTU_name_LIST = OTU_MATRIX_reverse.keys()
	
	biomarker_discovery_table_header_list = ['OTU name', 'metastats pValue', 'Kruskal-wallis pValue', 'indicator pValue', 'lefse pValue']
	significant_OTU_LIST = []
	significant_OTU_LIST.extend(significant_OTUs_DICT['metastats'].keys())
	significant_OTU_LIST.extend(significant_OTUs_DICT['kruskal_wallis'].keys())
	significant_OTU_LIST.extend(significant_OTUs_DICT['indicator'].keys())
	significant_OTU_LIST.extend(significant_OTUs_DICT['lefse'].keys())
	significant_OTU_LIST = list(set(significant_OTU_LIST))

	table_head_string = ''
	table_head_string += '<thead>\n'
	table_head_string += '\t<tr>\n'
	for each_head in biomarker_discovery_table_header_list:
		table_head_string += '\t\t<th class="text-center">' + each_head + '</th>\n'
	table_head_string += '\t</tr>\n'
	table_head_string += '</thead>\n'

	table_body_string = ''
	table_body_string += '\t<tbody>\n'
	for each_OTU in significant_OTU_LIST:
		table_body_string += '\t\t<tr>\n'

		table_body_string += '\t\t\t<td class="text-center">' + str(OTU_MATRIX_reverse[each_OTU]['OTU_ID']) + '</td>\n'
		if each_OTU in significant_OTUs_DICT['metastats'].keys():
			table_body_string += '\t\t\t<td class="text-center">' + str(significant_OTUs_DICT['metastats'][each_OTU][1]) + '</td>\n'
		else:
			table_body_string += '\t\t\t<td class="text-center">' + str('NS') + '</td>\n'
		if each_OTU in significant_OTUs_DICT['kruskal_wallis'].keys():
			table_body_string += '\t\t\t<td class="text-center">' + str(significant_OTUs_DICT['kruskal_wallis'][each_OTU][1]) + '</td>\n'
		else:
			table_body_string += '\t\t\t<td class="text-center">' + str('NS') + '</td>\n'
		if each_OTU in significant_OTUs_DICT['indicator'].keys():
			table_body_string += '\t\t\t<td class="text-center">' + str(significant_OTUs_DICT['indicator'][each_OTU][1]) + '</td>\n'
		else:
			table_body_string += '\t\t\t<td class="text-center">' + str('NS') + '</td>\n'
		if each_OTU in significant_OTUs_DICT['lefse'].keys():
			table_body_string += '\t\t\t<td class="text-center">' + str(significant_OTUs_DICT['lefse'][each_OTU][1]) + '</td>\n'
		else:
			table_body_string += '\t\t\t<td class="text-center">' + str('NS') + '</td>\n'
		table_body_string += '\t\t</tr>\n'
	table_body_string += '\t</tbody>\n'

	biomarker_table_html_string = """
		<div id="BIOMARKER_TABLE" class="container-fluid">
			<div class="row">
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-12">
					<div class="panel panel-default">
						<div class="panel-heading" style="border:0px; background-color: #F4F9FD;">
							<h3 class="panel-title">STATISTICS TABLE</h3>
						</div>
						<div class="panel-body" align="center">
							<!-- Table -->
							<div class="table-responsive">
							<table id="biomarker_table" class="table table-striped table-bordered table-hover small" cellspacing="0" width="100%" style="font-family:'Avenir';">
							""" + table_head_string + """
							""" + table_body_string + """
							</table>

						</div>
						</div>
						<div class="panel-footer">
						</div>
					</div>
				</div>
			</div>
		</div>
	"""
	return biomarker_table_html_string


def plotly_Relative_abundance_Barplot(shared_file_PATH, design_file_PATH, OTU_metadata_file_PATH):

	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)

	condensed_shared_DICT = {}
	condensed_shared_DICT = condense_shared_file_by_design(shared_file_PATH, design_MATRIX_forward[selected_design])

	#print condensed_shared_dict
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)

	design_name_LIST = condensed_shared_DICT.keys()
	OTU_name_list = condensed_shared_DICT[design_name_LIST[0]].keys()
	
	OTU_name_metadata_LIST = []
	for each_OTU in OTU_name_list:
		OTU_name_metadata_LIST.append(OTU_MATRIX_reverse[each_OTU]['OTU_ID'])

	longest_OTU_name = get_longest_length_string_list(OTU_name_metadata_LIST)
	#print longest_OTU_name
	
	Relative_abundance_barplot_TRACE_list = []
	Relative_abundance_barplot_LAYOUT = {}
	Relative_abundance_barplot_FIGURE_objects_dict = {}

	for each_design in design_name_LIST:
		Relative_abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
			x=OTU_name_metadata_LIST,
			y=condensed_shared_DICT[each_design].values(),
			name=each_design,
			hoverinfo='all',
		)
		)
	Relative_abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
		barmode='group',
		height=600,
		autosize=True,
		title='Relative abundance',
		showlegend=True,
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
			#color='#7f7f7f'
		),
		font=dict(family='Avenir', size=10),
		xaxis=dict(
			autorange=True,
		),
		yaxis=dict(
			title='Relative abundance',
			titlefont=dict(family='Avenir', size=12),
		),
		margin=dict(
			l=100,
			r=100,
			b=longest_OTU_name * 3,
			t=100,
		)
	)

	Relative_abundance_barplot_FIGURE_objects_dict['Relative_abundance'] = PLOTLY_GO.Figure(data=Relative_abundance_barplot_TRACE_list, layout=Relative_abundance_barplot_LAYOUT)

	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('RELATIVE_ABUNDANCE', 'BARPLOT', Relative_abundance_barplot_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Relative_Abundance_Piechart(shared_file_PATH, design_file_PATH, OTU_metadata_file_PATH):
	# ############################### CONDENSE SHARED FILE
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)

	condensed_shared_DICT = {}
	condensed_shared_DICT = condense_shared_file_by_design(shared_file_PATH, design_MATRIX_forward[selected_design])
	#print condensed_shared_dict
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	design_name_LIST = condensed_shared_DICT.keys()
	OTU_name_list = condensed_shared_DICT[design_name_LIST[0]].keys()

	OTU_name_metadata_LIST = []
	for each_OTU in OTU_name_list:
		OTU_name_metadata_LIST.append(OTU_MATRIX_reverse[each_OTU]['OTU_ID'])

	longest_OTU_name = get_longest_length_string_list(OTU_name_metadata_LIST)
	# ###################################  PLOTLY DATA CREATE
	Relative_Abundance_Piechart_TRACE_objects_dict = {}
	#Relative_Abundance_Piechart_LAYOUT_objects_dict = {}
	Relative_Abundance_Piechart_FIGURE_objects_dict = {}
	Relative_Abundance_Piechart_ANNOTATIONS_objects_dict = {}

	Domain_step = 1.0 / len(design_name_LIST)
	Domain_tolerance = 0.01
	#Domain_step = 0.2
	Domain_min = 0
	Domain_max = Domain_step
	Domain_list = [Domain_min, Domain_max]
	Domain_mid = (Domain_max + Domain_min) / 2

	for each_design in design_name_LIST:
		Relative_Abundance_Piechart_TRACE_objects_dict[each_design] = {}
		#Relative_Abundance_Piechart_LAYOUT_objects_dict[each_design] = {}
		Relative_Abundance_Piechart_ANNOTATIONS_objects_dict[each_design] = {}
		Relative_Abundance_Piechart_TRACE_objects_dict[each_design] = PLOTLY_GO.Pie(
			labels=OTU_name_metadata_LIST,
			values=condensed_shared_DICT[each_design].values(),
			name=each_design,
			type='pie',
			hoverinfo="label+percent+name",
			textinfo="all",
			hole=0.4,
			textposition="inside",
			domain=dict(
				x=Domain_list,
				y=[0.0, 1.0]
			)
		)
		Relative_Abundance_Piechart_ANNOTATIONS_objects_dict[each_design] = dict(
			showarrow=False,
			text=each_design,
			xanchor="center",
			xref="paper",
			yanchor="bottom",
			yref="paper",
			y=0.5,
			x=Domain_mid
		)

		Domain_tolerance = 0.01
		Domain_step = 1.0 / len(design_name_LIST)
		Domain_min = Domain_max + Domain_tolerance
		Domain_max = Domain_max + Domain_step + Domain_tolerance
		Domain_list = [Domain_min, Domain_max]
		Domain_mid = (Domain_min + Domain_max) / 2

	Relative_Abundance_Piechart_LAYOUT = PLOTLY_GO.Layout(
		title='Bacterial Taxonomic Composition',
		autosize=True,
		height=700,
		titlefont=dict(
			family='Avenir',
			size=16
		),
		showlegend=True,
		legend=dict(
			font=dict(family='Avenir', size=12),
			orientation="v",
			traceorder="normal"
		),
		font=dict(family='Avenir', size=15),
		margin=dict(
			l=100,
			r=100,
			b=100,
			t=100,
		),
		annotations=PLOTLY_GO.Annotations(Relative_Abundance_Piechart_ANNOTATIONS_objects_dict.values())
	)
	Relative_Abundance_Piechart_FIGURE_objects_dict['pie'] = PLOTLY_GO.Figure(data=Relative_Abundance_Piechart_TRACE_objects_dict.values(), layout=Relative_Abundance_Piechart_LAYOUT)
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('Relative_Abundance', 'PIECHART', Relative_Abundance_Piechart_FIGURE_objects_dict, mode_bar=False)
	return (plotly_html_string, plotly_script)


def plotly_Natural_Abundance_Barplot_multi(multi_sheets_shared_file_PATH, design_file_PATH, OTU_metadata_file_PATH):
	
	# ++++++++++++++++++++++++++++++ DESIGN MATRIX REVERSE AND FORWARD
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	#Forward: {'phinchID': {'24': ['SAMPLE_ALIAS_025'], '25': ['SAMPLE_ALIAS_026'], '26': ['SAMPLE_ALIAS_027'],
	#Reverse: {'SAMPLE_ALIAS_024': {'phinchID': '23', 'temp': '20.8 C', 'collection_date': '2012-12-12T17:00:00+08:00'
	# ----------------------------- END OF DESIGN MATRIX REVERSE AND FORWARD

	#Fisrst lets grab all OTU names
	OTU_name_list
	
	#STEP1: LOAD THE DICT
	multi_shared_file_DICT, multi_shared_file_LIST = multi_shared_file_to_dict_vertical(multi_sheets_shared_file_PATH)
	Natural_Abundance_barplot_FIGURE_objects_dict = {}
	visibility_flag_dict = {}
	for each_lineage in multi_shared_file_LIST:
		#
		OTU_name_list = []
		for each_head in multi_shared_file_DICT[each_lineage][1]:
			if each_head.lower() in ['label', 'groups', 'numotus', 'group']:
				continue
			else:
				OTU_name_list.append(each_head)
		OTU_color_dict = color_lover_generate(OTU_name_list)
		#
		sample_name_LIST = multi_shared_file_DICT[each_lineage][0]['Groups']
		print sample_name_LIST
		design_sample_name_LIST = []
		for each_sample in sample_name_LIST:
			design_sample_name_LIST.append(design_MATRIX_reverse[each_sample]['SAMPLE_ID'])
		# ######################## PLOTLY DATA CREATE
		Natural_Abundance_barplot_TRACE_list = []
		
		Natural_Abundance_barplot_LAYOUT = {}
		for each_OTU in OTU_name_list:
			Natural_Abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
					x=design_sample_name_LIST,
					y=multi_shared_file_DICT[each_lineage][0][each_OTU],
					name=each_OTU,
					hoverinfo='y+name',
					marker=dict(
						color=OTU_color_dict[each_OTU],
					)
				)
			)
		Natural_Abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
			barmode='stack',
			height=600,
			autosize=True,
			showlegend=True,
			legend=dict(
				font=dict(family='Avenir', size=10),
				orientation="v",
				traceorder="normal"
			),

			title='Bacterial Taxonomic Composition',
			hovermode='closest',
			titlefont=dict(
				family='Avenir',
				size=16
			),
			font=dict(family='Avenir', size=10),
			yaxis=dict(
				showgrid=True,
				title='Natural Abundance Value',
				titlefont=dict(
					family='Avenir',
					size=15,
				)
			),
			margin=dict(
				l=100,
				r=100,
				b=200,
				t=100,
			),

		)
		Natural_Abundance_barplot_FIGURE_objects_dict[each_lineage] = PLOTLY_GO.Figure(data=Natural_Abundance_barplot_TRACE_list, layout=Natural_Abundance_barplot_LAYOUT)
	#print Natural_Abundance_barplot_FIGURE_objects_dict.keys()
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('NATURAL_ABUNDANCE', 'STACKED_BARPLOT', Natural_Abundance_barplot_FIGURE_objects_dict, mode_bar=False)

	return (plotly_html_string, plotly_script)
	"""




	#LOAD UP DICT
	

	OTU_percentile_DICT, OTU_threshold_DICT = shared_file_to_OTU_percentile_dict_converter(shared_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_name_list = []
	for each_head in shared_file_vertical_LIST:
		if each_head.lower() in ['label', 'groups', 'numotus', 'group']:
			continue
		else:
			OTU_name_list.append(each_head)
	OTU_color_dict = color_lover_generate(OTU_name_list)

	#shared_file_horizontal_DICT, shared_file_horizontal_LIST = any_file_to_dict_converter_horizontal(shared_file_PATH, 1)
	sample_name_LIST = shared_file_vertical_DICT['Groups']
	design_sample_name_LIST = []
	for each_sample in sample_name_LIST:
		design_sample_name_LIST.append(design_MATRIX_reverse[each_sample]['SAMPLE_ID'])

	#print OTU_MATRIX_reverse
	#print OTU_MATRIX_forward
	# ######################## VISIBILE
	visibility_list_length = len(OTU_name_list)
	visibility_flag_dict = {}
	visibility_flag_dict['All'] = [True] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'] = [False] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')] = [False] * visibility_list_length
	visibility_flag_dict[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')] = [False] * visibility_list_length
	visibility_flag_dict['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')] = [False] * visibility_list_length
	# ######################## PLOTLY DATA CREATE
	Natural_Abundance_barplot_TRACE_list = []
	Natural_Abundance_barplot_FIGURE_objects_dict = {}
	Natural_Abundance_barplot_LAYOUT = {}

	for each_OTU in OTU_name_list:
		Natural_Abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
				x=design_sample_name_LIST,
				y=shared_file_vertical_DICT[each_OTU],
				name=OTU_MATRIX_reverse[each_OTU]['OTU_ID'],
				hoverinfo='y+name',
				marker=dict(
					color=OTU_color_dict[each_OTU],
				)
			)
		)
		each_OTU_index = OTU_name_list.index(each_OTU)
		# #############################################
		if each_OTU in OTU_percentile_DICT['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')]:
			visibility_flag_dict['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')]:
			visibility_flag_dict[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')]:
			visibility_flag_dict[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')][each_OTU_index] = True
		
		if each_OTU in OTU_percentile_DICT[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance']:
			visibility_flag_dict[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'][each_OTU_index] = True
		# #############################################
	VISIBLE_flag_list = []
	for each_percentage_dict in visibility_flag_dict:
			VISIBLE_flag_list.append(
				dict(
					args=[
						dict(
							#title=title_dictionary[each_index],
							visible=visibility_flag_dict[each_percentage_dict]
						)
					],
					label=each_percentage_dict,
					method='restyle',
				),
			)

	Natural_Abundance_barplot_LAYOUT = PLOTLY_GO.Layout(
		barmode='stack',
		height=600,
		autosize=True,
		showlegend=True,
		legend=dict(
			font=dict(family='Avenir', size=10),
			orientation="v",
			traceorder="normal"
		),

		title='Bacterial Taxonomic Composition',
		hovermode='closest',
		titlefont=dict(
			family='Avenir',
			size=16
		),
		font=dict(family='Avenir', size=10),
		yaxis=dict(
			showgrid=True,
			title='Natural Abundance Value',
			titlefont=dict(
				family='Avenir',
				size=15,
			)
		),
		margin=dict(
			l=100,
			r=100,
			b=200,
			t=100,
		),
		updatemenus=list([
			dict(
				direction="down",
				y=1,
				x=-0.05,
				font=dict(family='Avenir', size=12),
				buttons=VISIBLE_flag_list,
				type="buttons",
				active=0,
				visible=True
			)
		]
		)

	)
	Natural_Abundance_barplot_FIGURE_objects_dict['natural_abundance'] = PLOTLY_GO.Figure(data=Natural_Abundance_barplot_TRACE_list, layout=Natural_Abundance_barplot_LAYOUT)
	
	# ###################################  PLOTLY PLOT CREATE ###############################
	plotly_html_string, plotly_script = plotly_html_maker('NATURAL_ABUNDANCE', 'STACKED_BARPLOT', Natural_Abundance_barplot_FIGURE_objects_dict, mode_bar=False)

	return (plotly_html_string, plotly_script)
	"""
# ################################### VALIDATOR ############################## #


def design_file_validator(design_file_PATH, design_file_VALIDATED_PATH):
	design_file_DELIMITER_list = [',', '\t']
	design_file_DELIMITER = sniff_delimiter(design_file_PATH, design_file_DELIMITER_list)
	if design_file_DELIMITER not in ['\t']:
		print "design file delimiter: ", design_file_DELIMITER, " is not standard, It Should be Tab delimited"
		report("design file delimiter: " + design_file_DELIMITER + " is not standard, It Should be Tab delimited")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)

	#print design_file_OBJECT_dict
	print "Design file validation is in progress:"
	report("Design file validation is in progress:")

	design_file_HANDLE = open(design_file_PATH, 'rU')
	for each_line in design_file_HANDLE:
		line = each_line.split(design_file_DELIMITER)
		if len(line) < 2:
			print "Design file validation FAILED"
			error("Design file validation FAILED")
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
	design_file_HANDLE.close()

	print "Design file validation PASSED!!!"
	report("Design file validation PASSED!!!")
	copy_file(design_file_PATH, design_file_VALIDATED_PATH)
	return True


def taxonomy_reference_validator(taxonomy_reference_file_PATH, validated_taxonomy_reference_file_PATH):
	#reading taxonomy file and fix the taxonomy level
	return 1
# ################################### CONVERTER ############################## #


def biom_to_tsv_converter(processors, outputdir, stderr, stdout, run_pid, biom_exec_path, biom_file_PATH, tsv_file_PATH):

	space = ' '
	biom_string = ''
	biom_string += 'nohup' + space
	biom_string += biom_exec_path + space
	biom_string += 'convert -i' + space + biom_file_PATH + space + '-o' + space + tsv_file_PATH + space + '--to-tsv --header-key taxonomy --table-type="OTU table"'
	print "EXECUTING: ", biom_string
	report("EXECUTING: " + biom_string)
	exec_dict = {}
	exec_dict = execute([biom_string])
	if exec_dict["exitCode"] != 0:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	else:
		return True


def biom_to_shared_converter_old(biom_file_PATH, prefix_NAME, biom_file_TSV_FORMAT_PATH, biom_file_SHARED_FORMAT_PATH, processors, outputdir):
	#STEP1:Convert biom to TSV
	print "biom_to_tsv_converter:"
	biom_exec_path = 'biom'
	flag, stderr = execute_functions(biom_to_tsv_converter, processors, outputdir, 'multi', 'mafft', biom_exec_path, biom_file_PATH, biom_file_TSV_FORMAT_PATH)
	if flag is False:
		print "Execution of biom_to_tsv_converter failed!!!"
		error("Execution of biom_to_tsv_converter failed!!!")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		report("biom_to_tsv_converter executed successfully!!!")
		print "biom_to_tsv_converter executed successfully!!!"

	#STEP2: LOADING TABLE INTO PANDAS DATAFRAME
	biom_pandas_DATAFRAME = pandas.read_table(biom_file_TSV_FORMAT_PATH, index_col=0, skiprows=1, low_memory=False, encoding='utf-8')

	OTU_name_list = []
	OTU_name_list = biom_pandas_DATAFRAME.index.values.tolist()

	Sample_name_list = []
	Sample_name_list = biom_pandas_DATAFRAME.columns[:-1].values.tolist()

	Taxonomy_name_list = []
	Taxonomy_name_list = biom_pandas_DATAFRAME.taxonomy.values.tolist()

	#Merging OTU and TAX
	OTU_TAX_name_list = []
	for OTU, TAX in itertools.izip(OTU_name_list, Taxonomy_name_list):
		OTU_TAX_name = ''
		if pandas.isnull(OTU) is False:
			OTU_TAX_name += slugify(OTU)
			#OTU_TAX_name += OTU
		if pandas.isnull(TAX) is False:
			OTU_TAX_name += ';' + slugify(TAX)
			#OTU_TAX_name += ';' + TAX
			if OTU_TAX_name[-1] != ';':
				OTU_TAX_name += ';'
		if OTU_TAX_name == '':
			OTU_TAX_name += 'NO_OTU_' + str(numpy.random.randint(1, 100000))
		OTU_TAX_name_list.append(OTU_TAX_name)

	#CREATING SHARED STRING
	shared_file_string = 'label' + '\t' + 'Groups' + '\t' + 'numOtus' + '\t' + list_to_string(OTU_TAX_name_list, '\t') + '\n'
	numOtus_value = str(len(OTU_name_list))

	for each_Sample in Sample_name_list:
		sample_value_list = []
		sample_value_list = biom_pandas_DATAFRAME[each_Sample].values.tolist()

		sample_value_list = map(int, sample_value_list)
		sample_value_list = map(str, sample_value_list)
		#print map(str, sample_value_list)
		#print len(sample_value_list)
		shared_file_string += 'BIOM_SLICER' + '\t' + slugify(each_Sample) + '\t' + numOtus_value + '\t' + list_to_string(sample_value_list, '\t') + '\n'

	#shared_file_string = shared_file_string[:-1]
	write_string_down(shared_file_string, biom_file_SHARED_FORMAT_PATH)
	return True


def design_to_control_converter(design_file_PATH, control_file_PATH):
	string = ''
	f = open(design_file_PATH, 'rU')
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		string += line[0] + '\n'
	write_string_down(string, control_file_PATH)
	return True


def shared_dict_to_file_converter(shared_file_PATH, shared_HEADER_list, length_corrected_HEADER_list, shared_CONTENT_dict, prefix_NAME):
	shared_file_CONTENT_string = ''
	shared_file_CONTENT_string = list_to_string(length_corrected_HEADER_list) + '\n'
	shared_file_ROW_length = len(shared_CONTENT_dict['label'])
	numOTUS_value = str(len(shared_HEADER_list) - 3)
	for each_ROW in range(shared_file_ROW_length):
		shared_line = ''
		for each_HEAD in shared_HEADER_list:
			if each_HEAD not in shared_CONTENT_dict:
				continue
			elif each_HEAD == 'label':
				shared_line += prefix_NAME + '\t'
			elif each_HEAD == 'numOtus':
				shared_line += numOTUS_value + '\t'
			else:
				shared_line += shared_CONTENT_dict[each_HEAD][each_ROW] + '\t'
		shared_line = shared_line[:-1]
		shared_line += '\n'
		shared_file_CONTENT_string += shared_line
	write_string_down(shared_file_CONTENT_string, shared_file_PATH)
	return True


def filter_list_to_control_file_converter(filter_list, control_file):
	string = ''
	string = '\n'.join(filter_list)
	write_string_down(string, control_file)
	return True


def list_string_to_float_converter(list_string):
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


def shared_file_to_dict_converter(shared_file_PATH):
	#Read mothur shared file into dict and header list
	shared_DATAFRAME = pandas.read_table(shared_file_PATH, low_memory=False, encoding='utf-8', skip_blank_lines=True)

	header_list = []
	header_list = shared_DATAFRAME.columns.values.tolist()
	shared_DICT = {}
	for each_column in header_list:
		each_column_value_list = shared_DATAFRAME[each_column].values.tolist()
		shared_DICT[str(each_column)] = map(str, each_column_value_list)
	return shared_DICT


def design_file_to_dict_converter(design_file_PATH):
	#Reading
	design_DATAFRAME = get_pandas_DATAFRAME(design_file_PATH)
	header_list = []
	header_list = design_DATAFRAME.columns.values.tolist()
	design_DICT = {}
	for each_column in header_list:
		each_column_value_list = design_DATAFRAME[each_column].values.tolist()
		design_DICT[str(each_column)] = map(str, each_column_value_list)
	header_LIST = map(str, header_list)
	return (design_DICT, header_LIST)


def any_file_to_dict_converter_vertical(file_PATH, index_row=None):
	any_DATAFRAME = get_pandas_DATAFRAME(file_PATH)
	any_vertical_LIST = []
	any_vertical_LIST = any_DATAFRAME.columns.values.tolist()
	any_vertical_DICT = {}
	for each_column in any_vertical_LIST:
		each_column_value_list = any_DATAFRAME[each_column].values.tolist()
		any_vertical_DICT[str(each_column)] = map(str, each_column_value_list)
	any_vertical_LIST = map(str, any_vertical_LIST)
	return (any_vertical_DICT, any_vertical_LIST)


def any_file_to_dict_converter_horizontal(file_PATH, index_col=None):
	any_DATAFRAME = get_pandas_DATAFRAME(file_PATH)
	if index_col is None:
		index_col = 0
	any_horizontal_DICT = {}
	any_horizontal_LIST = []
	for index, row in any_DATAFRAME.iterrows():
		row = map(str, row.tolist())
		any_horizontal_DICT[row[index_col]] = row
		any_horizontal_LIST.append(row[index_col])
	return (any_horizontal_DICT, any_horizontal_LIST)


def biom_to_shared_converter(biom_file_PATH, shared_file_PATH, TAXONOMY_EXIST, sample_metadata_file_PATH, OTU_metadata_file_PATH, multi_sheets_shared_file_PATH, multi_sheets_OTU_metadata_file_PATH):
	#1: load biom table
	biom_table_HANDLE = biom.load_table(biom_file_PATH)

	#2: Read meta data from BIOM file
	sample_name_LIST = map(str, biom_table_HANDLE.ids(axis='sample'))
	OTU_name_LIST = map(str, biom_table_HANDLE.ids(axis='observation'))
	#3: OTU metadata extraction
	OTU_metadata_string = ''
	OTU_alias_LIST = []
	if biom_table_HANDLE.metadata(axis='observation') is not None:
		metadata_OTU_OBJECT = biom_table_HANDLE.metadata(axis='observation')
		OTU_metadata_string = 'OTU_ALIAS' + '\t' + 'OTU_ID' + '\t' + list_to_string(metadata_OTU_OBJECT[0].keys()) + '\n'
		for each_OTU in OTU_name_LIST:
			each_OTU_index = OTU_name_LIST.index(each_OTU)
			OTU_alias_name = 'OTU_ALIAS_' + str(each_OTU_index + 1).zfill(3)
			OTU_alias_LIST.append(OTU_alias_name)
			OTU_metadata_string += OTU_alias_name + '\t' + each_OTU + '\t'
			for each_metadata_OBJECT_key, each_metadata_OBJECT_value in metadata_OTU_OBJECT[each_OTU_index].iteritems():
				if each_metadata_OBJECT_key.lower() == 'taxonomy':
					TAXONOMY_EXIST = True
					corrected_taxonomy_string = remove_ambiguity_from_taxonomy(each_metadata_OBJECT_value)
					OTU_metadata_string += corrected_taxonomy_string + '\t'
				elif type(each_metadata_OBJECT_value) is list:
					OTU_metadata_string += list_to_string(each_metadata_OBJECT_value, ';') + '\t'
				else:
					OTU_metadata_string += str(each_metadata_OBJECT_value) + '\t'
			OTU_metadata_string = OTU_metadata_string[:-1]
			OTU_metadata_string += '\n'
	else:
		OTU_metadata_string = 'OTU_ALIAS' + '\t' + 'OTU_ID' + '\n'
		for each_OTU in OTU_name_LIST:
			each_OTU_index = OTU_name_LIST.index(each_OTU)
			OTU_alias_name = 'OTU_ALIAS_' + str(each_OTU_index + 1).zfill(3)
			OTU_alias_LIST.append(OTU_alias_name)
			OTU_metadata_string += OTU_alias_name + '\t' + each_OTU + '\n'
	write_string_down(OTU_metadata_string, OTU_metadata_file_PATH)

	#4: SAMPLE metadata extraction
	sample_metadata_string = ''
	sample_alias_DICT = {}
	if biom_table_HANDLE.metadata(axis='sample') is not None:
		metadata_sample_OBJECT = biom_table_HANDLE.metadata(axis='sample')
		sample_metadata_string = 'SAMPLE_ID' + '\t' + 'SAMPLE_ALIAS' + '\t' + list_to_string(metadata_sample_OBJECT[0].keys()) + '\n'
		for each_sample in sample_name_LIST:
			each_sample_index = sample_name_LIST.index(each_sample)
			sample_alias_name = 'SAMPLE_ALIAS_' + str(each_sample_index + 1).zfill(3)
			sample_alias_DICT[each_sample] = sample_alias_name
			sample_metadata_string += each_sample + '\t' + sample_alias_name + '\t' + list_to_string(metadata_sample_OBJECT[each_sample_index].values()) + '\n'
	else:
		sample_metadata_string = 'SAMPLE_ID' + '\t' + 'SAMPLE_ALIAS' + '\t' + 'DEFAULT_DESIGN' + '\n'
		for each_sample in sample_name_LIST:
			each_sample_index = sample_name_LIST.index(each_sample)
			sample_alias_name = 'SAMPLE_ALIAS_' + str(each_sample_index + 1).zfill(3)
			sample_alias_DICT[each_sample] = sample_alias_name
			sample_metadata_string += each_sample + '\t' + sample_alias_name + '\t' + 'BIOM_SLICER' + '\n'
	write_string_down(sample_metadata_string, sample_metadata_file_PATH)

	#5: generate shared table
	shared_file_string = ''
	shared_file_string += 'label' + '\t' + 'Groups' + '\t' + 'numOtus' + '\t' + list_to_string(OTU_alias_LIST, '\t') + '\n'
	numOtus_value = str(len(OTU_alias_LIST))
	for each_sample in sample_name_LIST:
		sample_value_list = biom_table_HANDLE.data(each_sample, axis='sample').tolist()
		sample_value_list = map(int, sample_value_list)
		shared_file_string += 'BIOM_SLICER' + '\t' + sample_alias_DICT[each_sample] + '\t' + numOtus_value + '\t' + list_to_string(sample_value_list, '\t') + '\n'
	write_string_down(shared_file_string, shared_file_PATH)
	#MULTI SHETTS
	if TAXONOMY_EXIST is True:
		multi_sheets_shared_file_function(shared_file_PATH, OTU_metadata_file_PATH, multi_sheets_shared_file_PATH, multi_sheets_OTU_metadata_file_PATH)
	return True


def biom_to_pandas_converter(biom_file_PATH):
	#1: load biom table
	biom_table_HANDLE = biom.load_table(biom_file_PATH)
	
	#2: Read meta data from BIOM file
	sample_ID_LIST = map(str, biom_table_HANDLE.ids(axis='sample'))
	OTU_ID_LIST = map(str, biom_table_HANDLE.ids(axis='observation'))
	#print OTU_ID_LIST
	#['0', '1', '2', '3', '4', '5', '6', '7'
	#print sample_ID_LIST
	#['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593', 'PC.355', 'PC.607', 'PC.634']

	#3: CREATE DATA DICT
	BIOM_data_DICT = {}
	for each_sample in sample_ID_LIST:



	sys.exit(2)



def design_file_to_design_MATRIX_converter(design_file_PATH):
	design_Vertical_DICT, design_Vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	design_MATRIX_forward = {}
	design_MATRIX_reverse = {}
	sample_name_LIST = design_Vertical_DICT[design_Vertical_LIST[0]]
	
	for each_design in design_Vertical_LIST[1:]:
		design_MATRIX_forward[each_design] = {}
		for each_sample in sample_name_LIST:
			each_sample_index = design_Vertical_DICT[design_Vertical_LIST[0]].index(each_sample)
			group_name = design_Vertical_DICT[each_design][each_sample_index]
			if group_name not in design_MATRIX_forward[each_design]:
				design_MATRIX_forward[each_design][group_name] = []
				design_MATRIX_forward[each_design][group_name].append(each_sample)
			else:
				design_MATRIX_forward[each_design][group_name].append(each_sample)
	#print design_MATRIX_forward
	for each_sample in sample_name_LIST:
		each_sample_index = design_Vertical_DICT[design_Vertical_LIST[0]].index(each_sample)
		design_MATRIX_reverse[each_sample] = {}
		for each_design in design_Vertical_LIST[1:]:
			design_MATRIX_reverse[each_sample][each_design] = design_Vertical_DICT[each_design][each_sample_index]
	#print design_MATRIX_reverse
	return(design_MATRIX_forward, design_MATRIX_reverse)


def shared_file_to_OTU_percentile_dict_converter(shared_file_PATH, mode=None):

	OTU_value_DICT = {}
	shared_file_Vertical_DICT, shared_file_Vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	for each_OTU in shared_file_Vertical_LIST[3:]:
		OTU_value_LIST = []
		OTU_value_LIST = map(int, shared_file_Vertical_DICT[each_OTU])
		OTU_total_value = sum(OTU_value_LIST)
		OTU_value_DICT[each_OTU] = OTU_total_value
	OTU_threshold_DICT = {}
	OTU_threshold_DICT['25%'] = numpy_percentile(OTU_value_DICT.values(), 25)
	OTU_threshold_DICT['50%'] = numpy_percentile(OTU_value_DICT.values(), 50)
	OTU_threshold_DICT['75%'] = numpy_percentile(OTU_value_DICT.values(), 75)
	
	OTU_percentile_DICT = {}
	OTU_percentile_DICT['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')] = []
	OTU_percentile_DICT[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')] = []
	OTU_percentile_DICT[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')] = []
	OTU_percentile_DICT[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'] = []

	for each_OTU in OTU_value_DICT:
		if OTU_value_DICT[each_OTU] < OTU_threshold_DICT['25%']:
			OTU_percentile_DICT['OTU Abundance < ' + format(OTU_threshold_DICT['25%'], ',')].append(each_OTU)
		
		if OTU_value_DICT[each_OTU] > OTU_threshold_DICT['25%'] and OTU_value_DICT[each_OTU] <= OTU_threshold_DICT['50%']:
			OTU_percentile_DICT[format(OTU_threshold_DICT['25%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['50%'], ',')].append(each_OTU)

		if OTU_value_DICT[each_OTU] > OTU_threshold_DICT['50%'] and OTU_value_DICT[each_OTU] <= OTU_threshold_DICT['75%']:
			OTU_percentile_DICT[format(OTU_threshold_DICT['50%'], ',') + ' < OTU Abundance <= ' + format(OTU_threshold_DICT['75%'], ',')].append(each_OTU)

		if OTU_value_DICT[each_OTU] > OTU_threshold_DICT['75%']:
			OTU_percentile_DICT[format(OTU_threshold_DICT['75%'], ',') + '< OTU Abundance'].append(each_OTU)

	return (OTU_percentile_DICT, OTU_threshold_DICT)


def multi_sheets_to_single_converter(multi_sheets_shared_file_PATH, taxonomy_level, single_sheets_shared_file_PATH):
	#Fist Load the specified_sheets
	shared_file_data_frame = pandas.read_excel(multi_sheets_shared_file_PATH, sheetname=taxonomy_level)
	shared_file_data_frame.to_csv(single_sheets_shared_file_PATH, sep='\t')
	return True


def multi_shared_file_to_dict_vertical(multi_sheets_shared_file_PATH):
	multi_shared_file_DICT = {}
	multi_shared_file_LIST = []
	lineage_Name_LIST = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	multi_shared_file_LIST = lineage_Name_LIST
	for each_lineage in lineage_Name_LIST:
		shared_file_data_frame = pandas.read_excel(multi_sheets_shared_file_PATH, sheetname=each_lineage)
		any_vertical_LIST = []
		any_vertical_LIST = shared_file_data_frame.columns.values.tolist()
		any_vertical_DICT = {}
		for each_column in any_vertical_LIST:
			each_column_value_list = shared_file_data_frame[each_column].values.tolist()
			any_vertical_DICT[str(each_column)] = map(str, each_column_value_list)
		any_vertical_LIST = map(str, any_vertical_LIST)
		multi_shared_file_DICT[each_lineage] = (any_vertical_DICT, any_vertical_LIST)
	return (multi_shared_file_DICT, multi_shared_file_LIST)
# ################################### UTILITIES ############################## #


def check_it_and_remove_it(filename, noreport=False):
	try:
		os.remove(filename)
		if noreport is False:
			pass
	except OSError:
		pass


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def slugify(target_word):
	text = target_word.replace(' ', r'_').replace('\\', r'_').replace('`', r'_').replace('*', r'_').replace('{', r'_').replace('}', r'_').replace('[', r'_').replace(']', r'_').replace('(', r'_').replace(')', r'_').replace('>', r'_').replace('#', r'_').replace('+', r'_').replace('-', r'_').replace('.', r'_').replace('!', r'_').replace('$', r'_').replace("'", r'_').replace('"', r'_').replace('\/', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_')
	return text


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(map(str, query_list))


def write_string_down(new_string, file_name):
	f = open(file_name, 'w')
	f.write(new_string)
	f.close()
	return True


def sniff_delimiter(template_file_PATH, delimiter_list=None):
	if delimiter_list is None:
		delimiter_list = ['\t']
	template = open(template_file_PATH, 'rU')
	try:
		dialect = csv.Sniffer().sniff(template.readline(), delimiter_list)
	except csv.Error:
		print "The assigned delimiter:", delimiter_list, " is not correct for this file: ", template_file_PATH
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
			# Looking for exactly the same extension
				if ext in ext_list:
					container.append(os.path.join(root, currentFile))
			elif mode == 'multiple':
				for each_ext in ext_list:
					if ext in each_ext:
						container.append(os.path.join(root, currentFile))
			elif mode == 'partial':
				# when this extension is part of the actual extension of files
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


def percentage(part, whole):

	return 100 * float(part) / float(whole)


def zip_it(target_file_NAME, zip_file_NAME):
	zip_file_HANDLE = zipfile.ZipFile(zip_file_NAME, 'w')
	zip_file_HANDLE.write(target_file_NAME)
	zip_file_HANDLE.close()
	return True


def color_lover_generate(target_list, color_scale=None):

	dict_of_symbol_color = {}
	for each_target in target_list:
		dict_of_symbol_color[each_target] = random_color()
	return dict_of_symbol_color


def random_color():
	golden_ratio_conjugate = 0.618033988749895
	randvalue = random.randint(1, 256)
	randvalue += golden_ratio_conjugate
	randvalue %= 1
	#hsv_to_rgb(randvalue, 0.5, 0.95)
	color = "#%06x" % random.randint(1, 0xFFFFFE)
	#print color
	return color


def sort_dictionary_by_value(target_dict, mode=None):
	sorted_dict = {}
	if mode is None:
		sorted_dict = collections.OrderedDict(sorted(target_dict.items(), key=operator.itemgetter(1), reverse=True))
	elif mode is 'Desc':
		sorted_dict = collections.OrderedDict(sorted(target_dict.items(), key=operator.itemgetter(1), reverse=True))
	elif mode is 'Asc':
		sorted_dict = collections.OrderedDict(sorted(target_dict.items(), key=operator.itemgetter(1), reverse=False))
	return sorted_dict.keys()


def remove_extension_files(outputdir, extension):
	extension_list = []
	extension_list.append(extension)
	scanned_container = []
	flag = scandirs(outputdir, scanned_container, extension_list, 'partial')
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


def numpy_percentile(freq_list, threshold):
	numpy_array = numpy.array(freq_list)
	return numpy.percentile(numpy_array, threshold)


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


def remove_line_from_file(the_file, line_number, new_file):
	string = ''
	f = open(the_file, 'rU')
	for idx, i in enumerate(f):
		if idx < line_number:
			continue
		string += i
	write_string_down(string, new_file)
	return True


def median(lst):
	sortedLst = sorted(lst)
	lstLen = len(lst)
	index = (lstLen - 1) // 2

	if (lstLen % 2):
		return sortedLst[index]
	else:
		return (sortedLst[index] + sortedLst[index + 1]) / 2.0


def remove_value_from_list(the_list, the_value):
	new_list = []
	for i in the_list:
		if i == the_value:
			continue
		new_list.append(i)
	return new_list


def numpy_variance(freq_list):
	numpy_array = numpy.array(freq_list)
	return numpy.var(numpy_array, dtype=numpy.float64)


def get_extension(file_PATH):

	return file_PATH.split('.')[-1].lower()


def add_extension(file_PATH, extension):

	os.rename(file_PATH, file_PATH + extension)
	return file_PATH + extension


def get_pandas_DATAFRAME(file_PATH):
	extension = get_extension(file_PATH)
	if extension in ['txt', 'tsv', 'csv']:
		return pandas.read_table(file_PATH, low_memory=False, encoding='utf-8', skip_blank_lines=True)
	elif extension in ['xlsx', 'xls']:
		return pandas.read_excel(file_PATH)
	else:
		print "Unknow extension"
		error("Unknow extension")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)


def get_intersection(list_A, list_B):

	return list(set(list_A).intersection(list_B))


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


def get_longest_length_string_list(string_LIST):

	return len(max(string_LIST, key=len))


def percentage_value(whole, threshold):
	if whole == 0:
		return 0
	return float(threshold) * float(whole) / 100.0
# ################################### SPECIFIC ############################### #


def mothur_shared_parser(shared_file):
	rarefile = open(shared_file, 'rU')
	header_list = []
	header_list = rarefile.readline().rstrip().split('\t')
	#print header_list
	rarefile.close()
	rarefile = open(shared_file, 'rU')
	reader = csv.DictReader(rarefile, header_list, delimiter="\t")
	columns_dict = {}
	for row in reader:
		for key, value in row.items():
			if key == value:
				continue
			elif key in columns_dict:
				columns_dict[key].append(value)
			else:
				columns_dict[key] = [value]
	return (header_list, columns_dict)


def match_shared_with_design(shared_file_PATH, design_file_PATH, control_file_PATH, new_shared_file_PATH, new_design_file_PATH, prefix_NAME, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	# we update shared file based on this control file
	f = open(control_file_PATH, 'rU')
	control_STRING = ''
	control_data = []
	for i in f:
		i = i.rstrip()
		control_data.append(i)
		control_STRING += i + '-'
	control_STRING = control_STRING[:-1]
	f.close()

	flag, stderr = execute_functions(mothur_get_groups, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, control_STRING)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container_list = []
		extension_list = ['.pick' + extension]
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container_list[0], new_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(new_shared_file_PATH)
	# ####################MATCH DESIGN FILE
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_shared_parser(new_shared_file_PATH)

	sample_name_list = columns_dict['Groups']
	new_design_string = 'SAMPLE_ID\tTREATMENT\n'
	removed_sample = []
	file_HANDLE = open(design_file_PATH, 'rU')
	for each_line in file_HANDLE:
		line = each_line.split('\t')
		if line[0] in sample_name_list:
			new_design_string += line[0] + '\t' + line[1].rstrip() + '\n'
		else:
			removed_sample.append(line[0])
	print "The following samples has been removed from your design file"
	print removed_sample
	output_HANDLE = open(new_design_file_PATH, 'w')
	output_HANDLE.write(new_design_string)
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


def sort_shared_file_function_old(shared_file_PATH, sorted_shared_file_PATH, prefix_NAME):
	shared_file_HANDLE = open(shared_file_PATH, 'rU')
	shared_file_HEADER_list = []
	shared_file_HEADER_list = shared_file_HANDLE.readline().rstrip().split('\t')
	shared_file_OBJECT = csv.DictReader(shared_file_HANDLE, shared_file_HEADER_list, delimiter='\t')
	shared_file_OBJECT_dict = {}
	for each_OBJECT in shared_file_OBJECT:
		for key, value in each_OBJECT.items():
			if key == value:
				continue
			elif key in shared_file_OBJECT_dict:
				shared_file_OBJECT_dict[key].append(value)
			else:
				shared_file_OBJECT_dict[key] = [value]
	shared_file_SUMMED_dict = {}
	for each_HEAD in shared_file_HEADER_list:
		if each_HEAD.lower() in ['label', 'numotus', 'groups']:
			continue
		else:
			each_HEAD_TOTAL_value = sum(map(float, shared_file_OBJECT_dict[each_HEAD]))
			shared_file_SUMMED_dict[each_HEAD] = each_HEAD_TOTAL_value
	shared_file_SORTED_HEADER_list = sort_dictionary_by_value(shared_file_SUMMED_dict, mode='Desc')
	shared_file_NEW_HEADER_list = ['label']
	shared_file_NEW_HEADER_list.append('Groups')

	shared_file_NEW_HEADER_list.append('numOtus')

	shared_file_NEW_HEADER_list.extend(shared_file_SORTED_HEADER_list)
	if len(shared_file_NEW_HEADER_list) > 100:
		shared_file_NEW_HEADER_list = shared_file_NEW_HEADER_list[:100]
	length_corrected_HEADER_list = []
	for each_OTU_name in shared_file_NEW_HEADER_list:
		length_corrected_HEADER_list.append(OTU_name_reduction(each_OTU_name))

	#print shared_file_NEW_HEADER_list
	flag = shared_dict_to_file_converter(sorted_shared_file_PATH, shared_file_NEW_HEADER_list, length_corrected_HEADER_list, shared_file_OBJECT_dict, prefix_NAME)
	if flag is True:
		pass
	return True


def sort_shared_file_function(shared_file_PATH, design_file_PATH, sorted_shared_file_PATH, sorted_design_file_PATH):
	shared_file_DICT_vertical_DICT, shared_file_DICT_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	shared_file_SUMMED_dict = {}
	for each_HEAD in shared_file_DICT_vertical_LIST:
		if each_HEAD.lower() in ['label', 'numotus', 'groups']:
			continue
		else:
			each_HEAD_TOTAL_value = sum(map(float, shared_file_DICT_vertical_DICT[each_HEAD]))
			shared_file_SUMMED_dict[each_HEAD] = each_HEAD_TOTAL_value
	sorted_OTU_header_LIST = sort_dictionary_by_value(shared_file_SUMMED_dict, mode='Desc')
	
	shared_file_string = ''
	shared_file_string += 'label' + '\t' + 'Groups' + '\t' + 'numOtus' + '\t' + list_to_string(sorted_OTU_header_LIST, '\t') + '\n'
	numOtus_value = str(len(sorted_OTU_header_LIST))
	
	for each_sample in shared_file_DICT_vertical_DICT['Groups']:
		shared_file_string += 'BIOM_SLICER' + '\t' + each_sample + '\t' + numOtus_value + '\t'
		each_sample_index = shared_file_DICT_vertical_DICT['Groups'].index(each_sample)
		for each_OTU in sorted_OTU_header_LIST:
			shared_file_string += shared_file_DICT_vertical_DICT[each_OTU][each_sample_index] + '\t'
		shared_file_string = shared_file_string[:-1]
		shared_file_string += '\n'

	write_string_down(shared_file_string, sorted_shared_file_PATH)
	copy_file(design_file_PATH, sorted_design_file_PATH)
	return True


def OTU_name_reduction(OTU_name):

	new_name = ''
	if len(OTU_name) > 60:
		if ';' in OTU_name:
			OTU_list = OTU_name.split(';')
		new_name += OTU_list[-3] + ';' + OTU_list[-2] + ';' + OTU_list[-1]
		return new_name
	else:
		return OTU_name


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


def shared_file_group_name_fixer(shared_file_PATH):
	f = open(shared_file_PATH, 'rU')
	shared_string = ''
	for i in f:
		line = i.split('\t')
		if line[0] == 'label':
			if line[1] == 'Groups':
				shared_string += i
			else:
				shared_string += line[0] + '\t' + 'Groups' + '\t' + '\t'.join(line[2:])
		else:
			shared_string += i
	write_string_down(shared_string, shared_file_PATH)
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


def get_low_abundance_sample_list(shared_file, mothur_exec_path, processors, outputdir):
	# ################First calculating number of sequences in each sample
	flag, stderr = execute_functions(mothur_nseqs, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file)
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
	flag = remove_extension_files(outputdir, '.rabund')
	nseqs_file = scanned_container[0]
	# ################PARSING NSEQS RESULT
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_result_parser(nseqs_file)
	check_it_and_remove_it(nseqs_file)
	integer_nseqs_list = map(float, columns_dict['nseqs'])
	min_threshold = 1
	#max_threshold = 99
	min_percentile_value = numpy_percentile(integer_nseqs_list, min_threshold)
	#max_percentile_value = numpy_percentile(integer_nseqs_list, max_threshold)

	sample_freq_dict = {}
	for each_sample in columns_dict['group']:
		each_sample_index = columns_dict['group'].index(each_sample)
		sample_freq_dict[each_sample] = float(columns_dict['nseqs'][each_sample_index])
	low_abundance_sample_list = []
	for each_sample in sample_freq_dict:
		if sample_freq_dict[each_sample] < min_percentile_value:
			low_abundance_sample_list.append(each_sample)
	return low_abundance_sample_list


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
		flag = scandirs(outputdir, scanned_container, extension_list, 'partial')
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
	#FIXER
	shared_file_group_name_fixer(new_shared_file)
	return True


def get_low_abundance_otu_list(shared_file):
	condensed_otu_dict = condense_shared_file_by_itself(shared_file)
	OTU_freq_list = condensed_otu_dict.values()
	min_threshold = 10
	min_percentile_value = numpy_percentile(OTU_freq_list, min_threshold)
	low_abundance_OTU_list = ['unclassified']
	for each_OTU in condensed_otu_dict:
		if condensed_otu_dict[each_OTU] < min_percentile_value:
			low_abundance_OTU_list.append(each_OTU)
	return low_abundance_OTU_list


def condense_shared_file_by_itself(shared_file):
	# a dict of sum of the count of each otu
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)

	sample_name_list = shared_dict.keys()
	OTU_list = shared_dict[sample_name_list[0]].keys()

	condensed_shared_dict = {}

	for each_OTU in OTU_list:
		OTU_total_count_list = []
		for each_sample in sample_name_list:
			OTU_total_count_list.append(float(shared_dict[each_sample][each_OTU]))
		condensed_shared_dict[each_OTU] = sum(OTU_total_count_list)
	return condensed_shared_dict


def parse_shared_file(shared_file):
	header = []
	columns = {}
	header, columns = mothur_result_parser(shared_file)
	shared_dict = {}
	for idx, each_sample in enumerate(columns['Groups']):
		if columns['Groups'][idx] not in shared_dict:
			shared_dict[each_sample] = {}
	for head in header:
		if head in ['label', 'numOtus', 'Groups']:
			continue
		for idx, each_value in enumerate(columns[head]):
			sample_name = columns['Groups'][idx]
			shared_dict[sample_name][head] = each_value

	return shared_dict


def parse_pcoa_axis_file(axis_file):
	header = []
	columns = {}
	header, columns = mothur_result_parser(axis_file)
	shared_dict = {}
	for idx, each_sample in enumerate(columns['group']):
		if columns['group'][idx] not in shared_dict:
			shared_dict[each_sample] = {}
	for head in header:
		if head in ['label', 'numOtus', 'group']:
			continue
		for idx, each_value in enumerate(columns[head]):
			sample_name = columns['group'][idx]
			shared_dict[sample_name][head] = each_value

	return shared_dict


def remove_OTU_shared_design(shared_file, design_file, otu_control_file, new_shared_file, new_design_file, name, mothur_exec_path, processors, outputdir):
	# First we are reading sample listed in control file in string and list
	f = open(otu_control_file, 'rU')
	otu_control_list = []
	for i in f:
		i = i.rstrip()
		otu_control_list.append(i)
	f.close()

	# #####################################################
	remove_otu_list = list(set(otu_control_list))
	#print "Removing Following OTUs:", remove_otu_list
	header = []
	columns = {}
	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	for head in header:
		if head.lower() in ['groups', 'numotus', 'label']:
			continue
		elif head not in remove_otu_list:
			OTU_dict[head] = columns[head]

	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Groups', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'

	for i in range(len(columns['label'])):

		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Groups'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	#FIXER
	shared_file_group_name_fixer(new_shared_file)
	copy_file(design_file, new_design_file)
	return True


def biomarker_discovery(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	significant_OTUs_dict = {}
	significant_OTUs_list = []
	#kruskal_name = 'kruskal_wallis_' + temp_name
	significant_OTUs_dict['kruskal_wallis'], kw_list = kw_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir)
	significant_OTUs_list.extend(kw_list)

	#////////////////////////////////////
	#lefse_name = 'lefse_' + temp_name
	#significant_OTUs_dict['lefse'], lefse_list = lefse_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir)
	#significant_OTUs_list.extend(lefse_list)

	#////////////////////////////////////
	#indicator_name = 'indicator_' + temp_name
	#significant_OTUs_dict['indicator'], indicator_list = indicator_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir)
	#significant_OTUs_list.extend(indicator_list)

	#////////////////////////////////////
	#metastats_name = 'metastats_' + temp_name
	significant_OTUs_dict['metastats'], metastats_list = metastats_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir)
	significant_OTUs_list.extend(metastats_list)

	if len(significant_OTUs_list) > 0:
		ranked_shared = outputdir + name + '_ranked_shared_file.txt'
		flag = keep_OTUs_in_shared(shared_file, significant_OTUs_list, ranked_shared)
		if flag is False:
			print "Warning"
	else:
		print "No Significant Otus Detected"
		return shared_file

	return ranked_shared


def keep_OTUs_in_shared(shared_file, otu_control_list, new_shared_file):
	removing_otu_list = list(set(otu_control_list))
	#print "Removing Following OTUs:", removing_otu_list
	header = []
	columns = {}

	header, columns = mothur_shared_parser(shared_file)
	OTU_dict = {}
	for head in header:
		if head.lower() in ['groups', 'numotus', 'label']:
			continue
		if head not in removing_otu_list:
			continue
		else:
			OTU_dict[head] = columns[head]

	numotus = str(len(OTU_dict.keys()))
	new_shared_string = ''
	new_shared_string = list_to_string(['label', 'Groups', 'numOtus', ] + OTU_dict.keys())
	new_shared_string += '\n'

	for i in range(len(columns['label'])):

		new_shared_string += columns['label'][i]
		new_shared_string += '\t'
		new_shared_string += columns['Groups'][i]
		new_shared_string += '\t'
		new_shared_string += numotus
		new_shared_string += '\t'
		for each_otu in OTU_dict.keys():
			new_shared_string += str(OTU_dict[each_otu][i])
			new_shared_string += '\t'
		new_shared_string += '\n'
	write_string_down(new_shared_string, new_shared_file)
	#FIXER
	shared_file_group_name_fixer(new_shared_file)
	return True


def kw_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	kw_dict = {}
	kw_list = []
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
			return (kw_dict, kw_list)
	kw_dict = {}
	kw_list = []
	header = []
	columns = {}
	header, columns = mothur_result_parser(kw_container[0])
	this_list = list_string_to_float_converter(columns[header[2]])
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


def metastats_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	metastats_dict = {}
	metastats_list = []
	path, absname, ext = split_file_name(shared_file)
	flag, stderr = execute_functions(mothur_metastats, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file, design_file)
	if flag is False:
		print "Execution of mothur_kruskal_wallis failed!!!"
	else:
		metastats_container = []
		extension_list = ['.metastats.logfile']
		flag = scandirs(outputdir, metastats_container, extension_list, 'partial')
		if flag is False:
			print "This extension is not availble: ", extension_list

		else:
			for each_log in metastats_container:
				check_it_and_remove_it(each_log)

		metastats_container = []
		extension_list = ['.metastats']
		flag = scandirs(outputdir, metastats_container, extension_list, 'partial')
		if flag is False:
			print "This extension is not availble: ", extension_list
			return (metastats_dict, metastats_list)
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
		this_list = list_string_to_float_converter(columns[header[7]])
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


def indicator_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	indicator_dict = {}
	indicator_list = []
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
			return (indicator_dict, indicator_list)
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


def lefse_analyser(shared_file, design_file, name, mothur_exec_path, processors, outputdir):
	lefse_dict = {}
	lefse_list = []
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
			return (lefse_dict, lefse_list)

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


def relative_abundance_shared_design_maker(shared_file, design_file, name, new_shared_file, new_design_file, mothur_exec_path, processors, outputdir):
	path, absname, ext = split_file_name(shared_file)
	#method = 'totalgroup'
	flag, stderr = execute_functions(mothur_relative_abundance, processors, outputdir, 'multi', 'mothur', mothur_exec_path, shared_file)
	if flag is False:
		print "Execution of mothur_normalized failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file)
		scanned_container = []
		extension_list = ['.relabund']
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
	#FIXER
	shared_file_group_name_fixer(new_shared_file)
	# we filter
	f = open(new_shared_file, 'rU')
	group_list = []
	for i in f:
		i = i.rstrip()
		line = i.split('\t')
		if line[1].lower() == 'groups':
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


def condense_shared_file_by_design_old(shared_file, design_file):
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)

	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file)
	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)

	reverse_design_dict = design_MATRIX_forward[selected_design]
	
	sample_name_list = shared_dict.keys()
	OTU_list = shared_dict[sample_name_list[0]].keys()

	condensed_shared_dict = {}
	for each_design in reverse_design_dict:
		condensed_shared_dict[each_design] = {}
		sample_list = reverse_design_dict[each_design]
		for each_OTU in OTU_list:
			OTU_total_count_list = []
			for each_sample in sample_list:
				OTU_total_count_list.append(float(shared_dict[each_sample][each_OTU]))
			condensed_shared_dict[each_design][each_OTU] = []
			condensed_shared_dict[each_design][each_OTU] = OTU_total_count_list
	#print condensed_shared_dict[reverse_design_dict.keys()[0]]
	return condensed_shared_dict


def distance_matrix_analyser(shared_file, calc, matrix_name, mothur_exec_path, absname, processors, outputdir):
	#calc = calc.split(':')[0]
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


def design_metadata_unifier(design_file_PATH, sample_metadata_file_PATH, unified_design_file_PATH):

	match_design_flag = True
	sample_vertical_DICT, sample_vertical_LIST = any_file_to_dict_converter_vertical(sample_metadata_file_PATH)
	if design_file_PATH is not None and isFileExist(design_file_PATH) is True:
		design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
		#Step1: Check the validity of design file
		intersection_DICT = {}
		length_LIST = []
		for each_design in design_vertical_LIST:
			intersection_DICT[each_design] = get_intersection(sample_vertical_DICT[sample_vertical_LIST[0]], design_vertical_DICT[each_design])
			length_LIST.append(len(intersection_DICT[each_design]))
		if not any(length_LIST):
			match_design_flag = False
			report("Design file does not match with the Biom file, Ignoring design file")
			print "Design file does not match with the Biom file, Ignoring design file"
	else:

		match_design_flag = False

	unified_design_file_string = ''
	if match_design_flag is True:
		#retreive max target
		max_value = 0
		max_target = ''
		for each_design in design_vertical_LIST:
			if len(intersection_DICT[each_design]) > max_value:
				max_value = len(intersection_DICT[each_design])
				max_target = each_design
				max_target_index = design_vertical_LIST.index(max_target)
		#what is max_target? is the column name from design_file which infered to matched with the first column in sample_metadata
		#LOAD design data
		sample_design_horizontal_DICT = {}
		sample_design_horizontal_header_LIST = []
		sample_design_horizontal_DICT, sample_design_horizontal_header_LIST = any_file_to_dict_converter_horizontal(design_file_PATH, max_target_index)
		#print sample_design_horizontal_DICT
		#Load metadata data
		sample_metadata_horizontal_DICT = {}
		sample_metadata_horizontal_header_LIST = []
		sample_metadata_horizontal_DICT, sample_metadata_horizontal_header_LIST = any_file_to_dict_converter_horizontal(sample_metadata_file_PATH)
		#print sample_metadata_horizontal_DICT
		#print design_vertical_LIST
		design_vertical_LIST.pop(max_target_index)
		#print sample_vertical_LIST
		sample_header_list = []
		sample_header_list.append(sample_vertical_LIST[1])
		sample_header_list.append(sample_vertical_LIST[0])
		sample_header_list.extend(sample_vertical_LIST[2:])
		unified_design_file_string += list_to_string(sample_header_list, '\t') + '\t' + list_to_string(design_vertical_LIST, '\t') + '\n'
		for each_sample in sample_design_horizontal_header_LIST:
			sample_design_row = sample_design_horizontal_DICT[each_sample]
			sample_design_row.pop(max_target_index)
			sample_value_list = []
			sample_value_list.append(sample_metadata_horizontal_DICT[each_sample][1])
			sample_value_list.append(sample_metadata_horizontal_DICT[each_sample][0])
			sample_value_list.extend(sample_metadata_horizontal_DICT[each_sample][2:])
			#sample_metadata_row = sample_metadata_horizontal_DICT[each_sample][1:]
			unified_design_file_string += list_to_string(sample_value_list, '\t') + '\t' + list_to_string(sample_design_row, '\t') + '\n'
	elif match_design_flag is False:
		sample_metadata_horizontal_DICT = {}
		sample_metadata_horizontal_header_LIST = []
		sample_metadata_horizontal_DICT, sample_metadata_horizontal_header_LIST = any_file_to_dict_converter_horizontal(sample_metadata_file_PATH)
		unified_header_list = []
		unified_header_list.append(sample_vertical_LIST[1])
		unified_header_list.append(sample_vertical_LIST[0])
		unified_header_list.extend(sample_vertical_LIST[2:])
		unified_design_file_string += list_to_string(unified_header_list, '\t') + '\n'
		for each_sample in sample_metadata_horizontal_header_LIST:
			unified_value_list = []
			unified_value_list.append(sample_metadata_horizontal_DICT[each_sample][1])
			unified_value_list.append(sample_metadata_horizontal_DICT[each_sample][0])
			unified_value_list.extend(sample_metadata_horizontal_DICT[each_sample][2:])
			#sample_metadata_row = sample_metadata_horizontal_DICT[each_sample][1:]
			unified_design_file_string += list_to_string(unified_value_list, '\t') + '\n'
	write_string_down(unified_design_file_string, unified_design_file_PATH)
	#print unified_design_file_string
	return True


def update_shared_file_using_design_function(shared_file_PATH, design_file_PATH, updated_shared_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	control_string = '-'.join(design_horizontal_LIST)

	flag, stderr = execute_functions(mothur_get_groups, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, control_string)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container_list = []
		extension_list = ['.pick' + extension]
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container_list[0], updated_shared_file_PATH)
	shared_file_group_name_fixer(updated_shared_file_PATH)
	return True


def match_shared_design_function(shared_file_PATH, design_file_PATH, matched_shared_file_PATH, matched_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#1 match shared_file with design files
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	control_string = '-'.join(design_horizontal_LIST)
	flag, stderr = execute_functions(mothur_get_groups, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, control_string)
	if flag is False:
		print "Execution of mothur_get_groups failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container_list = []
		extension_list = ['.pick' + extension]
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
		print "Scanning.."
		if flag is False:
			print "Failed :("
			print "This extension is not available: ", extension_list
			sys.exit(2)
		else:
			print "NICE :) Found it"
			counter = 1
			for file in scanned_container_list:
				print "File#", str(counter), ":", file
				counter += 1
	os.rename(scanned_container_list[0], matched_shared_file_PATH)
	shared_file_group_name_fixer(matched_shared_file_PATH)

	#2 matche design file with shared file path
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(matched_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, matched_design_file_PATH)
	return True


def taxonomy_inference_engine(raw_taxonomy_file_PATH, processed_taxonomy_file_PATH):
	
	return 1


def core_microbiome_function(shared_file_PATH, design_file_PATH, coremicrobiomed_shared_file_PATH, coremicrobiomed_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#RUN GET CORE MICROBIOME
	min_abundance_ratio = '1'
	flag, stderr = execute_functions(mothur_get_coremicrobiome, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, min_abundance_ratio)
	if flag is False:
		print "Execution of mothur_get_coremicrobiome failed!!!"
	path, absname, extension = split_file_name(shared_file_PATH)
	scanned_container_list = []
	extension_list = ['.core.microbiomelist']
	flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
	print "Scanning.."
	if flag is False:
		print "Failed :("
		print "This extension is not available: ", extension_list
		sys.exit(2)
	else:
		print "NICE :) Found it"
		counter = 1
		for file in scanned_container_list:
			print "File#", str(counter), ":", file
			counter += 1
	#PARSE RESULT FILE
	core_microbiome_file_PATH = add_extension(scanned_container_list[0], '.txt')
	core_microbiome_horizontal_DICT, core_microbiome_horizontal_LIST = any_file_to_dict_converter_horizontal(core_microbiome_file_PATH)
	Selected_OTU_LIST = []
	Selected_OTU_file_PATH = outputdir_PATH + 'BIOM_SLICER_SELECTED_OTU_file.txt'
	Selected_OTU_LIST = core_microbiome_horizontal_DICT['1'][1].split(',')
	Selected_OTU_string = list_to_string(Selected_OTU_LIST, '\n')
	write_string_down(Selected_OTU_string, Selected_OTU_file_PATH)
	#Filter in shared_file with selected otu
	flag, stderr = execute_functions(mothur_get_otus, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, Selected_OTU_file_PATH)
	if flag is False:
		print "Execution of mothur_get_coremicrobiome failed!!!"
	path, absname, extension = split_file_name(shared_file_PATH)
	scanned_container_list = []
	extension_list = ['.pick' + extension]
	flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
	print "Scanning.."
	if flag is False:
		print "Failed :("
		print "This extension is not available: ", extension_list
		sys.exit(2)
	else:
		print "NICE :) Found it"
		counter = 1
		for file in scanned_container_list:
			print "File#", str(counter), ":", file
			counter += 1
	os.rename(scanned_container_list[0], coremicrobiomed_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(coremicrobiomed_shared_file_PATH)
	#MATCH DESIGN FILE
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(coremicrobiomed_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, coremicrobiomed_design_file_PATH)
	return True


def which_design(design_MATRIX_forward):
	unique_LIST = []
	selected_design = ''
	for each_design_group in design_MATRIX_forward.keys():
		unique_LIST = list(set(design_MATRIX_forward[each_design_group]))
		if len(unique_LIST) < 6 and len(unique_LIST) > 1:
			selected_design = each_design_group
	if selected_design == '':
		return 'DEFAULT_DESIGN'
	else:
		return selected_design


def sample_abundance_scanning_function(shared_file_PATH, design_file_PATH, abundance_scanned_shared_file_PATH, abundance_scanned_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#Check the abundance of each sample
	shared_file_horizontal_DICT, shared_file_horizontal_LIST = any_file_to_dict_converter_horizontal(shared_file_PATH, 1)
	total_Abundance_DICT = {}
	for each_sample in shared_file_horizontal_DICT:
		total_Abundance_DICT[each_sample] = sum(map(int, shared_file_horizontal_DICT[each_sample][3:]))
	#print total_Abundance_DICT
	min_threshold = 10
	max_threshold = 90
	min_percentile_value = numpy_percentile(total_Abundance_DICT.values(), min_threshold)
	max_percentile_value = numpy_percentile(total_Abundance_DICT.values(), max_threshold)
	low_abundance_sample_LIST = []
	high_abundance_sample_LIST = []
	normal_abundance_sample_LIST = []
	for each_sample in total_Abundance_DICT:
		if total_Abundance_DICT[each_sample] < min_percentile_value:
			low_abundance_sample_LIST.append(each_sample)
		elif total_Abundance_DICT[each_sample] > max_percentile_value:
			high_abundance_sample_LIST.append(each_sample)
		else:
			normal_abundance_sample_LIST.append(each_sample)
	removing_samples_LIST = []
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)

	if len(low_abundance_sample_LIST) == 0:
		print "NO LOW ABUNDANCE SAMPLES DETECTED"
	else:
		print "LOW ABUNDANCE SAMPLES: "
		for each_sample in low_abundance_sample_LIST:
			print design_MATRIX_reverse[each_sample]['SAMPLE_ID']
		removing_samples_LIST.extend(low_abundance_sample_LIST)
	if len(high_abundance_sample_LIST) == 0:
		print "NO HIGH ABUNDANCE SAMPLES DETECTED"
	else:
		print "HIGH ABUNDANCE SAMPLES: "
		for each_sample in high_abundance_sample_LIST:
			print design_MATRIX_reverse[each_sample]['SAMPLE_ID']
		removing_samples_LIST.extend(high_abundance_sample_LIST)
	if len(removing_samples_LIST) == 0:
		print "No abnormally abundant sample detected"
		copy_file(shared_file_PATH, abundance_scanned_shared_file_PATH)
		copy_file(design_file_PATH, abundance_scanned_design_file_PATH)
	else:
		removing_samples_string = ''
		removing_samples_string = list_to_string(removing_samples_LIST, '-')
		flag, stderr = execute_functions(mothur_remove_groups, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, removing_samples_string)
		if flag is False:
			print "Execution of mothur_get_groups failed!!!"
		else:
			path, absname, extension = split_file_name(shared_file_PATH)
			scanned_container = []
			extension_list = ['.pick' + extension]
			flag = scandirs(outputdir_PATH, scanned_container, extension_list, 'partial')
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
		os.rename(scanned_container[0], abundance_scanned_shared_file_PATH)
		#FIX GROUP NAME
		shared_file_group_name_fixer(abundance_scanned_shared_file_PATH)
		#MATCH DESIGN FILE
		design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
		shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(abundance_scanned_shared_file_PATH)
		design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
		
		matched_design_file_string = ''
		matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

		for each_sample in shared_file_vertical_DICT['Groups']:
			if each_sample in design_horizontal_LIST:
				matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
		write_string_down(matched_design_file_string, abundance_scanned_design_file_PATH)

	return True


def relative_abundance_function(shared_file_PATH, design_file_PATH, relative_abundance_shared_file_PATH, relative_abundance_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#1 Relative abundance calculation
	flag, stderr = execute_functions(mothur_relative_abundance, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH)
	if flag is False:
		print "Execution of mothur_normalized failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container = []
		extension_list = ['.relabund']
		flag = scandirs(outputdir_PATH, scanned_container, extension_list)
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
	os.rename(scanned_container[0], relative_abundance_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(relative_abundance_shared_file_PATH)
	#Match shared file
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(relative_abundance_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, relative_abundance_design_file_PATH)
	return True


def normalize_abundance_function(shared_file_PATH, design_file_PATH, normalize_abundance_shared_file_PATH, normalize_abundance_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#1 Normalize
	normalize_method = 'totalgroup'
	flag, stderr = execute_functions(mothur_normalize, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, normalize_method)
	if flag is False:
		print "Execution of mothur_normalized failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container = []
		extension_list = ['.norm.shared']
		flag = scandirs(outputdir_PATH, scanned_container, extension_list, 'partial')
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
	os.rename(scanned_container[0], normalize_abundance_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(normalize_abundance_shared_file_PATH)
	#Match shared file
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(normalize_abundance_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, normalize_abundance_design_file_PATH)
	return True


def OTU_abundance_scanning_function(shared_file_PATH, design_file_PATH, OTU_abundance_scanned_shared_file_PATH, OTU_abundance_scanned_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#1 heuristic method to find the nseqs value
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_total_count_LIST = []
	for each_OTU in shared_file_vertical_LIST[3:]:
		OTU_total_count_LIST.append(sum(map(int, shared_file_vertical_DICT[each_OTU])))
	min_threshold = 10
	#print OTU_total_count_LIST
	rare_value = numpy_percentile(OTU_total_count_LIST, min_threshold)
	#sys.exit(2)
	#2 remove rare
	nseqs = str(int(rare_value))
	flag, stderr = execute_functions(mothur_remove_rare, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, nseqs)
	if flag is False:
		print "Execution of mothur_normalized failed!!!"
	else:
		path, absname, extension = split_file_name(shared_file_PATH)
		scanned_container = []
		extension_list = ['.pick' + extension]
		flag = scandirs(outputdir_PATH, scanned_container, extension_list, 'partial')
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
	os.rename(scanned_container[0], OTU_abundance_scanned_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(OTU_abundance_scanned_shared_file_PATH)
	#Match shared file
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(OTU_abundance_scanned_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, OTU_abundance_scanned_design_file_PATH)
	return True


def biomarker_discovery_functions(shared_file_PATH, design_file_PATH, BIOMARKERED_shared_file_PATH, BIOMARKERED_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):
	#DESIGN PROCESS
	selcted_design_file_PATH = outputdir_PATH + 'selected_design_file.txt'
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	design_file_vertical_DICT, design_file_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	selected_design = which_design(design_MATRIX_forward)
	selected_design_string = 'SAMPLE_ALIAS\tDESIGN\n'
	for each_sample in design_file_vertical_DICT['SAMPLE_ALIAS']:
		each_sample_index = design_file_vertical_DICT['SAMPLE_ALIAS'].index(each_sample)
		selected_design_string += each_sample + '\t' + design_file_vertical_DICT[selected_design][each_sample_index] + '\n'
	write_string_down(selected_design_string, selcted_design_file_PATH)

	#Forward: {'phinchID': {'24': ['SAMPLE_ALIAS_025'], '25': ['SAMPLE_ALIAS_026'], '26': ['SAMPLE_ALIAS_027'],
	#Reverse: {'SAMPLE_ALIAS_024': {'phinchID': '23', 'temp': '20.8 C', 'collection_date': '2012-12-12T17:00:00+08:00'
	
	significant_OTUs_DICT = {}
	significant_OTUs_LIST = []
	#KRUSKAL-WALLIS
	significant_OTUs_DICT['kruskal_wallis'], kw_list = kw_analyser(shared_file_PATH, selcted_design_file_PATH, 'BIOM_SLICER', mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	significant_OTUs_LIST.extend(kw_list)
	#METASTATS
	significant_OTUs_DICT['metastats'], metastats_list = metastats_analyser(shared_file_PATH, selcted_design_file_PATH, 'BIOM_SLICER', mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	significant_OTUs_LIST.extend(metastats_list)
	#INDICATOR
	significant_OTUs_DICT['indicator'], indicator_list = indicator_analyser(shared_file_PATH, selcted_design_file_PATH, 'BIOM_SLICER', mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	significant_OTUs_LIST.extend(indicator_list)
	#LEFSE
	significant_OTUs_DICT['lefse'], lefse_list = lefse_analyser(shared_file_PATH, selcted_design_file_PATH, 'BIOM_SLICER', mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	significant_OTUs_LIST.extend(lefse_list)
	#print significant_OTUs_DICT

	#biomarker_table_html_string = biomarker_discovery_html_table_function(significant_OTUs_DICT, OTU_metadata_file_PATH)
	#return biomarker_table_html_string
	biomarker_OTU_LIST = []
	biomarker_OTU_LIST.extend(kw_list)
	biomarker_OTU_LIST.extend(metastats_list)
	biomarker_OTU_LIST.extend(indicator_list)
	biomarker_OTU_LIST.extend(lefse_list)
	biomarker_OTU_LIST = list(set(biomarker_OTU_LIST))
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	OTU_name_LIST = shared_file_vertical_LIST[3:]
	threshold = 10
	otu_percentile = percentage_value(len(OTU_name_LIST), threshold)
	if len(biomarker_OTU_LIST) < otu_percentile:
		print "Biomarker discovery result are too stringent, will skip thus step"
		copy_file(shared_file_PATH, BIOMARKERED_shared_file_PATH)
		copy_file(design_file_PATH, BIOMARKERED_design_file_PATH)
		return True

	#print biomarker_OTU_LIST

	biomarker_OTU_file_PATH = outputdir_PATH + 'BIOM_SLICER_BIOMARKER_OTU_file.txt'
	biomarker_OTU_string = list_to_string(biomarker_OTU_LIST, '\n')
	write_string_down(biomarker_OTU_string, biomarker_OTU_file_PATH)
	#Filter in shared_file with selected otu
	flag, stderr = execute_functions(mothur_get_otus, processors_COUNT, outputdir_PATH, 'multi', 'mothur', mothur_exec_PATH, shared_file_PATH, biomarker_OTU_file_PATH)
	if flag is False:
		print "Execution of mothur_get_coremicrobiome failed!!!"
	path, absname, extension = split_file_name(shared_file_PATH)
	scanned_container_list = []
	extension_list = ['.pick' + extension]
	flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'partial')
	print "Scanning.."
	if flag is False:
		print "Failed :("
		print "This extension is not available: ", extension_list
		sys.exit(2)
	else:
		print "NICE :) Found it"
		counter = 1
		for file in scanned_container_list:
			print "File#", str(counter), ":", file
			counter += 1
	os.rename(scanned_container_list[0], BIOMARKERED_shared_file_PATH)
	#FIXER
	shared_file_group_name_fixer(BIOMARKERED_shared_file_PATH)
	#MATCH DESIGN FILE
	design_horizontal_DICT, design_horizontal_LIST = any_file_to_dict_converter_horizontal(design_file_PATH)
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(BIOMARKERED_shared_file_PATH)
	design_vertical_DICT, design_vertical_LIST = any_file_to_dict_converter_vertical(design_file_PATH)
	
	matched_design_file_string = ''
	matched_design_file_string = list_to_string(design_vertical_LIST, '\t') + '\n'

	for each_sample in shared_file_vertical_DICT['Groups']:
		if each_sample in design_horizontal_LIST:
			matched_design_file_string += list_to_string(design_horizontal_DICT[each_sample], '\t') + '\n'
	write_string_down(matched_design_file_string, BIOMARKERED_design_file_PATH)
	return True


def condense_shared_file_by_design(shared_file_PATH, design_DICT):
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	
	OTU_list = shared_file_vertical_LIST[3:]
	condensed_shared_DICT = {}
	for each_design in design_DICT:
		condensed_shared_DICT[each_design] = {}
		sample_list = design_DICT[each_design]
		for each_OTU in OTU_list:
			OTU_total_count_list = []
			for each_sample in sample_list:
				each_sample_index = shared_file_vertical_DICT['Groups'].index(each_sample)
				OTU_total_count_list.append(float(shared_file_vertical_DICT[each_OTU][each_sample_index]))
			condensed_shared_DICT[each_design][each_OTU] = sum(OTU_total_count_list)
	return condensed_shared_DICT


def condense_shared_file_by_taxonomy(shared_file_PATH, OTU_metadata_file_PATH):
	#STEP1: LOADING TAXONOMY DICT
	lineage_LIST = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	taxonomy_DICT = {}
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	#print OTU_MATRIX_forward['taxonomy'].keys()
	taxonomy_LIST = []

	for each_lineage in lineage_LIST:
		lineage_index = lineage_LIST.index(each_lineage)
		taxonomy_LIST = []
		for each_taxonomy in OTU_MATRIX_forward['taxonomy'].keys():
			each_taxonomy_LIST = each_taxonomy.split(';')
			taxonomy_LIST.append(each_taxonomy_LIST[lineage_index])
		taxonomy_LIST = list(set(taxonomy_LIST))
		taxonomy_DICT[each_lineage] = taxonomy_LIST
	#STEP2: PARSE EACH TAXONOMY in sample
	OTU_taxonomic_DICT = {}
	for each_lineage in lineage_LIST:
		each_lineage_index = lineage_LIST.index(each_lineage)
		OTU_taxonomic_DICT[each_lineage] = {}
		for each_total_taxonomy in OTU_MATRIX_forward['taxonomy']:
			each_total_taxonomy_LIST = each_total_taxonomy.split(';')
			if each_total_taxonomy_LIST[each_lineage_index] in taxonomy_DICT[each_lineage]:
				OTU_taxonomic_DICT[each_lineage][each_total_taxonomy_LIST[each_lineage_index]] = OTU_MATRIX_forward['taxonomy'][each_total_taxonomy]
	#print OTU_taxonomic_DICT['Phylum']['p_TM7']
	#print OTU_taxonomic_DICT['Family']['f_F16']
	#'p_TM7': ['OTU_ALIAS_032'], 'p_Deferribacteres': ['OTU_ALIAS_047', 'OTU_ALIAS_116'], 'p_Proteobacteria': ['OTU_ALIAS_305'], '*k_Bacteria': ['OTU_ALIAS_025', 'OTU_ALIAS_037', 'OTU_ALIAS_061', 'OTU_ALIAS_110', 'OTU_ALIAS_119',
	#STEP3: SUMS UP
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	lineage_DICT = {}

	for each_lineage in lineage_LIST:
		lineage_DICT[each_lineage] = {}
		for each_entity in OTU_taxonomic_DICT[each_lineage]:
			list_of_lists = []
			if len(OTU_taxonomic_DICT[each_lineage][each_entity]) == 1:
				lineage_DICT[each_lineage][each_entity] = shared_file_vertical_DICT[OTU_taxonomic_DICT[each_lineage][each_entity][0]]
			else:
				for each_OTU in OTU_taxonomic_DICT[each_lineage][each_entity]:
					list_of_lists.append(map(int, shared_file_vertical_DICT[each_OTU]))
				lineage_DICT[each_lineage][each_entity] = map(str, [sum(x) for x in zip(*list_of_lists)])
	print lineage_DICT


def multi_sheets_shared_file_function(shared_file_PATH, OTU_metadata_file_PATH, multi_sheets_shared_file_PATH, multi_sheets_OTU_metadata_file_PATH):
	#STEP1: create new OTU metadata file path including lineage information
	lineage_Name_LIST = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	lineage_LIST = []
	lineage_DICT = {}
	for each_lineage in lineage_Name_LIST:
		lineage_index = lineage_Name_LIST.index(each_lineage)
		lineage_LIST = []
		for each_taxonomy in OTU_MATRIX_forward['taxonomy'].keys():
			each_taxonomy_LIST = each_taxonomy.split(';')
			lineage_LIST.append(each_taxonomy_LIST[lineage_index])
		lineage_LIST = list(set(lineage_LIST))
		lineage_DICT[each_lineage] = lineage_LIST
	#print lineage_DICT
	#'Kingdom': ['k_Bacteria'], 'Family': ['f_Bacteroidaceae', '*o_RF39', 'f_Streptococcaceae', 'f_Ruminococcaceae', 'f_Rikenellaceae', 'f_Helicobacteraceae', 'f_F16', '*o_Bacteroidales', 'f_Clostridiales_Family_XIII_Incertae_Sedis'

	#STEP2: create multi_sheets_OTU_metadata_file
	OTU_ALIAS_taxonomic_name_DICT = {}
	lineage_info_string = 'OTU_ALIAS\ttaxonomy_level\ttaxonomy\n'
	OTU_counter = 1
	for each_lineage in lineage_Name_LIST:
		OTU_ALIAS_taxonomic_name_DICT[each_lineage] = {}
		for each_taxonomy in lineage_DICT[each_lineage]:
			OTU_ALIAS_name = 'OTU_ALIAS_' + str(OTU_counter).zfill(3)
			lineage_info_string += OTU_ALIAS_name + '\t' + each_lineage + '\t' + each_taxonomy + '\n'
			OTU_ALIAS_taxonomic_name_DICT[each_lineage][OTU_ALIAS_name] = each_taxonomy
			OTU_counter += 1
	write_string_down(lineage_info_string, multi_sheets_OTU_metadata_file_PATH)
	
	#STEP3: Create taxonomy_DICT of each OTU
	OTU_taxonomic_DICT = {}
	for each_lineage in lineage_Name_LIST:
		each_lineage_index = lineage_Name_LIST.index(each_lineage)
		OTU_taxonomic_DICT[each_lineage] = {}
		
		for each_total_taxonomy in OTU_MATRIX_forward['taxonomy']:
			each_total_taxonomy_LIST = each_total_taxonomy.split(';')
			
			if each_total_taxonomy_LIST[each_lineage_index] not in OTU_taxonomic_DICT[each_lineage]:
				OTU_taxonomic_DICT[each_lineage][each_total_taxonomy_LIST[each_lineage_index]] = OTU_MATRIX_forward['taxonomy'][each_total_taxonomy]
			else:
				OTU_taxonomic_DICT[each_lineage][each_total_taxonomy_LIST[each_lineage_index]].extend(OTU_MATRIX_forward['taxonomy'][each_total_taxonomy])
			OTU_taxonomic_DICT[each_lineage][each_total_taxonomy_LIST[each_lineage_index]] = list(set(OTU_taxonomic_DICT[each_lineage][each_total_taxonomy_LIST[each_lineage_index]]))

	#print len(OTU_taxonomic_DICT['Kingdom']['k_Bacteria'])
	#'p_TM7': ['OTU_ALIAS_032'], 'p_Deferribacteres': ['OTU_ALIAS_047', 'OTU_ALIAS_116'], 'p_Proteobacteria': ['OTU_ALIAS_305'], '*k_Bacteria': ['OTU_ALIAS_025', 'OTU_ALIAS_037', 'OTU_ALIAS_061', 'OTU_ALIAS_110', 'OTU_ALIAS_119',
	
	#STEP4: COLLAPSE OTU VALUE OF EACH SAMPLE
	OTU_collapsed_value_DICT = {}
	shared_file_vertical_DICT, shared_file_vertical_LIST = any_file_to_dict_converter_vertical(shared_file_PATH)
	for each_lineage in lineage_Name_LIST:
		OTU_collapsed_value_DICT[each_lineage] = {}
		for each_entity in OTU_taxonomic_DICT[each_lineage]:
			list_of_lists = []
			if len(OTU_taxonomic_DICT[each_lineage][each_entity]) == 1:
				OTU_collapsed_value_DICT[each_lineage][each_entity] = shared_file_vertical_DICT[OTU_taxonomic_DICT[each_lineage][each_entity][0]]
			else:
				for each_OTU in OTU_taxonomic_DICT[each_lineage][each_entity]:
					list_of_lists.append(map(int, shared_file_vertical_DICT[each_OTU]))
				OTU_collapsed_value_DICT[each_lineage][each_entity] = map(str, [sum(x) for x in zip(*list_of_lists)])
	#print OTU_collapsed_value_DICT
	#{'Kingdom': {'k_Bacteria': ['2', '1', '0', '5', '2', '22', '0', '23', '4']}, 'Family': {'f_Bacteroidaceae': ['25',

	#STEP4: CREATE MULTI SHEETS EXCEL FILE
	OTU_value_DICT = {}
	OTU_Name_LIST = ['label', 'Groups', 'numOtus']
	excel_writer = pandas.ExcelWriter(multi_sheets_shared_file_PATH)
	for each_lineage in lineage_Name_LIST:
		OTU_value_DICT = {}
		OTU_Name_LIST = ['label', 'Groups', 'numOtus']
		numOtus_value = len(OTU_collapsed_value_DICT[each_lineage].keys())
		OTU_Name_LIST.extend(OTU_collapsed_value_DICT[each_lineage].keys())
		OTU_value_DICT['label'] = shared_file_vertical_DICT['label']
		OTU_value_DICT['Groups'] = shared_file_vertical_DICT['Groups']
		OTU_value_DICT['numOtus'] = [numOtus_value] * len(shared_file_vertical_DICT['Groups'])
		OTU_value_DICT.update(OTU_collapsed_value_DICT[each_lineage])
		lineage_data_frame = pandas.DataFrame.from_dict(OTU_value_DICT)
		
		lineage_data_frame.to_excel(excel_writer, sheet_name=each_lineage, columns=OTU_Name_LIST, encoding='string', index=False)
	excel_writer.save()
	sys.exit(2)
	return True


def complexity_reduction_factory(shared_file_PATH, design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH):

	#STEP1: SAMPLE ABUNDANCE SCANNING
	# ++++++++++++++++++++++++++++++ SAMPLE ABUNDANCE SCANNING FUNCTION
	print "SAMPLE ABUNDANCE SCANNING FUNCTION is in progress"
	report("SAMPLE ABUNDANCE SCANNING FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	sample_abundance_scanned_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_SAMPLE_ABUNDANCE_SCANNED_SHARED_file_STEP2.txt'
	sample_abundance_scanned_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_SAMPLE_ABUNDANCE_SCANNED_DESIGN_file_STEP2.txt'
	flag = sample_abundance_scanning_function(shared_file_PATH, design_file_PATH, sample_abundance_scanned_shared_file_PATH, sample_abundance_scanned_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	if flag is True:
		pass
	else:
		print "ABORT"
		sys.exit(2)
	print "SAMPLE ABUNDANCE SCANNING FUNCTION PASSED!!!"
	report("SAMPLE ABUNDANCE SCANNING FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF ABUNDANCE SCANNING FUNCTION

	#STEP2: OTU ABUNDANCE SCANNING
	# ++++++++++++++++++++++++++++++ OTU ABUNDANCE SCANNING FUNCTION
	print "OTU ABUNDANCE SCANNING FUNCTION is in progress"
	report("OTU ABUNDANCE SCANNING FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	OTU_abundance_scanned_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_OTU_ABUNDANCE_SCANNED_SHARED_file_STEP2.txt'
	OTU_abundance_scanned_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_OTU_ABUNDANCE_SCANNED_DESIGN_file_STEP2.txt'
	flag = OTU_abundance_scanning_function(sample_abundance_scanned_shared_file_PATH, sample_abundance_scanned_design_file_PATH, OTU_abundance_scanned_shared_file_PATH, OTU_abundance_scanned_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	if flag is True:
		pass
	else:
		print "ABORT"
		sys.exit(2)
	print "OTU ABUNDANCE SCANNING FUNCTION PASSED!!!"
	report("OTU ABUNDANCE SCANNING FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF ABUNDANCE SCANNING FUNCTION

	#STEP2: COREMICROBIOME CALCULATION
	# ++++++++++++++++++++++++++++++ CORE MICROBIOME FUNCTION
	print "CORE MICROBIOME FUNCTION is in progress"
	report("CORE MICROBIOME FUNCTION file is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	coremicrobiomed_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_COREMICROBIOMED_shared_file_STEP2.txt'
	coremicrobiomed_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_COREMICROBIOMED_design_file_STEP2.txt'
	flag = core_microbiome_function(OTU_abundance_scanned_shared_file_PATH, OTU_abundance_scanned_design_file_PATH, coremicrobiomed_shared_file_PATH, coremicrobiomed_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	if flag is True:
		pass
	else:
		print "ABORT"
		sys.exit(2)
	print "CORE MICROBIOME FUNCTION PASSED!!!"
	report("CORE MICROBIOME FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF CORE MICROBIOME FUNCTION

	#STEP4: BIOMARKER DISCOVERY
	# +++++++++++++++++++++++++++++ BIOMARKER DISCOVERY
	print "BIOMARKER DISCOVERY FUNCTION is in progress"
	report("BIOMARKER DISCOVERY file is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	biomarkered_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_BIOMARKERED_shared_file_STEP2.txt'
	biomarkered_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_BIOMARKERED_design_file_STEP2.txt'

	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(design_file_PATH)
	selected_design = ''
	selected_design = which_design(design_MATRIX_forward)
	if selected_design != 'DEFAULT_DESIGN':
		flag = biomarker_discovery_functions(coremicrobiomed_shared_file_PATH, coremicrobiomed_design_file_PATH, biomarkered_shared_file_PATH, biomarkered_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
		if flag is True:
			pass
		else:
			print "ABORT"
			sys.exit(2)
	else:
		copy_file(coremicrobiomed_shared_file_PATH, biomarkered_shared_file_PATH)
		copy_file(coremicrobiomed_design_file_PATH, biomarkered_design_file_PATH)
	print "BIOMARKER DISCOVERY FUNCTION PASSED!!!"
	report("BIOMARKER DISCOVERY FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF BIOMARKER DISCOVERY

	#STEP3: ASSIGNING NATURAL ABUNDANCE FILES
	# +++++++++++++++++++++++++++++ ASSIGNING NATURAL ABUNDANCE FILES
	NATURAL_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_NATURAL_ABUNDANCE_shared_file_STEP2.txt'
	NATURAL_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_NATURAL_ABUNDANCE_design_file_STEP2.txt'
	copy_file(biomarkered_shared_file_PATH, NATURAL_abundance_shared_file_PATH)
	copy_file(biomarkered_design_file_PATH, NATURAL_abundance_design_file_PATH)
	# ----------------------------- END OF ASSIGNING NATURAL ABUNDANCE FILES

	#STEP3: NORMALIZE ABUNDANCE
	# +++++++++++++++++++++++++++++ NORMALIZE ABUNDANCE FUNCTION
	print "NORMALIZE ABUNDANCE FUNCTION is in progress"
	report("NORMALIZE ABUNDANCE FUNCTION file is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	NORMALIZED_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_NORMALIZED_ABUNDANCE_shared_file_STEP2.txt'
	NORMALIZED_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_NORMALIZED_ABUNDANCE_design_file_STEP2.txt'
	flag = normalize_abundance_function(NATURAL_abundance_shared_file_PATH, NATURAL_abundance_design_file_PATH, NORMALIZED_abundance_shared_file_PATH, NORMALIZED_abundance_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	if flag is True:
		pass
	else:
		print "ABORT"
		sys.exit(2)
	print "RELATIVE ABUNDANCE FUNCTION PASSED!!!"
	report("RELATIVE ABUNDANCE FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF RELATIVE ABUNDANCE FUNCTION

	#STEP5: RELATIVE ABUNDANCE
	# +++++++++++++++++++++++++++++ RELATIVE ABUNDANCE FUNCTION
	print "RELATIVE ABUNDANCE FUNCTION is in progress"
	report("RELATIVE ABUNDANCE FUNCTION file is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	RELATIVE_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_RELATIVE_ABUNDANCE_shared_file_STEP2.txt'
	RELATIVE_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_RELATIVE_ABUNDANCE_design_file_STEP2.txt'
	flag = relative_abundance_function(NATURAL_abundance_shared_file_PATH, NATURAL_abundance_design_file_PATH, RELATIVE_abundance_shared_file_PATH, RELATIVE_abundance_design_file_PATH, mothur_exec_PATH, processors_COUNT, outputdir_PATH)
	if flag is True:
		pass
	else:
		print "ABORT"
		sys.exit(2)
	print "RELATIVE ABUNDANCE FUNCTION PASSED!!!"
	report("RELATIVE ABUNDANCE FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF RELATIVE ABUNDANCE FUNCTION

	# ++++++++++++++++++++++++++++++ SORT NATURAL ABUNDANCE SHARED FILE
	NATURAL_sorted_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_NATURAL_SORTED_shared_file_STEP2.txt'
	NATURAL_sorted_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_NATURAL_SORTED_design_file_STEP2.txt'
	sort_shared_file_function(NATURAL_abundance_shared_file_PATH, NATURAL_abundance_design_file_PATH, NATURAL_sorted_abundance_shared_file_PATH, NATURAL_sorted_abundance_design_file_PATH)
	# ------------------------------ END OF SORT NATURAL ABUNDANCE SHARED FILE

	# ++++++++++++++++++++++++++++++ SORT NORMALIZED ABUNDANCE SHARED FILE
	NORMALIZED_sorted_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_NORMALIZED_SORTED_shared_file_STEP2.txt'
	NORMALIZED_sorted_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_NORMALIZED_SORTED_design_file_STEP2.txt'
	sort_shared_file_function(NORMALIZED_abundance_shared_file_PATH, NORMALIZED_abundance_design_file_PATH, NORMALIZED_sorted_abundance_shared_file_PATH, NORMALIZED_sorted_abundance_design_file_PATH)
	# ------------------------------ END OF SORT NATURAL ABUNDANCE SHARED FILE

	# ++++++++++++++++++++++++++++++ SORT RELATIVE ABUNDANCE SHARED FILE
	RELATIVE_sorted_abundance_shared_file_PATH = outputdir_PATH + 'BIOM_SLICER_RELATIVE_SORTED_shared_file_STEP2.txt'
	RELATIVE_sorted_abundance_design_file_PATH = outputdir_PATH + 'BIOM_SLICER_RELATIVE_SORTED_design_file_STEP2.txt'
	sort_shared_file_function(RELATIVE_abundance_shared_file_PATH, RELATIVE_abundance_design_file_PATH, RELATIVE_sorted_abundance_shared_file_PATH, RELATIVE_sorted_abundance_design_file_PATH)
	# ------------------------------ END OF SORT RELATIVE ABUNDANCE SHARED FILE
	
	return(NATURAL_sorted_abundance_shared_file_PATH, NATURAL_sorted_abundance_design_file_PATH, NORMALIZED_sorted_abundance_shared_file_PATH, NORMALIZED_sorted_abundance_design_file_PATH, RELATIVE_sorted_abundance_shared_file_PATH, RELATIVE_sorted_abundance_design_file_PATH)


def remove_ambiguity_from_taxonomy(taxonomy_LIST):
	#Default
	taxonomy_identifier = ['k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']
	#taxonomy_LIST = taxonomy_string.split(';')
	#First clean up
	
	current_TAX_LIST = []
	for each_tax_level in taxonomy_LIST:
		sluggified_tax_level = slugify(str(each_tax_level))
		if sluggified_tax_level.lower() in ['root', 'unknown', 'domain', ' ', 'k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']:
			pass
		else:
			current_TAX_LIST.append(sluggified_tax_level)

	last_index = 0
	tax_string = ''

	for each_level in taxonomy_identifier:
		each_level_index = taxonomy_identifier.index(each_level)
		try:
			current_TAX_LIST[each_level_index]
			last_index = each_level_index
			last_level = each_level
			if current_TAX_LIST[last_index][0:2] != each_level:
				tax_string += last_level + current_TAX_LIST[last_index] + ';'
			else:
				tax_string += current_TAX_LIST[last_index] + ';'
		except IndexError:
			if current_TAX_LIST[last_index][0:2] not in taxonomy_identifier:
				tax_string += '*' + last_level + current_TAX_LIST[last_index] + ';'
			else:
				tax_string += '*' + current_TAX_LIST[last_index] + ';'
	return tax_string.replace('_;', ';')
# ################################### MOTHUR FUNCTIONS ####################### #


def test_mothur(processors, outputdir, stderr, stdout, run_pid, mothur_exec_PATH):
	# test mothur to see if it is working and grab the version of mothur by scanning its output log
	mothur_input_dictionary = {}
	command = 'get.current'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_PATH
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
	parameter_list.append(',' + space + 'iters=100, calc=chao, groupmode=T')
	parameter_list.append(',' + space + 'freq=' + str(freq_value))
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
	parameter_list.append(',' + space + 'size=100000, iters=10000, calc=nseqs-simpson-invsimpson-shannon-sobs-ace')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_nseqs(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file):
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
	parameter_list.append(',' + space + 'calc=nseqs')
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


def mothur_relative_abundance(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file):
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


def mothur_corr_axes(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, axis_file, shared_file=None, design_file=None):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'corr.axes'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	if shared_file is None:
		parameter_list.append('metadata=' + design_file)
	elif design_file is None:
		parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'method=spearman')
	parameter_list.append(',' + space + 'axes=' + axis_file)
	parameter_list.append(',' + space + 'numaxes=3')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_coremicrobiome(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, min_abundance_ratio):
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.coremicrobiome'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'abundance=' + min_abundance_ratio)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_get_otus(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, accnos_file):
	space = ' '
	mothur_input_dictionary = {}
	command = 'get.otus'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'accnos=' + accnos_file)
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
	#parameter_list.append(',' + space + 'norm=5000')
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True


def mothur_remove_rare(processors, outputdir, stderr, stdout, run_pid, mothur_exec_path, shared_file, nseqs):
	# grab path of paired fastq files and save it into container
	space = ' '
	mothur_input_dictionary = {}
	command = 'remove.rare'
	mothur_input_dictionary['command'] = command
	mothur_input_dictionary['mothur_exec_path'] = mothur_exec_path
	mothur_input_dictionary['processors'] = processors
	mothur_input_dictionary['outputdir'] = outputdir
	mothur_input_dictionary['nohup_in'] = 'nohup'
	mothur_input_dictionary['nohup_out'] = '> ' + stdout + ' 2> ' + stderr + ' & echo $! > ' + run_pid
	parameter_list = []
	parameter_list.append('shared=' + shared_file)
	parameter_list.append(',' + space + 'nseqs=' + nseqs)
	mothur_input_dictionary['parameters'] = parameter_list
	make_file_object = mothur_process(mothur_input_dictionary)
	exec_dict = {}
	flag, exec_dict = make_file_object.execute_mothur_command()
	if flag is False:
		print "[FATAL-ERROR]: Commandline can not be executed!!!"
		print "ABORTING!!!"
		sys.exit(2)
	return True
# ################################### MAIN FUNCTION ########################## #


def main(argv):

	report_string = ''
	# ++++++++++++++++++++++++++++++ PARSE INPUT ARGUMENTS
	parser = argparse.ArgumentParser()
	main_file = parser.add_argument_group('Main file parameters')
	main_file.add_argument("--biom", help="Universal microbiom abundance matrix file(http://biom-format.org)", action='store')
	main_file.add_argument("--design", help="Design file: Tab delimited file to assign samples to a specific treatments, or other categories.", action='store')
	args = parser.parse_args()
	# ------------------------------ END OF PARSE INPUT ARGUMENTS

	# ++++++++++++++++++++++++++++++ BEURACRATICS PROCEDURES
	report_string += "######################################################################################################################################\n"
	print "######################################################################################################################################"
	report_string += "biom SLICER 1.0 EXECUTION HAS INITIATED" + '\n'
	print "biom SLICER 1.0 EXECUTION HAS INITIATED"
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

	# ++++++++++++++++++++++++++++++ OUTPUT DIRECTORY CHECKING
	args.outputdir = DEFAULT_OUTPUTDIR
	global report_file
	report_file = args.outputdir + "biom_slicer_report.txt"
	check_it_and_remove_it(report_file, True)
	report(report_string)
	# ------------------------------ END OF OUTPUT DIRECTORY CHECKING

	# ++++++++++++++++++++++++++++++ PROCESSORS CHECKING
	args.processors = DEFAULT_PROCESSORS
	# ------------------------------ END OF PROCESSORS CHECKING

	# ++++++++++++++++++++++++++++++ PREFIX NAME CHECKING
	args.prefix = DEFAULT_PREFIX
	# ------------------------------ END OF PREFIX NAME CHECKING

	# ++++++++++++++++++++++++++++++ EXECUTIVE DIRECTORY CHECKING
	args.execdir = DEFAULT_EXECDIR
	# ------------------------------ END OF EXECUTIVE DIRECTORY CHECKING

	# ++++++++++++++++++++++++++++++ CHECKING EXECUTIVES
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY[--execdir]"
	report("VERIFYING THE SANITY/VERSION OF EXECUTABLES IN EXECUTIVE DIRECTORY[--execdir]")
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
		error("Python version is older than 2.7")
		print "Python version is older than 2.7"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	report("1: MOTHUR")
	print "1: MOTHUR"
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	mothur_exec_PATH = args.execdir + 'mothur'
	if isFileExist(mothur_exec_PATH) is False:
		error("Your mothur file path has Access/Exist issue")
		print "Your mothur file path has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	else:
		print "mothur execution file path is: ", mothur_exec_PATH
		report("mothur execution file path is: " + mothur_exec_PATH)
		#report("Testing mothur executables: ")
		print "Testing mothur executable: "
		report("Testing mothur executable: ")
		flag, stderr = execute_functions(test_mothur, args.processors, args.outputdir, 'multi', 'mothur', mothur_exec_PATH)
		if flag is False:
			error("[" + FAILED_MARK + "]")
			print "[", FAILED_MARK, "]"
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
			report("Mothur executables responding successfully!!!")
			print "Mothur executables responding successfully!!!"
			report("version of mothur executables: " + target_lines[0])
			print "version of mothur executables:", target_lines[0]
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	report("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# ------------------------------ END OF CHECKING EXECUTIVES
	# ------------------------------ END OF BEURACRATICS PROCEDURES

	# ++++++++++++++++++++++++++++++ CHECKING INPUTS
	print "\n###################################################################"
	report("\n###################################################################")
	print "VERIFYING THE SANITY/VALIDITY OF INPUT FILES"
	report("VERIFYING THE SANITY/VALIDITY OF INPUT FILES")
	print "###################################################################\n"
	report("###################################################################\n")

	# ++++++++++++++++++++++++++++++ CHECKING biom FILE
	if args.biom is None:
		args.biom = DEFAULT_TESTDIR + "ZAC_biom.txt"
	if isFileExist(args.biom) is False:
		error("[--biom]: biom file has Access/Exist issue")
		print "[--biom]: biom file has Access/Exist issue"
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	report("biom_to_shared_converter: ")
	print "biom_to_shared_converter: "
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))

	biom_to_excel_converter(args.biom)

	shared_file_PATH = args.outputdir + args.prefix + '_shared_file_STEP1.txt'
	sample_metadata_file_PATH = args.outputdir + args.prefix + '_sample_metadata_file_STEP1.txt'
	OTU_metadata_file_PATH = args.outputdir + args.prefix + '_OTU_metadata_file_STEP1.txt'
	multi_sheets_shared_file_PATH = args.outputdir + args.prefix + '_multi_sheets_shared_file_STEP1.xlsx'
	multi_sheets_OTU_metadata_file_PATH = args.outputdir + args.prefix + 'multi_sheets_OTU_metadata_file_STEP1.txt'
	flag = biom_to_shared_converter(args.biom, shared_file_PATH, TAXONOMY_EXIST, sample_metadata_file_PATH, OTU_metadata_file_PATH, multi_sheets_shared_file_PATH, multi_sheets_OTU_metadata_file_PATH)
	if flag is False:
		print "ABORTING!!!"
		error("ABORTING!!!")
	else:
		print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
		report("biom_to_shared_converter executed successfully!!!")
		print "biom_to_shared_converter executed successfully!!!"
		
		args.shared = shared_file_PATH
	# ------------------------------ END OF CHECKING biom FILE

	# ++++++++++++++++++++++++++++++ CHECKING DESIGN FILE
	print "Design file checking and validation:"
	report("Design file checking and validation:")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	unified_design_file_PATH = args.outputdir + args.prefix + '_unified_design_file_STEP1.txt'
	design_metadata_unifier(args.design, sample_metadata_file_PATH, unified_design_file_PATH)
	args.design = unified_design_file_PATH
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "Design file checking and validation ended successfully!!!"
	report("Design file checking and validation ended successfully!!!")
	# ------------------------------ END OF CHECKING DESIGN FILE
	# ------------------------------ END OF CHECKING INPUTS
	print "###################################################################\n"
	report("###################################################################\n")

	# ++++++++++++++++++++++++++++++ MATCH SHARED FILE WITH DESIGN FILE
	print "Matching shared file with design file is in progress"
	report("Matching shared file with design file is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	matched_shared_file_PATH = args.outputdir + args.prefix + '_MATCHED_shared_file_STEP1.txt'
	matched_design_file_PATH = args.outputdir + args.prefix + '_MATCHED_design_file_STEP1.txt'
	flag = match_shared_design_function(args.shared, args.design, matched_shared_file_PATH, matched_design_file_PATH, mothur_exec_PATH, args.processors, args.outputdir)
	if flag is True:
		args.shared = matched_shared_file_PATH
		args.design = matched_design_file_PATH
	print "Matching shared file with design file PASSED!!!"
	report("Matching shared file with design file PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "###################################################################\n"
	report("###################################################################\n")
	# ----------------------------- END OF MATCH SHARED FILE WITH DESIGN FILE
	FINAL_STRING = ''
	FINAL_STRING = html_visualizer(FINAL_STRING)

	# ++++++++++++++++++++++++++++++ DESIGN MATRIX REVERSE AND FORWARD
	design_MATRIX_forward, design_MATRIX_reverse = design_file_to_design_MATRIX_converter(args.design)
	OTU_MATRIX_forward, OTU_MATRIX_reverse = design_file_to_design_MATRIX_converter(OTU_metadata_file_PATH)
	#Forward: {'phinchID': {'24': ['SAMPLE_ALIAS_025'], '25': ['SAMPLE_ALIAS_026'], '26': ['SAMPLE_ALIAS_027'],
	#Reverse: {'SAMPLE_ALIAS_024': {'phinchID': '23', 'temp': '20.8 C', 'collection_date': '2012-12-12T17:00:00+08:00'
	# ----------------------------- END OF DESIGN MATRIX REVERSE AND FORWARD
	
	# +++++++++++++++++++++++++++++ BRIEF STATISTICS TABLE FUNCTION
	print "STATISTICS TABLE FUNCTION is in progress"
	report("STATISTICS TABLE FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	#brief_statistics_table_html_string = ''
	#brief_statistics_table_html_string = brief_statistics_table_function(args.shared, design_MATRIX_reverse)
	statistics_table_html_string = statistics_table_function(args.shared, design_MATRIX_forward, design_MATRIX_reverse)
	FINAL_STRING = html_visualizer(FINAL_STRING, statistics_table_html_string)
	print "STATISTICS TABLE FUNCTION PASSED!!!"
	report("STATISTICS TABLE FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF BRIEF STATISTICS TABLE

	# +++++++++++++++++++++++++++++ RAREFACTION CURVE LINECHART FUNCTION
	print "RAREFACTION CURVE LINECHART FUNCTION is in progress"
	report("RAREFACTION CURVE LINECHART FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_rarefaction_html_string = ''
	plotly_rarefaction_javascript_string = ''
	plotly_rarefaction_html_string, plotly_rarefaction_javascript_string = plotly_Rarefaction_Curve_Linechart_DROPDOWN(args.shared, design_MATRIX_forward, design_MATRIX_reverse, mothur_exec_PATH, args.processors, args.outputdir)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_rarefaction_html_string, plotly_rarefaction_javascript_string)
	print "RAREFACTION CURVE LINECHART FUNCTION PASSED!!!"
	report("RAREFACTION CURVE LINECHART FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF RAREFACTION CURVE LINECHART FUNCTION

	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# ++++++++++++++++++++++++++++++ COMPLEXITY REDUCTION FACTORY
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	NATURAL_shared_file_PATH, NATURAL_design_file_PATH, NORMALIZED_shared_file_PATH, NORMALIZED_design_file_PATH, RELATIVE_shared_file_PATH, RELATIVE_design_file_PATH = complexity_reduction_factory(args.shared, args.design, mothur_exec_PATH, args.processors, args.outputdir)
	# ------------------------------------------------------------
	# ----------------------------- END OF COMPLEXITY REDUCTION FACTORY
	# ------------------------------------------------------------
	
	# +++++++++++++++++++++++++++++ NATURAL ABUNDANCE BARPLOT FUNCTION
	print "NATURAL ABUNDANCE BARPLOT FUNCTION is in progress"
	report("NATURAL ABUNDANCE BARPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_natural_abundance_barplot_html_string = ''
	plotly_natural_abundance_barplot_javascript_string = ''
	plotly_natural_abundance_barplot_html_string, plotly_natural_abundance_barplot_javascript_string = plotly_Natural_Abundance_Barplot_multi(multi_sheets_shared_file_PATH, args.design, OTU_metadata_file_PATH)
	#plotly_natural_abundance_barplot_html_string, plotly_natural_abundance_barplot_javascript_string = plotly_Natural_Abundance_Barplot(NATURAL_shared_file_PATH, NATURAL_design_file_PATH, OTU_metadata_file_PATH)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_natural_abundance_barplot_html_string, plotly_natural_abundance_barplot_javascript_string)
	print "NATURAL ABUNDANCE BARPLOT FUNCTION PASSED!!!"
	report("NATURAL ABUNDANCE BARPLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF NATURAL ABUNDANCE BARPLOT FUNCTION
	
	# +++++++++++++++++++++++++++++ ALPHA DIVERSITY BOXPLOT FUNCTION
	print "ALPHA DIVERSITY BOXPLOT FUNCTION is in progress"
	report("ALPHA DIVERSITY BOXPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_alpha_diversity_boxplot_html_string = ''
	plotly_alpha_diversity_boxplot_javascript_string = ''
	plotly_alpha_diversity_boxplot_html_string, plotly_alpha_diversity_boxplot_javascript_string = plotly_Alpha_Diversity_Index_Boxplot_DROPDOWN(NORMALIZED_shared_file_PATH, NORMALIZED_design_file_PATH, mothur_exec_PATH, args.processors, args.outputdir)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_alpha_diversity_boxplot_html_string, plotly_alpha_diversity_boxplot_javascript_string)
	print "ALPHA DIVERSITY BOXPLOT FUNCTION PASSED!!!"
	report("ALPHA DIVERSITY BOXPLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF ALPHA DIVERSITY BOXPLOT FUNCTION

	# +++++++++++++++++++++++++++++ RELATIVE ABUNDANCE BARPLOT FUNCTION
	print "RELATIVE ABUNDANCE BARPLOT FUNCTION is in progress"
	report("RELATIVE ABUNDANCE BARPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_Relative_abundance_Barplot_html_string = ''
	plotly_Relative_abundance_Barplot_javascript_string = ''
	plotly_Relative_abundance_Barplot_html_string, plotly_Relative_abundance_Barplot_javascript_string = plotly_Relative_abundance_Barplot(RELATIVE_shared_file_PATH, RELATIVE_design_file_PATH, OTU_metadata_file_PATH)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Relative_abundance_Barplot_html_string, plotly_Relative_abundance_Barplot_javascript_string)
	print "RELATIVE ABUNDANCE BARPLOT FUNCTION PASSED!!!"
	report("RELATIVE ABUNDANCE BARPLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF RELATIVE ABUNDANCE BARPLOT FUNCTION
	
	# +++++++++++++++++++++++++++++ RELATIVE ABUNDANCE PIECHART FUNCTION
	print "RELATIVE ABUNDANCE PIECHART FUNCTION is in progress"
	report("RELATIVE ABUNDANCE PIECHART FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_Relative_Abundance_Piechart_html_string = ''
	plotly_Relative_Abundance_Piechart_javascript_string = ''
	plotly_Relative_Abundance_Piechart_html_string, plotly_Relative_Abundance_Piechart_javascript_string = plotly_Relative_Abundance_Piechart(RELATIVE_shared_file_PATH, RELATIVE_design_file_PATH, OTU_metadata_file_PATH)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Relative_Abundance_Piechart_html_string, plotly_Relative_Abundance_Piechart_javascript_string)
	print "RELATIVE ABUNDANCE PIECHART FUNCTION PASSED!!!"
	report("RELATIVE ABUNDANCE PIECHART FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF RELATIVE ABUNDANCE PIECHART FUNCTION
	"""
	# +++++++++++++++++++++++++++++ PCOA FUNCTION
	if NO_BETA is False:
		print "PCOA FUNCTION is in progress"
		report("PCOA FUNCTION is in progress")
		print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
		plotly_Principal_Coordinate_Analysis_html_string = ''
		plotly_Principal_Coordinate_Analysis_javascript_string = ''
		plotly_Principal_Coordinate_Analysis_html_string, plotly_Principal_Coordinate_Analysis_javascript_string = plotly_Principal_Coordinate_Analysis(args.shared, args.design, mothur_exec_PATH, args.prefix, args.processors, args.outputdir)
		FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Principal_Coordinate_Analysis_html_string, plotly_Principal_Coordinate_Analysis_javascript_string)
		print "PCOA FUNCTION PASSED!!!"
		report("PCOA FUNCTION PASSED!!!")
		print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF PCOA FUNCTION

	# +++++++++++++++++++++++++++++ DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION
	print "DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION is in progress"
	report("DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_Dominant_Entities_Distribution_Boxplot_html_string = ''
	plotly_Dominant_Entities_Distribution_Boxplot_javascript_string = ''
	plotly_Dominant_Entities_Distribution_Boxplot_html_string, plotly_Dominant_Entities_Distribution_Boxplot_javascript_string = plotly_Dominant_Entities_Distribution_Boxplot(relative_shared_file, relative_design_file)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Dominant_Entities_Distribution_Boxplot_html_string, plotly_Dominant_Entities_Distribution_Boxplot_javascript_string)
	print "DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION PASSED!!!"
	report("DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF DOMINANT ENTITIES DISTRIBUTION BOX PLOT FUNCTION
	"""
	# +++++++++++++++++++++++++++++ FINALIZING
	FINAL_STRING = html_visualizer(FINAL_STRING)
	flag = remove_extension_files(CURRENT_PATH, '.logfile')
	write_string_down(FINAL_STRING, 'biom_slicer_result.html')
	zip_file_NAME = 'biom_slicer_result.zip'
	zip_it('biom_slicer_result.html', zip_file_NAME)
	print "BIOM SLICER EXECUTION COMPLETED."
	report("BIOM SLICER EXECUTION COMPLETED.")
	report("Completion time: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "Completion time: ", time.strftime("%Y-%m-%d %H:%M:%S")
	# ----------------------------- END OF FINALIZING


# ################################### FINITO ################################# #
if __name__ == "__main__": main(sys.argv[1:])
