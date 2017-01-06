# ################################### INFO ##################################### #
# Microbiome Slicer 1.0
# BY: Amir Shams
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
												$('#brief_statistics_table').DataTable( {
													"order": [[ 1, "asc" ]],
													"paging": true,
													"info": false,
													"ordering": true,
													"searching": false,
													"scrollY": "600px",
													"lengthMenu": [[25, 50, -1], [25, 50, "All"]]
												} );
												
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


def plotly_Natural_Abundance_Barplot_NORMAL(shared_file, design_file):
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


def plotly_Rarefaction_Curve_Linechart_DROPDOWN(shared_file, design_file, mothur_exec_path, processors, outputdir):
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


def plotly_Alpha_Diversity_Index_Boxplot_DROPDOWN(shared_file, design_file, mothur_exec_path, processors, outputdir):
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


def brief_statistics_table_function(shared_file_PATH, design_file_PATH):
	headers_list = []
	columns_dict = {}
	headers_list, columns_dict = mothur_shared_parser(shared_file_PATH)
	design_dictionary = design_dict_maker(design_file_PATH)
	biom_matrix_header_list = ['Groups', 'Sample Name', 'Total number', 'Classified number', 'Unclassified number', 'Classified percentage', 'Unclassified percentage']
	no_unclassified_flag = True
	if 'unclassified' in headers_list:
		no_unclassified_flag = False
	else:
		no_unclassified_flag = True
	thead_string = '								<thead>\n								<tr>\n'
	for thead in biom_matrix_header_list:
		thead_string += '									<th class="text-center">' + thead + '</th>\n'
	thead_string += '								</tr>\n								</thead>\n'
	tbody_string = '								<tbody>\n'
	bacterial_list = []
	
	for each_head in headers_list:
		if each_head.lower() in ['label', 'numotus', 'groups']:
			continue
		else:
			bacterial_list.append(each_head)
	biom_matrix_length = len(columns_dict['numOtus'])

	for each_row in range(biom_matrix_length):
		each_row_list = []
		for each_bacteria in bacterial_list:
			each_row_list.append(columns_dict[each_bacteria][each_row])
		each_row_list_int = list(map(int, each_row_list))
		total_read_count = int(sum(each_row_list_int))
		#we test for being outliers

		tbody_string += '									<tr>\n'
		tbody_string += '									<td class="text-center">' + str(design_dictionary[columns_dict['Groups'][each_row]]) + '</td>\n'
		tbody_string += '									<td class="text-center">' + str(columns_dict['Groups'][each_row]) + '</td>\n'
		
		tbody_string += '									<td class="text-center">' + str(total_read_count) + '</td>\n'
		if no_unclassified_flag is False:
			unclassified_read_count = int(columns_dict['unclassified'][each_row])
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


def plotly_Bacterial_relative_abundance_Barplot(shared_file, design_file):
	condensed_shared_dict = {}
	condensed_shared_dict = condense_shared_file_by_design(shared_file, design_file)
	
	design_name_list = condensed_shared_dict.keys()
	OTU_name_list = condensed_shared_dict[design_name_list[0]].keys()
	#print sample_name_list
	#print OTU_name_list
	#print shared_dict[sample_name_list[0]]
	
	Bacterial_relative_abundance_barplot_TRACE_list = []
	Bacterial_relative_abundance_barplot_LAYOUT = {}
	Bacterial_relative_abundance_barplot_FIGURE_objects_dict = {}

	for each_design in design_name_list:
		Total_value_of_each_OTU_list = []
		for each_OTU in OTU_name_list:
			Total_value_of_each_OTU_list.append(sum(condensed_shared_dict[each_design][each_OTU]))
		Bacterial_relative_abundance_barplot_TRACE_list.append(PLOTLY_GO.Bar(
			x=OTU_name_list,
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
			title='Normalized abundance',
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


def plotly_Relative_Abundance_Piechart(shared_file, design_file):
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


def biom_to_shared_converter(biom_file_PATH, prefix_NAME, biom_file_TSV_FORMAT_PATH, biom_file_SHARED_FORMAT_PATH, processors, outputdir):
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
		shared_file_string += 'biom_SLICER' + '\t' + slugify(each_Sample) + '\t' + numOtus_value + '\t' + list_to_string(sample_value_list, '\t') + '\n'
		
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
	text = target_word.replace(' ', r'_').replace('\\', r'_').replace('`', r'_').replace('*', r'_').replace('{', r'_').replace('}', r'_').replace('[', r'_').replace(']', r'_').replace('(', r'_').replace(')', r'_').replace('>', r'_').replace('#', r'_').replace('+', r'_').replace('-', r'_').replace('.', r'_').replace('!', r'_').replace('$', r'_').replace("'", r'_').replace('"', r'_').replace('\/', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace('__', r'_').replace(';_', r';')
	return text


def list_to_string(query_list, delimiter=None):
	# convert list to string by specified delimiter
	if delimiter is None:
		delimiter = '\t'
	return delimiter.join(query_list)


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
			elif mode == 'partial':
				for each_ext in ext_list:
					if ext in each_ext:
						container.append(os.path.join(root, currentFile))
			elif mode == 'ex_partial':
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
# ################################### SPECIFIC ############################### #


def make_default_design(shared_file, prefix, DEFAULT_DESIGN):
	header_list = []
	columns_dict = {}
	sample_list = []
	header_list, columns_dict = mothur_shared_parser(shared_file)
	design_string = 'sample\tdesign\n'
	if 'Groups' in columns_dict:
		sample_list = columns_dict['Groups']
	elif 'Group' in columns_dict:
		sample_list = columns_dict['Group']
	for i in sample_list:
		design_string += i + '\t' + prefix + '\n'
	write_string_down(design_string, DEFAULT_DESIGN)
	return True


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
		flag = scandirs(outputdir_PATH, scanned_container_list, extension_list, 'ex_partial')
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


def sort_shared_file_function(shared_file_PATH, sorted_shared_file_PATH, prefix_NAME):
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


def OTU_name_reduction(OTU_name):
	new_name = ''
	if len(OTU_name) > 100:
		new_name += OTU_name[:50]
		new_name += '...'
		new_name += OTU_name[-50:]
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


def condense_shared_file_by_design(shared_file, design_file):
	shared_dict = {}
	shared_dict = parse_shared_file(shared_file)

	reverse_design_dict = {}
	reverse_design_dict = design_dict_maker(design_file, 'reverse')

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
	biom_file_TSV_FORMAT_PATH = args.outputdir + args.prefix + '_biom_file_TSV_FORMAT_STEP1.txt'
	biom_file_SHARED_FORMAT_PATH = args.outputdir + args.prefix + '_biom_file_SHARED_FORMAT_STEP1.txt'
	flag = biom_to_shared_converter(args.biom, args.prefix, biom_file_TSV_FORMAT_PATH, biom_file_SHARED_FORMAT_PATH, args.processors, args.outputdir)
	if flag is False:
		print "ABORTING!!!"
		error("ABORTING!!!")
	else:
		print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
		report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
		report("biom_to_shared_converter executed successfully!!!")
		print "biom_to_shared_converter executed successfully!!!"
		args.shared = biom_file_SHARED_FORMAT_PATH

	# ------------------------------ END OF CHECKING biom FILE

	# ++++++++++++++++++++++++++++++ CHECKING DESIGN FILE
	print "Design file checking and validation:"
	report("Design file checking and validation:")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	if args.design is None:
		# DESIGN FILE MAKING
		DEFAULT_DESIGN = args.outputdir + 'biom_slicer_default_design.txt'
		flag = make_default_design(args.shared, args.prefix, DEFAULT_DESIGN)
		if flag is False:
			error("[--design]: Something is wrong during design file making")
			print "[--design]: Something is wrong during design file making"
			print "ABORTING!!!"
			error("ABORTING!!!")
			sys.exit(2)
		else:
			args.design = DEFAULT_DESIGN
			NO_BETA = True
	elif isFileExist(args.design) is False:
		print "[--design]: Design file has Access/Exist issue"
		error("[--design]: Design file has Access/Exist issue")
		print "ABORTING!!!"
		error("ABORTING!!!")
		sys.exit(2)
	elif isFileExist(args.design) is True:
		NO_BETA = False
	design_file_VALIDATED_PATH = args.outputdir + args.prefix + '_VALIDATED_design_file_STEP1.txt'
	design_file_validator(args.design, design_file_VALIDATED_PATH)
	args.design = design_file_VALIDATED_PATH
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
	control_from_design_file_PATH = args.outputdir + args.prefix + '_Matched_CONTROL_from_DESIGN_STEP2.txt'
	flag = design_to_control_converter(args.design, control_from_design_file_PATH)
	if flag is True:
		pass
	matched_shared_file_PATH = args.outputdir + args.prefix + '_Matched_shared_file_STEP2.txt'
	matched_design_file_PATH = args.outputdir + args.prefix + '_Matched_design_file_STEP2.txt'
	flag = match_shared_with_design(args.shared, args.design, control_from_design_file_PATH, matched_shared_file_PATH, matched_design_file_PATH, args.prefix, mothur_exec_PATH, args.processors, args.outputdir)
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
	
	# +++++++++++++++++++++++++++++ BRIEF STATISTICS TABLE FUNCTION
	print "BRIEF STATISTICS TABLE FUNCTION is in progress"
	report("BRIEF STATISTICS TABLE FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	brief_statistics_table_html_string = ''
	brief_statistics_table_html_string = brief_statistics_table_function(args.shared, args.design)
	FINAL_STRING = html_visualizer(FINAL_STRING, brief_statistics_table_html_string)
	print "BRIEF STATISTICS TABLE FUNCTION PASSED!!!"
	report("BRIEF STATISTICS TABLE FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF BRIEF STATISTICS TABLE

	# ++++++++++++++++++++++++++++++ SORT SHARED FILE & FILTER
	sorted_shared_file_PATH = args.outputdir + args.prefix + '_SORTED_shared_file_STEP3.txt'
	sort_shared_file_function(args.shared, sorted_shared_file_PATH, args.prefix)
	args.shared = sorted_shared_file_PATH
	# ------------------------------ END OF SORT SHARED FILE & FILTER

	# +++++++++++++++++++++++++++++ NATURAL ABUNDANCE BARPLOT FUNCTION
	print "NATURAL ABUNDANCE BARPLOT FUNCTION is in progress"
	report("NATURAL ABUNDANCE BARPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_natural_abundance_barplot_html_string = ''
	plotly_natural_abundance_barplot_javascript_string = ''
	plotly_natural_abundance_barplot_html_string, plotly_natural_abundance_barplot_javascript_string = plotly_Natural_Abundance_Barplot_NORMAL(args.shared, args.design)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_natural_abundance_barplot_html_string, plotly_natural_abundance_barplot_javascript_string)
	print "NATURAL ABUNDANCE BARPLOT FUNCTION PASSED!!!"
	report("NATURAL ABUNDANCE BARPLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF NATURAL ABUNDANCE BARPLOT FUNCTION

	# +++++++++++++++++++++++++++++ RAREFACTION CURVE LINECHART FUNCTION
	print "RAREFACTION CURVE LINECHART FUNCTION is in progress"
	report("RAREFACTION CURVE LINECHART FUNCTIONis in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_rarefaction_html_string = ''
	plotly_rarefaction_javascript_string = ''
	plotly_rarefaction_html_string, plotly_rarefaction_javascript_string = plotly_Rarefaction_Curve_Linechart_DROPDOWN(args.shared, args.design, mothur_exec_PATH, args.processors, args.outputdir)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_rarefaction_html_string, plotly_rarefaction_javascript_string)
	print "RAREFACTION CURVE LINECHART FUNCTION PASSED!!!"
	report("RAREFACTION CURVE LINECHART FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF RAREFACTION CURVE LINECHART FUNCTION

	# +++++++++++++++++++++++++++++ ALPHA DIVERSITY BOXPLOT FUNCTION
	print "ALPHA DIVERSITY BOXPLOT FUNCTION is in progress"
	report("ALPHA DIVERSITY BOXPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_alpha_diversity_boxplot_html_string = ''
	plotly_alpha_diversity_boxplot_javascript_string = ''
	plotly_alpha_diversity_boxplot_html_string, plotly_alpha_diversity_boxplot_javascript_string = plotly_Alpha_Diversity_Index_Boxplot_DROPDOWN(args.shared, args.design, mothur_exec_PATH, args.processors, args.outputdir)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_alpha_diversity_boxplot_html_string, plotly_alpha_diversity_boxplot_javascript_string)
	print "ALPHA DIVERSITY BOXPLOT FUNCTION PASSED!!!"
	report("ALPHA DIVERSITY BOXPLOT FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF ALPHA DIVERSITY BOXPLOT FUNCTION

	# +++++++++++++++++++++++++++++ REMOVING LOW ABUNDANCE SAMPLES FUNCTION
	low_abundance_sample_list = []
	low_abundance_sample_list = get_low_abundance_sample_list(args.shared, mothur_exec_PATH, args.processors, args.outputdir)
	if len(low_abundance_sample_list) == 0:
		print "LOW ABUNDANCE SAMPLE COUNT:", len(low_abundance_sample_list)
	else:
		print "LOW ABUNDANCE SAMPLE COUNT:", len(low_abundance_sample_list)
		print "REMOVING LOW ABUNDANCE SAMPLES!!"
		print low_abundance_sample_list
		filter_list_control_file = args.outputdir + args.prefix + '_SAMPLE_filter_list_control_file.txt'
		filter_list_to_control_file_converter(low_abundance_sample_list, filter_list_control_file)
		samples_removed_shared_file = args.outputdir + args.prefix + '_SAMPLES_REMOVED_shared_file_STEP3.txt'
		samples_removed_design_file = args.outputdir + args.prefix + '_SAMPLES_REMOVED_design_file_STEP3.txt'
		flag = remove_sample_shared_design(args.shared, args.design, filter_list_control_file, samples_removed_shared_file, samples_removed_design_file, args.prefix, mothur_exec_PATH, args.processors, args.outputdir)
		if flag is True:
			print "Shared file initial update is successfull."
			print "# ##########################"
			print "Shared file replaced with Step2 Shared file(SAMPLES_REMOVED_shared_file_STEP3.txt)"
			args.shared = samples_removed_shared_file
			print "Design file replaced with Step2 Design file(SAMPLES_REMOVED_shared_file_STEP3.txt)"
			args.design = samples_removed_design_file
			print "# ##########################"
		else:
			print "Something is wrong with update_shared_file."
			sys.exit(2)
	print "SAMPLE ABUNDANCE EVALUATION STEP PASSED!!"
	# ----------------------------- END OF REMOVING LOW ABUNDANCE SAMPLES FUNCTION

	# +++++++++++++++++++++++++++++ REMOVING LOW ABUNDANCE OTUS FUNCTION
	low_abundance_OTU_list = []
	low_abundance_OTU_list = get_low_abundance_otu_list(args.shared)
	if len(low_abundance_OTU_list) == 0:
		print "LOW ABUNDANCE OTU COUNT:", len(low_abundance_OTU_list)
	else:
		print "LOW ABUNDANCE OTU COUNT:", len(low_abundance_OTU_list)
		print "REMOVING LOW ABUNDANCE OTUS!!"
		print low_abundance_OTU_list
		filter_list_control_file = args.outputdir + args.prefix + '_OTU_filter_list_control_file.txt'
		filter_list_to_control_file_converter(low_abundance_OTU_list, filter_list_control_file)
		
		otu_remove_shared_file = args.outputdir + args.prefix + '_OTUs_REMOVED_shared_file_STEP4.txt'
		otu_remove_design_file = args.outputdir + args.prefix + '_OTUs_REMOVED_design_file_STEP4.txt'
		print "remove_otu_file is set: so I am going to remove otu only listed."
		flag = remove_OTU_shared_design(args.shared, args.design, filter_list_control_file, otu_remove_shared_file, otu_remove_design_file, args.prefix, mothur_exec_PATH, args.processors, args.outputdir)
		if flag is True:
			print "remove_otu_shared_design is successfull."
			print "# ##########################"
			print "Shared file replaced with Step4 Shared file(OTUs_REMOVED_shared_file_STEP4.txt)"
			args.shared = otu_remove_shared_file
			print "Design file replaced with Step4 Design file(OTUs_REMOVED_design_file_STEP4.txt)"
			args.design = otu_remove_design_file
			print "# ##########################"
		else:
			print "Something is wrong with remove_otu_shared_design."
			sys.exit(2)
	print "OTU ABUNDANCE EVALUATION STEP PASSED!!"
	# ----------------------------- END OF REMOVING LOW ABUNDANCE OTUS FUNCTION

	# +++++++++++++++++++++++++++++ biomARKER DISCOVERY FUNCTION
	if NO_BETA is False:
		ranked_shared = biomarker_discovery(args.shared, args.design, args.prefix, mothur_exec_PATH, args.processors, args.outputdir)
		args.shared = ranked_shared
	# ----------------------------- END OF biomARKER DISCOVERY FUNCTION

	# +++++++++++++++++++++++++++++ RELATIVE ABUNDANCE CALCULATION FUNCTION
	relative_shared_file = args.outputdir + args.prefix + '_RELATIVE_shared_file_STEP7.txt'
	relative_design_file = args.outputdir + args.prefix + '_RELATIVE_design_file_STEP7.txt'
	flag = relative_abundance_shared_design_maker(args.shared, args.design, args.prefix, relative_shared_file, relative_design_file, mothur_exec_PATH, args.processors, args.outputdir)
	# ----------------------------- END OF RELATIVE ABUNDANCE CALCULATION FUNCTION

	# +++++++++++++++++++++++++++++ RELATIVE ABUNDANCE BARPLOT FUNCTION
	print "RELATIVE ABUNDANCE BARPLOT FUNCTION is in progress"
	report("RELATIVE ABUNDANCE BARPLOT FUNCTION is in progress")
	print "Execution started at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution started at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	plotly_Bacterial_relative_abundance_Barplot_html_string = ''
	plotly_Bacterial_relative_abundance_Barplot_javascript_string = ''
	plotly_Bacterial_relative_abundance_Barplot_html_string, plotly_Bacterial_relative_abundance_Barplot_javascript_string = plotly_Bacterial_relative_abundance_Barplot(relative_shared_file, relative_design_file)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Bacterial_relative_abundance_Barplot_html_string, plotly_Bacterial_relative_abundance_Barplot_javascript_string)
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
	plotly_Relative_Abundance_Piechart_html_string, plotly_Relative_Abundance_Piechart_javascript_string = plotly_Relative_Abundance_Piechart(relative_shared_file, relative_design_file)
	FINAL_STRING = html_visualizer(FINAL_STRING, plotly_Relative_Abundance_Piechart_html_string, plotly_Relative_Abundance_Piechart_javascript_string)
	print "RELATIVE ABUNDANCE PIECHART FUNCTION PASSED!!!"
	report("RELATIVE ABUNDANCE PIECHART FUNCTION PASSED!!!")
	print "Execution completed at ", time.strftime("%Y-%m-%d %H:%M:%S")
	report("Execution completed at " + time.strftime("%Y-%m-%d %H:%M:%S"))
	# ----------------------------- END OF RELATIVE ABUNDANCE PIECHART FUNCTION

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
	
	# +++++++++++++++++++++++++++++ FINALIZING
	FINAL_STRING = html_visualizer(FINAL_STRING)
	flag = remove_extension_files(CURRENT_PATH, '.logfile')
	write_string_down(FINAL_STRING, 'biom_slicer_result.html')
	zip_file_NAME = 'biom_slicer_result.zip'
	zip_it('biom_slicer_result.html', zip_file_NAME)
	print "biom SLICER EXECUTION COMPLETED."
	report("biom SLICER EXECUTION COMPLETED.")
	report("Completion time: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	print "Completion time: ", time.strftime("%Y-%m-%d %H:%M:%S")
	# ----------------------------- END OF FINALIZING


# ################################### FINITO ################################# #
if __name__ == "__main__": main(sys.argv[1:])
