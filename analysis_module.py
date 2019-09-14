##########################################
# 
# Written by Chris Karwin, 2/11/15, UCI
# Module name: analysis
# Classes: 
#	1) Analysis: superclass
#	2) Iteration(Analysis)
#	3) PyLikelihood(Iteration,Analysis)
# 	4) Plots(Analysis)
#
#	Purpose: Perform a binned analysis of the Fermi-LAT data.
#
###########################################

###########################################
#
# Index of functions:
#
#	Analysis Class(superclass):
#		-gtselect()
#		-maketime()
#		-cmap()
#		-ccube()
#		-expCube()	
#		-expMap()
#		-srcMaps(xml_file, outfile)
#		-run_gtlike(xml_file, outfile)
#		-model_map(src_file, xml_file, outfile, out_type)
#		-Residual_map(cmap, model, outfile, operation)
#		-Residual_Map2(cmap, model, ccube, modelcube)
#		-counts_map_small
#		-Make_TS(model_file_xml, xref, yref, nxpix, nypix, bin_size, outfile)
#		-gtobsSim(infile, srclist, seed, name)
#
#	Iteration(Analysis):
#		-ff_source(source, ff)
#		-ff_by_name(root, name_list, ff)
#		-iterate(starting_xml_infile, src_file, mass) 
#
#	PyLikelihood(Iteration, Analysis)
#		-run_gtlike(src_file, xml_infile, xml_outfile, mass, pass_count)
#		
#	Plots(Analysis, Iteration)
#		-plot_main(x_list, y_list, savefig, Title, x_label, y_label)
#		-Counts_spectra_plots()
#		-calculate_flux(model_file, exp_file, name, low, high)
#		-pylikelihood_flux(src, xml, source_name_list)
#		-pylikelihood_counts(input_list, source_plot_list)
#		-pyLikelihood_fractional_residuals(src_file,xml_file)		
#
############################################


############################################
# Methodology:
#
# For any given analysis do the following:
#	
#	1) Define a new directory
#	2) Make a parameter card
#		-This file defines all the neccessary parameters for the analysis
#		-The parameter card must follow the same format as that in parameters.txt
#	3) Make a TS file
#		-This file contains the point sources to be iterated over
#		-The TS file must follow the same format as TS.txt
#	4) Make a client code
#		-This is a python script that specifies the neccessary functions to be called
#		-The file client.py can be used as a template
#	5) If one is scanning over a given parameter, then make a loop file		
#		-This file is used to submit multiple jobs
#		-The file loop.py can be used as a template
#
###########################################

###########################################
#
# Analysis class: superclass
# Note: Instances of the Analysis class take parameters.txt and TS.txt as an input 
#

#Imports:
import gt_apps as my_apps
from GtApp import GtApp
import pyfits
import xml.etree.ElementTree as ET
import pyLikelihood
from BinnedAnalysis import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from pylab import *
import time
import matplotlib.gridspec as gridspec
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.table import Table

class Analysis:
	
# constructor method initializes parameters passed from parameters.txt
	
	def __init__(self, parameters, TS):
		
		self.file_parameters = parameters
		list = []
		with open(self.file_parameters,'r') as f:
			for line in f:
				split = line.split()
				list.append(split[2])
		self.name = list[0] #name of analysis
		self.nxpix = int(list[1]) #number of x pixels
		self.nypix = int(list[2]) #number of y pixels
		self.xref = float(list[3]) #center of roi (galactic longitude, l)
		self.yref = float(list[4]) #center of roi (galactic latitude, b)
		self.binsz = float(list[5]) #degrees per pixel
		self.coordsys = list[6] #coordinate system (Galactic coordinates)
		self.emin = int(list[7]) #min energy of analysis
		self.emax = int(list[8]) #max energy of analysis
		self.enumbins = int(list[9]) #number of energy bins
		self.proj = list[10] #projection method
		self.scfile = list[11] #spacecraft file 
		self.infile = list[12] #input data file
		self.irfs = list[13] #instrument response functions to use		
 		self.evclass = int(list[14]) #source/clean, etc. 
		self.evtype = int(list[15]) #back/front conversion, etc. 
		self.reduced_x = int(list[16]) #number of reduced x pixels, in accordance with ccube requirements
		self.reduced_y = int(list[17]) #number of reduced y pixels, in accordance with ccube requirements
		self.path = list[18] #path to data files
		self.ROI_RA = list[19] #RA for ROI center; for data from NASA server use 'INDEF'
		self.ROI_DEC = list[20] #DEC for ROI center; for data from NASA server use 'INDEF'
		self.ROI_radius = list[21] #ROI radius; for data from NASA server use 'INDEF'
		self.t_min = list[22] #start time for observations in MET; for data from NASA server use 'INDEF'	
		self.t_max = list[23] #end time for observations in MET; for data from NASA server use 'INDEF'
		self.z_max = list[24] #max zenith angle
		self.file_TS = TS
		with open(self.file_TS,'r') as f:
			self.TS = eval(f.read())
		

	#run gtselect (this essentially makes the zenith cut on the initial data)
	
	def gtselect(self):
		
		my_apps.filter['evclass'] = self.evclass
		my_apps.filter['ra'] = self.ROI_RA 
		my_apps.filter['dec'] = self.ROI_DEC 
		my_apps.filter['rad'] = self.ROI_radius  
		my_apps.filter['emin'] = self.emin
		my_apps.filter['emax'] = self.emax
		my_apps.filter['zmax'] = self.z_max
		my_apps.filter['tmin'] = self.t_min
		my_apps.filter['tmax'] = self.t_max 
		my_apps.filter['infile'] = '%s_binned_events.txt' % self.name
		my_apps.filter['evtype'] = self.evtype
		my_apps.filter['outfile'] =  '%s_binned_filtered.fits' % self.name
		my_apps.filter.run()
		return 'outfile'
	
	#maketime

	def maketime(self):

		my_apps.maketime['scfile'] = self.scfile
		my_apps.maketime['roicut'] = 'no'
		my_apps.maketime['evfile'] = '%s_binned_filtered.fits' % self.name
		my_apps.maketime['outfile'] = '%s_binned_gti.fits' % self.name
		my_apps.maketime['filter'] = 'DATA_QUAL==1 && (LAT_CONFIG==1)' 
		my_apps.maketime.run()
		return 'outfile'

	#make a 2D counts map

	def cmap(self):
		
		my_apps.counts_map['algorithm'] = 'CMAP'
		my_apps.counts_map['evfile'] = '%s_binned_gti.fits' % self.name
		my_apps.counts_map['outfile'] = '%s_cmap.fits' % self.name
		my_apps.counts_map['nxpix'] = self.nxpix
		my_apps.counts_map['nypix'] = self.nypix
		my_apps.counts_map['binsz'] = self.binsz
		my_apps.counts_map['coordsys'] = self.coordsys
		my_apps.counts_map['xref'] = self.xref
		my_apps.counts_map['yref'] = self.yref	
		my_apps.counts_map['proj'] = self.proj
		my_apps.counts_map.run()
		return 'outfile'

	#make a 3D counts map
	#It is very important to change the number of pixels in this step!
		
	def ccube(self):
		
		my_apps.counts_map['algorithm'] = 'CCUBE'
		my_apps.counts_map['evfile'] = '%s_binned_gti.fits' % self.name
		my_apps.counts_map['outfile'] = '%s_ccube.fits' % self.name
		my_apps.counts_map['nxpix'] = self.reduced_x
		my_apps.counts_map['nypix'] = self.reduced_y
		my_apps.counts_map['binsz'] = self.binsz
		my_apps.counts_map['coordsys'] = self.coordsys
		my_apps.counts_map['xref'] = self.xref
		my_apps.counts_map['yref'] = self.yref
		my_apps.counts_map['proj'] = self.proj
		my_apps.counts_map['ebinalg'] = 'LOG'
		my_apps.counts_map['emin'] = self.emin
		my_apps.counts_map['emax'] = self.emax
		my_apps.counts_map['enumbins'] = self.enumbins
		my_apps.counts_map.run()
		return 'outfile'

	#generate livetime cube

	def expCube(self):

		my_apps.expCube['evfile'] = '%s_binned_gti.fits' % self.name
		my_apps.expCube['scfile'] =  self.scfile
		my_apps.expCube['outfile'] = '%s_binned_ltcube.fits' % self.name
		my_apps.expCube['dcostheta'] = 0.025
		my_apps.expCube['binsz'] = 1
		my_apps.expCube.run()
		return 'outfile'

	
	#generate exposure map(here we are just doing all-sky)

	def expMap(self):
		
		expMapBinned = GtApp('gtexpcube2' , 'Likelihood')
		expMapBinned.command()

		expMapBinned['infile'] = '%s_binned_ltcube.fits' % self.name
		expMapBinned['cmap'] = 'none'
		expMapBinned['outfile'] = '%s_binned_expMap.fits' % self.name
		expMapBinned['irfs'] = self.irfs
		expMapBinned['nxpix'] = 1800 #all sky
		expMapBinned['nypix'] = 900  #all sky	
		expMapBinned['binsz'] = self.binsz
		expMapBinned['coordsys'] = self.coordsys
		expMapBinned['xref'] = self.xref
		expMapBinned['yref'] =  self.yref
		expMapBinned['emin'] = self.emin
		expMapBinned['emax'] = self.emax
		expMapBinned['proj'] = 'STG'
		expMapBinned['enumbins'] = self.enumbins
		expMapBinned.run()
		return 'outfile'

	#generate source map

	def srcMaps(self, xml_file, outfile):

		srcMapsBinned = GtApp('gtsrcmaps' , 'Likelihood') 
		srcMapsBinned.command()
	
		srcMapsBinned['expcube'] = '%s%s_binned_ltcube.fits' % (self.path,self.name)
		srcMapsBinned['cmap'] = '%s%s_ccube.fits' % (self.path,self.name)
		srcMapsBinned['srcmdl'] = xml_file
		srcMapsBinned['bexpmap'] = '%s%s_binned_expMap.fits' % (self.path,self.name)
		srcMapsBinned['outfile'] = outfile
		srcMapsBinned['irfs'] = self.irfs
		srcMapsBinned['ptsrc'] = 'yes'
		
		srcMapsBinned.run()
		return 'outfile'
	
	#run gtlike

	def run_gtlike(self, xml_file, outfile):

		likeBinned = GtApp('gtlike' , 'Likelihood')
		likeBinned.command()

		likeBinned['sfile'] = outfile
		likeBinned['statistic'] = 'BINNED'
		likeBinned['cmap'] = '%s_binned_srcmaps.fits' % self.name
		likeBinned['bexpmap'] = '%s_binned_expMap.fits' % self.name
		likeBinned['expcube'] = '%s_binned_ltcube.fits' % self.name
		likeBinned['srcmdl'] = xml_file
		likeBinned['irfs'] = self.irfs
		likeBinned['optimizer'] = 'NEWMINUIT'
		likeBinned['results'] = '%s_results.dat' % self.name
		likeBinned['specfile'] = '%s_counts-spectra.fits' % self.name
		likeBinned.run()
		return 		

	#make model map

	def model_map(self, src_file, xml_file, outfile, outtype): #outtype is either 'cmap' or 'ccube'
		
		my_apps.model_map['srcmaps'] = src_file
		my_apps.model_map['srcmdl'] = xml_file
		my_apps.model_map['outfile'] = outfile
		my_apps.model_map['irfs'] = self.irfs
		my_apps.model_map['expcube'] = '%s%s_binned_ltcube.fits' % (self.path,self.name)
		my_apps.model_map['bexpmap'] = '%s%s_binned_expMap.fits' % (self.path,self.name)
		my_apps.model_map['outtype'] = outtype
		my_apps.model_map.run()
		return 'outfile'		
	
	#produce total residual map
	
	def Residual_map(self, cmap, model, outfile, operation):
		
			
		hdulist_A = pyfits.open(cmap)
		cmap = hdulist_A[0].data
		hdulist_B = pyfits.open(model)
		model = hdulist_B[0].data
	

		for x in range(0,self.reduced_x):
			for y in range(0,self.reduced_y):
				if operation == 'minus':
					model[y][x] = cmap[y][x] - model[y][x] 
					#model[y][x] = (cmap[y][x] - model[y][x])/math.sqrt(model[y][x]) #for template map
				if operation == 'plus':
					model[y][x] = cmap[y][x] + model[y][x]
		hdulist_B.writeto(outfile)
	
	#produce redisual maps in 3 energy bins
	#0-3 corresponds to 1.0 - 1.6 GeV
	#3-10 corresponds to 1.6 - 10 GeV
	#10-100 corresponds to 10 - 100 GeV
	
	#For M31:
	#0-6 corresponds to 300 MeV - 1.2 GeV
	#6-15 corresponds to 1.2 GeV - 9.5 GeV
	#15-30 corresponds to 9.5 GeV - 300 GeV

	def Residual_Map2(self, cmap, model, ccube, modelcube):

		
		hdulist_A = pyfits.open(cmap)
		cmap = hdulist_A[0].data
		hdulist_B = pyfits.open(model)
		model = hdulist_B[0].data
		hdulist_C = pyfits.open(ccube)
		ccube = hdulist_C[0].data
		hdulist_D = pyfits.open(modelcube)
		modelcube = hdulist_D[0].data
		
		energy = hdulist_C[1].data
				
		M31_2 = [[0,10],[10,30]]
		M31_bins = [[0,10],[10,18],[18,30]] #300 MeV - 3 GeV, 3 GeV - 18.9 GeV, 18.9 GeV - 300 GeV
		GC_bins = [[0,10],[10,20]]
		GC_bins_3 = [[0,5],[5,13],[13,20]] #1-3.2, 3.2-20.0,20.0-100	
		HE_defecit = [[0,13],[13,17],[17,20]]
		All = [[0,30]]
		for sublist in GC_bins_3:
		
			for x in range(0,self.reduced_x):
				for y in range(0,self.reduced_y):
					cmap[x][y] = 0
		
			for E in range(sublist[0],sublist[1]):
				for x in range(0,self.reduced_x):
					for y in range(0,self.reduced_y):
						cmap[x][y] =  cmap[x][y] + ccube[E][x][y]
			for x in range(0,self.reduced_x):
				for y in range(0,self.reduced_y):
					model[x][y] = 0


			for E in range(sublist[0],sublist[1]):
				for x in range(0,self.reduced_x):
					for y in range(0,self.reduced_y):
						model[x][y] = model[x][y] + modelcube[E][x][y]

			for x in range(0,self.reduced_x):
				for y in range(0,self.reduced_y):
					#model[x][y] = (cmap[x][y] - model[x][y]) #normal spatial residuals
					model[x][y] = (cmap[x][y] - model[x][y])/model[x][y] #for template map
			#name = self.name + '_' + str(sublist[0]) + '-' + str(sublist[1]) + '_residual.fits'
			name = self.name + '_' + str(sublist[0]) + '-' + str(sublist[1]) + '_fractional_residual.fits'#for template map
			hdulist_B.writeto(name)
	
			

	#match the size of the cmap to the model map for comparison in residual map

	def counts_map_small(self):
		
		my_apps.counts_map['algorithm'] = 'CMAP'
		my_apps.counts_map['evfile'] = '%s_binned_gti.fits' % self.name
		my_apps.counts_map['outfile'] = '%s_cmap_small.fits' % self.name
		my_apps.counts_map['nxpix'] = self.reduced_x
		my_apps.counts_map['nypix'] = self.reduced_y
		my_apps.counts_map['binsz'] = self.binsz
		my_apps.counts_map['coordsys'] = self.coordsys
		my_apps.counts_map['xref'] = self.xref
		my_apps.counts_map['yref'] = self.yref
		my_apps.counts_map['proj'] = self.proj
		my_apps.counts_map.run()
		return 'outfile'

	#make a TS map

	def Make_TS(self, model_file_xml, xref, yref, nxpix, nypix, bin_size, outfile):

		my_apps.TsMap['statistic'] = "Binned"
		my_apps.TsMap['scfile'] = self.scfile
		my_apps.TsMap['evfile'] = "%s%s_binned_gti.fits" % (self.path, self.name)
		my_apps.TsMap['bexpmap'] = "%s%s_binned_expMap.fits" % (self.path, self.name)
		my_apps.TsMap['expcube'] = "%s%s_binned_ltcube.fits" % (self.path, self.name)
		my_apps.TsMap['srcmdl'] = model_file_xml
		my_apps.TsMap['cmap'] = "%s%s_ccube.fits" % (self.path, self.name)
		my_apps.TsMap['irfs'] = self.irfs
		my_apps.TsMap['optimizer'] = "NEWMINUIT"
		my_apps.TsMap['outfile'] = outfile
		my_apps.TsMap['nxpix'] = nxpix
		my_apps.TsMap['nypix'] = nypix
		my_apps.TsMap['binsz'] = bin_size
		my_apps.TsMap['coordsys'] = 'GAL'
		my_apps.TsMap['xref'] = xref
		my_apps.TsMap['yref'] = yref
		my_apps.TsMap['proj'] = self.proj
		my_apps.TsMap.run()
		return 'outfile'
	
	def gtobsSim(self,infile,srclist,seed,name):
		
		simulate = GtApp('gtobssim','Likelihood')
		simulate.command()	
		
		simulate['scfile'] = '/pub/ckarwin/GC_Analysis/Data/photon_file_IG-ft2.fits'
		simulate['infile'] = 'input_library.xml'
		simulate['srclist'] = 'source_list.txt'
		simulate['evroot'] = name
		simulate['simtime'] = 410227199
		simulate['tstart'] = 239557417
		simulate['use_ac'] = 'yes'
		simulate['ra'] = 266.4
		simulate['dec'] = -28.9
		simulate['radius'] = 11
		simulate['irfs'] = 'P7REP_CLEAN_V15::FRONT'
		simulate['seed'] = seed
		simulate.run()

#################################################################
#
# Iteration Class(Analysis)
# Customization: Define gtlike to run in an iterative manner
#

class Iteration(Analysis):

	def ff_source(self, source, ff):
		
		ff = str(ff)
		spectrum = source[0]
		if spectrum.attrib['type'] == 'PowerLaw':
			spectrum[0].attrib['free'] = ff
			spectrum[1].attrib['free'] = ff

		if spectrum.attrib['type'] == 'LogParabola':
			spectrum[0].attrib['free'] = ff
			spectrum[1].attrib['free'] = ff
			spectrum[2].attrib['free'] = ff
			spectrum[3].attrib['free'] = '0'
			return
		if spectrum.attrib['type'] == 'PowerLaw2':
			spectrum[0].attrib['free'] = ff
			spectrum[1].attrib['free'] = ff
			return
		if spectrum.attrib['type'] == 'ConstantValue':
			return
		if spectrum.attrib['type'] == 'PLSuperExpCutoff':
			spectrum[0].attrib['free'] = ff
			spectrum[1].attrib['free'] = ff
			spectrum[2].attrib['free'] = '0'
			spectrum[3].attrib['free'] = ff
			spectrum[4].attrib['free'] = '0'
                        return
	
	def ff_by_name(self, root, name_list, ff):
		for source in root:
			if source.attrib['name'] in name_list:
			 	self.ff_source(source, str(ff))
			else: self.ff_source(source, str(int(not bool(ff))))
		return	


	def iterate(self, starting_xml_file, src_file, mass):
        	
		readfile = starting_xml_file
                pass_count = 0
		for sublist in self.TS:	
			pass_count += 1 #iteration number
                	tree = ET.parse(readfile)
                	root = tree.getroot()
                	self.ff_by_name(root, sublist, 1) #frees point sources in sublist, and fixes all other point sources
		
			#The following components must be reset to their original value following each iteration:
			for source in root:
				if source.attrib['name'] == 'H2pi0R1':
					source[0][0].attrib['value'] = '1'
			
				if source.attrib['name'] == "HIpi0R1":
					source[0][0].attrib['value'] = '1'
			
				if source.attrib['name'] == "ICR1":
					source[0][0].attrib['value'] = '1'
			
				if source.attrib['name'] == "W28":
					source[0][0].attrib['value'] = "1.095921774"
					source[0][1].attrib['value'] = "2.748630021"
					source[0][2].attrib['value'] = "0.346455084"
			
				if source.attrib['name'] == "Brem":
					source[0][0].attrib['value'] = '1'

				if source.attrib['name'] == "HIIpi0R1to17":
					source[0][0].attrib['value'] = '1'
				
				if source.attrib['name'] == "H2pi0R14to17":
					source[0][0].attrib['value'] = '1'


			        if source.attrib['name'] == "H2pi0R13" or source.attrib['name'] == "H2pi0A5" or source.attrib['name'] == "H2pi0A6toA8" \
				or source.attrib['name'] == "H2pi0A6" or source.attrib['name'] == "H2pi0A7" or source.attrib['name'] == "H2pi0A8":
					source[0][0].attrib['value'] = '1'	
					
				if source.attrib['name'] == "HIpi0R13" or source.attrib['name'] == "HIpi0A5" or source.attrib['name'] == "HIpi0A6" \
				or source.attrib['name'] == "HIpi0A7" or source.attrib['name'] == "HIpi0A8":
					source[0][0].attrib['value'] = '1'
				
					
				#if source.attrib['name'] == "HIpi0A5":
				#	source[0][0].attrib['value'] = '0.65'

	
				if source.attrib['name'] == "HIpi0R14to17" or source.attrib['name'] == "HIpi0R14to17_BPL_constrained":
					source[0][0].attrib['value'] = '1'

			        if source.attrib['name'] == "ICR13" or source.attrib['name'] == "ICA5" or source.attrib['name'] == "ICA6" \
				or source.attrib['name'] == "ICA7" or source.attrib['name'] == "ICA8" or source.attrib['name'] == "ICA6-A7":
			       		source[0][0].attrib['value'] = '1'
				
				#if source.attrib['name'] == "ICR14to17":
				#	source[0][0].attrib['value'] = '1'
				
				if source.attrib['name'] == "EGBfree" or source.attrib['name'] == "EGBfree_fixed":				
					source[0][0].attrib['value'] = '1.0'

				if source.attrib['name'] == "H2pi0R13_Index" or source.attrib['name'] == "Brem_Index" or source.attrib['name'] == "HIIpi0R1to17_Index":
					source[0][0].attrib['value'] = '1'	
				        source[0][1].attrib['value'] = '2'
					source[0][0].attrib['free'] = '1'
					source[0][1].attrib['free'] = '1'

				if source.attrib['name'] == "HIpi0R13_Index" or source.attrib['name'] == "HIpi0A5_Index" or source.attrib['name'] == "HIpi0A6_Index" \
					or  source.attrib['name'] == "HIpi0A7_Index" or source.attrib['name'] == "HIpi0A8_Index":
					source[0][0].attrib['value'] = '1'
					#source[0][1].attrib['value'] = '0'
					source[0][0].attrib['free'] = '1'
					#source[0][1].attrib['free'] = '1'

				if source.attrib['name'] == "HIpi0R14to17_Index":
					source[0][0].attrib['value'] = '1'
                                        #source[0][1].attrib['value'] = '0'                       
					source[0][0].attrib['free'] = '1'
					#source[0][1].attrib['free'] = '1'

				if source.attrib['name'] == "ICR13_Index" or source.attrib['name'] == "ICA5_Index" or source.attrib['name'] == "ICA6-A7_Index" or source.attrib['name'] == "ICA8_Index":
					source[0][0].attrib['value'] = '1'
                                        source[0][1].attrib['value'] = '0'                       
					source[0][0].attrib['free'] = '1'
					source[0][1].attrib['free'] = '1'

				if source.attrib['name'] == "ICR14to17_Index":
					source[0][0].attrib['value'] = '1'
                                        source[0][1].attrib['value'] = '2'                       
					source[0][0].attrib['free'] = '1'
					source[0][1].attrib['free'] = '1'
				
				if source.attrib['name'] == "H2pi0R13_BPL" or source.attrib['name'] == "HIpi0R13_BPL"\
				or source.attrib['name'] == "HIpi0R14to17_BPL" or source.attrib['name'] == "ICR13_BPL" \
				or source.attrib['name'] == "ICR14to17_BPL" or source.attrib['name'] =="HII_CR_BPL" \
				or source.attrib['name'] == "HIIpi0R1to17_BPL" or source.attrib['name'] == "BPL_DM_R1":
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '0'
					source[0][2].attrib['value'] = '0'
					source[0][3].attrib['value'] = '2'


		         	if source.attrib['name'] == "Disk Component(NFW_A)":
					source[0][0].attrib['value'] = "1" 
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(NFW)":
					source[0][0].attrib['value'] = "1" 
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(Gas)":
					source[0][0].attrib['value'] = "1"
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(Gauss1)":
					source[0][0].attrib['value'] = "1"
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(Gauss2)":
					source[0][0].attrib['value'] = "1"
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(Gauss5)":
					source[0][0].attrib['value'] = "1"
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == "Disk Component(Gauss10)":
					source[0][0].attrib['value'] = "1"
					source[0][1].attrib['value'] = "-2."
					source[0][3].attrib['value'] = "3000."
			
				if source.attrib['name'] == 'up_type':
					source[0][0].attrib['value'] = '1'
				if source.attrib['name'] == 'down_type':
					source[0][0].attrib['value'] = '1'
				if source.attrib['name'] == 'leptons':
					source[0][0].attrib['value'] = '1'

			
				if source.attrib['name'] == 'ZZp':
					source[0][0].attrib['value'] = '1'
				if source.attrib['name'] == 'WWp':
					source[0][0].attrib['value'] = '1'
				if source.attrib['name'] == 'Zgp':
					source[0][0].attrib['value'] = '1'
				if source.attrib['name'] == 'hh':
					source[0][0].attrib['value'] = '1'

				if source.attrib['name'] == 'gll_iem_v06':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '0'
					source[0][0].attrib['free'] = '1'
					source[0][1].attrib['free'] = '0'
				
				if source.attrib['name'] == 'iso_P8R2_SOURCE_V6_v06':
					source[0][0].attrib['value'] = '1'
			
				if source.attrib['name'] == 'iso_P8R2_CLEAN_V6_v06':
					source[0][0].attrib['value'] = '1'
					
				if source.attrib['name'] == 'iso_P8R2_ULTRACLEANVETO_V6_v06':
                                        source[0][0].attrib['value'] = '1'
	
			#	if source.attrib['name'] == '3FGL J0042.5+4117':
			#		source[0][0].attrib['free'] = '1'
				#	source[0][1].attrib['free'] = '1'
				#	source[0][0].attrib['value'] = '1'
				#	source[0][1].attrib['value'] = '2'
			
			        if source.attrib['name'] == 'DM_R1':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
			    	        source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'

			     	if source.attrib['name'] == 'DM_R2':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
					
				if source.attrib['name'] == 'DM_R3' or source.attrib['name'] == 'DM_3_4_N' or source.attrib['name'] == 'DM_R3S'\
					or source.attrib['name'] == "DM_3_4_Q1" or source.attrib['name'] == "DM_3_4_Q2"\
					or source.attrib['name'] == "DM_3_4_Q3" or source.attrib['name'] == "DM_3_4_Q4":
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
						
				if source.attrib['name'] == 'DM_R4' or source.attrib['name'] == 'DM_3_4_S' or source.attrib['name'] == 'DM_R4S':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
						
				if source.attrib['name'] == 'DM_R5' or source.attrib['name'] == 'DM_R5N' or source.attrib['name'] == 'DM_R5S':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
						
				if source.attrib['name'] == 'DM_R6' or source.attrib['name'] == 'DM_R6N' or source.attrib['name'] == 'DM_R6S':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
						
				if source.attrib['name'] == 'DM_R7':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][3].attrib['value'] = '10000'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					source[0][3].attrib['free'] = '1'
				
				if source.attrib['name'] == 'PL_DM_R1':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'

				if source.attrib['name'] == 'PL_DM_R2':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
					
				if source.attrib['name'] == 'PL_DM_R3':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
						
				if source.attrib['name'] == 'PL_DM_R4':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
						
				if source.attrib['name'] == 'PL_DM_R5':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
						
				if source.attrib['name'] == 'PL_DM_R6':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
						
				if source.attrib['name'] == 'PL_DM_R3_Q1':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '2'	
					source[0][0].attrib['free'] = '1'		
					source[0][1].attrib['free'] = '1'
			
				if source.attrib['name'] == 'Ring_1' or source.attrib['name'] == 'Ring_2' \
				or source.attrib['name'] == 'Ring_3' or source.attrib['name'] == 'Ring_4' \
				or source.attrib['name'] == 'Ring_5' or source.attrib['name'] == 'Ring_6' \
				or source.attrib['name'] == 'Ring_3_Q1' or source.attrib['name'] == 'Ring_3_Q2' \
				or source.attrib['name'] == 'Ring_3_Q3' or source.attrib['name'] == 'Ring_3_Q4'\
		   		or source.attrib['name'] == 'Ring_4_Q1' or source.attrib['name'] == 'Ring_4_Q2'\
				or source.attrib['name'] == 'Ring_4_Q3' or source.attrib['name'] == 'Ring_4_Q4':
					source[0][0].attrib['value'] = '30'
					source[0][1].attrib['value'] = '2'
					source[0][2].attrib['value'] = '-3'						
					source[0][3].attrib['value'] = '-3'
					source[0][4].attrib['value'] = '-3'
					source[0][5].attrib['value'] = '-3'
					source[0][6].attrib['value'] = '-3'
					source[0][7].attrib['value'] = '-3'
					source[0][8].attrib['value'] = '-3'
					source[0][9].attrib['value'] = '-3'
					source[0][10].attrib['value'] = '-3'
					source[0][11].attrib['value'] = '-3'
					source[0][12].attrib['value'] = '-3'
					source[0][13].attrib['value'] = '-3'
					source[0][14].attrib['value'] = '-3'
					source[0][15].attrib['value'] = '-3'			

				if source.attrib['name'] == 'Modified_Isotropic':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '-3'
					source[0][2].attrib['value'] = '-3'						
					source[0][3].attrib['value'] = '-3'
					source[0][4].attrib['value'] = '-3'
					source[0][5].attrib['value'] = '-3'
					source[0][6].attrib['value'] = '-3'
					source[0][7].attrib['value'] = '-3'
					source[0][8].attrib['value'] = '-3'
					source[0][9].attrib['value'] = '-3'
  			     		source[0][10].attrib['value'] = '-3'	
						
				if source.attrib['name'] == 'Limited_Ring_1' or source.attrib['name'] == 'Limited_Ring_2' \
				or source.attrib['name'] == 'Limited_Ring_3' or source.attrib['name'] == 'Limited_Ring_4':
					source[0][0].attrib['value'] = '1'
					source[0][1].attrib['value'] = '-3'
					source[0][2].attrib['value'] = '-3'						
					source[0][3].attrib['value'] = '-3'
					source[0][4].attrib['value'] = '-3'
					source[0][5].attrib['value'] = '-3'
					source[0][6].attrib['value'] = '-3'
					source[0][7].attrib['value'] = '-3'
					source[0][8].attrib['value'] = '-3'
					source[0][9].attrib['value'] = '-3'
					source[0][10].attrib['value'] = '-3'
					source[0][11].attrib['value'] = '-3'
	
				if source.attrib['name'] == 'P_Index_IEM':
					source[0][0].attrib['free'] = '1'
					source[0][0].attrib['value'] = '1'
				

				if source.attrib['name'] == 'Ring_1_fixed_index' or source.attrib['name'] == 'Ring_2_fixed_index' \
				or source.attrib['name'] == 'Ring_3_fixed_index' or source.attrib['name'] == 'Ring_4_fixed_index':
					source[0][0].attrib['value'] = '1'
					source[0][0].attrib['free'] = '1'
					
			
			xml_infile = self.name +  '_pass' + str(pass_count) + '_in.xml'
                	tree.write(xml_infile)
			xml_outfile = self.name + '_pass' + str(pass_count) + '_out.xml'
			self.run_gtlike(src_file, xml_infile, xml_outfile, mass, pass_count)
			readfile = xml_outfile #this becomes infile for next iteration

#########################################################
#
# PyLikelihood class
# Purpose: Run likelihood using pyLikelihood

class pyLikelihood(Iteration, Analysis):

	def run_gtlike(self, src_file, xml_infile, xml_outfile, mass, pass_count):
        	
		
		my_expCube = self.path + self.name + '_binned_ltcube.fits'
		my_ExpMap = self.path + self.name + '_binned_expMap.fits'
		
		obs = BinnedObs(srcMaps=src_file, expCube=my_expCube, binnedExpMap=my_ExpMap,irfs=self.irfs)
		like = BinnedAnalysis(obs, xml_infile, optimizer='NEWMINUIT')
		likeobj = pyLike.NewMinuit(like.logLike)
		z = like.fit(verbosity=0,covar=True,optObject=likeobj)
		like.logLike.writeXml(xml_outfile)
		
		covariance_name = 'covariance_' + str(pass_count) + '.txt'
		covariance = like.covariance	
		f = open(covariance_name, 'w')
		f.write(str(covariance))
		f.close()
	
		#Uncomment to record PS flux after each iteratin:	
#		PS_flux_list = []
#		for each in self.TS[pass_count - 1]:		
#			source_name = each
#			print source_name
#			flux = like.flux(source_name, 1000,100000,energyFlux=False)#units of ph/cm^2/s
#			flux_error = like.fluxError(source_name,1000,100000,energyFlux=False)
#			PS_flux_list.append([source_name,flux,flux_error])
#		
#		PS_flux_file = 'PS_flux_pass_' + str(pass_count) + '.txt'
#		f = open(PS_flux_file,'w')
#		f.write(str(PS_flux_list)) 
#		f.close()


		return 
	
					

###########################################################
#
# Plots class
# Purpose: Plot residuals and counts
#

class Plots(Analysis, Iteration):

	def plot_main(self, x_list, y_list, savefig, Title, x_label, y_label):
		
		fontsize = 16
		axwidth = 3
		fontsize = 18
		plt.rc('axes',linewidth=axwidth)
		font = {'family': 'monospace','weight':'bold','size':14}
		plt.rc('font',**font)	
		fig = plt.figure(figsize=(13,13))
		ax = fig.add_subplot(1,1,1)
	
		plt.plot(x_list,y_list,ls='',marker='o')
		ax.set_title(Title,y=1.04,fontsize=30,fontweight='bold')
		plt.xlabel(x_label,fontsize=25,fontweight='bold')
		plt.ylabel(y_label,fontsize=25,fontweight='bold')

		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		plt.savefig(savefig)
		plt.show()
	
	def Counts_spectra_plots(self):


		#The first extension of counts_spectra.fits contains the counts spectrum data, 
		#which is a binnary table hdu object. The columns are the sources, and the rows are the energy bins. Here we define the column names:
		
		hdulist = pyfits.open(self.file)
		cols = hdulist[1].columns
		sources = cols.names

		counts_data = Table.read(self.file, format = 'fits', hdu = 1) 

		#make an array for the observed counts
		obs_counts = np.array(counts_data[cols.names[0]])
		obs_err = obs_counts ** (0.5)

		#make an array for the model counts by summing over all the sources
		model_counts = np.zeros(20)
		for s in sources[1:]:
			model_counts = model_counts + np.array(counts_data[s])

		#calculate the residual plot:
		Residual = (obs_counts - model_counts)/model_counts
		res_err = obs_err/model_counts

		#You will plot residuals versus energy. You must determine which point within each energy bin to use. Here we choose the geometric mean:
		bin_data = Table.read(self.file, format = 'fits', hdu = 3)
		bin_min = np.array(bin_data['E_MIN'])
		bin_max = np.array(bin_data['E_MAX'])
		bins = (bin_min * bin_max)**(0.5)

		#plot the residuals versus energy
		plt.semilogx(bins, Residual, marker = 'o', ls = '', color = 'Black') 
		plt.errorbar(bins, Residual, yerr = res_err, ls = '', color = 'Black')
		plt.axhline(y = 0, ls = 'dotted', color = 'Black')
		plt.title('%s Total Residual' % self.name )
		plt.xlabel('Energy (Mev)')
		plt.axis([1000,100000,-0.4,0.6])
		plt.ylabel('(Counts - Model) / Model')
		plt.savefig('Residual.png')
		plt.close()

		#Now to plot counts for individual sources

		iso_counts = np.array(counts_data['EGBfree'])
		H2_counts = np.array(counts_data['H2pi0R1'])
		HI_counts = np.array(counts_data['HIpi0R1'])
		IC_counts = np.array(counts_data['ICR1']) 

		# get model and "observed" counts for point source at galactic center
		gc_source = 'PGW_0108'
		gc_mod = np.array(counts_data[gc_source])
		gc_obs = obs_counts * (gc_mod/model_counts)
		gc_err = obs_counts ** (0.5)

		# plot count and model arrays, and error bars
		plt.errorbar(bins, obs_counts, yerr = obs_err, label = 'Total Observed', marker = '.', ls = '')
		plt.loglog(bins, model_counts, label = 'Total Model', ls = 'solid')
		plt.loglog(bins, iso_counts, label = 'EGB', ls = 'dotted')
		plt.loglog(bins, H2_counts, label = 'H2, R1', ls = 'dashed')
		plt.loglog(bins, HI_counts, label = 'HI, R1', ls = 'dashed')
		plt.loglog(bins, IC_counts, label = 'IC, R1', ls = 'dashed')
		plt.errorbar(bins, gc_obs, yerr = gc_err, label = 'GC Observed', marker = '.', ls = '')
		plt.loglog(bins, gc_mod, label = 'GC Model', ls = 'dashdot')
		plt.title('%s Total Counts' % self.name)
		plt.xlabel('Energy (MeV)')
		plt.ylabel('Total counts')
		plt.legend(bbox_to_anchor=(0, 0), loc=3, fontsize=12)
		plt.savefig('counts.png')
		plt.close()
	
			

	def calculate_flux(self,model_file,exp_file,name,low,high):
			

			
			print 'Calculating flux for ' + model_file + '...'			
			
			w_model = WCS(model_file)
			w_exposure = WCS(exp_file)

			hdulist_A = pyfits.open(model_file)
			tbdata_A = hdulist_A[0].data 
			size = tbdata_A.shape
			len_x = size[2]
			len_y = size[1]
			len_E = size[0]
			energy_A = hdulist_A[2].data #for model cube
			if 'ccube' in model_file:
				energy_A = hdulist_A[1].data #for ccube
			hdulist_B = pyfits.open(exp_file)
			tbdata_B = hdulist_B[0].data
			flux = np.zeros(len_E)
			counts_flux = np.zeros(len_E)
			raw_counts = np.zeros(len_E)
			energy_list_plot = np.zeros(len_E)
			dN_dE = np.zeros(len_E)
			x_error = []
			for E in range(0,len_E):
				print 'Summing over x and y for energy bin: ' + str(E) 
  				total_counts = 0
				total_time = 0
				for x in range(0,len_x):
					for y in range(0,len_y):
						#if low < (x-69)**2 + (y-69)**2 <= high:
						#if tbdata_A[E][y][x] > 0:	
							model_x, model_y = w_model.wcs_pix2world(x,y,0,0)[:-1]
							exposure_x, exposure_y = w_exposure.wcs_world2pix(model_x,model_y,0,0)[:-1]
							
							delta_E = (energy_A[E]['E_MAX']/1000)-(energy_A[E]['E_MIN']/1000) #for model cube
							energy = ((energy_A[E]['E_MAX']/1000)*(energy_A[E]['E_MIN']/1000))**0.5 #for model cube
							#delta_E = (energy_A[E][2]/1000)-(energy_A[E][1]/1000) #for ccube
							#energy = ((energy_A[E][2]/1000)*(energy_A[E][1]/1000))**0.5 #for ccube
							this_energy = energy_A[E]['E_MIN']/1000 + delta_E #for writing ascii file
							energy_list_plot[E] = this_energy 
							counts = tbdata_A[E][y][x]
							time = (tbdata_B[E][int(exposure_y)][int(exposure_x)]*tbdata_B[E+1][int(exposure_y)][int(exposure_x)])**0.5
							flux[E] = flux[E] + (counts/delta_E)*(energy**2)/time
							dN_dE[E] = dN_dE[E] + (counts/delta_E)/time
							counts_flux[E] = counts_flux[E] + counts/time
							raw_counts[E] = raw_counts[E] + counts
				#print 
				#print 'Flux for this energy bin [MeV / cm^2 / s]'
				#print flux[E]	
				#print		
				x_error.append(delta_E)
			print
			print 'counts flux:'
			print counts_flux
			print

			print 
			print 'raw counts:'
			print raw_counts
			print 

			print
			print  'flux list:'
			print flux
			print
			print 'energy error:'
			print x_error
			print
			
			total_flux = 0
			for each in flux:
				total_flux = total_flux + each
		
			print
			print 'total flux: '
			print total_flux
			print		

			total_counts_flux = 0
			for each in counts_flux:
				total_counts_flux = total_counts_flux + each
		
			print
			print 'total counts flux: '
			print total_counts_flux
			print		


			camma_list = []
			for each in flux:
				camma_list.append(each)
			f = open('ultimate_results.txt','a')
			f.write('\n\n')
			f.write(str(model_file) + ' = ' + str(camma_list))
			f.close()

	
			#use atropy to write an ascii file, i.e. .dat file:	
			
			line_1 = "Total counts flux [ph/cm^2/s]: " + str(total_counts_flux) 
			line_2 = "Total energy flux [MeV/cm^2/s]: " + str(total_flux) 

			data = Table({"Energy [MeV]":energy_list_plot,"energy_flux [MeV/cm^2/s]":flux,"counts_flux [ph/cm^2/s]":counts_flux, "differential_flux [ph/cm^2/s/MeV]":dN_dE},names=["Energy [MeV]","energy_flux [MeV/cm^2/s]","counts_flux [ph/cm^2/s]", "differential_flux [ph/cm^2/s/MeV]"])			
			data.meta['comments'] = [line_1,line_2]	
			file_name = name + '_flux.dat'
			ascii.write(data, file_name, delimiter="\t")
			


			return total_counts_flux	 
		
	def pylikelihood_flux(self, src, xml, source_name_list):		
		
		axwidth = 3
		axlength = 12
		fontsize=20
		plt.rc('axes',linewidth=axwidth)
		font = {'family': 'monospace','weight':'bold','size':14}
		plt.rc('font',**font)
		fig = plt.figure(figsize=(13,8))
		#gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
		#gs.update(hspace=0.0)
		#ax1  = plt.subplot(gs[0])
		ax1 = fig.add_subplot(1,1,1)

		my_expCube = self.path + self.name + '_binned_ltcube.fits'
		my_binnedExpMap = self.path + self.name + '_binned_expMap.fits'
		
	
		obs = BinnedObs(src, expCube=my_expCube, \
		binnedExpMap=my_binnedExpMap,irfs=self.irfs)
		like = BinnedAnalysis(obs, xml, optimizer='NEWMINUIT')
		z = like.fit(verbosity=0,covar=True)	
		
		print
		print 'covariance: ' + str(like.covariance)
		print
		covariance = like.covariance	
		f = open('pylikelihood_cov.txt', 'w')
		f.write(str(covariance))
		f.close()

		print
		print like.model
		print

		print
		print '-log(L): ' + str(z) 
		print
		
		f = open('likelihood.txt','a')
		f.write(str(z))
		f.close()

		like.logLike.writeXml('pylikelihood_fit.xml')
		#extract the energy intervals:		
		ccube = self.path + self.name + '_ccube.fits'
		hdulist = pyfits.open(ccube)
		my_energy = hdulist[1].data
		
		results_list = []
		for source_name in source_name_list:
			
			flux_list = []
			flux_error_list = []
			Normalized_flux_error_list = []
			energy_list = []
			energy_error_list = []
			counts_flux  = []	
			
			##################
			# Sanity Checks:
			print '##################'
			print "Here are some sanity checks calculated with like.flux"
			print "Use these values if working with energy flux"
			print 
			print 'source name:'
			print source_name
			print "total energy flux b/n 1 GeV - 100 GeV"			
			true_energy_flux = like.flux(source_name,1000,100000,energyFlux=True)
			true_energy_flux_error = like.fluxError(source_name,1000,100000,energyFlux=True)
			print true_energy_flux  
			print true_energy_flux_error
			print 
			print "total counts flux b/n 1 GeV - 100 GeV"
			print like.flux(source_name,1000,100000,energyFlux=False)
                        print like.fluxError(source_name,1000,100000,energyFlux=False)
			print 
			#print "total energy flux b/n 100 MeV - 100 GeV"                 
                        #print like.flux(source_name,100,100000,energyFlux=True)
			#print like.fluxError(source_name,100,100000,energyFlux=True)
			#print
                        #print "total counts flux b/n 100 MeV = 100 GeV"
                        #print like.flux(source_name,100,100000,energyFlux=False)
			#print like.fluxError(source_name,100,100000,energyFlux=False) 
			#print '###################'

			for each in my_energy:
				
				E_low = each[1] / 1000.0 #the energy is given in kev, but it needs to be in MeV
				E_high = each[2] / 1000.0
				delta_E = E_high - E_low
				energy = (E_low*E_high)**0.5
				Energy = (E_low + E_high) / 2.0
			
				
				counts =  like.flux(source_name,E_low,E_high,energyFlux=False) #this gives integrated value of [photons / (cm^2*s)]  
				flux = counts/delta_E*(energy**2)
				flux_error = like.fluxError(source_name,E_low,E_high,energyFlux=False)
				Normalized_flux_error = flux_error/delta_E*(energy**2)
				
				energy_list.append(float('{:.3f}'.format(Energy)))
				energy_error_list.append(float('{:.3f}'.format(delta_E / 2.0)))	
				flux_list.append(float('{:.3e}'.format(flux)))
				Normalized_flux_error_list.append(float('{:.3e}'.format(Normalized_flux_error)))
				counts_flux.append(float('{:.3e}'.format(counts)))
				flux_error_list.append(float('{:.3e}'.format(flux_error)))
		
			total_counts_error = 0
			total_counts_flux = 0
			total_flux = 0
			total_flux_error = 0
			for each in flux_list:
				total_flux = total_flux + each			
			for each in Normalized_flux_error_list:
				total_flux_error = total_flux_error + each
			for each in counts_flux:
				total_counts_flux = total_counts_flux + each
			for each in flux_error_list:
				total_counts_error = total_counts_error + each
			print 
			print str(source_name) + ':' 
			print
			print 'energy_list:' + str(energy_list) 
			print
			print 'energy_error_list: ' + str(energy_error_list)
			print
			print 'total flux: ' + str(true_energy_flux) + ' MeV/cm^2/s'
			print
			print 'total flux error: ' + str(true_energy_flux_error) + ' MeV/cm^2/s'
			print
			luminosity = true_energy_flux * 1/ (6.2415*10**5)
			luminosity_error = true_energy_flux_error * 1/ (6.2415*10**5)
			print 'total luminosity: ' + str(luminosity) + ' erg/cm^2/s'
			print
			print 'total luminosity error: ' + str(luminosity_error) + ' erg/cm^2/s'
			print
			print 'flux list: ' + str(flux_list) 
			print 
			print 'flux error list: ' + str(Normalized_flux_error_list)
			print
			print 'total counts flux: ' +  str(total_counts_flux) + ' ph/cm^2/s'
			print
			print 'total counts flux error: ' + str(total_counts_error) + ' ph/cm^2/s'
			print
			print
                        print 'total counts intensity: ' +  str(total_counts_flux/0.2352) + ' ph/cm^2/s/sr'
                        print
                        print 'total counts intensity error: ' + str(total_counts_error/0.2352) + ' ph/cm^2/s/sr'
                        print
			results_list.append([source_name, total_flux, total_flux_error])
			ax1.loglog(energy_list,flux_list, label = source_name)
			#ax1.errorbar(energy_list,flux_list,yerr=Normalized_flux_error_list,color='black',lw=3,markeredgewidth=3,capsize=6)	


			##############################################################
			# use pandas to create data frame and save output as .csv file
		 
			#use atropy to write an ascii file, i.e. .dat file:	
			
			line_1 = "Total counts flux [ph/cm^2/s]: " + str(total_counts_flux) + ' +/- ' + str(total_counts_error)
			line_2 = "Total energy flux [MeV/cm^2/s]: " + str(true_energy_flux) + ' +/- ' + str(true_energy_flux_error)
			line_3 = "Total luminosity [erg/cm^2/s]: " + str(luminosity) + ' +/- ' + str(luminosity_error)
			data = Table({"Energy [MeV]":energy_list,"energy_flux [MeV/cm^2/s]":flux_list,"energy_flux_error [MeV/cm^2/s]":Normalized_flux_error_list,"counts_flux [ph/cm^2/s]":counts_flux,"counts_flux_error [ph/cm^2/s]":flux_error_list},names=["Energy [MeV]","energy_flux [MeV/cm^2/s]", "energy_flux_error [MeV/cm^2/s]","counts_flux [ph/cm^2/s]", "counts_flux_error [ph/cm^2/s]"])			
			data.meta['comments'] = [line_1,line_2,line_3]	
			file_name = str(source_name) + '_flux.dat'
			ascii.write(data, file_name, delimiter="\t")
			
			#use pandas to write a dataframe, i.e. .csv file:
			#data = {"Energy [MeV]":energy_list,"Total energy flux":flux_list,"Total energy flux error":total_flux_error}
			#df = pd.DataFrame(data)
			#file_name = str(source_name) + '_flux.csv'
			#df.to_csv(file_name, sep='/t')
		
		plt.legend(loc=1, frameon=False,ncol=3)	
#		plt.axhline(0.0,ls=':')
		plt.xlabel('Energy (MeV)',fontsize=25,fontweight='bold')
		plt.ylabel('E$\mathregular{^{2}}$ dN/dE [MeV cm$\mathregular{^{-2}}$s$\mathregular{^{-1}}$]',fontsize=25,fontweight='bold')
		#plt.title('M31 ($\mathregular{30^\circ}$ ROI)',fontsize=32,fontweight='bold',y=1.04)
#		#plt.xlim((1000,100000))
		plt.ylim((1e-7,5e-4))
		
		for tick in ax1.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for tick in ax1.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for line in ax1.get_xticklines() + ax1.get_yticklines():
			line.set_markersize(18)
			line.set_markeredgewidth(3)
		for tick in ax1.yaxis.get_minor_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		ax1.tick_params(which='minor',width=axwidth,length=axlength)	
		#plt.savefig('flux.png',bbox_inches='tight')
		plt.show()
		plt.close()
		

		
	def pyLikelihood_counts(self, input_list, source_plot_list):

	 	
		axwidth = 3
		axlength = 12
		fontsize=20
		plt.rc('axes',linewidth=axwidth)
		font = {'family': 'monospace','weight':'bold','size':14}
		plt.rc('font',**font)
		fig = plt.figure(figsize=(13,8))
		ax = plt.gca()

		my_expCube = self.path + self.name + '_binned_ltcube.fits'
		my_binnedExpMap = self.path + self.name + '_binned_expMap.fits'

		for sublist in input_list:
			src = sublist[0]
			xml = sublist[1]
			
			obs = BinnedObs(src, expCube=my_expCube, \
			binnedExpMap=my_binnedExpMap,irfs=self.irfs)
			like = BinnedAnalysis(obs, xml, optimizer='NEWMINUIT')

			E = (like.energies[:-1] + like.energies[1:])/2

			#sum_model = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
			#dark_matter = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
			
			for sourceName in source_plot_list:

				plt.loglog(E,like._srcCnts(sourceName), label = sourceName)
			#for sourceName in like.sourceNames():
			#	sum_model = sum_model + like._srcCnts(sourceName)
		#	plt.loglog(E,sum_model, label = 'sum_model')
			
		#	print 
		#	print 'Sum Model:'
		#	print sum_model
		#	print


			#Fermi_counts = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
			#for sourceName in like.sourceNames():
			# 	if '3FGL' in sourceName:
		#			Fermi_counts = Fermi_counts + like._srcCnts(sourceName)
		#	print 
		#	print 'Fermi Counts:'
		#	print Fermi_counts
		#	print
		
		#	NPS_counts = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
	#		for sourceName in like.sourceNames():
	#		 	if 'PS' in sourceName:
	#				NPS_counts = NPS_counts + like._srcCnts(sourceName)
	#		print 
	#		print 'NPS Counts:'
	#		print NPS_counts
	#		print
			

	
			for sourceName in source_plot_list:
				counts_list = like._srcCnts(sourceName)
				total = 0 
				for each in counts_list:
					total = total + each
				print
				print 'name: ' + str(sourceName)
				print 'counts: ' + str(total)
				print 'counts list:' 
				print counts_list
				print 
			#for sourceName in ['up_type','down_type','leptons']:
			#	dark_matter = dark_matter + like._srcCnts(sourceName)
			#plt.loglog(E,dark_matter, label = 'total_dark_matter')
		#plt.loglog(E,like.nobs)	
		#plt.errorbar(E,like.nobs,yerr=np.sqrt(like.nobs), fmt='o',label='Counts') 
		plt.legend(ncol=3, loc=2, frameon=False)
		#plt.axvline(x=50118.72,ls='--')
		plt.ylim((0.1,5000))	
		plt.xlim((1e3,1e5))
		plt.xlabel('Energy (MeV)',fontsize=25,fontweight='bold')
		plt.ylabel('Events',fontsize=25,fontweight='bold')
		plt.title('Counts Spectrum (standard 3FGL)',fontsize=32,fontweight='bold')
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for line in ax.get_xticklines() + ax.get_yticklines():
			line.set_markersize(18)
			line.set_markeredgewidth(3)
		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		ax.tick_params(which='minor',width=axwidth,length=axlength)	
		plt.savefig('%s_counts.png' % self.name)
		plt.show()
		plt.close()	
	

	def pyLikelihood_fractional_residuals(self,src_file,xml_file):
		


		axwidth = 3
		axlength = 12
		fontsize=20
		plt.rc('axes',linewidth=axwidth)
		font = {'family': 'monospace','weight':'bold','size':14}
		plt.rc('font',**font)
		fig = plt.figure(figsize=(13,8))
		ax = plt.gca()
		
		my_expCube = "ltcube_00.fits"
		my_binnedExpMap = "bexpmap_00.fits"
		#my_expCube = self.path + self.name + '_binned_ltcube.fits'
		#my_binnedExpMap = self.path + self.name + '_binned_expMap.fits'		
		src = src_file
		xml = xml_file
			

		obs = BinnedObs(src, expCube=my_expCube, \
		binnedExpMap=my_binnedExpMap,irfs=self.irfs)
		like = BinnedAnalysis(obs, xml, optimizer='NEWMINUIT')
		E = (like.energies[:-1] + like.energies[1:])/2

		sum_model = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
		for sourceName in like.sourceNames():
			#print sourceName
			sum_model = sum_model + like._srcCnts(sourceName)	
		
		
		resid = (like.nobs - sum_model)/sum_model
		resid_err = (np.sqrt(like.nobs)/sum_model)
		plt.semilogx(E, resid, marker = 'o', ms=8,markeredgecolor = 'none',ls = '', color = 'black')
		plt.errorbar(E, resid, yerr = resid_err, ls = '', color = 'black', label='_nolegend_')
		
		
		print
		print 'Energy:'
		print E
		print 
		print 'Residuals:'
		print resid
		print 
		print 'Residual Error:'
		print resid_err
		print

	
	
		#use atropy to write an ascii file, i.e. .dat file:	
				
		data = Table({"Energy [MeV]":E,"Residuals":resid,"Residual Error":resid_err},names=["Energy [MeV]", "Residuals", "Residual Error"])				
		file_name = 'fractional_count_residuals.dat'
		ascii.write(data, file_name, delimiter="\t")
			
		#plt.legend(bbox_to_anchor=(0.02, 0.95), loc=2, frameon=False)	
		plt.axhline(0.0,ls=':')
		plt.xlabel('Energy (MeV)',fontsize=25,fontweight='bold')
		plt.ylabel('(Data-Model)/Model',fontsize=25,fontweight='bold')
		#plt.title('Pseudoscalar Interactions, NFW $\mathregular{\gamma = 1.2}$',fontsize=32,fontweight='bold',y=1.04)
		#plt.xlim((1e3,1e5))
		plt.ylim((-0.1,0.1))
		#plt.axhline(0.05,ls='--')
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		for line in ax.get_xticklines() + ax.get_yticklines():
			line.set_markersize(18)
			line.set_markeredgewidth(3)
		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontsize(fontsize)
			tick.label1.set_fontweight('bold')
		ax.tick_params(which='minor',width=axwidth,length=axlength)	
		plt.savefig('%s_Residuals.png' % self.name )
		plt.show()
		plt.close()
		
