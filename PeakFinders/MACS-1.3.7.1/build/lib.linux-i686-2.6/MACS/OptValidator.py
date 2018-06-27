# Time-stamp: <2009-10-16 14:02:04 Tao Liu>

"""Module Description

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import os
import re
import logging
from math import log
from MACS.IO import BEDParser, ELANDResultParser, ELANDMultiParser, PairEndELANDMultiParser, SAMParser, BAMParser, BowtieParser
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def opt_validate ( optparser ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()

    # treatment file
    if not options.tfile:       # only required argument
        optparser.print_help()
        sys.exit(1)

    # format

    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.build = ELANDResultParser().build_fwtrack
    elif options.format == "BED":
        options.build = BEDParser().build_fwtrack
    elif options.format == "ELANDMULTI":
        options.build = ELANDMultiParser().build_fwtrack
    elif options.format == "ELANDMULTIPET":
        options.build = PairEndELANDMultiParser().build_fwtrack
    elif options.format == "SAM":
        options.build = SAMParser().build_fwtrack
    elif options.format == "BAM":
        options.build = BAMParser().build_fwtrack
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.build = BowtieParser().build_fwtrack
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # for ELANDMULTIPET format
    if options.format == "ELANDMULTIPET":
        # treatment files
        fs = options.tfile.split(',')
        if len(fs) != 2:
            logging.error("Only two filenames are acceptable! But you provided %d:'%s'" % (len(fs),options.tfile))
            sys.exit(1)
        options.tfile = fs
        if not os.path.isfile(options.tfile[0]) or not os.path.isfile(options.tfile[1]):
            logging.error("No such file: %s or %s!" % (options.tfile[0],options.tfile[1]))
            sys.exit(1)
        # input files
        if options.cfile:
            fs = options.cfile.split(',')
            if len(fs) != 2:
                logging.error("Only two filenames are acceptable! But you provided %d:'%s'" % (len(fs),options.cfile))
                sys.exit(1)
            options.cfile = fs
            if not os.path.isfile(options.cfile[0]) or not os.path.isfile(options.cfile[1]):
                logging.error("No such file: %s or %s!" % (options.cfile[0],options.cfile[1]))
                sys.exit(1)
    else:
        if not os.path.isfile (options.tfile):
            logging.error("No such file: %s!" % options.tfile)
            sys.exit(1)

        # input file
        if options.cfile and not os.path.isfile (options.cfile):
            logging.error("No such file: %s!" % options.cfile)
            sys.exit(1)

    # lambda set
    lambdaset_hit =  re.match("^(\d+)\,(\d+)\,(\d+)$",options.lambdaset)
    if lambdaset_hit:
        lambdaset = sorted(map(int,lambdaset_hit.groups()))
    else:
        logging.error("Lambda set string must be like \"1000,5000,10000\", which means three integers seperated by commas.")
        sys.exit(1)
    options.lambdaset = lambdaset # overwrite it


    # shiftsize>0
    if options.shiftsize <=0 :
        logging.error("--shiftsize must > 0!")
        sys.exit(1)

    # -10*log10 pvalue
    options.log_pvalue = log(options.pvalue,10)*-10
    
    # uppercase the format string 
    options.format = options.format.upper()
	
    # output filenames
    options.peakxls = options.name+"_peaks.xls"
    options.peakbed = options.name+"_peaks.bed"
    options.zwig_tr = options.name+"_treat_afterfiting"
    options.zwig_ctl= options.name+"_control_afterfiting"
    options.negxls  = options.name+"_negative_peaks.xls"
    options.diagxls = options.name+"_diag.xls"
    options.modelR  = options.name+"_model.r"

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
	
    options.error   = logging.critical		# function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # options argument text
    options.argtxt = "\n".join((
            "# ARGUMENTS LIST:",\
            "# name = %s" % (options.name),\
            "# format = %s" % (options.format),\
            "# ChIP-seq file = %s" % (options.tfile),\
            "# control file = %s" % (options.cfile),\
            "# effective genome size = %.2e" % (options.gsize),\
            "# tag size = %d" % (options.tsize),\
            "# band width = %d" % (options.bw),\
            "# model fold = %s" % (options.mfold),\
            "# pvalue cutoff = %.2e" % (options.pvalue),\
            "# Ranges for calculating regional lambda are : peak_region," + ",".join(map(str,options.lambdaset))))
    
    # wig file?
    subdir = options.name+"_MACS_wiggle"
    if options.store_wig:
        # check subdir
        if os.path.exists(subdir):
            options.error("./%s exists! Unable to create directory to store wiggle file!!" % (subdir))
            sys.exit(1)
        else:
            options.wig_dir_tr = subdir+"/treat/"
            options.wig_dir_ctl= subdir+"/control/"

    return options
