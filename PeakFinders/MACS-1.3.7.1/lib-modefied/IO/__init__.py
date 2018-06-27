# Time-stamp: <2009-10-07 16:46:20 Tao Liu>

"""Module for FWTrackI, BEDParser and ELANDParser classes

Copyright (c) 2007 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import re
import sys
import logging
import struct
from array import array
from random import sample as random_sample
from Constants import *
# ------------------------------------
# constants
# ------------------------------------
__version__ = "TabIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "FWTrackI, BEDParser and ELANDParser classes"

# ------------------------------------
# Misc functions
# ------------------------------------


# ------------------------------------
# Classes
# ------------------------------------

class PeakIO:
    """IO for peak information.

    """
    def __init__ (self):
        self.peaks = {}
    
    def add (self, chromosome, start, end, summit=None, 
             peak_height=None, number_tags=None, 
             pvalue=None, fold_enrichment=None, fdr=None):
        """items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       peak_height,number_tags,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, pvalue_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low and p[7]<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
    
    def ave_fold_enrichment (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        t = 0
        for chrom in chrs:
            x += len(peaks[chrom])
            for p in peaks[chrom]:
                t+=p[7]
        return t/x

    def max_fold_enrichment (self):
        """Return the maximum fc.
        
        """
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            if peaks[chrom]:
                m = max([i[7] for i in peaks[chrom]])
                if m>x:
                    x=m
        return x
        
    
    def tobed (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (chrom,peak[0],peak[1])
        return text

    def towig (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (peak[0],peak[1])
        return text
        
    def init_from_dict (self, data):
        """Initialize the data from a dictionary. Improper assignment
        will damage the data structure.
        
        """
        self.peaks = {}
        chrs = data.keys()
        chrs.sort()
        for chrom in chrs:
            self.peaks[chrom]=[]
            a = self.peaks[chrom].append
            for i in data[chrom]:
                a(i)

    def merge_overlap ( self ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in xrange(len(peaks_chr)):
                if not prev_peak:
                    prev_peak = peaks_chr[i]
                    continue
                else:
                    if peaks_chr[i][0] <= prev_peak[1]:
                        s_new_peak = prev_peak[0]
                        e_new_peak = peaks_chr[i][1]
                        l_new_peak = e_new_peak-s_new_peak
                        if peaks_chr[i][4] > prev_peak[4]:
                            summit_new_peak = peaks_chr[i][3]
                            h_new_peak = peaks_chr[i][4]
                        else:
                            summit_new_peak = prev_peak[3]
                            h_new_peak = prev_peak[4]
                        prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
                    else:
                        n_append(prev_peak)
                        prev_peak = peaks_chr[i]
            if prev_peak:
                n_append(prev_peak)
        del peaks
        self.peaks = new_peaks
        return True

    def overlap_with_other_peaks (self, peaks2, cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialzed as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = peaks1.keys()
        chrs2 = peaks2.keys()
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(peaks1[k])
            rl2_k = iter(peaks2[k])
            tmp_n = False
            try:
                r1 = rl1_k.next()
                r2 = rl2_k.next()
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = rl1_k.next()
                        tmp_n = False
                    else:
                        r2 = rl2_k.next()
            except StopIteration:
                continue
        return total_num
        

class FWTrackI:
    """Fixed Width Ranges along the whole genome (commonly with the
    same annotation type), which are stored in a dict.

    Ranges are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.

    Example:
       >>> tabfile = TabFile("tag.bed",format="bed",mode="r")
       >>> track = FWTrackI()
       >>> for (chromosome,rg) in tabfile:
       ...    track.add_range(chromosome,rg)
       >>> track.get_ranges_by_chr["chr1"] # all ranges in chr1 
    """
    def __init__ (self,fw=0,anno=""):
        """fw is the fixed-width for all ranges
        """
        self.fw = fw
        self.__ranges = {}
        self.__comments = {}    # store the comments if available for each range
        self.__well_merged = False
        self.total = 0					# total tags
        self.total_unique = 0		# total unique tags
        self.annotation = anno   # need to be figured out

    def add_loc (self, chromosome, fiveendpos, strand):
        """Add a range to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, neg for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__ranges.has_key(chromosome):
            self.__ranges[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
            self.__comments[chromosome] = [array(UBYTE2,[]),array(UBYTE2,[])] # for (+,-)
        self.__ranges[chromosome][strand].append(fiveendpos)
        self.__comments[chromosome][strand].append(1)
        self.total+=1

    def get_ranges_by_chr (self, chromosome):
        """Return array of locations by chromosome.

        Not recommanded! Use generate_rangeI_by_chr() instead.
        """
        if self.__ranges.has_key(chromosome):
            return self.__ranges[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_comments_by_chr (self, chromosome):
        """Return array of comments by chromosome.

        """
        if self.__comments.has_key(chromosome):
            return self.__comments[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__ranges.keys()
        l.sort()
        return l

    def length (self):
        """Total sequenced length = total number of tags * width of tag		
        """
        return self.total*self.fw

    def sort (self):
        """Naive sorting for tags. After sorting, comments are massed
        up.

        
        """
        for k in self.__ranges.keys():
            self.__ranges[k][0] = sorted(self.__ranges[k][0])
            self.__ranges[k][1] = sorted(self.__ranges[k][1])

    def merge_overlap (self):
        """merge the SAME ranges. Record the duplicate number in self.__comments{}
        
        *Note: different with the merge_overlap() in TrackI class,
        which merges the overlapped ranges.
        """
        self.total = 0
        self.total_unique = 0
        for k in self.__ranges.keys(): # for each chromosome
            # + strand
            plus = sorted(self.__ranges[k][0])
            if len(plus) <1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
                new_plus = []
                new_plus_c = []
            else:
                plus_c = self.__comments[k][0]
                (new_plus,new_plus_c) = (array(BYTE4,[plus[0]]),array(UBYTE2,[1]))
            
                pappend = new_plus.append
                pcappend = new_plus_c.append
                n = 0                # the position in new list
                for p in plus[1:]:
                    if p == new_plus[n]:
                        try:
                            new_plus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 65535 + strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_plus_c[n]=65535

                    else:
                        pappend(p)
                        pcappend(1)
                        n += 1
                self.total_unique +=  len(new_plus)
                self.total += sum(new_plus_c)
            # - strand
            minus = sorted(self.__ranges[k][1])
            if len(minus) <1:
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
                new_minus = []
                new_minus_c = []
            else:
                minus_c = self.__comments[k][0]
                (new_minus,new_minus_c) = (array(BYTE4,[minus[0]]),array(UBYTE2,[1]))
            
                mappend = new_minus.append
                mcappend = new_minus_c.append
                n = 0                # the position in new list
                for p in minus[1:]:
                    if p == new_minus[n]:
                        try:
                            new_minus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 65535 - strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_minus_c[n]=65535
                    else:
                        mappend(p)
                        mcappend(1)
                        n += 1
                self.total_unique +=  len(new_minus)
                self.total += sum(new_minus_c)

            self.__ranges[k]=[new_plus,new_minus]
            self.__comments[k]=[new_plus_c,new_minus_c]
            self.__well_merged = True
		
    def merge_plus_minus_ranges_w_duplicates (self):
        """Merge minus strand ranges to plus strand. The duplications
        on a single strand is erased. But if the same position is on
        both pos and neg strand, keep them both.
        
        Side effect: Reset the comments. self.total_unique is set to None.
        """
        self.total_unique = None
        self.total = 0
        for chrom in self.__ranges.keys():
            (plus_tags,minus_tags) = self.__ranges[chrom]
            new_plus_tags = array(BYTE4,[])
            #reset comments
            self.__comments[chrom][0] = array(UBYTE2,[])
            self.__comments[chrom][1] = array(UBYTE2,[])
            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    new_plus_tags.append(plus_tags[ip])
                    ip += 1
                else:
                    new_plus_tags.append(minus_tags[im])
                    im += 1
            if im < lenm:
                # add rest of minus tags
                new_plus_tags.extend(minus_tags[im:])
            if ip < lenp:
                # add rest of plus tags
                new_plus_tags.extend(plus_tags[ip:])

            self.__ranges[chrom] = [new_plus_tags,[]]
            self.total += len(new_plus_tags)

    def sample (self, percent):
        """Sample the tags for a given percentage.
        
        Side effect: self.total_unique is set to None, and comments are unset.
        """
        self.total = 0
        self.total_unique = None
        for key in self.__ranges.keys():
            num = int(len(self.__ranges[key][0])*percent)
            self.__ranges[key][0]=array(BYTE4,sorted(random_sample(self.__ranges[key][0],num)))
            num = int(len(self.__ranges[key][1])*percent)
            self.__ranges[key][1]=array(BYTE4,sorted(random_sample(self.__ranges[key][1],num)))
            self.total += len(self.__ranges[key][0]) + len(self.__ranges[key][1])
            self.__comments[key] = [[],[]]
            
    def __str__ (self):
        return self.__to_wiggle()
        
    def __to_wiggle (self):
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for key in self.__ranges.keys():
            t += "variableStep chrom=%s span=%d strand=0\n" % (key,self.fw)
            for i in self.__ranges[key][0]:
                t += str(i)+"\n"
            t += "variableStep chrom=%s span=%d strand=1\n" % (key,self.fw)
            for i in self.__ranges[key][1]:
                t += str(i)+"\n"
        return t


class TrackI:
    """Variable Width Ranges along the whole genome (commonly with the
    same annotation type), which are stored in a dict.

    Ranges are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    def __init__ (self,anno=""):
        """
        """
        self.__ranges = {}
        self.__comments = {}    # store the comments if available for each range
        self.__well_sorted = False
        self.__well_merged = False
        self.total = 0
        self.total_w_duplicates = 0
        self.annotation = anno   # need to be figured out

    def add_range (self, chromosome, range):
        """Add a range to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        range -- RangeI objectleverage
        """
        if not isinstance(range,RangeI):
            raise Exception("Must add a RangeI object!")
        if not self.__ranges.has_key(chromosome):
            self.__ranges[chromosome] = [[],[]] # for (+strand, -strand)
            self.__comments[chromosome] = [[],[]] # for (+,-)
        #print "add: "+str(range)
        if range.strand == 1: 
            self.__ranges[chromosome][0].append(range.start)
            self.__comments[chromosome][0].append(1)
        elif range.strand == -1: 
            self.__ranges[chromosome][1].append(range.end)
            self.__comments[chromosome][1].append(1)
        self.total+=1
        self.__well_sorted = False
        self.__well_merged = False
        
    def add_loc (self, chromosome, fiveendpos, strand):
        """Add a range to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, neg for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__ranges.has_key(chromosome):
            self.__ranges[chromosome] = [[],[]] # for (+strand, -strand)
            self.__comments[chromosome] = [[],[]] # for (+,-)
        #print "add: "+str(range)
        self.__ranges[chromosome][strand].append(fiveendpos)
        self.__comments[chromosome][strand].append(1)
        self.total+=1
        self.__well_sorted = False
        self.__well_merged = False

    def get_ranges_by_chr (self, chromosome):
        """Return array of locations by chromosome.

        Not recommanded! Use generate_rangeI_by_chr() instead.
        """
        if self.__ranges.has_key(chromosome):
            return self.__ranges[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_comments_by_chr (self, chromosome):
        """Return array of comments by chromosome.

        """
        if self.__comments.has_key(chromosome):
            return self.__comments[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))


    def generate_rangeI_by_chr_strand (self, chromosome,strand):
        """A generator to return RangeI objects on given chromosome and given strand.

        strand: 0 for plus, 1 for minus strand

        Recommended function.
        """
        if self.__ranges.has_key(chromosome):
            for i in self.__ranges[chromosome][strand]:
                yield RangeI(start=i,end=i+self.fw,strand=strand) # buggy!!
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__ranges.keys()
        l.sort()
        return l

    def length (self):
        return self.total*self.fw

    def intersect(self, frs2, verbose=False):
        """Calculate how many ranges in frs2(a TrackI object) are intersected by self. Self and frs2 must each be sorted and got rid of overlaps by merge_overlap() function.

        If verbose is set to True, return a tuple (number_of_hits, detail_intersection_info),
        else only return the number.
        """
        total_num = 0
        self_chrs = self.get_chr_names()
        frs2_chrs = frs2.get_chr_names()
        v=[]
        for k in self_chrs:
            if not frs2_chrs.count(k):
                continue
            self_rl_k = self.get_ranges_by_chr(k) 
            frs2_rl_k = frs2.get_ranges_by_chr(k)
            rself = self_rl_k.next()
            if verbose: t = k+":"+str(rself)
            rfrs2 = frs2_rl_k.next()
            try:
                while (True):
                    if rself.end < rfrs2.end:
                        # rself not include rfrs2
                        rself = self_rl_k.next()
                        if verbose:
                            if t.find("\t") != -1: 
                                v.append(t)
                            t=k+":"+str(rself)
                        #print "self"+str(rself)
                    else:       # here rself.end must > rfrs2.end
                        if rself.start <= rfrs2.start: 
                            # rself include rfrs2
                            total_num += 1
                            if verbose: t+="\t"+str(rfrs2)
                        rfrs2 = frs2_rl_k.next()
                        #print "frs2"+str(rfrs2)
            except StopIteration:
                if verbose:
                    if t.find("\t") != -1: 
                        v.append(t)
                continue
        if verbose:
            return (total_num,v)
        else:
            return total_num


    def merge_overlap (self):
       """merge the SAME ranges. Record the duplicate number in self.__comments{}

       *Note: different with the merge_overlap() in TrackI class, which merges the overlapped ranges. 

       """
       if not self.__well_sorted:
           self.sort()
       self.total = 0
       self.total_w_duplicates = 0
       for k in self.__ranges.keys(): # for each chromosome
           #logging.debug(" ... mergeing chr: %s" % (k))
           (plus,minus) = self.__ranges[k]
           (plus_c,minus_c) = self.__comments[k]
           if len(plus)<1:
               (new_plus,new_plus_c) = (plus,plus_c)
           else:
               (new_plus,new_plus_c) = ([plus[0]],[1])
               # + strand
               n = 0                # the position in new list
               for p in plus[1:]:
                   if p == new_plus[n]:
                       new_plus_c[n]+=1
                   else:
                       new_plus.append(p)
                       new_plus_c.append(1)
                       n += 1
               self.total +=  len(new_plus)
               self.total_w_duplicates += sum(new_plus_c)
           # - strand
           if len(minus)<1:
               (new_plus,new_plus_c) = (minus,minus_c)
           else:
               (new_minus,new_minus_c) = ([minus[0]],[1])
               n = 0                # the position in new list
               for p in minus[1:]:
                   if p == new_minus[n]:
                       new_minus_c[n]+=1
                   else:
                       new_minus.append(p)
                       new_minus_c.append(1)
                       n += 1
               self.total +=  len(new_minus)
               self.total_w_duplicates += sum(new_minus_c)
       self.__ranges[k]=[new_plus,new_minus]
       self.__comments[k]=[new_plus_c,new_minus_c]
       self.__well_merged = True

    def merge_plus_minus_ranges_w_duplicates (self):
        """Merge minus strand ranges to plus strand. Reset the comments. keep duplicates.

        """
        self.total = 0
        self.total_w_duplicates = None
        for chrom in self.__ranges.keys():
            (plus_tags,minus_tags) = self.__ranges[chrom]
            new_plus_tags = []
            #reset comments
            self.__comments[chrom][0] = []
            self.__comments[chrom][1] = []
            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    new_plus_tags.append(plus_tags[ip])
                    ip += 1
                else:
                    new_plus_tags.append(minus_tags[im])
                    im += 1
            if im < lenm:
                # add rest of minus tags
                new_plus_tags.extend(minus_tags[im:])
            if ip < lenp:
                # add rest of plus tags
                new_plus_tags.extend(plus_tags[ip:])

            self.__ranges[chrom] = [new_plus_tags,None]
            self.total += len(new_plus_tags)

    def status (self):
        """Return the information whether or not the track is sorted and/or merged.

        Return value: tuple (boolean,boolean) means (sorted?,merged?)
        """
        return (self.__well_sorted,self.__well_merged)

    def sort (self):
        """Sort all ranges in list. Recommanded after you import all ranges.
        """
        if self.__well_sorted:
            pass
        for key in self.__ranges.keys():
            self.__ranges[key][0].sort(lambda x,y: cmp(x,y))
            self.__ranges[key][1].sort(lambda x,y: cmp(x,y))
        self.__well_sorted = True


    def __str__ (self):
        return self.__to_wiggle()
        
    def __to_wiggle (self):
        t = "track type=wiggle_0 name=\"variableStep\" description=\"%s\"\n" % (self.annotation)
        for key in self.__ranges.keys():
            t += "variableStep chrom=%s span=%d\n" % (key,self.fw)
            for i in self.__ranges[key][0]:
                t += str(i)+"\n"
            for i in self.__ranges[key][1]:
                t += str(i)+"\n"
        return t

class WigTrackI:
    """Designed only for wig files generated by MACS/pMA2C/MAT(future
    version). The limitation is 'span' parameter for every track must
    be the same.
    
    """
    def __init__ (self):
        self.__data = {}
        self.span=0
        self.maxvalue =-10000
        self.minvalue = 10000

    def add_loc (self,chromosome,pos,value):
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
        self.__data[chromosome][0].append(pos)
        self.__data[chromosome][1].append(value)
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value
        
    def sort (self):
        """Naive sorting for tags. After sorting, counts are massed
        up.

        Note: counts are massed up, so they will be set to 1 automatically.
        """
        for k in self.__data.keys():
            self.__data[k] = sorted(self.__data[k])
            

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([pos],[value])
        """
        if self.__data.has_key(chromosome):
            return self.__data[chromosome]
        else:
            return None
            #raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__data.keys())
        return l

    def write_wig (self, fhd, name, shift=0):
        """Write all data to fhd in Wiggle Format.

        shift will be used to shift the coordinates. default: 0
        """
        chrs = self.get_chr_names()
        fhd.write("track type=wiggle_0 name=\"%s\"\n" % (name))
        for chrom in chrs:
            fhd.write("variableStep chrom=%s span=%d\n" % (chrom,self.span))
            (p,s) = self.__data[chrom]
            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                fhd.write("%d\t%.4f\n" % (pos+shift,score))                

    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score > cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def filter_score_below (self, cutoff=0):
        """Keep points below a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score < cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def write_gff (self, fhd, shift=0, source=".", feature="."):
        """Write all data to fhd in GFF format.

        shift will be used to shift the coordinates. default: 0
        """
        assert isinstance(fhd,file)
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                fhd.write(
                    "\t".join( (chrom,source,feature,
                                str(pi-shift),str(pi-shift+self.span-1),
                                str(si),'+','.'
                                ) )+"\n"
                    )

    def write_bed (self, fhd):
        """Write all data to fhd in BED format.
        
        """
        pass

    def remove_redundant (self):
        """Remove redundant position, keep the highest score.
        
        """
        chrs = set(self.__data.keys())
        ndata = {}
        for chrom in chrs:
            ndata[chrom] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
            nd_p_append = ndata[chrom][0].append
            nd_s_append = ndata[chrom][1].append
            (p,s) = self.__data[chrom]
            prev_p = None
            prev_s = None
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if not prev_p:
                    prev_p = pi
                    prev_s = si
                else:
                    if pi == prev_p:
                        if si>prev_s:
                            prev_s = si
                    else:
                       nd_p_append (prev_p)
                       nd_s_append (prev_s)
            nd_p_append (prev_p)
            nd_s_append (prev_s)
        del self.__data
        self.__data = ndata

    def find_peaks (self, bw=None):
        """A naive peak finding algorithm to find all local maximum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            increase = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si > prev_s:
                    # increase
                    increase = True
                elif si < prev_s:
                    # decrease
                    if increase:
                        # prev_p is a summit
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s > prev_peak_s:
                                # new summit is high
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new summit
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    increase = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def find_valleys (self, bw=None):
        """A naive peak finding algorithm to find all local minimum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            decrease = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si < prev_s:
                    # decrease
                    decrease = True
                elif si > prev_s:
                    # increase
                    if decrease:
                        # prev_p is a valley
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s < prev_peak_s:
                                # new summit is lower
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new valley
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    decrease = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def total (self):
        t = 0
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            t += len(p)
        return t


class BEDParser:
    """File Parser Class for tabular File.

    """
    def __init__ (self):
        """All parameters except for format are the same as python's
        builtin file class.

        Format parameter can be "bed","gff" or other particular string
        for other format of file. For example, if your file is '\\t'
        delimited, and the first column of the file is chromosome
        name, the second column is start position, third column is end
        position, fourth column is the strand, the coordinates are
        0-indexed and the range is closed, then you should write the
        format string as "123401\\t" (six numbers and the delimiter).

        Note: Use the fifth and sixth number in format string to
        indicate whether the coordinates are 0-indexed (0) or
        1-indexed (1) and whether the range is opened (0) or closed
        (1) (i.e. whether the end position is included in the range
        or not)
        """
        self.__format__ = "123600\t"
        self.__delimiter__ = self.__format__[6:]
        self.__chr_col__ = int(self.__format__[0])-1
        self.__start_col__ = int(self.__format__[1])-1
        self.__end_col__ = int(self.__format__[2])-1
        self.__strand_col__ = int(self.__format__[3])-1
        self.__start_shift__ = int(self.__format__[4])
        self.__end_shift__ = int(self.__format__[5])

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note: All ranges will be merged (exclude the same
        range) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackI object.

        if do_merge is False, it will not merge the same range after
        the track is built.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline.startswith("track"): return ("track line",None,None) # track line is skipped
        if thisline.startswith("#"): return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()

        chromname = thisfields[0]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

        if len(thisfields) < 6 : # default pos strand if no strand
                                 # info can be found
            return (chromname,
                    int(thisfields[1]),
                    0)
        else:
            if thisfields[5] == "+":
                return (chromname,
                        int(thisfields[1]),
                        0)
            elif thisfields[5] == "-":
                return (chromname,
                        int(thisfields[2]),
                        1)
            else:
                raise self.StrandFormatError(thisline,thisfields[5])

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)


class ELANDResultParser:
    """File Parser Class for tabular File.

    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)


        thisfields = thisline.split()
        thistaglength = len(thisfields[1])

        chromname = thisfields[6]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass


        if thisfields[2] == "U0" or thisfields[2]=="U1" or thisfields[2]=="U2":
            strand = thisfields[8]
            if strand == "F":
                return (chromname,
                        int(thisfields[7])-1,
                        0)
            elif strand == "R":
                return (chromname,
                        int(thisfields[7])+thistaglength-1,
                        1)
            else:
                raise self.StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

class ELANDMultiParser:
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()
        thistagname = thisfields[0]         # name of tag
        thistaglength = len(thisfields[1]) # length of tag

        if len(thisfields) < 4:
            return (None,None,None)
        else:
            thistaghits = sum(map(int,thisfields[2].split(':')))
            if thistaghits > 1:
                # multiple hits
                return (None,None,None)
            else:
                (chromname,pos) = thisfields[3].split(':')

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                strand  = pos[-2]
                if strand == "F":
                    return (chromname,
                            int(pos[:-2])-1,
                            0)
                elif strand == "R":
                    return (chromname,
                            int(pos[:-2])+thistaglength-1,
                            1)
                else:
                    raise self.StrandFormatError(thisline,strand)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

class PairEndELANDMultiParser:
    """File Parser Class for two ELAND multi Files for Pair-End sequencing.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, lfhd, rfhd, dist=200):
        """Build FWTrackI from all lines, return a FWTrackI object.

        lfhd: the filehandler for left tag file
        rfhd: the filehandler for right tag file
        dist: the best distance between two tags in a pair

        The score system for pairing two tags:

        score = abs(abs(rtag-ltag)-200)+error4lefttag+error4righttag

        the smaller score the better pairing. If the score for a
        pairing is bigger than 200, this pairing will be discarded.

        Note only the best pair is kept. If there are over two best
        pairings, this pair of left and right tags will be discarded.

        Note, the orders in left tag file and right tag file must
        match, i.e., the Nth left tag must has the same name as the
        Nth right tag.

        Note, remove comment lines beforehand.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        lnext = lfhd.next
        rnext = rfhd.next
        self.dist = dist
        try:
            while 1:
                lline = lnext()
                rline = rnext()
                (chromname,fpos,strand) = self.__fw_parse_line(lline,rline)

                i+=1
                if i == 1000000:
                    m += 1
                    logging.info(" %d" % (m*1000000))
                    i=0
                if not fpos or not chromname:
                    continue

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                fwtrack.add_loc(chromname,fpos,strand)

        except StopIteration:
            pass
        return fwtrack
    
    def __fw_parse_line (self, leftline, rightline ):
        # >HWI-EAS275_5:4:100:340:1199/1	GTGCTGGTGGAGAGGGCAAACCACATTGACATGCT	2:1:0	chrI.fa:15061365F0,15068562F0,chrIV.fa:4783988R1
        # >HWI-EAS275_5:4:100:340:1199/2	GGTGGTGTGTCCCCCTCTCCACCAGCACTGCGGCT	3:0:0	chrI.fa:15061451R0,15068648R0,15071742R0

        leftfields = leftline.split()
        lefttaglength = len(leftfields[1]) # length of tag
        rightfields = rightline.split()
        righttaglength = len(rightfields[1]) # length of tag

        if len(rightfields) < 4 or len(leftfields) < 4:
            # one of the tag cann't be mapped to genome
            return (None,None,None)
        else:
            lefthits = self.__parse_line_to_dict(leftfields[3])
            righthits = self.__parse_line_to_dict(rightfields[3])            
            parings = []

            for seqname in lefthits.keys():
                if not righthits.has_key(seqname):
                    continue
                else:
                    leftpses = lefthits[seqname] # pse=position+strand+error
                    rightpses = righthits[seqname]
                    for (lp,ls,le) in leftpses:
                        for (rp,rs,re) in rightpses:
                            # try to pair them
                            if ls == 'F':
                                if rs == 'R':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),0) )
                                        parings.append((score,seqname,lp,0))
                                else:
                                    # strands don't match
                                    continue
                            else:
                                if rs == 'F':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),1) )
                                        parings.append((score,seqname,lp,1))
                                else:
                                    # strands don't match
                                    continue
            if not parings:
                return (None,None,None)
            parings.sort()
            if len(parings)>1 and parings[0][0] == parings[1][0]:
                # >2 best paring, reject!
                return (None,None,None)
            else:
                return parings[0][1:]                                
                    
    def __parse_line_to_dict ( self, linestr ):
        items = linestr.split(',')
        hits = {}
        for item in items:
            if item.find(':') != -1:
                # a seqname section
                (n,pse) = item.split(":") # pse=position+strand+error
                try:
                    n = n[:n.rindex(".fa")]
                except ValueError:
                    pass

                hits[n]=[]
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
            else:
                # only pse section
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
        return hits

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

### Contributed by dawe
class SAMParser:
    """File Parser Class for SAM File.

    Each line of the output file contains at least: 
    1. Sequence name 
    2. Bitwise flag
    3. Reference name
    4. 1-based leftmost position fo clipped alignment
    5. Mapping quality
    6. CIGAR string
    7. Mate Reference Name
    8. 1-based leftmost Mate Position
    9. Inferred insert size
    10. Query sequence on the same strand as the reference
    11. Query quality
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        if thisline[0] == "@": return ("comment line",None,None) # header line started with '@' is skipped
        thisfields = thisline.split()
        thistagname = thisfields[0]         # name of tag
        thisref = thisfields[2]
        bwflag = int(thisfields[1])
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None, None, None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (None, None, None)   # not a proper pair
            if bwflag & 8:
                return (None, None, None)   # the mate is unmapped
            p1pos = int(thisfields[3]) - 1
            p2pos = int(thisfields[7]) - 1
            if p1pos > p2pos:
                # this pair is the farthest one, skip it
                return (None, None, None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = int(thisfields[3]) - 1 + len(thisfields[9])	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0
            thisstart = int(thisfields[3]) - 1	

        try:
            thisref = thisref[:thisref.rindex(".fa")]
        except ValueError:
            pass

        return (thisref, thisstart, thisstrand)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

class BAMParser:
    """File Parser Class for BAM File.

    File is gzip-compatible and binary.
    Information available is the same that is in SAM format.
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        references = []
        # move to pos 4, there starts something
        fhd.seek(4)
        header_len =  struct.unpack('<i', fhd.read(4))[0]
        fhd.seek(header_len + fhd.tell())
        # get the number of chromosome
        nc = struct.unpack('<i', fhd.read(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fhd.read(4))[0]
            references.append(fhd.read(nlength)[:-1])
            # jump over chromosome size, we don't need it
            fhd.seek(fhd.tell() + 4)
        
        while 1:
            try:
                entrylength = struct.unpack('<i', fhd.read(4))[0]
            except struct.error:
                break
            (chrid,fpos,strand) = self.__fw_binary_parse(fhd.read(entrylength))                    
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if fpos >= 0:
                fwtrack.add_loc(references[chrid],fpos,strand)
        fhd.close()
        return fwtrack
    
    def __fw_binary_parse (self, data ):
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return (None,-1,None)

        thisref = struct.unpack('<i', data[0:4])[0]
        thisstart = struct.unpack('<i', data[4:8])[0]
        (cigar, bwflag) = struct.unpack('<hh', data[12:16])
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None, -1, None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (None, -1, None)   # not a proper pair
            if bwflag & 8:
                return (None, -1, None)   # the mate is unmapped
            p1pos = thisstart
            p2pos = struct.unpack('<i', data[24:28])[0]
            if p1pos > p2pos:
                # this pair is the farthest one, skip it
                return (None, -1, None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = thisstart + struct.unpack('<i', data[16:20])[0]	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0

        return (thisref, thisstart, thisstrand)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

### End ###

class BowtieParser:
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note: All ranges will be merged (exclude the same
        range) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackI object.

        if do_merge is False, it will not merge the same range after
        the track is built.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        """
        The following definition comes from bowtie website:
        
        The bowtie aligner outputs each alignment on a separate
        line. Each line is a collection of 8 fields separated by tabs;
        from left to right, the fields are:

        1. Name of read that aligned

        2. Orientation of read in the alignment, - for reverse
        complement, + otherwise

        3. Name of reference sequence where alignment occurs, or
        ordinal ID if no name was provided

        4. 0-based offset into the forward reference strand where
        leftmost character of the alignment occurs

        5. Read sequence (reverse-complemented if orientation is -)

        6. ASCII-encoded read qualities (reversed if orientation is
        -). The encoded quality values are on the Phred scale and the
        encoding is ASCII-offset by 33 (ASCII char !).

        7. Number of other instances where the same read aligns
        against the same reference characters as were aligned against
        in this alignment. This is not the number of other places the
        read aligns with the same number of mismatches. The number in
        this column is generally not a good proxy for that number
        (e.g., the number in this column may be '0' while the number
        of other alignments with the same number of mismatches might
        be large). This column was previously described as"Reserved".

        8. Comma-separated list of mismatch descriptors. If there are
        no mismatches in the alignment, this field is empty. A single
        descriptor has the format offset:reference-base>read-base. The
        offset is expressed as a 0-based offset from the high-quality
        (5') end of the read.

        """
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline.startswith("#"): return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()

        chromname = thisfields[2]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

            if thisfields[1] == "+":
                return (chromname,
                        int(thisfields[3]),
                        0)
            elif thisfields[1] == "-":
                return (chromname,
                        int(thisfields[3])+len(thisfields[4]),
                        1)
            else:
                raise self.StrandFormatError(thisline,thisfields[1])

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

