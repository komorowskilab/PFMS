#################################################
#########                              ##########
#########      README for SeqSite      ##########
#########                              ##########
#################################################

SeqSite is for pinpointing transcription factor binding sites 
from ChIP-seq data. It can locate closely spaced binding sites
and isolated binding sites in detected binding regions at a 
high resolution. 


##############
  1. INSTALL
##############

1.1 download the most updated SeqSite source code package from 
http://bioinfo.au.tsinghua.edu.cn/seqsite/

1.2 Untar the source code package, and 'cd' into the dir 
SeqSite_1.1.2_src.
  
  tar zxvf SeqSite_1.1.2_src.tar.gz
  cd SeqSite_1.1.2_src

1.3 Type 'make' to generate the binary file SeqSite.

  make

1.4 Type 'make install' to copy the executable file SeqSite 
to your directory for binary files: ~/bin 

The default installing path is ~/bin. 
please specify BIN_DIR in Makefile, if you want to install 
SeqSite to anywhere else.

  make install

1.5 Type 'SeqSite' to get the information how to run it.

  SeqSite -h


#######################
  2. RUN SeqSite NOW!
#######################

2.1 Input files for SeqSite

ChIP-seq tags: a BED file with 4 fields required: chrId, start, end, and strand
Control tags:  a BED file with 4 fields required: chrId, start, end, and strand

It is recommended to run SeqSite with control data, although it is 
not required.

Users can download PERL script provided on our website 
http://bioinfo.au.tsinghua.edu.cn/seqsite/
to convert other formats to BED. 


2.2 Usage

SeqSite [options] <input.bed> <output.bar> <output.bed>
        input.bed    ChIP-seq data in BED format (4 fields required: chrId, start, end, and strand)
        output.bar   BAR file containing binding sites identified
        output.bed   BED file containing binding regions detected
Options: (* advanced)
        -c <string>  control data in BED format (4 fields required) (default: not use)
        -g <int>     effective genome size (default: 2.4e+9 for the human genome)
        -d <int>     * tag clustering distance (default: 30)
        -n <int>     * min tag count in a tag cluster (default: 10)
        -S           * filter single-strand tag clusters (default: not filter)
        -l <double>  * average DNA fragment length (default: estimate from data)
        -t <int>     * top <int>% tag clusters for frag. length estimating (default: 5)
        -p <double>  p-value cutoff for binding region detection (default: 1e-3)
        -f <double>  FDR for binding region detection (default: 0.1)
        -s <int>     * arm length for smoothing tag signal (default: 20)
        -k <int>     * kernel density bandwidth for smoothing tag signal (default: use -s)
        -w <int>     * experimental motif width (default: 20)
        -F           * filter out the duplicate reads (default: FALSE)
        -q           quiet: no screen display (default: show progress)
Help Options:
        -h           show this help message
        -v           show version information
        -a           about SeqSite


2.3 Output files of SeqSite

2.3.1 BED file for binding regions

Each column of the BED file represents:
chr#, start, end, read-count|fold-change|p-value|q-value, score, strand(+)

2.3.2 BAR file for binding sites

Each column of the BAR file represnets:
chr#, position, p-value, fold-change, q-value, R-square, slope(normalized)


################
  3. EXAMPLES
################

We provide an example for a quick start.

3.1 Please download the following data files first:
GABP ChIP-seq data: http://bioinfo.au.tsinghua.edu.cn/seqsite/files/GABP.bed.gz
Control data:       http://bioinfo.au.tsinghua.edu.cn/seqsite/files/RX_noIP.bed.gz

3.2 Unzip the files
  
  gunzip GABP.bed.gz
  gunzip RX_noIP.bed.gz

3.3 Run SeqSite 

  SeqSite -c ./RX_noIP.bed ./GABP.bed GABP.SeqSite.BS.bar GABP.SeqSite.BR.bed


##########################
  4. BUGS and QUESTIONS
##########################

Please return any bug reports and questions to 
Xi Wang (wang-xi05@mails.tsinghua.edu.cn).


