# PFMS 
Peak-finder meta server for ChIP-seq data analaysis: http://bioinf.icm.uu.se/~pfms/

developed at LCB - Uppsala University by Marcin Kruczyk and Husen M. Umer.


     *************Download the distributed package (PFMS-1.3.zip)**********
--Extract

unzip -q PFMS-1.3.zip

cd PFMS-1.3/


     *************Installation****************

--Section 3.2.1 in the manual (for users with root access)

sudo python setup.py install

PFMetaserver -i Input_file.bed -o test1 -bed -sissr -findpeaks -erange -max_cpu_use 2 -min_cpu 3 -quantile 80 -minFN



--Section 3.2.2 in the manual  (for normal users)

python setup.py install -normal

python PFMetaserver.py -i Input_file.bed -o test1 -bed -sissr -findpeaks -erange -max_cpu_use 2 -min_cpu 3 -quantile 80 -minFN


     *************Test Installation****************************

python tests.py [install-normal | install-super]


     *************Test PFMS****************************

Test usage:

python tests.py input_file_path.bed output_label control_file_path.bed no-hpeak [no-macs] ... [no-hpeak] [no-findpeaks] [-super]


Or mannualy run these tests and check the output directory


--Tests:

PFMetaserver -i Input.bed -o test1-bed6 -macs -sissr -seqsite -findpeaks -cisgenome -erange

PFMetaserver -i Input.bed -control Treat.bed -o test2-bedctr6 -macs -sissr -seqsite -findpeaks -cisgenome -erange

PFMetaserver -i Input2.bed -o test3-wig4 -macs -findpeaks -cisgenome -seqsite

PFMetaserver -i Input.bed -control Treat.bed -o test4-wigctr3 -macs -findpeaks -cisgenome



     *************Remove PFMS****************************

--remove the installed version:
python setup.py remove
