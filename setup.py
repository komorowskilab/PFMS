#! /usr/bin/python
'''
@intention: to build distribution package of PFMS with 'python setup.py sdist' command
	- and installing PFMS by users with 'python setup.py install [-normal]' command
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 18 Mar 2011
@version: 1.1

'''
import distutils
from distutils.core import setup
import os
import sys
from tests import test_install
from tests import test_normal_install

def main():

    distutils.debug = "true"
    distutils.DISTUTILS_DEBUG = "true"
    setup_args = dict(name='PFMS',
                      version='1.1',
                      description='PeakFinder Ranking',
                      author='Marcin Kruczyk & Husen Umer',
                      author_email='marcin.kruczy@icm.uu.se; husen.umer@icm.uu.se',
                      url='http://www.lcb.uu.se/',
                      license='GNU',
                      long_description='Find the best pest genome peaks',
                      platforms=['Linux', 'MACOS', 'Windows'],
                      packages=['src'],
                      package_dir={'src': 'src'},
                      scripts=['src/PFMetaserver'],
                      data_files=[('', ['README.txt'])],

                      package_data={'src': ['../PeakFinders/pfms.conf',
                      '../PeakFinders/sissrs_v1.4/sissrs.pl',
                      '../PeakFinders/SeparateReads.jar',
                      '../PeakFinders/SeqSite1.0/*',
                      '../PeakFinders/findpeaks/*',
                      '../PeakFinders/Erange/commoncode/rnapath/*',
                      '../PeakFinders/Erange/commoncode/*.*',
                      '../PeakFinders/cisGenome-2.0/bin/*',
                      '../PeakFinders/cisGenome-2.0/datatable/*',
                      '../PeakFinders/cisGenome-2.0/src/*',
                      '../PeakFinders/MACS-1.3.7.1/*.py',
                      '../PeakFinders/MACS-1.3.7.1/*.cfg',
                      '../PeakFinders/MACS-1.3.7.1/*.in',
                      '../PeakFinders/MACS-1.3.7.1/PKG-INFO',
                      '../PeakFinders/MACS-1.3.7.1/00README',
                      '../PeakFinders/MACS-1.3.7.1/COPYING',
                      '../PeakFinders/MACS-1.3.7.1/INSTALL',
                      '../PeakFinders/MACS-1.3.7.1/TODO',
                      '../PeakFinders/MACS-1.3.7.1/ChangeLog',
                      '../PeakFinders/MACS-1.3.7.1/bin/*',
                      '../PeakFinders/MACS-1.3.7.1/lib/*.*',
                      '../PeakFinders/MACS-1.3.7.1/lib/IO/*',
                      '../PeakFinders/MACS-1.3.7.1/MACS.egg-info/*'
                      '../PeakFinders/HPeak/HPeak-1.1/*.*',
                      '../PeakFinders/HPeak/HPeak-1.1/data/*'
                      ],
                      },

                      classifiers=[
                      'Development Status :: Under Developement',
                      'Environment :: Console',
                      'Intended Audience :: Chip-Seq Users, Bioinformaticians',
                      'License :: Artistic License',
                      'Operating System :: MacOS :: MacOS X',
                      'Operating System :: Microsoft :: Windows',
                      'Operating System :: Linux',
                      'Programming Language :: Python',
                      ]
                      )
    setup( ** setup_args)

def install_BEDTools(BEDTools_path, root_path):
    os.chdir(BEDTools_path)
    os.system("make")
    os.chdir(root_path)

def copy_BEDTools_bin(BEDTools_path, destination_path, root_path):
    os.chdir(BEDTools_path)
    os.system("sudo cp bin/* "+ destination_path)
    os.chdir(root_path)

def install_macs(macs_path, sudo, root_path):
    os.chdir(macs_path)
    os.system(sudo + " chmod +x macs")
    os.chdir(root_path)
    
def install_cisgenome(cisgenome_path, root_path):
    os.chdir(cisgenome_path)
    os.system("make bin")
    os.system("make clean")
    os.chdir(root_path)

def install_seqsite(seqsite_path, root_path):
    os.chdir(seqsite_path)
    os.system("make")
    os.chdir(root_path)

def install_hpeak(hpeak_path, root_path):
    os.chdir(hpeak_path)
    os.system("g++ -o chiphmm chiphmm.cpp")
    os.system("g++ -o hmmminus chiphmmminus.cpp")
    os.chdir(root_path)

if __name__ == '__main__':
    if not float(sys.version[:3]) >= 2.4:
        sys.stderr.write("Error: Python 2.4 or greater is required! python 2.6.2 is recommended!\n")
        sys.exit(1)
    #install the pfs if arg[1]==install
    if 'install' in sys.argv:
        if '-normal' in sys.argv:
            print "Installing the Peakfinders"
            install_BEDTools("./PeakFinders/BEDTools-Version-2.16.2/", os.getcwd()) 
            
            install_cisgenome("./PeakFinders/cisGenome-2.0/src/", os.getcwd())
            install_hpeak("./PeakFinders/HPeak/HPeak-1.1/", os.getcwd())
            install_seqsite("./PeakFinders/SeqSite1.0/", os.getcwd())
            install_macs("./PeakFinders/MACS-1.3.7.1/bin/", "", os.getcwd())
            
            test_normal_install()

        else:
            main()
            sudo = "sudo " #on windows set it to null
            if sys.platform.lower().startswith("win"):
                sudo = ''
            #if the PeakFinders directory is not exist then create on
            if not os.path.isdir(sys.prefix + "/PeakFinders"):
                os.system(sudo + "mkdir " + sys.prefix + "/PeakFinders")
            os.system(sudo + "cp -f -r PeakFinders/* " + sys.prefix + "/PeakFinders/")
            print "Installing the Peakfinders"
            install_BEDTools("./PeakFinders/BEDTools-Version-2.16.2/", os.getcwd()) 
            copy_BEDTools_bin("PeakFinders/BEDTools-Version-2.16.2/", "/usr/local/bin/", os.getcwd())
            
            install_cisgenome(sys.prefix + "/PeakFinders/cisGenome-2.0/src/", "./")
            install_hpeak(sys.prefix + "/PeakFinders/HPeak/HPeak-1.1/", "./")
            install_seqsite(sys.prefix + "/PeakFinders/SeqSite1.0/", "./")
            install_macs(sys.prefix + "/PeakFinders/MACS-1.3.7.1/bin/", sudo, "./")
            test_install()

    #Remove the installed files
    elif "remove" in sys.argv:
        if os.access(sys.prefix + "/PeakFinders/", os.F_OK):
            try:
                os.removedirs(get_python_lib() + "/src")
            except:
                os.system("rm -r " + sys.prefix + "/PeakFinders/")

        if os.access(sys.prefix + "/README.txt", os.F_OK):
            os.remove(sys.prefix + "/README.txt")

        if os.access("/usr/bin/" + "/PFMetaserver", os.F_OK):
            os.remove("/usr/bin/" + "/PFMetaserver")
        if os.access("/usr/local/bin/" + "/PFMetaserver", os.F_OK):
            os.remove("/usr/local/bin/" + "/PFMetaserver")

        from distutils.sysconfig import get_python_lib

        if os.access(get_python_lib() + "/src", os.F_OK):
            try:
                os.removedirs(get_python_lib() + "/src")
            except:
                os.system("rm -r " + get_python_lib() + "/src")
        if os.access(sys.prefix + "/local" + get_python_lib()[4:] + "/src", os.F_OK):
            try:
                os.removedirs(sys.prefix + "/local" + get_python_lib()[4:] + "/src")
            except:
                os.system("rm -r " + sys.prefix + "/local" + get_python_lib()[4:] + "/src")

        if os.access(get_python_lib() + "/PeakFinders", os.F_OK):
            try:
                os.removedirs(get_python_lib() + "/PeakFinders")
            except:
                os.system("rm -r " + get_python_lib() + "/PeakFinders")
        if os.access(sys.prefix + "/local" + get_python_lib()[4:] + "/PeakFinders", os.F_OK):
            try:
                os.removedirs(sys.prefix + "/local" + get_python_lib()[4:] + "/PeakFinders")
            except:
                os.system("rm -r " + sys.prefix + "/local" + get_python_lib()[4:] + "/PeakFinders")

        if os.access(get_python_lib() + "/PFMS-1.0.egg-info", os.F_OK):
            os.remove(get_python_lib() + "/PFMS-1.0.egg-info")
        if os.access(sys.prefix + "/local" + get_python_lib()[4:] + "/PFMS-1.0.egg-info", os.F_OK):
            os.remove(sys.prefix + "/local" + get_python_lib()[4:] + "/PFMS-1.0.egg-info")

    else:
        main()
