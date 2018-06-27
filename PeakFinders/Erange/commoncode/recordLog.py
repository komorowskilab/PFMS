#
#  recordLog.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 12/14/08.
#

from commoncode import writeLog
import sys

if '-verbose' in sys.argv or len(sys.argv) < 4:
    print '%s: version 1.0' % sys.argv[0]
    
if len(sys.argv) < 4:
    print 'usage: python %s logFile messenger message [-verbose]' % sys.argv[0]
    sys.exit(1)

writeLog(sys.argv[1], sys.argv[2], sys.argv[3])