"""
   python debug.py --proc=?               (gdb mode)
   python debug.py --idb --proc=?         (idb mode)
   python debug.py --pid=?               (gdb mode)
   python debug.py --idb --pid=?         (idb mode)
"""
import os, getopt, sys

try:
    opts, args = getopt.getopt(sys.argv[1:], "ip:",["idb","proc=","pid="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    print __doc__
    sys.exit(2)
idb = False
for o, a in opts:
    if o in ("-i", "--idb"):
        idb = True
    elif o in ("-p", "--proc"):
        try: 
           iproc = int(a)
        except:
           iproc = 0
        host=(os.uname())[1]
        os.system("top -b -u fuentes  -n 1 > dbgtmp%s%d" % (host,iproc))
        os.system("grep -i imaging dbgtmp%s%d > dbg%s%d" % (host,iproc,host,iproc))
        os.system("rm dbgtmp%s%d " % (host,iproc) )
        #import data
        topinfo=open("dbg%s%d" % (host,iproc) ,"r")
        data = [line.split() for line in topinfo]
        topinfo.close()
        pids = range(len(data))
        for i in range(len(data)):
            pids[i]=int(data[i][0])
        pids.sort()
        PID=pids[iproc]
    elif o in ("--pid"):
        PID = int(a)
    else:
        assert False, "unhandled option"


if(idb):
  os.system( "idb -pid %d %s/exec/imaging_%s " %  
             (PID,os.environ.get("WORK"),os.environ.get("PETSC_ARCH")) )
else:
  os.system( "gdb %s/exec/imaging_%s %d" %  
             (os.environ.get("WORK"),os.environ.get("PETSC_ARCH"),PID) )
#os.system("rm dbg%s%d " % (host,iproc))
