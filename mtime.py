import time
"""
this module is used to calculate the time used.
"""

class mtime(object):

    """Docstring for mtime. """

    def __init__(self):
        """ """
        self.time={}

    def start(self, arg1="all"):
        """TODO: Docstring for start.

        :arg1: time name
        :returns: current time

        """
        self.time[arg1] = time.time()

    def end(self, arg1="all",iodata=0,Ncal=0):
        """TODO: Docstring for end.

        :arg1: time name
        :returns: this time name all used times.

        """
        self.time[arg1] = time.time() - self.time[arg1]
        print(20*"=","time")
        print("%s time  used is : %fs"%(arg1, self.time[arg1]))
        if Ncal != 0:
            print("%s Gflops is : %f"%(arg1,Ncal/(1024**3)/self.time[arg1]))
        if iodata != 0:
            print("%s io is %f GB/s"%(arg1,iodata/(1024**3)/self.time[arg1]))
        print(20*"=")
