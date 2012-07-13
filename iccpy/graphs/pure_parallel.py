import sys
from subprocess import call, Popen
from time import sleep

def sample_func(i):
    """ a sample function """
    print i,'says hello'
    


    
def parallel_run(func, call_count, delay=1):
    """ calls func(i) for i in 0,...,call_count-1, each from a different process """
    
    if len(sys.argv)==1:
        # in the master iteration
        code_name = sys.argv[0]
        for i in range(call_count):
            iter = '%d'%i
            Popen(' '.join(['python', code_name, iter]), shell=True)
            sleep(delay)
    elif len(sys.argv)==2:
        # we are one of the slaves, call job i
        i = int(sys.argv[1])
        func(i)

        


if __name__=='__main__':
    parallel_run(sample_func, 10)
