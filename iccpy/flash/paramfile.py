import sys

def update_params(par_file, newparams, buf=sys.stdout):
    """ take parameter file (par_file) and make the corresponding text updated with the (param, newvalue) pairs in the newparams dict """
    f = open(par_file, 'r')

    upperparams = dict((p.upper(), p) for p in newparams)
    
    for line in f:
        # go through each line, test if we should use the old parameter or new

        sline = line.lstrip() # remove leading whitespace

        if len(sline)==0:
            # whitespace
            buf.write(line)
            continue        

        if sline[0]=='#':
            # just a comment
            buf.write(line)
            continue

        
        if '=' not in line:
            # not a parameter!
            raise Exception('Not a valid flash file with line \n%s'%line)

        param = sline.split('=')[0].rstrip()

        if param.upper() not in upperparams:
            # not in our new list
            buf.write(line)
            continue            

        newval = newparams[upperparams[param.upper()]]

        l = '%s   =  %s # UPDATED PARAMETER\n'%(param, newval)
        buf.write(l)
        
    f.close()

def params_to_dict(par_file):
    
    f = open(par_file, 'r')
    res = {}
    for line in f:
        # go through each line, test if we should use the old parameter or new

        sline = line.lstrip() # remove leading whitespace

        if len(sline)==0:
            # whitespace
            continue        

        if sline[0]=='#':
            # just a comment
            continue

        
        if '=' not in line:
            # not a parameter!
            raise Exception('Not a valid flash file with line \n%s'%line)

        param, val = sline.split('=')
        if '#' in val:
            val = val[:val.find('#')]
        
        res[param.rstrip()] = val
        
    f.close()
    return res


if __name__=='__main__':
    par_file = '/gpfs/data/rmdq85/disk_fback/matrix_Sig_fgas_extraKS/L5_Sigma8_fgas5_KSexact/flash.par'
    params = {'lrefine_max':'2', 'lrefine_min':'0'}
    update_params(par_file, params)
