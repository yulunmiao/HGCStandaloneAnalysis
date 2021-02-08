import os, errno
import abc

def create_dir(d):
    try:        
        os.makedirs(d)              
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

class SubmitBase(object):
    def __init__(self, outDir, eosDirOut, bfield, params, eosDirIn=''):
        #variables
        self.outDir = outDir
        self.eosDirIn = eosDirIn
        self.eosDirOut = eosDirOut
        self.p = params
        self.bfield = bfield

        self._sanity_checks()

        self.gran_tag = '$(GRAN)'
        self.en_tag   = '$(ENERGY)'
        self.eta_tag  = '$(ETA)'
        self.run_tag  = '$(Step)'
        self._not_implemented_error = 'Please implement the method in your derived class.'

        #lambda functions
        self.clean_tag = lambda t: t.strip('$').strip('(').strip(')')
        self.shellify_tag = lambda t: t.replace('(', '{').replace(')','}')
        self.remove_dot = lambda t: t.replace('.', '')

        #other operations
        create_dir(self.outDir)

    def _sanity_checks(self):
        if self.bfield not in ('BON', 'BOFF'):
            raise ValueError('[SubmitBase::_sanity_checks] The magnetic filed must be either ON or OFF.')
        if self.p.nRuns < 1 or self.p.nRuns > 1e5:
            raise ValueError('[SubmitBase::_sanity_checks] Please check the number of runs.')
        
    @abc.abstractmethod
    def write_shell_script_file(self):
        raise NotImplementedError(self._not_implemented_error)

    @abc.abstractmethod
    def write_condor_submission_file(self):
        raise NotImplementedError(self._not_implemented_error)
    
    def _unique_name(self, elems):
        """
        Args: -elems: list or tuple of 2-element tuples (key-value pairs)
        Note: using a dict instead does not preserve element order.
        """
        for i,t in enumerate(elems):
            if i==0: #first
                s = str(t[1])
            elif i == len(elems)-1: #last
                s += '.' + str(t[1])
            else:
                s += '_' + t[0] + str(t[1])
        return s
