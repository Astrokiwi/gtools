import os
import io
import sys

class suppressed:
    """Suppresses the output of function calls to a module or object. Currently only works for function calls - this breaks
    access to variables etc.
    
    Example usage:
    
    import noisy_module_with_too_much_output
    quiet_module = output_suppressor.suppressed(noisy_module_with_too_much_output)
    quiet_module.some_function_in_original_module(*args,**kwargs)
    """
    
    def __init__(self,x):
        self._wrapped_class = x
    
    def __getattr__(self,name):
        return self._silenced(self._wrapped_class.__dict__[name])

    def _silenced(self,f):
        def g(*args,**kwargs):
            """Adapted from https://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python"""
            # open 2 fds
            null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
            # save the current file descriptors to a tuple
            save = os.dup(1), os.dup(2)
            # put /dev/null fds on 1 and 2
            os.dup2(null_fds[0], 1)
            os.dup2(null_fds[1], 2)

            output = f(*args,**kwargs)

            # restore file descriptors so I can print the results
            os.dup2(save[0], 1)
            os.dup2(save[1], 2)
            # close the temporary fds
            os.close(null_fds[0])
            os.close(null_fds[1])

            return output
        return g
