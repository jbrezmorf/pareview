from contextlib import contextmanager
import sys

sys.setrecursionlimit(15)

class XX:
    '''
    Running guard
    '''
    @contextmanager
    def running_guard(self, persistent_object):
        if not hasattr(persistent_object,'Running'):
            persistent_object.Running=False

        try:
            if not persistent_object.Running:
                persistent_object.Running=True
                yield
            else:
                print "Recursion. Skip."
        finally:
            persistent_object.Running=False

    def g(self, func):
        func()
    
    def f(self):
        global glob
        print "before with"
        with self.running_guard(glob):
            print "in with"
            self.g(self.f)


glob=XX()

XX().f()

            