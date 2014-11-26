def no_recursion(func):
    def inner(*args, **kwargs):
        if not hasattr(func, 'callcount'):
            func.callcount = 0
        if func.callcount >= 1:
            print "recursion detected %s calls deep. exiting." % func.callcount
            return None
        else:
            func.callcount+=1
            x=func(*args, **kwargs)
            func.callcount-=1
            return x        
    return inner  
  
  
  
@no_recursion  
def f(obj):
    return obj()
  
def g():
    f(None)

f(g)    
    
