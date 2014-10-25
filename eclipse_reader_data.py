import threading

if (not hasattr(self,"info_done_event")):
    self.info_done_event=threading.Event()

self.info_done_event.wait(1)    
# main RequestData script  
self.code.RequestData(self)