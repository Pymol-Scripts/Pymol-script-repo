from time import time
import warnings

class Event:
    """Base class for ViewerFramework events.
"""

    def __init__(self, *args, **kw):
        """  """
        self.timestamp = time()
        self.args = args
        self.kw = kw


class EventHandler:
    """This mix-in class adds methods for registening functions called
listeners to be called upon a particular Event.
"""

    def __init__(self):
        self.eventListeners = {}


    def registerListener(self, event, function):
        """
        registers a function to be called for a given event.
        event has to be a class subclassing VFEvent

        None <- registerListener(event, function)

        arguments:
            event: event class
            function: callable object that will be called with the
                      event instance as an argument.
        """
        assert issubclass(event, Event)
        assert callable(function)

        if not self.eventListeners.has_key(event):
            self.eventListeners[event] = [function]
        else:
            if function in self.eventListeners[event]:
                warnings.warn('function %s already registered for event %s'%(
                    function,event))
            else:
                self.eventListeners[event].append(function)


    def dispatchEvent(self, event):
        """call all registered listeners for this event type.
arguments:
    event: instance of an event
"""
        assert isinstance(event, Event)
        if self.eventListeners.has_key(event.__class__):
            for func in self.eventListeners[event.__class__]:
                func(event)
 
