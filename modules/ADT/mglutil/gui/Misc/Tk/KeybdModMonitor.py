import Tkinter

from mglutil.util.callback import CallbackManager

class KeyboardModifierMonitor:

    def __init__(self):
        self.kbdModifier = {
            'Shift_L':0,
            'Alt_L':0,
            'Control_L':0,
            'Shift_R':0,
            'Alt_R':0,
            'Control_R':0,
            }
	self.keybdModifierCallbacksDown = {
	    'Shift_L':CallbackManager(),
            'Alt_L':CallbackManager(),
            'Control_L':CallbackManager(),
            'Shift_R':CallbackManager(),
            'Alt_R':CallbackManager(),
            'Control_R':CallbackManager(),
	}
	self.keybdModifierCallbacksUp = {
	    'Shift_L':CallbackManager(),
            'Alt_L':CallbackManager(),
            'Control_L':CallbackManager(),
            'Shift_R':CallbackManager(),
            'Alt_R':CallbackManager(),
            'Control_R':CallbackManager(),
	}

 
    def modifierDown(self, event):
        """track changes in SHIFT, CONTROL, ALT kyes positions"""
        if event.keysym in ['Shift_L', 'Shift_R', 'Control_L', 'Control_R',
                            'Alt_L', 'Alt_R']:
            self.kbdModifier[event.keysym] = 1
            # grab all event to make sure get the key release event even
            # if the mouse is outside the application

            # we have problems with this because when we release we loose
            # the grab. As a consequence, if SHIFT+buttton1 was used to start
            # a rubberband and SHIFT is released BEFORE the button,
            # We loose motion and release event and the line stops moving
            # and is never deleted :(

            # this was an attempt to have tha canvas set the grab. But then
            # modifier event are not caught !
            
            #self.oldgrab = self.master.grab_current()
            #print 'setting global grab', self.oldgrab
            #self.master.grab_set_global()
	    self.keybdModifierCallbacksDown[event.keysym].CallCallbacks(event)
	    

    def modifierUp(self, event):
        """track changes in SHIFT, CONTROL, ALT kyes positions"""
        if event.keysym in ['Shift_L', 'Shift_R', 'Control_L', 'Control_R',
                            'Alt_L', 'Alt_R']:
            self.kbdModifier[event.keysym] = 0
            # release the grab. Release must be done on button release event
            # this is the Problem. if SHFT is released before button we loose
            # button motion and button release events after that.
            # Seems that a solution to this would require this object to also
            # monitor mouse buttons and release the grab after the last release
            # of either the button or the modifier.
            self.master.grab_release()

            #if self.oldgrab:
            #    self.oldgrab.grab.set()
	    self.keybdModifierCallbacksUp[event.keysym].CallCallbacks(event)


    def isShift(self):
        return self.kbdModifier['Shift_L'] or self.kbdModifier['Shift_R']


    def isControl(self):
        return self.kbdModifier['Control_L'] or self.kbdModifier['Control_R']


    def isAlt(self):
        return self.kbdModifier['Alt_L'] or self.kbdModifier['Alt_R']


    def getModifier(self):
        if self.kbdModifier['Shift_L']: return 'Shift_L'
        elif self.kbdModifier['Control_L']: return 'Control_L'
        elif self.kbdModifier['Alt_L']: return 'Alt_L'
        elif self.kbdModifier['Shift_R']: return 'Shift_R'
        elif self.kbdModifier['Control_R']: return 'Control_R'
        elif self.kbdModifier['Alt_R']: return 'Alt_R'
        else: return 'None'
