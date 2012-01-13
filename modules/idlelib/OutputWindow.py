import re
import tkMessageBox

import utils
import IOBinding
from editorpage import OutputPage
from EditorWindow import EditorWindow

class OutputWindow(EditorWindow):
    """An editor window that can serve as an output file.

    Also the future base class for the Python shell window.
    This class has no input facilities.
    """

    def __init__(self, *args):
        EditorWindow.__init__(self, start_page=OutputPage, *args)
        self.remove_tab_controls()

    # Act as output file

    def write(self, s, tags=(), mark="insert", text=None):
        if text is None:
            page = self.current_page
            if not page: # tab was destroyed
                return
            text = page.text
        # Tk assumes that byte strings are Latin-1;
        # we assume that they are in the locale's encoding
        if isinstance(s, str):
            try:
                s = unicode(s, IOBinding.encoding)
            except UnicodeError:
                # some other encoding; let Tcl deal with it
                pass
        text.insert(mark, s, tags)
        text.see(mark)
        text.update()

    def writelines(self, l):
        map(self.write, l, text=self.current_page.text)

    def flush(self):
        pass
