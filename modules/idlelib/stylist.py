from configHandler import idleConf

TTK = idleConf.GetOption('main', 'General', 'use-ttk', type='int')

class PoorManStyle(object):
    def __init__(self, parent, styles=None, cfgstyles=None):
        self.parent = parent
        self.cfgstyles = cfgstyles
        self.styles = styles

    def configure(self, style, lookup=None, background=None):
        if style not in self.cfgstyles: # passed wrong style probably
            return

        widget = getattr(self.parent, self.cfgstyles[style])
        if lookup:
            return widget.cget('bg')

        widget.configure(bg=background)

    def style_it(self, w, style):
        if TTK:
            w['style'] = style
            return

        if not style in self.styles: # may not need to be styled
            return

        w.configure(**self.styles[style])
