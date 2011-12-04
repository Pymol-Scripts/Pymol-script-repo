# Exceptions used on both versions of tabbedpages
class InvalidNameError(Exception): pass
class AlreadyExistsError(Exception): pass

def get_tabbedpage():
    """Returns the TabbedPageSet available for use."""
    try:
        from tabbedpages_new import TabbedPageSet
    except ImportError:
        from tabbedpages_old import TabbedPageSet

    return TabbedPageSet
