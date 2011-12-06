#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
from compiler.ast import Module
import StringIO, copy, getopt
import os, sys, unittest, urlparse, signal, time, warnings, subprocess
from ConfigParser import ConfigParser, NoSectionError, NoOptionError
from ZSI.wstools.TimeoutSocket import TimeoutError
from ZSI.generate import commands

"""Global Variables:
    CONFIG_FILE -- configuration file 
    CONFIG_PARSER -- ConfigParser instance
    DOCUMENT -- test section variable, specifying document style.
    LITERAL -- test section variable, specifying literal encodings.
    BROKE -- test section variable, specifying broken test.
    TESTS -- test section variable, whitespace separated list of modules.
    SECTION_CONFIGURATION -- configuration section, turn on/off debuggging.
    TRACEFILE -- file class instance.
    TOPDIR -- current working directory
    MODULEDIR  -- stubs directory 
    PORT -- port of local container
    HOST -- address of local container
    SECTION_SERVERS -- services to be tested, values are paths to executables.
"""
CONFIG_FILE = 'config.txt'
CONFIG_PARSER = ConfigParser()
DOCUMENT = 'document'
LITERAL = 'literal'
BROKE = 'broke'
TESTS = 'tests'
SECTION_CONFIGURATION = 'configuration'
SECTION_DISPATCH = 'dispatch'
TRACEFILE = sys.stdout
TOPDIR = os.getcwd()
MODULEDIR = os.path.join(TOPDIR, 'stubs')
SECTION_SERVERS = 'servers'

CONFIG_PARSER.read(CONFIG_FILE)

DEBUG = CONFIG_PARSER.getboolean(SECTION_CONFIGURATION, 'debug')
SKIP = CONFIG_PARSER.getboolean(SECTION_CONFIGURATION, 'skip')
TWISTED = CONFIG_PARSER.getboolean(SECTION_CONFIGURATION, 'twisted')
LAZY = CONFIG_PARSER.getboolean(SECTION_CONFIGURATION, 'lazy')
OUTPUT = CONFIG_PARSER.get(SECTION_CONFIGURATION, 'output') or sys.stdout

if DEBUG:
    from ZSI.wstools.logging import setBasicLoggerDEBUG
    setBasicLoggerDEBUG()

sys.path.append('%s/%s' %(os.getcwd(), 'stubs'))
ENVIRON = copy.copy(os.environ)
ENVIRON['PYTHONPATH'] = ENVIRON.get('PYTHONPATH', '') + ':' + MODULEDIR


def _SimpleMain():
    """Gets tests to run from configuration file.
    """
    unittest.TestProgram(defaultTest="all")
main = _SimpleMain


def _TwistedMain():
    """Gets tests to run from configuration file.
    """
    from twisted.internet import reactor
    reactor.callWhenRunning(_TwistedTestProgram, defaultTest="all")
    reactor.run(installSignalHandlers=0)
if TWISTED: main = _TwistedMain


def _LaunchContainer(cmd):
    '''
    Parameters:
        cmd -- executable, sets up a ServiceContainer or ?
    '''
    host = CONFIG_PARSER.get(SECTION_DISPATCH, 'host')
    port = CONFIG_PARSER.get(SECTION_DISPATCH, 'port')
    try:
        process = subprocess.Popen(['python', cmd, port], env=ENVIRON)
    except:
        print >>sys.stderr, 'error executing: %s' %cmd
        raise
    time.sleep(3)
    return process


class _TwistedTestProgram(unittest.TestProgram):

    def runTests(self):
        from twisted.internet import reactor
        if self.testRunner is None:
            self.testRunner = unittest.TextTestRunner(verbosity=self.verbosity)

        result = self.testRunner.run(self.test)
        reactor.stop()
        return result.wasSuccessful()



class ConfigException(Exception):
    """Exception thrown when configuration settings arent correct.
    """
    pass

class TestException(Exception):
    """Exception thrown when test case isn't correctly set up.
    """
    pass


class ServiceTestCase(unittest.TestCase):
    """Conventions for method names:
    test_net*
    -- network tests
    
    test_local*
    -- local tests
    
    test_dispatch*
    -- tests that use the a spawned local container
    
    class attributes: Edit/Override these in the inheriting class as needed
        out -- file descriptor to write output to
        name -- configuration item, must be set in class.
        url_section -- configuration section, maps a test module 
           name to an URL.
        client_file_name --
        types_file_name --
        server_file_name --
    """
    out = OUTPUT
    name = None
    url_section = 'WSDL'
    client_file_name = None
    types_file_name = None
    server_file_name = None
    
    def __init__(self, methodName):
        """
        parameters:
           methodName -- 
        instance variables:
            client_module
            types_module
            server_module
            processID
            done

        """
        self.methodName = methodName
        self.url = None
        self.wsdl2py_args = []
        self.wsdl2dispatch_args = []
        self.portkwargs = {}
        self.client_module = self.types_module = self.server_module = None
        self.done = False

        if TWISTED:
            self.wsdl2py_args.append('--twisted')

        if LAZY:
            self.wsdl2py_args.append('--lazy')

        unittest.TestCase.__init__(self, methodName)

    write = lambda self, arg: self.out.write(arg)

    if sys.version_info[:2] >= (2,5):
        _exc_info = unittest.TestCase._exc_info
    else:
        _exc_info = unittest.TestCase._TestCase__exc_info

    def __call__(self, *args, **kwds):
        self.run(*args, **kwds)

    def run(self, result=None):
        if result is None: result = self.defaultTestResult()
        result.startTest(self)
        testMethod = getattr(self, self.methodName)
        try:
            try:
                self.setUp()
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())
                return

            ok = False
            try:
                t1 = time.time()
                pyobj = testMethod()
                t2 = time.time()
                ok = True
            except self.failureException:
                result.addFailure(self, self._exc_info())
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())

            try:
                self.tearDown()
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())
                ok = False
            if ok: 
                result.addSuccess(self)
                print>>self
                print>>self, "|"+"-"*60
                print>>self, "|  TestCase: %s" %self.methodName
                print>>self, "|"+"-"*20
                print>>self, "|  run time:   %s ms" %((t2-t1)*1000)
                print>>self, "|  return  :   %s" %pyobj
                print>>self, "|"+"-"*60

        finally:
            result.stopTest(self)




    def getPortKWArgs(self):
        kw = {}
        if CONFIG_PARSER.getboolean(SECTION_CONFIGURATION, 'tracefile'):
            kw['tracefile'] = TRACEFILE
        
        kw.update(self.portkwargs)
        return kw
    
    def _setUpDispatch(self):
        """Set this test up as a dispatch test.
        url -- 
        """
        host = CONFIG_PARSER.get(SECTION_DISPATCH, 'host')
        port = CONFIG_PARSER.get(SECTION_DISPATCH, 'port')
        path = CONFIG_PARSER.get(SECTION_DISPATCH, 'path')
        
        scheme = 'http'
        netloc = '%s:%s' %(host, port)
        params = query = fragment = None
        
        self.portkwargs['url'] = \
            urlparse.urlunparse((scheme,netloc,path,params,query,fragment))
        
    _wsdl = {}
    def _generate(self):
        """call the wsdl2py and wsdl2dispatch scripts and
        automatically add the "-f" or "-u" argument.  Other args
        can be appended via the "wsdl2py_args" and "wsdl2dispatch_args"
        instance attributes.
        """
        url = self.url
        if os.path.isfile(url):
            url = os.path.abspath(url)

        if SKIP:
            ServiceTestCase._wsdl[url] = True
            return
        
        ServiceTestCase._wsdl[url] = False
        try:
            os.mkdir(MODULEDIR)
        except OSError, ex:
            pass

        os.chdir(MODULEDIR)
        if MODULEDIR not in sys.path:
            sys.path.append(MODULEDIR)
 
        try:
            commands.wsdl2py([url] + self.wsdl2py_args)
            ServiceTestCase._wsdl[url] = True
        finally:
            os.chdir(TOPDIR)
            
    _process = None
    _lastToDispatch = None
    def setUp(self):
        """Generate types and services modules once, then make them
        available thru the *_module attributes if the *_file_name 
        attributes were specified.
        """
        section = self.url_section
        name = self.name
        if not section or not name:
            raise TestException, 'section(%s) or name(%s) not defined' %(
                section, name)
            
        if not CONFIG_PARSER.has_section(section):
            raise TestException,\
                'No such section(%s) in configuration file(%s)' %(
                self.url_section, CONFIG_FILE)

        self.url = CONFIG_PARSER.get(section, name)
        
        status = ServiceTestCase._wsdl.get(self.url)
        if status is False:
            self.fail('generation failed for "%s"' %self.url)
            
        if status is None:
            self._generate()
            
        # Check for files
        tfn = self.types_file_name
        cfn = self.client_file_name
        sfn = self.server_file_name
    
        files = filter(lambda f: f is not None, [cfn, tfn,sfn])
        if None is cfn is tfn is sfn:
            return
        
        for n,m in map(lambda i: (i,__import__(i.split('.py')[0])), files):
            if tfn is not None and tfn == n:
                self.types_module = m
            elif cfn is not None and cfn == n:
                self.client_module = m
            elif sfn is not None and sfn == n:
                self.server_module = m
            else: 
                self.fail('Unexpected module %s' %n)

        # DISPATCH PORTION OF SETUP
        if not self.methodName.startswith('test_dispatch'):
            return
        
        self._setUpDispatch()
        if ServiceTestCase._process is not None:
            return

        try:
            expath = CONFIG_PARSER.get(SECTION_DISPATCH, name)
        except (NoSectionError, NoOptionError), ex:
            self.fail('section dispatch has no item "%s"' %name)

        if ServiceTestCase._lastToDispatch == expath:
            return
        
        if ServiceTestCase._lastToDispatch is not None:
           ServiceTestCase.CleanUp()
        
        ServiceTestCase._lastToDispatch = expath
        ServiceTestCase._process = \
            _LaunchContainer(os.path.join(os.path.abspath(TOPDIR), 
                                                          *expath.split('/')))
       
            
    def CleanUp(cls):
        """call this when dispatch server is no longer needed,
        maybe another needs to be started.  Assumption that
        a single "Suite" uses the same server, once all the
        tests are run in that suite do a cleanup.
        """
        if cls._process is None:
            return
        os.kill(cls._process.pid, signal.SIGKILL)
        cls._process = None
    CleanUp = classmethod(CleanUp)
         
         
class ServiceTestSuite(unittest.TestSuite):
    """A test suite is a composite test consisting of a number of TestCases.

    For use, create an instance of TestSuite, then add test case instances.
    When all tests have been added, the suite can be passed to a test
    runner, such as TextTestRunner. It will run the individual test cases
    in the order in which they were added, aggregating the results. When
    subclassing, do not forget to call the base class constructor.
    """
    def __init__(self, tests=()):
        unittest.TestSuite.__init__(self, tests)

    def __call__(self, result):
        # for python2.4
        return self.run(result)

    def addTest(self, test):
        unittest.TestSuite.addTest(self, test)

    def run(self, result):
        for test in self._tests:
            if result.shouldStop:
                break
            test(result)
            
        ServiceTestCase.CleanUp()
        return result
    

