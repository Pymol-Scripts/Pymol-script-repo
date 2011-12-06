import string, types, Numeric


class ConstrainedParameterSet:

    def __init__(self):
        #conDict records constraints between parameters
        #key is parm1Name,  name of parameter1 , 
        #value is list of triples: (parm2Name, func, args) 
        #changes in parm1 cause parm2 to be updated by func
        #eg: to enforce 
        #self.center = 2*self.offset
        #so that changes in self.offset force changes in self.center
        #self.conDict['offset'] = [('center', Numeric.multiply, ('offset, 2.0'))]
        #to enforce the reciprocal constraint 
        #so that changes in self.center force changes in self.offset
        #self.conDict['center'] = [('offset', Numeric.multiply, ('center,0.5'))]
        #suitable functions are Numeric.divide and Numeric.multiply
        #DO SOMETHING SO THIS DOESN'T GET into an endless loop
        #
        self.conDict = {}
        #rangeDict provides methods of checking validity 
        #of a given value for key
        #possible types methods include:
        #list of discrete values, tuple defining valid range, a type
        #or a function which returns 1 if value is valid or 0 if not
        self.rangeDict = {}
        # values can be integers, floats, strings....???
        self.typesList = [type(1), type(1.0), type('a')]


    def tie(self, parm1Name, parm2Name, func, args):
        #eg:
        # changes in self.center force  changes in self.offset
        # self.tie('center','offset',Numeric.multiply,'center,2.0')
        if not self.conDict.has_key(parm1Name):
            self.conDict[parm1Name] = []
        self.conDict[parm1Name].append((parm2Name, func, args))
        #is this at all clear?
        #cD = self.conDict
        #cD[parm1Name] = cD.get(parm1Name, []).append((parm2Name, func, args))


    def updateConstraints(self, parmName):
        #called when parm changes to update other linked parms
        conList = self.conDict.get(parmName, None)
        if not conList:
            # do nothing + return
            print 'no constraints on ', parmName
            return 
        #conList has tuples (func, args)
        #eg: sample value in conList
        # (parm2Name, Numeric.multiply, (parmName, 0.5))
        for parm2Name, func, args in conList:
            #FIX THIS:
            #to update self.parm2Name:
            #   need to get (self.center, 0.5) from args='center, 0.5'
            setattr(self,parm2Name, apply(func, eval('self.'+args)))


    def untie(self, parm1Name, parm2Name, func, args):
        #eg:
        #g.untie('center','offset',Numeric.multiply,'center,2.0')
        conList = self.conDict.get(parm1Name, None)
        if not conList:
            print 'no constraints on ', parm1Name
            return "ERROR"
        if (parm2Name, func, args) not in conList:
            print '(%s,%s,%s) not in %s constraints'%(parm2Name, func, args, parm1Name)
            return "ERROR"
        self.conDict[parm1Name].remove((parm2Name, func, args))

 
    def setRange(self, parm, range):
        #range can be list, interval tuple, type or validation func
        #FIX THIS: do some range validation
        self.rangeDict[parm] = range


    def validateValue(self, parm, value):
        rangeD = self.rangeDict
        if not rangeD.has_key(parm):
            #nothing specified for this parm
            return value
        range = rangeD[parm]
        if type(range)==types.ListType:
            if value in range:
                return value
            else:
                return "ERROR: value not in range list"
        elif type(range)==types.TupleType:
            if value>=range[0]and value<=range[1]:
                return value
            else:
                return "ERROR: value not in range interval"
        elif range in self.typesList:
            if type(value)==range:
                return value
            else:
                return "ERROR: value not specified type"
        else:
            #only thing left is validation function
            ok = apply(range, value)
            if ok:
                return value
            else:
                return "ERROR: value failed validation func"


    def update(self, parm, value):
        check = self.validateValue(parm, value)
        if check!=value:
            print 'failed validation:\n', check
            return "ERROR"
        self.updateConstraints(parm)

            
    def fix(self, parmName):
        #????????????????????????????
        #this method makes self.parmName constant
        #by removing any constraints which force it to change
        for k, v in self.conDict.items():
            for triple in v:
                if triple[0]==parmName:
                    self.conDict[k].remove(triple)




class GeneralRegularGridDescriptor(ConstrainedParameterSet):

    keywords = ['center',
                'offset',
                'length',
                'nbGridPoints',
                'gridSpacing'
                ]

    def __init__(self, **kw):

        ConstrainedParameterSet.__init__(self)

        for k in self.keywords:
            setattr(self, k, kw.get(k, Numeric.array((0.,0.,0.))))

        consDict = kw.get('consDict', None)
        if consDict:
            for k, v in consDict.items():
                self.tie(k, v[0], v[1], v[2])

        rangeDict = kw.get('rangeDict', None)
        if rangeDict:
            for k, v in rangeDict.items():
                self.setRange(k, v)

        fixed = kw.get('fixed', None)
        if fixed:
            #fixed should be a list of parameters to fix
            #same problem of type(parm)...a string???
            for p in fixed:
                self.fix(p)
