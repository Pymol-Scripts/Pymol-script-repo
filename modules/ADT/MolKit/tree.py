#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/tree.py,v 1.70.4.2 2011/06/10 21:36:12 sargis Exp $
#
# $Id: tree.py,v 1.70.4.2 2011/06/10 21:36:12 sargis Exp $
#

"""
This module implements the classes TreeNode and TreeNodeSet.
TreeNodes are objects that can be used to build hierachical structures where
each node is a list of its children. Each Child knows its Parent.
An optional 'elementType' argument allows to specify the type of elements
a node can hold. When specified it is used during the 'adoption' process of a
child to assert the child's type validity.
It is also possible to specify the class type of a list of specialized
TreeNodes. I.e. if a node called Residue is derived from TreeNode, one can
specify that a list of Residues should be a ResidueSet (probably derived from
TreeNodeSet).
Also it is possible to specify the type of a set of children of a given node.
i.e. for a Residue node we can specify the a subset of its children will be
an AtomSet (derived from TreeNodeSet).
Finally, the 'children' member can be aliased to a new name. I.e. the children
of a Residue element res can be addressed as res.children or res.atoms if
'atoms' is used as an alias for children in the Residue object.

TreeNodeSets is a class which represents lists of TreeNodes. It inherits
from ListSet which enables Boolean operations among TreeNodeSets.
A 'get' method is provided to filter objects out from the set using a user
specified criterion which is either a function to be applied to each object of the set
or a string which forms the basis for a subselect.
"""

import string, types, re
from MolKit.listSet import ListSet

verbose = False


def evalString(str):
    if len(str)==0:
        return
    try:
        function = eval("%s"%str)
    except:
        #try:
            obj = compile(str, '<string>', 'exec')
            exec(obj)
            function = eval(obj.co_names[0])
        #except:
        #    raise ValueError
    return function


class TreeNodeSetSelector:
    """
object used to select subsets of objects in a TreeNodeSet.
--
Its select method is called from TreeNodeSet.get method to
select and return a subset of nodes based on 'selectionString' which can be 
either a function or a string.  This method always returns a TreeNodeSet,
possibly empty, of the same type of 'nodes' and  a message string used to indicate 
portions of the selectionString which failed to select any nodes.
--
Here, the type of selectionString can be either
    - function: the function is called for each object in the set. All objects
which the function return True are selected and returned. NOTE: to set the
stringRepr of the subset returned using a lambda expression, enclose
the lambda expression in quotation marks: "lambda x: x.element=='Zn'"
    - string: the selectionString can be a single 'item' or a list of 'items'
Each item in turn can be either a function or a string....
--
Order of string processing:
1: the selectionString is checked to see if it is an integer. If
so, it is used as an index into nodes.
2: the selectionString is checked to see if it is a comma-separated list of items.
The item or items are processed sequentially using a 'TreeNodeSetSelector' of the 
appropriate level to select subsets of nodes which are added together and finally returned.
3. '$' matches the last node in the set: mols.chains.residues.get('$') returns
the last residue in the set.
4. '#' is used for relative indices: for example, mols.chains.residues.get('#1')
returns the first residue in each chain.
5. ranges of nodes are supported and specified as 'exp1-exp2'. A range returns
nodes between the first match of exp1 and the last match of exp2, inclusive.
6. In general, strings are matched to the 'name' attribute of the node using
the python re module.  MOST standard regular expression syntax is supported: 
for example, nodes.get('[OCN]') would select all atoms whose names are O, C or N.
7. THERE ARE CURRENTLY TWO EXCEPTIONS TO STRAIGHTFORWARD USE OF THE RE MODULE:
     -> we support ranges specified using '-' which is replaced by ':'
     -> we support using wildcard '*' by itself. Here "*" is ALWAYS replaced by '.*'
--
***********************************************************************************
***********************************************************************************
To use standard python re syntax, please use the objectsFromStringRE method.
***********************************************************************************
***********************************************************************************
--
TreeNodeSetSelectors derived from this class support a few other types of
items such as NamedSets (eg for ChainSetSelector: proteic, dna) and matchSequence
(for one-letter residue names)
"""

    def __init__(self):
        self.level = TreeNodeSet


    def select(self, nodes, selectionString, sets=None, caseSensitive=True,
               escapeCharacters=False):
        """
        select and return a subset of nodes based on 'selectionString' which can be 
either a function or a string.  This method always returns a TreeNodeSet,
possibly empty, of the same type of 'nodes' and  a message string used to indicate 
portions of the selectionString which failed to select any nodes.
--
Optional input parameters include:
    - sets: an instance of MolKit/Sets class [which is derived from python dict
class] whose keys are strings and whose values are sets.  These sets can be
selected by 'names' matching the keys in this dictionary.
    - caseSensitive: by default, regular expression matching using the python re
module is caseSensitive.  This is changed by setting caseSensitive to False.
[use other than caseSensitive is currently experimental only] 
    - escapeCharacters: possibly supports using '\\' [currently experimental only] 
"""
        assert nodes.__class__==self.level

        self.sets = sets
        self.caseSensitive = caseSensitive
        self.escapeCharacters = escapeCharacters
        
        # if selectionString is empty, return
        zeroSelect = [] # list of items in selection string that select nothing
        if selectionString=='':
            return nodes, zeroSelect
        
        try:
            int(selectionString)
            int_string = True
        except:
            int_string = False 
        # split selection string on commas unless it is callable
        if  callable(selectionString):
            selectionString = str(selectionString)
            selList = [selectionString]
        elif  selectionString in nodes.name and int_string is True:
            #bug discovered from 11206900816.pqr 
            #molecule name which is a number...
            ind = nodes.name.index(selectionString)
            return nodes[ind:ind+1], zeroSelect
        else:
            selList = string.split(selectionString, ',')

        # build a set of selected nodes by selecting w/ each item in list
        selNodes = self.level([])
        #selNodes = nodes.__class__([])
        for item in selList:
            newNodes = self.processListItem(nodes, item, sets)
            if newNodes:
                selNodes = selNodes | newNodes
            else:
                # FIXME remove item from selectionString for stringRepr
                zeroSelect.append(item)
        selNodes.setStringRepr(selectionString)
        return selNodes, zeroSelect

    

    def processListItem(self, nodes, item, sets=None):
        """describe what this function does
$ to select last object
callables that can be objects or strings
ranges specified using -
relative numbers, i.e. the position in the list (1-based)
sets

"""
        #classifies each item in selList and calls appropriate function

        #FIX THIS: what about an id ==' '
        #item = item.strip()
                
        # handle '$' to select last item in set
        if item=='$':
            newNodes = nodes[-1]
            #return nodes.__class__(newNodes)
            return self.level([newNodes])

        #first detect callable objects:
        if callable(item):
            #print "is callable"
            result = filter(item, nodes.data)
            if len(result)==len(nodes.data):
                return nodes
            else: 
                return self.level(result)
        else:
            try:
                func = evalString(item)
                if callable(func):
                    result = filter(func, nodes.data)
                    if len(result)==len(nodes.data):
                        return nodes
                    else:
                        #return nodes.__class__(result)
                        return self.level(result)
            except:
                pass

        # check for ranges that do not contain [ because these are regexp
        if string.find(item, '-')!=-1 and string.find(item, '[')==-1:        
            #call range w/ nodes here:
            newNodes = self.getRange( nodes, item )
            return newNodes

        # handle sets
        if sets and item in sets.keys():
            return sets[item]

        # handle relative numbers
        if item is not None and len(item) and item[0]=='#':
            item = item[1:]
            newNodes = self.getRelativeIndex(nodes, item)
            return newNodes

        # handle numbers
        try:
            item = int(item)
            if len(nodes.data) > item:
                return self.level([nodes.data[item]])
                #?????
                #return self.level([nodes.data[item-1]])
            else:
                return self.level([])
        except:
            pass

        #  else try a regexp
        newNodes = self.regexp(nodes, item)
        return newNodes


    def getRange(self, nodes, item):
        #general range
        if len(nodes)<2:
            return None
        levItList=string.split(item, '-')
        firstNodes = self.processListItem(nodes, levItList[0].strip())
        lastNodes = self.processListItem(nodes, levItList[1].strip())
        if firstNodes and lastNodes:        
            return self.rangeMatch(nodes, firstNodes[0], lastNodes[-1])
        else:
            return None

    
    def rangeMatch(self, nodes, fr, to):
        ##given the two nodes, get a range
        indexfro = nodes.data.index(fr)
        indextoo = nodes.data.index(to)
        #assert indexfro<=indextoo
        if indexfro<=indextoo:
            #return nodes.ReturnType(nodes.data[indexfro:indextoo+1])
            return self.level(nodes.data[indexfro:indextoo+1])
        else: return None


    def regexp(self, nodes, item):
        #use self.procFunction to build a regexp
        #Chain strings match to 'id', so use 'objectsFromStringField'
        if not self.caseSensitive:
            if self.escapeCharacters:
                item = self.processStringcIWEC(item)
                #newNodes = self.processStringcIWEC(item)
            else:
                item = self.processStringcI(item)
        else:
            item = self.processStringcS(item)
        return self.level(nodes.objectsFromString(item))
    

    def processStringcIWEC(self,someString):
        # COMMENT
        import string
        strList = string.split(someString, ',')
        retExp = ''
        specialList = ['?','*','.','$','#', ':', '-']
        numbList = ['0','1','2','3','4','5','6','7','8','9']
        for i in range(len(strList)):
            item = strList[i]
            newExp = ''
            escape = 0
            startbrace = 0
            closebrace = 0
            ctr = 0
            if item[0]=='\\': 
                newExp=newExp+item[1]
                ctr = 2
            for c in item[ctr:]:
                if c == '\\': 
                    escape = 1
                    ###why continue here???
                    #continue 
                if c == '[':
                    startbrace = startbrace + 1
                if c == ']':
                    startbrace = startbrace - 1
                    closebrace = closebrace + 1
                if escape or c in specialList or c in numbList:
                    newExp = newExp+c
                    escape = 0
                else:
                    if startbrace:
                        if c =='^'or c =='[': 
                            newExp = newExp + c
                        else:
                            newExp = newExp + string.upper(c)+ string.lower(c)
                    elif closebrace:
                        newExp = newExp + c
                        closebrace = closebrace -1
                    else:
                        newExp = newExp + '['+string.upper(c)+ string.lower(c)+']'
            retExp = retExp + newExp
            if i < len(strList)-1:
                retExp = retExp +','
        someString = retExp
        return  self.processStringcI(someString)


    def processStringcI(self, someString):
        import string
        someString = string.replace(someString, '?', '.')
        #if there are any commas, preceed them by (?i)
        someString = string.replace(someString, ',','(?i),')
        #in any case add one (?i) at the end
        someString = someString + '(?i)'
        return self.processStringcS(someString)


    def processStringcS(self, someString):
        import string
        #in all cases do these things:
        if type(someString)==types.StringType:
            #protect [A-Z]
            if someString.find(']')==-1:
                someString = string.replace(someString, '-', ':')
            #convert * to .* except for \*
            if someString.find('\*')>-1:
                #mask \* before replacing * with .*
                someString.replace('\*', '\\')
                someString.replace('*', '.*')
                #unmask \\ after replacing * with .*
                someString.replace('\\', '\*')
            else:
                someString = string.replace(someString, '*', '.*')
        return someString


    def getRelativeIndex(self, nodes, item):
        # index into list of siblings, ie nodes with same parent
        if nodes[0].parent==None:  #catch parentless nodes here
            try:
                index = int(item)
                if index<len(nodes):
                    return nodes[index:index+1]
                else:
                    print "invalid index:", item, " for ", len(nodes), " parentless nodes"
                    return  self.level([])
            except:
                print "invalid index:", item, " for parentless nodes"
                return  self.level([])

        number = int(item)
        if number==0:
            print "0 is not a valid relative index: valid relative indices start at 1"
            return self.level([])

        l=[]
        parentNodes = nodes.parent.uniq()
        if max(parentNodes)==None:  #what???, checking for  [None, None, None,...]
            if len(nodes.data)>=number:
                return self.level(nodes.data[number-1])
                #return nodes.__class__(nodes.data[number-1])
        else:
            for item in parentNodes:
                if len(item.children)>=number:
                    l.append(item.children[number-1])

        #return nodes.__class__(l)
        return self.level(l)
    
            

class TreeNodeSet(ListSet):
    """Class to represent a set of nodes from a a tree"""

## save 2 function calls
##      def __init__(self, objects=None, elementType=None):
##          """TreeNodeSet constructor"""

##          ListSet.__init__(self, objects, elementType)

    def __hash__(self):
        return id(self)

    
    def __repr__(self):
        if len(self.data):
            ob = self.data[0]
            if self.stringRepr and len(self.stringRepr)>30:
                strRepr = self.stringRepr[:30]+'...'
            else:
                strRepr = self.stringRepr
            return '<%s instance> holding %d %s, "%s"' %(
                self.__class__.__name__, len(self.data), ob.__class__.__name__,
                strRepr)
        else:
            return '<%s instance> empty'% ( self.__class__.__name__, )


    def getParentsNoStringRepr(self):
        """Returns the TreeNodeSet of all parents but does not build a string
repr.  This should be used when temporary sets not used as command arguments
are built.
"""
        parents = []
        for o in self.data:
            parents.append(o.parent)
        if len(parents)>0:
            return parents[0].setClass(parents)
        else:
            return TreeNodeSet([])

        
    def ReturnType(self, result):
        """Try to be clever about what we return"""
        if result is None or len(result)==0:
            return self.__class__([])
        if len(result) > 0:
            from molecule import BondSet
            if isinstance(result[0], TreeNode):
                return result[0].setClass( result )
            
            elif isinstance(result[0], BondSet):
                flat = []
                for i in result:
                    n = len(flat)
                    flat[n:n+len(i)] = i.data
                return BondSet(flat)
                
            elif isinstance(result[0], TreeNodeSet):
                # we want to flatten a list of TreeNodeSets into a single
                # TreeNodeSet
                if len(result)==1:
                    return result[0]
                # MS TOTAL performance killer but maintains stringrepr
                flat = result[0].copy()
                for set in result[1:]:
                    flat.extend(set)
                return flat
##                 flat = [] # efficient flattening
##                 for i in result:
##                     n = len(flat)
##                     flat[n:n+len(i)] = i.data
##                 return result[0].__class__( flat )
            elif hasattr(result[0], 'setClass') and \
                 issubclass(result[0].setClass, TreeNodeSet):
                # this is use to turn lists of Bonds into BondSets
                return result[0].setClass(result)
            else: # here we could build Numeric arrays for numeric values
                return result


    def __getattr__(self, member):
        #if len(self.data)==0: return self.ReturnType([])
        if member[:2]=='__':
            if self.__dict__.has_key(member):
                return self.__dict__[member]
            else:
                raise AttributeError('member %s not found'%member)

        res = self.ReturnType(ListSet.__getattr__(self, member))
        if len(res)==0:
            return res
        elif res[0]==None: 
            return res  #this returns a list [None,None,None...]
        # make the stringRepr more concise using the fact that member
        # is either parent or children and always selects the whole set
        if member=='parent':
            stringRepr=""
            if len(res.data)>1:
                prev = res[0]
                l = len(prev.children)
                stringRepr = (prev.full_name()+';')*l
                i = l
                while i < len(res):
                    o = res[i]
                    name = o.full_name()
                    l = len(o.children)
                    stringRepr += (name+';')*l
                    i += l
                    if l == 0: break                    
            res.stringRepr = stringRepr

        elif member=='children' or \
                 (self.data[0].__dict__.has_key(member) and \
                  self.data[0].children is self.data[0].__dict__[member]):
            ## if we take all children of all elements in the set we want to
            ## optimize the stringRepr
            ## this happens if the member asked for is 'children' or
            ## the attribute called member in the elements of the set
            ## is the same as obj.children for each obj in the set e.g.
            ## chainSet.residues is the same as chainSet.children
            res.stringRepr = None
            if self.stringRepr:
                res.stringRepr = self.extendStringRepr(self.stringRepr)
            elif verbose:
                import traceback
                traceback.print_stack()
                print 'TreeNodeSet getattr on sets with no stringRepr:', repr(self), member
                res.stringRepr = None
            
        return res


#    def objectsFromStringField(self, selectionString, field='name'):
#        """find and return the list of nodes in this TreeNodeSet whose
#attribute 'field' matches the regular expression '^selectionString$'
#"""
#        return self.objectsFromStringFieldRE('^%s$'%selectionString, field)
#
#
#    def objectsFromStringFieldRE(self, selectionString, field):
#        """find and return the list of nodes in this TreeNodeSet whose
#attribute 'field' matches the regular expression '^selectionString$'
#"""
#        s2 = '^'+ s + '$'
#        prog = re.compile(s2)
#        #prog = re.compile("^%s$"%s)
#        match = eval('filter(lambda x, prog=prog: prog.search(x.%s), self.data)'%f)
#        if len(match)>0:
#            return match
#        else:
#            return None
#

    def objectsFromString(self, selectionString, field='name'):
        """find and return the list of nodes in this set whose name match
the regular expression '^selectionString$'
"""
        #print "start oFS: selectionString=", selectionString
        selectionString = selectionString.replace("+", "\\+") #escape + to handle cases such as 'monomer + agonist.pdb'
        return self.objectsFromStringRE( '^'+selectionString+'$', field )


    def objectsFromStringRE(self, reg, field):
        """find the list of nodes in this set whose attribute 'field' matches
the regular expression 'reg'.  Return the list of matched nodes.
When nothing matches 'reg', return an empty list.
"""
        #print "start oFSRE: reg=", reg
        prog = re.compile("%s"%reg)
        result = []
        for node in self.data:
            if prog.search(getattr(node, field)):
                result.append(node)
        #if len(result):
        #    return result
        #else:
        #    return None
        return result


    def getSelector(self):
        if self.selector==None:
            self.selector = TreeNodeSetSelector()
        return self.selector


    def get(self, selectionString, selector=None, sets=None,
             caseSensitive=True, escapeCharacters=False,
             returnMsg=False):

        """
        select and return the elements of the set based on 'selectionString'. 
The selectionString can be a single 'item' or a list of 'items'
--
Each item can be either a function or a string:
    - function: the function is called for each object in the set. All objects
which the function return True are selected and returned. NOTE: to set the
stringRepr of the subset returned using a lambda expression, enclose
the lambda expression in quotation marks: "lambda x: x.element=='Zn'"
    - string: a 'TreeNodeSetSelector' of the appropriate level is used to select 
and return the elements of the set based on 'selectionString'. (See
TreeNodeSetSelector for more details).
--
optional input parameters include:
    - sets: an instance of MolKit/Sets class [which is derived from python dict
class] whose keys are strings and whose values are sets.  These sets can be
selected by 'names' matching the keys in this dictionary.
    - caseSensitive: by default, regular expression matching using the python re
module is caseSensitive.  This is changed by setting caseSensitive to False.
    - escapeCharacters: possibly supports using '\\' [currently: experimental ONLY] 
    - returnMsg: return a message which is a list of 'items' which do not select anything. 
By default, no message is returned. 
"""
        selector = self.getSelector()
        if type(selectionString) in types.StringTypes:
            result, msg = selector.select(self, selectionString, sets=sets,
                        caseSensitive=caseSensitive,
                        escapeCharacters=escapeCharacters)
            result = self.ReturnType(result)
            selectionStringRepr = '(%s\s\%s)'%(self.stringRepr, selectionString)
            result.setStringRepr(selectionStringRepr)
            if returnMsg:
                result = (result, msg)
            return result
        elif callable(selectionString):
            result = filter(selectionString, self.data)
            if len(result)==len(self.data):
                return self
            else:
                return self.ReturnType(result)
            # FIXME . it would be nice to save a lambda 
            # function here for instance
            #selectionStringRepr = result.full_name()
            #result.setStringRepr(selectionStringRepr)
            #return self.ReturnType(result)
        else:
            raise RuntimeError("argument has to be a function or a string")


    def truncateStringRepr(self, strRepr, num_levels=1, 
                ops=['+','-','s','^','&']):
        if strRepr is not None:
            #eg stringSel:A:/+/1crn: :
            sub_list = strRepr.split('/')
            for j in range(len(sub_list)):
                if sub_list[j] in ops:
                    continue
                for ctr in range(num_levels):
                    ind = sub_list[j].rfind(':')
                    sub_list[j] = sub_list[j][:ind]
            #finally, build back the overall strings with /+/ etc
            newStringRepr = sub_list[0]
            index = 1
            if len(sub_list)>1:
                if sub_list[0]=='':
                    newStringRepr = '/'+sub_list[1] + '/'
                    index = 2
                for item in sub_list[index:]:
                    newStringRepr += '/' + item
        return newStringRepr



    def extendStringRepr(self, strRepr, num_levels=1, 
                    ops=['+','-','s','^','&']):
        #print "strRepr=", strRepr
        newStringRepr = ''
        if strRepr is not None:
            #eg stringSel:A:/+/1crn: :
            sub_list = strRepr.split('/')
            #print "sub_list=", sub_list, " with len=", len(sub_list)
            for j in range(len(sub_list)):
                #print 'j=', j
                if sub_list[j] in ops:
                    #print "skipping ", sub_list[j]
                    continue
                sub_list[j] += ':' * num_levels
            #finally, build back the overall strings with /+/ etc
            newStringRepr = sub_list[0]
            index = 1
            #print "newStringRep=", newStringRepr
            if len(sub_list)>1:
                if sub_list[0]=='':
                    #print "in empty string if"
                    newStringRepr = '/'+sub_list[1] + '/'
                    index = 2
                for item in sub_list[index:]:
                    newStringRepr += '/' + item
        #print "esr: returning ", newStringRepr
        return newStringRepr


    def setLevel(self, what, uniq=1):
        """
        find nodes of the given type in the tree and update the stringRepr. 
        When change is to level above, always do a uniq of the result and 
        truncate the stringRepr before the last colon [:]. When change is to level 
        below, add a colon
        """
        if self.elementType == what:
            return self

        elif len(self) == 0:
            try:
                return what().setClass([])
            except:
                raise RuntimeError ("could not find level of type %s"%what)

        else:
            stringRepr = self.stringRepr
            levelBelow = self[0].isAbove(what)
            if levelBelow>1:
                exec("result=self"+".children"*levelBelow)
                if uniq:
                    result = result.uniq()
                if stringRepr is not None:
                    stringRepr = self.extendStringRepr(stringRepr, levelBelow)
                result.stringRepr = stringRepr
                return result

            elif levelBelow==1:
                result = self.children
                if uniq:
                    result=result.uniq()
                if stringRepr is not None:
                    stringRepr = self.extendStringRepr(stringRepr, levelBelow)
                result.stringRepr = stringRepr
                return result
            else:
                levelAbove = self[0].isBelow(what)
                if levelAbove>1:
                    exec("result=self"+".parent"*levelAbove)
                    if uniq:
                        result = result.uniq()
                    if stringRepr is not None:
                        result.stringRepr = self.truncateStringRepr(stringRepr, levelAbove)
                    return result

                elif levelAbove==1:
                    result = self.parent
                    if uniq:
                        result=result.uniq()
                    if stringRepr is not None:
                        result.stringRepr = self.truncateStringRepr(stringRepr, levelAbove)
                    return result

                else:
                    raise RuntimeError ("could not find level of type %s"%what)
                    


    def findType(self, what, uniq=0):
        """
        Find nodes of the given type in the tree. When searches above
        does always a uniq of the result.
        """
        if self.elementType == what:
            return self

        elif len(self) == 0:
            try:
                return what().setClass([])
            except:
                raise RuntimeError ("could not find level of type %s"%what)
        

        else:
            levelBelow = self[0].isAbove(what)
            if levelBelow>1:
                exec("result=self"+".children"*levelBelow)
                if uniq:
                    result = result.uniq()
                return result

            elif levelBelow==1:
                result = self.children
                if uniq:
                    result=result.uniq()
                return result
            else:
                levelAbove = self[0].isBelow(what)
                if levelAbove>1:
                    exec("result=self"+".parent"*levelAbove)
                    if uniq:
                        result = result.uniq()
                    return result

                elif levelAbove==1:
                    result = self.parent
                    if uniq:
                        result=result.uniq()
                    return result

                else:
                    raise RuntimeError ("could not find level of type %s"%what)
                    

    def findChildrenOfType(self, what):
        """for a set of nodes, go down the tree until we find Nodes of the
        given type and merge the result"""
        if len(self.data)==0: return self
        if isinstance(self.data[0], what): return self
        nodeSet = self.findType(what)
        for node in self.data[1:]:
            nodes = node.findType(what)
            nodeSet = nodeSet + nodes
        return nodeSet


    def findParentsOfType(self, what):
        """for a set of nodes, go up the tree until we find Nodes of the
        given type and merge the result"""
        if len(self.data)==0:
            return self
        if isinstance(self.data[0], what):
            return self
        nodeSet = self.findType(what)
        for node in self.data[1:]:
            nodes = node.findType(what)
            nodeSet = nodeSet + nodes
        return nodeSet
        

    def NodesFromName(self, nameStr):
        """
retrieves nodes from the tree using a name.
Name is a string as produced by
full_name() i.e. node names separated by ':' going from root to leaf
comma ',' separated list of names are allowed as well as
range specified using the '-' character   
"""
        from MolKit.stringSelector import StringSelector
        selector = StringSelector()
        result, msg = selector.select(self, nameStr)
        return result
##        result = TreeNodeSet()
##        # if the nameStr end with the ';' character we remove it, else it
##        # expands the empty string after the last ; to all molecules
##        if len(nameStr)>1:
##            if nameStr[-1]==';':
##                nameStr = nameStr[:-1]
##        for name in string.split(nameStr, ';'):
##            names = string.split(name, ':')
##            node = self
##            # narrow down selection using names at successive levels
##            for name in names:
##                nodes = node.get( name )
##                if nodes is not None and len(nodes)>0: 
##                    node = nodes.children
##                else:
##                    return result
##            if nodes:
##                result = result + nodes

##        result.setStringRepr(nameStr)
##        if len(result)>0:
##            return result[0].setClass(result)
##        else:
##            return result
    

    def full_name(self, useShortCut=1):
        """Build a string representation of the nodeset by concatenating
        names up to the root. The last level can be a list separated by comma
        If elements of set belong to more than 1 parent, a unique string is
        produced for each object and they are concatenated using a semi-colon
        ';' as a separator.
        
        ShortCut is a representation of self, a TreeNodeSet, which replaces an
        explicit naming of each node in each level by ':', a list slicing
        operator which returns the entire level. ShortCut can be used when
        all self.data has the same top and all children of self.data.top[0] 
        of every level down to and including self.elementType are in self.data. 
        useShortCut is a parameter which allows user to specify which type of 
        string to return. 

        for example, in the case of a TreeNodeSet containing all Residues in 
        1crn.pdb:
            ShortCut-> 
                    '1crn::'
            no ShortCut -> 
                    '1crn: :,CYS16,GLY42,CYS40,ALA9,CYS26,PRO36,ILE25,
            SER6,ARG17,ASP43,ALA45,CYS4,ILE35,ALA38,ALA24,PRO5,PRO41,THR28,
            GLY31,THR1,GLU23,PRO19,VAL8,ASN12,ILE33,ILE7,PHE13,TYR29,THR39,
            TYR44,ASN46,ASN14,THR2,CYS3,THR21,CYS32,ILE34,GLY37,LEU18,PRO22,
            ALA27,SER11,THR30,GLY20,ARG10,VAL15'
        """

        # empty set
        if len(self)==0: return ''
        if len(self.data)==0: return ''

        # all elements in the set have the same parent. If the number of
        # objects in the set is less then the number of children of top
        # object we return the full name of the parent followed by a comma
        # separated list of children names. Else we only return the full name
        # of the parent with one ':' for each level self.elementType is below
        # top
        if self.data[0].parent!=None:
            if useShortCut:
                a = self.data[0]
                all = TreeNodeSet([a.top]).findChildrenOfType(self.elementType) 
                if len(self)==len(all):
                    i=1
                    b = a
                    while b.parent != a.top:
                        i = i+1
                        b = b.parent
                    return a.top.name+':'*i

            parents = self.parent.uniq()
            if len(parents) == 1:
                name = parents[0].full_name()
                name = name + ":"+self.data[0].name
                for o in self.data[1:]: name = name+','+o.name
                return name

        #special case for molecules        
        if self.data[0].top==self.data[0]:
            name = self.data[0].name
            for o in self.data[1:]:
                name = name + ',' + o.name
            return name
            

        # else objects in set have no parent or different parents. We return
        # a list of full_name()
        name = self.data[0].full_name()
        for o in self.data[1:]:
            name = name + ';' + o.full_name()
        return name


    def sort(self, func = None):
        if len(self.data)==0: return
        
        if func is None:
            self.data.sort(self.data[0].compare)
        else:
            self.data.sort(func)


##      def compare(self, tn1, tn2):
##          """compare 2 tree nodes of a TreeNode Set by comparing there
##          uniqIndex"""
##          return cmp( tn1._uniqIndex, tn2._uniqIndex)


        
## class TreeNode(object):
class TreeNode:
    """Base class for node that can be used to build trees.
    Every node is by definition a list of its children.mol.
    a leaf is an empty list.
    a root has no parent.
    """
    _numberOfDeletedNodes  = 0

    def __str__(self):
        return repr(self)

##     def __del__(self):
##         self.__class__._numberOfDeletedNodes = self.__class__._numberOfDeletedNodes + 1

    
    def __init__(self, name='NoName', parent=None, elementType=None,
                 objects=None, childrenName=None, setClass=TreeNodeSet,
                 childrenSetClass=TreeNodeSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """TreeNode constructor.
        Arguments:
        optional objects (list)
        optional name (string)
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode)
                 represents the type of the children of that node
        optional childrenName (string) alias for children attribute
        optional setClass (class type) type of object used for selections
                 sets of such nodes
        optional childrenSetClass (class type) type of object used for
                 a selection of children of that node
        optional top (TreeNode) root of the tree
        """

##         assert TreeNodeSet in childrenSetClass.__bases__ or \
##                childrenSetClass is TreeNodeSet
        assert issubclass(childrenSetClass, TreeNodeSet) or \
               childrenSetClass is TreeNodeSet

##         self.children = TreeNodeSet(objects, elementType)
##         self.children.__class__ = childrenSetClass
        
        self.children = childrenSetClass(objects)
        self.setClass = setClass
        self.childrenSetClass = childrenSetClass
        self.childrenName = childrenName
        if childrenName:
            setattr(self, childrenName, self.children)
        self.name = name
        self.elementType = elementType
        self.parent = parent
        if parent is not None:
            parent.adopt(self, childIndex, assignUniqIndex)
#        if top:
#            assert isinstance(top, TreeNode)
        if top is None:
            self.top = self
        else:
            assert isinstance(top, TreeNode)
            self.top = top
        # WARNING this assumes the name of each child is unique
        #    also, this creates an additional reference to the object whcihc
        # should be handled at deletion time
        self.childByName = {} # {'name': childnode}


    def deleteSubTree(self):
        """ Function to actually delete all the reference to a TreeNode to
        Free the memory !!!!"""
        import sys
        # 1- Find the levels below the node to delete.
        #print '1', sys.getrefcount(self)
        levelsBelow = []
        level = self
        while 1:
            level = level.children
            levelsBelow.append(level[0].__class__)
            if not hasattr(level, 'children') or len(level.children)==0: break
            
        levelsToDelete = levelsBelow[:-1]
        # only the levels with children attribute are considered.
        # 2- Loop on the lowest levels and delete all the childrens
        #    then the nodes.
        levelsToDelete.reverse()
        if levelsToDelete:
            for lev in levelsToDelete:
                levelNodes = self.findType(lev)
                for i in xrange(len(levelNodes)):
                    while len(levelNodes[i].children)!=0:
                        levelNodes[i].children[0].__dict__.clear()
                        del(levelNodes[i].children[0])
                del(levelNodes)
        
        while len(self.children)!=0:
            self.children[0].__dict__.clear()
            del(self.children[0])
        return levelsBelow
        

    def adopt(self, child, index=None, assignUniqIndex=1, setChildrenTop=0):
        """Have a parent node adopt a child node"""

        if self.elementType is not None:
            assert isinstance(child, self.elementType)
        child.parent = self
##          if len(self.children)==0:
##              child._uniqIndex = 0
##          else:
##              child._uniqIndex = self.children[-1]._uniqIndex+1
        child._uniqIndex = len(self.children)
        if index is None:
            self.children.append(child)
        else:
            self.children.insert(index, child)
            if assignUniqIndex: self.assignUniqIndex()

        #also correct self.childrenName if it exists:
        if self.childrenName!=None:
            setattr(self, self.childrenName, self.children)
##            exec('self.'+self.childrenName+'=self.children')

        child.top = self.top
        if setChildrenTop:
            #possibly set top of any children of child to self.top
            grandChildren=child.children
            while len(grandChildren):
                grandChildren.top = self.top
                grandChildren = grandChildren.children

        self.childByName[child.name] = child


    def makeNameUniq(self, aliasList=None):
        """ None <- makeNameUniq() 
            checks that each name in self.children is unique.
            concatenates _uniqIndex to name when identical names are found"""

        for i in range(len(self.children)): 
            child=self.children[i]
            for j in range(i+1, len(self.children)):
                item = self.children[j] 
                if child.name == item.name:
                    item.name=child.name +'_'+ str(item._uniqIndex)
                    for alias in aliasList:
                        exec('item.'+alias+'=item.name')


    def assignUniqIndex(self):
        i = 0
        for c in self.children:
            c._uniqIndex = i
            i = i+1


    def remove(self, child, assignUniqIndex=1, cleanup=0):
        """remove a child"""

        assert child in self.children
        self.children.remove(child)
        #also correct self.childrenName if it exists:
        if self.childrenName!=None:
            setattr(self, self.childrenName, self.children)
#            exec('self.'+self.childrenName+'=self.children')
        #commented in next three lines 5/9
        if not len(self.children) and cleanup and self.parent:
            self.parent.remove(self, cleanup=cleanup)
            return
        if assignUniqIndex: self.assignUniqIndex()

    # FIXME added for testing under 2.2.1
    # should disappear in the future when we use Properties
    def __eq__(self, other):
        """compare 2 tree nodes by comparing their address in memory"""
        return id(self) == id(other)

    def __ne__(self, other):
        """compare 2 tree nodes by comparing their address in memory"""
        return id(self) != id(other)

    def __cmp__(self, other):
        """compare 2 tree nodes by comparing their address in memory"""
        return cmp( id(self), id(other))


    # this method is required because we implemented __cmp__
    # if missing, these objects are unhashable
    def __hash__ (self):
        """return a hash value for this object"""
        return id(self)


    def compare(self, one, other):
        """compare 2 tree nodes"""

        if other is None: return 1
        if id(one) == id(other): return 0
        
        # treat top level separately
        if id(one.top) < id(other.top): return -1
        elif id(one.top) > id(other.top): return 1

        # if we get here, we know one and other are in the same top
        # build a list of parents of one
        if one.parent:
            parent1 = [one]
            p = one.parent
            while p.parent:
                parent1.append(p)
                p = p.parent
        else:
            parent1 = []
            
        # build a list of parents of other
        if other.parent:
            parent2 = [other]
            p = other.parent
            while p.parent:
                parent2.append(p)
                p = p.parent
        else:
            parent2 = []

        # loop over parents from top to bottom (molecule->atoms)
        l1 = len(parent1)
        l2 = len(parent2)
        if l1 < l2: return -1
        elif l1 > l2: return 1
        else:
            for i in range(l1-1, -1, -1):
                p1 = parent1[i]._uniqIndex
                p2 = parent2[i]._uniqIndex
                if p1 < p2: return -1
                elif p1 > p2: return 1

    def __repr__(self):
        if hasattr(self, 'children') and  hasattr(self.children, 'data') and\
            len(self.children.data) > 0:
            return "<%s instance> %s with %d %s" %(self.__class__.__name__,
                      self.full_name(), len(self.children), self.elementType)
        else:
            return "<%s instance> %s" % (self.__class__.__name__,
                                      self.full_name())

    def getRoot(self):
        """returns the root of the tree this node belongs to"""
        n = self
        while hasattr(n, 'parent'): n = n.parent
        return n


    def getParentOfType(self, what):
        """return first parent of a givent type"""
        n = self
        while hasattr(n, 'parent') and not isinstance(n, what):
            n = n.parent
        if not hasattr(n, 'parent'):
            raise RuntimeError ("node %s has no parent of type %s" % (self,
                                                                      what))
        return n


    def findLevels(self, lastLevel = None):
        """goes down the tree until the children member is empty and return
        the class type of objects in children"""
        n = self
        cla = [self.__class__]
        while len(n.children) > 0:
            if lastLevel is not None and isinstance(n, lastLevel): break
            n = n.children[0]
            cla.append(n.__class__)
        return cla

    
    def findType(self, _what, uniq=0):
        """go down the tree until we find Nodes of the given type"""
        # Create the corresponding setClass...
        n = self.setClass([self])
        if n.elementType == _what: return n
        result = n.findType(_what, uniq=uniq)
        return result
    

    # FIXME .. should be more clever about going down
    # and try every member of type TreeNode
    def isAbove(self, Klass):
        """go down the tree until we find Nodes of the given type
        return the number of level above self at which we find Klass, else
        0 is retuned"""
        #assert type(Klass) is types.ClassType
        assert issubclass(Klass, TreeNode)
        n = self
        l = 1
        while (1):
            if len(n.children)==0: return 0
            if n.children[0].__class__ == Klass: return l
            n = n.children[0]
            l = l+1
            

    # FIXME .. should be more clever about going up
    # and try every member of type TreeNode
    def isBelow(self, Klass):
        """go up the tree until we find Nodes of the given type
        return the number of level above self at which we find Klass, else
        0 is retuned"""
        assert issubclass(Klass, TreeNode)
##         assert type(Klass) is types.ClassType
        n = self
        l = 1
        while (1):
            if n.parent is None:
                return 0
            if n.parent.__class__ == Klass: return l
            n = n.parent
            l = l+1

    def get(self, function):
        """select elements among the children of that node using a lambda
        function"""

        return self.children.get(function)


    def dump(self):
        """print out all members and their values"""
        for item in self.__dict__.items():
            if type(item[1])==types.ListType and len(repr(item[1]))>60:
                st = repr(item[1][0])
                s = "List of %d %s " % (len(item[1]), type(item[1][0]))
                s = s + st[:min(len(st), 56-len(s))] + ' ...'
            elif type(item[1])==types.TupleType and len(str(item[1]))>60:
                st = repr(item[1][0])
                s = "Tuple of %d %s " % (len(item[1]), type(item[1][0]))
                s = s + st[:min(len(st), 56-len(s))] + ' ...'
            else:
                s = repr(item[1])
            print "%-20s %-59s" % (item[0], s)


    def full_name(self):
        """Build the node's name by concatenating all names up to the root"""
        name = self.name
        a = self
        while a.parent is not None:
            a = a.parent
            name = a.name+':'+name
        return name
        

    def NodesFromName(self, name):
        """retrieves nodes using a name. Name is a string as produced by
        full_name() i.e. node names separated by ':' going from root to leaf
        """
        names = string.split(name, ':')
        if self.name!=names[0]:
            return None
        node = self
        for name in names[1:]:
            nodes = node.get( name )
            node = nodes.children
        return nodes

            
    def getPrevious(self):
        """ Get the previous TreeNode in a TreeNodeSet"""
        elementSet = self.parent.children
        if elementSet.index(self) != 0 :
            previous = elementSet[elementSet.index(self)-1]
        else:
            previous = None
        return previous


    def getNext(self):
        """Gets the next TreeNode in a TreeNodeSet"""
        elementSet = self.parent.children
        if elementSet.index(self) != len(elementSet)-1:
            next = elementSet[elementSet.index(self)+1]
        else:
            next = None
        return next
    

    def merge(self, right):
        """
        Merges two tree objects by creating a new object that is a 
        copy of the first object but has the children of both of the
        merged objects as its children.
        """
        # Right now it merges the tree given in argument into the left one.
        import copy
        assert self.__class__ == right.__class__
        # Create the new object as a copy of self
        # This is a shallow copy so self at the end will also be the merged
        # result. Unless we do a deep copy this is not useful.
        # the deepcopy doesn't work on a tree instance.
        
        #new = copy.copy(self)
        # the new object then adopts the children of the right tree.
        for item in right.children:
            #new.adopt(item)
            self.adopt(item)
        # The new tree is then returned.
        #return new


    def _copyNode(self, node, copyDict, nameExt):
        import copy
        if node not in copyDict.keys():
            newcopy = copy.copy(node)
            newcopy.name = node.name + nameExt
            #set children of newcopy to []
            #nb: never copy atoms
            newcopy.children = node.children.__class__([])
            copyDict[node] = newcopy
        else:
            newcopy = copyDict[node]
        if node.__class__!=self.__class__:
            if node.parent:
                self._copyNode(node.parent, copyDict, nameExt)
            else:
                newcopy.top = newcopy
        else:
            # for self level, do adoption of copy of self here
            #if split at residue level: chain adopts newcopy
            if node.parent and newcopy not in node.parent.children:
                node.parent.adopt(newcopy)
                #if node==self:
                #    print 'parent got newcopy'
            elif not node.parent: #at top or protein level
                newcopy.top = newcopy
                #if node==self:
                #    print 'top got newcopy'
        cDkeys = copyDict.keys()
        for item in node.children:
            if item in cDkeys:
                itemcopy = copyDict[item]
                if itemcopy not in newcopy.children:
                    newcopy.adopt(copyDict[item])
                itemcopy.top = newcopy.top
                #only remove the lowest level; back in split method


    def split(self, nodes, nameExt='_copy1'):
        copyDict = {}
        for parent in nodes.parent.uniq():
            self._copyNode(parent, copyDict, nameExt)
        #the last step is to reparent the nodes, the bottom level of newcopy
        for node in nodes:
            nodeParent = copyDict[node.parent]
            node.parent.remove(node, cleanup=1)
            nodeParent.adopt(node)
        #need to repair top link in children below nodes:
        top = nodes[0].top
        while len(nodes.children):
            nodes = nodes.children
            nodes.top = top
        return copyDict[self]

