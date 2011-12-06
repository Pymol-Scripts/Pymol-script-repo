#############################################################################
#
# Author: Ruth Huey, Michel Sanner
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################
#
# $Id: stringSelector.py,v 1.41 2008/11/21 20:02:01 vareille Exp $
#


from MolKit.tree import TreeNodeSet, evalString
from MolKit.molecule import MoleculeSet, AtomSet
from MolKit.protein import ProteinSet, ChainSet, ResidueSet
import types, string
from mglutil.util.misc import isInstance

class StringSelector:
    # add outer loop to split on ;
    # handle a:b:c:d;e:f;g;


    def select(self, nodes, selectionString, sets=None, caseSensitive=True,
               escapeCharacters=False):
        #print "in select with selectionString=", selectionString
        overallresults = None
        overallmsg = ""
        overall_class = None
        #eg: stringSel:A:PRO1:N;stringSel:B:PRO1:CA;
        if len(selectionString) and selectionString[-1]==';':
            selectionString = selectionString[:-1]
        #eg: stringSel:A:PRO1:N;stringSel:B:PRO1:CA
        ## mol1;mol2;mol3:::
        allSelectionStrings = selectionString.split(';')
        #print "allSelectionStrings=", allSelectionStrings
        all_tops = nodes.top.uniq()
        final_str_repr = ""  #keep track of string to assign at end
        for selString in allSelectionStrings:
           # print "processing selString=", selString
            lambda_index = selString.find('lambda ')
            if lambda_index == -1:
                setsStrings = selString.split(':')
            else:
                #split repair possible splits of lambda expressions
                #":::lambda x:x.top.name=='ind'"
                #lambda_index = 3
                setsStrings = selString[:lambda_index].split(':')
                #setsStrings = ['','','','']
                #replace the last one with join of split on lambda ':'
                lambda_exp_list = selString[lambda_index:].split(':')
                lambda_exp = lambda_exp_list[0]+':'+lambda_exp_list[1]
                setsStrings[-1] = lambda_exp
                if len(lambda_exp_list)>2:
                #there may be trailing stuff
                    setsStrings.extend(lambda_exp_list[2:]) 
                #setsStrings = selString.split(':')
            #print "setsStrings=", setsStrings
            len_setsStrings = len(setsStrings)
            results = all_tops
            msg = ''
            for i in range(len_setsStrings):
                #print "setsStrings[",i,"]=", setsStrings[i]
                if i==0 and len_setsStrings==1 and setsStrings[i]=='':
                    #special case of empty selection string
                    final_str_repr = ''
                    for m in all_tops:
                        final_str_repr += m.name + ','
                        final_str_repr = final_str_repr[:-1]
                elif setsStrings[i]!='':
                    #print "in elif loop"
                    these_sets = None
                    if sets is not None:
                        these_sets = sets.get(stype=results.__class__)
                    results, msgList = results.get(setsStrings[i], sets=these_sets,
                                         caseSensitive=caseSensitive,
                                         escapeCharacters=escapeCharacters,
                                         returnMsg=True)
                    #print "results=", results, "msgList=", msgList
                    for m in msgList:
                        setsStrings[i].replace(m, '')
                        msg = msg + str(m)
                    if len(results)==0:
                        #print "no results with ", setsStrings[i]
                        overallresults = None
                        break
                #check whether the rest of the setsStrings are empty, if so
                final_str_repr += setsStrings[i] + ':'
                results.stringRepr = final_str_repr
                if i<len(setsStrings)-1:
                    results = results.children
            final_str_repr = final_str_repr[:-1]
            results.stringRepr = final_str_repr
            if len(results):
                if overallresults is None:
                    overallresults = results
                    overall_class = results.__class__
                    overallmsg = msg
                else:
                    #levels are the same
                    if overall_class==results.__class__:
                        overallresults = overallresults + results
                        #does this happen automatically?
                        overallmsg = overallmsg + msg
                    else:
                        print "ERROR"
                        print "overall_class->", overall_class
                        print "results.__class__->", results.__class__
                        print "results=", results
                        raise RuntimeError
                if selString!=allSelectionStrings[-1]:
                    final_str_repr += ';'
                    final_str_repr.replace("/+/", ';')
                else:
                    overallresults.stringRepr =  final_str_repr

        if overallresults is None:
            overallmsg = selectionString
            if len(setsStrings)==1:
                overallresults = ProteinSet()
            elif len(setsStrings)==2:
                overallresults = ChainSet()
            elif len(setsStrings)==3:
                overallresults = ResidueSet()
            else:
                overallresults = AtomSet()

        return overallresults, overallmsg

            
class CompoundStringSelector:


    def check_for_lambda_expressions(self, selectionString):
        #special treatment for (expr\s\lambda x: len(x.bonds)==1)
        #check if lambda is present
        start_index = selectionString.find('lambda')
        if start_index==-1:
            #print "-1: c_f_l_e: returning ", selectionString
            return selectionString
        #convert  (expr\s\lambda x: len(x.bonds)==1) to 
        #..(expr\s\lambda x: len#x.bonds#==1)..
        ##??how to find the end of the lambda expression
        #case 1. lambda x: len#x.bonds#==1)\s\....
        end_index = selectionString[start_index:].find("\\s\\")
        #case 2. lambda x: len(x.bonds)==1)
        ss = selectionString[start_index:end_index]
        #case 2b. convert '(' and ')' to #
        ss1 = ss.replace('(', '##')
        ss2 = ss1.replace(')', '###')
        #print "c_f_l_e: returning ", ss2, ' and ',  True
        ss3 = selectionString[:start_index] + ss2 + ')'
        #print "final selectionstring=", ss3
        return ss3
         

    # use a StringSelector to select each chunk between operators
    def get_nodes(self, nodes, selectionString, sets=None):
        #print "\n\nin get_nodes with nodes=", nodes , " and selectionString=", selectionString
        #print "gn: selList = ", self.parse(selectionString)
        #(()/s/expr)...recurse innermost nested ()
        # which is between first_index and second_index
        msg = ''
        if selectionString.find('lambda')>-1 and selectionString.find("###")==-1:
            selectionString = self.check_for_lambda_expressions(selectionString)
        second_index = selectionString.find(')')
        first_index = selectionString.rfind('(')
        #li<0,ri<0;li>0,ri>0...ok
        #li<0,ri>0;li>0,ri<0...error
        if (second_index<0 and first_index>0) or (second_index>0 and first_index<0):
            raise ValueError, '%s badly nested selection string:'%selectionString
        elif second_index>-1 and first_index>-1:
            r_i = selectionString[:second_index].rfind('(')
            sub_sel_string = selectionString[r_i+1:second_index]
            these_nodes, msg = self.get_nodes(nodes, sub_sel_string, sets=sets)
            if first_index==0 and second_index==len(selectionString)-1:
                return these_nodes, msg
            else:  # use these_nodes to deal with the rest of it
                left_part = selectionString[:first_index]
                right_part = selectionString[second_index:]
                #print "pp: selectionString=", selectionString
                #print "pp: first_index=", first_index
                #print "pp: second_index=", second_index
                #print "pp: right_part=", right_part
                #print "pp: left_part=", left_part
                #print " calling pp with these_nodes=", len(these_nodes)
                #print " calling pp with right_part=", right_part
                #print " calling pp with left_part=", left_part
                return self.process_parts(these_nodes, right_part, left_part, sets=sets)
        else:   #case where second_index==-1 and first_index==-1:
            # there are NO parentheses here so just
            # do the subselect with this string and return nodes
            #print "gn: selectionString=", selectionString
            if selectionString.find("###")>0:
                #print "replacing parantheses"
                selectionString = selectionString.replace("###", ')')
                selectionString = selectionString.replace("##", '(')
                #print "Now selectionString is ", selectionString
            xx = selectionString.split('\\s\\')
            #print "xx=", xx
            selected, msg = self.process(nodes, xx[0], sets=sets)
            #print "selected=", selected
            if len(xx)==1:
                # no subset: finished processing
                return selected, msg 
            elif not len(selected):  #could be empty here
                return None, msg     #??proper return value??
            final = selected.get(xx[1])
            #print "returning len(final)=", len(final)
            return final, msg

    def do_op(self, nodes, selectionString, sets=None):
        #print " in do_op with ", nodes,  " and  ", selectionString
        msg_to_return = ""
        setsStrings = selectionString.split('/')
        #print "setsStrings=", setsStrings
        # create a StringSelector to handle each stringsel
        getSet = StringSelector()
        result, msg = getSet.select(nodes, setsStrings[0], sets=sets)
        #print "result =", result
        msg_to_return += msg
        result = result.copy()
        #print "do_op: setsStrings=", setsStrings
        if len(setsStrings)==1:
            #print "returning len(result)=", len(result)
            #print "returning msg_to_return=", msg_to_return
            return result, msg_to_return
        # initialize the stringRepr
        stringRepr = setsStrings[0]        
        for i in range(1, len(setsStrings), 2):
            op = setsStrings[i]
            arg = setsStrings[i+1]
            #check for subselect '(' here???
            tmp, msg = self.get_nodes(nodes, arg)
            #print "tmp=", tmp
            msg_to_return += msg
            if op=='|':
                result += tmp
                result = result.uniq()
            elif op=='-':
                result -= tmp
            elif op=='&':
                result &= tmp
            elif op=='+':
                result += tmp
            elif op=='^':
                result ^= tmp
            else:
                raise ValueError, '%s bad operation in selection string'%op
            #PROCESS the stringRepr
            # only add stringsel if something was selected
            if len(tmp):
                if len(stringRepr):
                    stringRepr += '/%s/'%op+arg
                elif op=='+':
                    stringRepr += arg
                else:
                    print "about to raise RuntimeError on ", op
                    raise RuntimeError('ERROR: selection string starting with operator which is not +, (%s)'%op)
            else: #add arg to msg_to_return
                msg_to_return += arg
        result.setStringRepr(stringRepr)
        return result, msg_to_return


    def process(self, nodes, criteria, sets=None):
        #print "in process with len(nodes)=", len(nodes), ' and criteria=', criteria
        #print "type(criteria)==", type(criteria)
        msg_to_return = ""
        if type(criteria)==types.StringType:
            selected, msg = self.do_op(nodes, criteria, sets=sets)
            #print " back in process with len(selected) = ", len(selected)
            msg_to_return += msg
        #elif type(criteria)==types.InstanceType:
        elif isInstance(criteria) is True:
            selected = criteria
        return selected, msg_to_return
            
        
    def process_parts(self, nodes, right_part, left_part, sets=None):
        #print "in pp:nodes=", nodes
        #print "in pp:right_part=", right_part
        #print "in pp:left_part=", left_part
        right_part = right_part.replace("###",')')
        right_part = right_part.replace("##",'(')
        msg_to_return = ""
        xx = right_part[1:].split('\\s\\')
        #print "xx=", xx
        if xx[0]==':':
            nodes = nodes.children
        elif xx[0]=='::':
            nodes = nodes.children.children
        elif xx[0]==':::':
            nodes = nodes.children.children.children
        result = nodes
        for item in xx[1:]:
            if item[-1]==')':
                item = item[:-1]
            result = result.get(item)
            #print "item=", item
            #print "result=", len(result)
            if not len(result):
                msg_to_return += item
        return result, msg_to_return


    def parse(self, selectionString):
        ##print "PARSE: ", selectionString
        sub_select_list = []
        op = ""
        last_piece = 0
        l_ct = 0
        r_ct = 0
        ss = selectionString
        for i in range(len(ss)):
            if ss[i]=='(':
                l_ct += 1
                if len(op):
                    op_list = op.split('/')
                    sub_select_list.append(op_list)
                    #print "op:", op
                    #print "now ss = ", ss[last_piece:]
                    op = ""
                    #print "set last_piece here"
                    last_piece = i
                elif last_piece<i and l_ct==0 and r_ct==0:
                    #check for trailing op
                    pp = ss[last_piece:i].split('/')
                    ##print "pp = ", pp
                    if len(pp)>1:
                       sub_select_list.append(pp[0])
                       sub_select_list.append(['' , pp[1]])
                    else:
                        ##print "adding ", ss[last_piece:i]
                        sub_select_list.append(ss[last_piece:i])
                    last_piece = i
                elif i>0 and last_piece==0 and r_ct==0 and l_ct==1:
                    #check for leading op
                    pp = ss[last_piece:i].split('/')
                    ##print "2: pp = ", pp
                    if len(pp)>1:
                       sub_select_list.append(pp[0])
                       sub_select_list.append(['' , pp[1]])
                    else:
                        ##print "adding ", ss[last_piece:i]
                        sub_select_list.append(ss[last_piece:i])
                    last_piece = i
            elif ss[i]==')':
                r_ct += 1
            #elif ss[i]==':' and l_ct==0 and l_ct==r_ct:
            #    sub_select_list.append(':')
            #    last_piece = i
            if l_ct>0 and l_ct==r_ct:
                #print "2:"
                sub_select_list.append(ss[last_piece:i+1])
                last_piece+=1
                r_ct = 0
                l_ct = 0
            elif last_piece>0 and r_ct==0 and l_ct==0:
                op+= ss[i]
        ##print "last_piece=", last_piece 
        if last_piece==0:
            sub_select_list.append(selectionString)
        ##print "sub_select_list=", sub_select_list 

        # parse selectionString into subselect (...) sections 
        # and ops
        # string/op/string
        # or
        # (....)* /op/(....)*
        # eg (....)/op/(....)/op/....
        # or ???
        return sub_select_list


    def select(self, nodes, selectionString, returnMsg=False, sets=None):
        #split string into list of lists and strings
        #lists result from string/op/string....
        #strings result from subselects (string\\s\\string)...
        token_list = self.parse(selectionString)
        if not len(token_list):
            return nodes, ''   #???
        # create a StringSelector to handle string tokens 
        getSet = StringSelector()
        ###initialize result from first item in token_list
        ##also initialize the return msg
        msg_to_return = ''
        first_item = token_list[0]   
        first_type = type(first_item)
        if first_type==types.StringType:  #a subselect '(.....)'
            result, msg = self.get_nodes(nodes, first_item, sets=sets)
        elif first_type==types.ListType:  # string/op/-> [string, op, ....]
            result, msg = getSet.select(nodes, first_item[0], sets=sets)
            op = first_item[1]
            token_list.insert(1, op)
        else:
            print " INVALID selection string ", selectionString
            return None, ''
        msg_to_return += msg
        result = result.copy()
        stringRepr = result.stringRepr

        #PROCESS THE REST OF token_list
        ##if len(token_list[1:]):
            ##print "for selectionString=", selectionString
            ##print "  the rest of the token_list=", token_list[1:]
        for i in range(1, len(token_list), 2):
            ##print "###########################"
            ##print "IN FINAL FOR LOOP", i, ":", token_list[i]
            ##print "###########################"
            op = token_list[i]
            #after the first, can only have ['',op,string] lists
            if type(op)==types.ListType:
                assert op[0]==''
                op = op[1]
            arg = token_list[i+1]
            nested = self.parse(arg)
            #???
            #print "nested=", nested
            if len(nested):
                tmp, msg = self.get_nodes(nodes, arg)
            msg_to_return += msg
            if op=='|':
                result += tmp
                result = result.uniq()
            elif op=='-':
                result -= tmp
            elif op=='&':
                result &= tmp
            elif op=='+':
                result += tmp
            elif op=='^':
                result ^= tmp
            else:
                raise ValueError, '%s bad operation in selection string'%op

            #PROCESS the stringRepr
            # only add stringsel if something was selected
            if len(tmp):
                if len(stringRepr):
                    stringRepr += '/%s/'%op+arg
                elif op=='+':
                    stringRepr += arg
                else:
                    print "about to raise RuntimeError on ", op
                    raise RuntimeError('ERROR: selection string starting with operator which is not +, (%s)'%op)
        if len(result):
            result.setStringRepr(stringRepr)
        #print "end of select: returnMsg=", returnMsg
        return result, msg_to_return


class MVStringSelector:


    def __init__(self, moleculeSet, selString, userPref = 'cS'):
        #print "moleculeSet=", moleculeSet;print "selString=", selString
        self.form = None
        self.userPref = userPref
        #make a list of the 4 lists
        self.selList = []
        for item in selString:
            itemList = string.split(item, ',')
            self.selList.append(itemList)
        self.procFunction = eval('self.processString%s' %userPref)
        #???should these be initialized here?
        self.moleculeSet=moleculeSet
        self.molSet=None
        self.chainSet=None
        self.resSet=None
        self.atomSet=None
        #self.sets__ = sets
        self.buildFDs()


    def buildFDs(self):
        self.molFD={}
        self.molFD['range']=self.getMolRange
        self.molFD['regexp']=self.getMolMatch
        self.molFD['relative']=self.getMolIndex
        self.molFD['index']=self.getMolIndex
        self.molFD['NamedResSet'] = self.doNothing
        self.molFD['NamedAtomSet'] = self.doNothing
        self.chainFD={}
        self.chainFD['range']=self.getChainRange
        self.chainFD['regexp']=self.getChainMatch
        self.chainFD['relative']=self.getChainRelIndex
        self.chainFD['index']=self.getChainIndex
        self.chainFD['NamedResSet'] = self.doNothing
        self.chainFD['NamedAtomSet'] = self.doNothing
        self.resFD={}
        self.resFD['range']=self.getResidueRange
        self.resFD['regexp']=self.getResidueMatch
        self.resFD['relative']=self.getResidueRelIndex
        self.resFD['index']=self.getResidueIndex
        self.resFD['NamedResSet'] = self.getNamedResSet
        self.resFD['NamedAtomSet'] = self.doNothing
        self.atomFD={}
        self.atomFD['range']=self.getAtomRange
        self.atomFD['regexp']=self.getAtomMatch
        self.atomFD['relative']=self.getAtomRelIndex
        self.atomFD['index']=self.getAtomIndex
        self.atomFD['NamedResSet'] = self.doNothing
        self.atomFD['NamedAtomSet'] = self.getNamedAtomSet


    def go(self):
        msgStr = None
        if self.moleculeSet: 
            self.molSet = self.getMolecules(self.moleculeSet, self.selList[0])
        else: 
            self.molSet = None
            return []
            #return [], msgStr

        #SPLIT here if mol.children==Atoms
        #possibly: self.Mols4levels=filter(lambda x: x.childrenSetClass == ChainSet, self.molSet)
        #then process the others separately....?????????
        #eg self.Mols2levels=filter(lambda x: x.childrenSetClass == AtomSet, self.molSet)
        #self.Mols2levels would get fed to self.getAtoms and two results ???merged???
        if not self.molSet:
            self.chainSet = None
            msgStr = str(self.selList[0])+ " selected no molecules"
            return []
            #return [], msgStr

        noChainMols=MoleculeSet(filter(lambda x: Chain not in x.levels, self.molSet))
        haveChainMols=MoleculeSet(filter(lambda x: Chain in x.levels, self.molSet))

        #build the residues belonging to molecules w/ no Chains (!WEIRD!)
        ncrs=ResidueSet()
        for item in noChainMols:
            if Residue in item.levels:
                 itemRes=item.allAtoms.parent.uniq()
                 ncrs= ncrs + itemRes

        if ncrs:
            noChainResSet=self.getResidues(ncrs,self.selList[2]) 
            if not noChainResSet: noChainResSet=ResidueSet()
        else:
            noChainResSet=ResidueSet()

        if len(haveChainMols):
            self.chainSet = self.getChains(haveChainMols.findType(Chain), self.selList[1])
        if self.chainSet:
            haveChainResSet=self.getResidues(self.chainSet.findType(Residue), self.selList[2])
            #also test noChainMols for residues
            if haveChainResSet:
                self.resSet=haveChainResSet+noChainResSet
            else:
                self.resSet = noChainResSet
        else: 
            self.resSet = noChainResSet
            ##don't return unless selList[2]!==['']
            if self.selList[1]!=['']:
                msgStr = str(self.selList[1])+ " selected no chains"
                return []
                #return [], msgStr

        #now: if self.selList for Chain and Residue level was empty, get the Atoms from noChains
        if self.selList[1]==[''] and self.selList[2]==['']:
            tla=AtomSet()
            for item in noChainMols:
                if Residue not in item.levels:
                    tla=tla + item.allAtoms
            twoLevelAtoms = tla
        else:
            twoLevelAtoms=AtomSet()
        if self.resSet:
            resAts=self.resSet.findType(Atom)
            if twoLevelAtoms:resAts=resAts+twoLevelAtoms
            self.atomSet = self.getAtoms(resAts, self.selList[3])
        else: 
            if self.selList[2]!=['']:
                msgStr = str(self.selList[2])+ " selected no residues"
                return []
                #return [], msgStr
            else: self.atomSet = self.getAtoms(twoLevelAtoms, self.selList[3])

        selNodes = self.atomSet
        #find correct levelType to return
        #need to split atomSet into two parts:
        if self.atomSet:
            haveChainAtoms = AtomSet(filter(lambda x: x.top!=x.parent,self.atomSet))
            haveNoChainAtoms=self.atomSet-haveChainAtoms
            if self.selList[3]==['']:
                #change atoms to residues
                if haveChainAtoms: selNodes = haveChainAtoms.parent.uniq()
                else: selNodes=ResidueSet()
                if self.selList[2]==['']:
                    #change residues to chains
                    if len(selNodes): 
                        selNodes = selNodes.parent.uniq()
                    if self.selList[1]==['']:
                        #change chains to molecules
                        if haveNoChainAtoms: noChainTops=haveNoChainAtoms.top.uniq()
                        else: noChainTops=ProteinSet()
                        if selNodes: selTops= selNodes.top.uniq()
                        else: selTops=ProteinSet()
                        selNodes = selTops+noChainTops
                        if self.selList[0]==['']:
                            #change molecules to molecules(?)in the case of no strs
                            if selNodes.__class__!=MoleculeSet:
                                selNodes = MoleculeSet(selNodes.top.uniq())
        else:
            msgStr = str(self.selList[3])+ " selected no atoms"
        for item in ['moleculeSet','molSet','chainSet','resSet','atomSet']:
             if hasattr(self, item): 
                 delattr(self, item)
#                 exec('del self.'+item)
        return selNodes
        #return selNodes, msgStr


    def getMolecules(self,nodes,selList):
        #for this function, nodes must be already be a MoleculeSet
        #molecule numbers are matched to self.vf.Mols indices by  getMolIndexMatch
        assert isinstance(nodes, MoleculeSet)
        if selList[0]!= '':
            selNodes= self.processList(nodes, selList,self.molFD)
        else: selNodes = nodes
        return selNodes


    def getChains(self,nodes,selList):
        #molecule numbers are matched to self.vf.Mols indices by  getMolIndexMatch
        assert isinstance(nodes, ChainSet)
        if selList[0]!= '':
            selNodes= self.processList(nodes, selList,self.chainFD)
        else: selNodes = nodes
        return selNodes


    def getResidues(self,nodes,selList):
        #for this function, nodes must be already be a ResidueSet
        #molecule numbers are matched to self.vf.Mols indices by  getMolIndexMatch
        assert isinstance(nodes, ResidueSet)
        if selList[0]!= '':
            selNodes= self.processList(nodes, selList,self.resFD)
            #print 'MVselector 3',sys.getrefcount(self.vf.Mols[0].chains[0].residues[8])
        else: selNodes = nodes
        return selNodes


    def getAtoms(self,nodes,selList):
        #for this function, nodes must be already be an AtomSet
        #molecule numbers are matched to self.vf.Mols indices by  getMolIndexMatch
        assert isinstance(nodes, AtomSet)
        if selList[0]!= '':
            selNodes= self.processList(nodes, selList,self.atomFD)
        else: selNodes = nodes
        return selNodes


    def processList(self, nodes, selList, FD):
        #first, check for empty string:
        if len(selList)==1 and selList[0]=='':
            return nodes        

        #otherwise:build a set of selected nodes by selecting w/ each item in list
        selNodes = None
        for item in selList:
            newNodes = self.processListItem(nodes, item, FD)
            if newNodes:
                if selNodes: selNodes = selNodes + newNodes
                else: selNodes = newNodes
        return selNodes
            

    def processListItem(self, nodes, item, FD):
        #classifies each item in selList and calls appropriate function

        if len(item)>1:
            while item[0]==' ':
                item=item[1:]
        #first detect 'lastitem' special char
        if item=='$':
            newNodes = nodes[-1]
            return nodes.ReturnType(newNodes)

        #first detect callable objects:
        try:
            func = evalString(item)
            newNodes = nodes.get(func)                
            return newNodes
        except:
            pass

        #detect ranges:
        if string.find(item, '-')!=-1 and string.find(item, '[')==-1:        
            #call range w/ nodes here:
            newNodes = FD['range'](item, nodes )
            return newNodes
        
        if item in residueList_.keys():
            newNodes = FD['NamedResSet'](item,nodes)
            return newNodes

        if item in atomList_.keys():
            newNodes = FD['NamedAtomSet'](item,nodes)
            return newNodes

        #next detect relative numbers
        if item[0]=='#':
            item = item[1:]
            newNodes = FD['relative'](item, nodes)
            return newNodes

        #next detect numbers
        try:
            item = int(item)
            newNodes = FD['index'](item, nodes)
            return newNodes        
        except:
            newNodes = FD['regexp'](item, nodes)
            return newNodes        


    #now all the functions:
    def rangeMatch(self, nodes,fr,to):
        ##given the two nodes, get a range
        indexfro = nodes.data.index(fr)
        indextoo = nodes.data.index(to)
        #assert indexfro<=indextoo
        if indexfro<=indextoo:
            return nodes.ReturnType(nodes.data[indexfro:indextoo+1])
        else: return None
        
        
    def doNothing(self, item, nodes):
        return None


    def stringMatch(self, nodes, regexp, field):
        return nodes.ReturnType(nodes.objectsFromString(regexp, field))
        #return nodes.ReturnType(nodes.objectsFromStringField(regexp, field))
                

    #FOR MOLECULES:
    def getMolRange(self, item, nodes):
        if len(nodes)<2:
            return None
        levItList=string.split(item, '-')
        firstNodes = self.processListItem(nodes, levItList[0], self.molFD)
        lastNodes = self.processListItem(nodes, levItList[1], self.molFD)
        if firstNodes and lastNodes:        
            return self.rangeMatch(nodes,firstNodes[0],lastNodes[-1])
        else:
            return None
        

    def getMolIndex(self,item,nodes):
        #have to make sure it's a valid number
        try:
            number = int(item)
        except:
            msgStr = item + " is invalid Molecule index"
            #self.vf.warningMsg(msgStr)   
            return None 
        #indices start w/ 1:
        return nodes[number-1:number]        


    def getMolMatch(self,item, nodes):
        #use self.procFunction to build a regexp
        #Mol strings match to 'name', so just use 'objectsFromString'
        reItem = self.procFunction(item)
        return nodes.ReturnType(nodes.objectsFromString(reItem))
        
        
    #FOR CHAINS:
    def getChainRange(self, item, nodes):
        if len(nodes)<2:
            return None
        levItList=string.split(item, '-')
        #Should chain range be on a per molecule basis????
        # this selects range of all chains 
        firstNodes = self.processListItem(nodes, levItList[0], self.chainFD)
        lastNodes = self.processListItem(nodes, levItList[1], self.chainFD)
        if firstNodes and lastNodes:
            return self.rangeMatch(nodes,firstNodes[0],lastNodes[-1])
        else:
            return None
        

    def getChainIndex(self, item, nodes):
        #this will get the item-th chain in ChainSet
        try:
            number = int(item)
        except:
            msgStr = item + " is invalid relative index"
            #self.vf.warningMsg(msgStr)   
            return None
        return nodes[number-1:number]        
        

    def getChainRelIndex(self,item,nodes):
        try:
            number = int(item)
        except:
            msgStr = str(item) + " is invalid relative index"
            #self.vf.warningMsg(msgStr)   
            return None
        #indices start w/ 1:
        parentNodes = nodes[0].parent.setClass(nodes.parent.uniq())
        l=[]
        for item in parentNodes:
           if len(item.children)>=number:
                l.append(item.children[number-1]) 
        return ChainSet(l)


    def getChainMatch(self,item, nodes):
        #use self.procFunction to build a regexp
        #Chain strings match to 'id', so use 'objectsFromStringField'
        reItem = self.procFunction(item)
        return nodes.ReturnType(nodes.objectsFromString(reItem,'id'))
        #return nodes.ReturnType(nodes.objectsFromStringField(reItem,'id'))
        

    #FOR RESIDUES:
    def getResidueRange(self, item, nodes):
        #this needs to be done on a PER CHAIN basis:
        if len(nodes) <2: return None
        levItList=string.split(item, '-')
        selNodes = None
        parentNodes = nodes[0].parent.setClass(nodes.parent.uniq())
        for par in parentNodes:
            nds = ResidueSet(filter(lambda x, par=par: x.parent==par, nodes))
            if len(nds)<2: continue
            firstNodes = self.processListItem(nds, levItList[0], self.resFD)
            lastNodes = self.processListItem(nds, levItList[1], self.resFD)
            if firstNodes and lastNodes: 
                newNodes = self.rangeMatch(nds,firstNodes[0],lastNodes[-1])
                if newNodes:
                    if selNodes: selNodes = selNodes + newNodes
                    else: selNodes = newNodes
            else: continue
        return selNodes


    def getResidueIndex(self, item, nodes):
        #residue indices are strings
        item = str(item)
        ans= filter(lambda x, item = item: x.number==item,nodes)
        return ResidueSet(ans)
        

    def getResidueRelIndex(self,item,nodes):
        try:
            number = int(item)
        except:
            msgStr = str(item) + " is invalid relative index"
            #self.vf.warningMsg(msgStr)   
            return None
        #indices start w/ 1:
        parentNodes = nodes[0].parent.setClass(nodes.parent.uniq())
        l=[]
        for par in parentNodes:
           if len(par.children)>=number:
                l.append(par.children[number-1]) 
        return ResidueSet(l)


    def getResidueMatch(self,item, nodes):
        #use self.procFunction to build a regexp
        #any item starting w/ a digit is matched to number field
        #Residue strings match to 'type', so use 'objectsFromStringField'
        reItem = self.procFunction(item)
        try:
            t = int(item[0])
            return nodes.ReturnType(nodes.objectsFromString(reItem, 'number'))
            #return nodes.ReturnType(nodes.objectsFromStringField(reItem, 'number'))
        except ValueError:
            return nodes.ReturnType(nodes.objectsFromString(reItem))
            #return nodes.ReturnType(nodes.objectsFromStringField(reItem, 'type'))


    def getNamedResSet(self, item,nodes):
        #here get all residues w/ name in  residueList_[item]
        rlist = residueList_[item]
        ans = filter(lambda x, rlist=rlist, nodes=nodes: x.type in rlist, nodes)
        return ResidueSet(ans)


    #FOR ATOMS:
    def getAtomRange(self, item, nodes):
        if len(nodes)<2:
            return None
        levItList=string.split(item, '-')
        if len(levItList)!=2: return None
        if levItList[0][0]=='#' or levItList[1][0]=='#':
            return self.getAtomRelRange(item,nodes)
        firstNodes = self.processListItem(nodes, levItList[0], self.chainFD)
        lastNodes = self.processListItem(nodes, levItList[1], self.chainFD)
        if firstNodes and lastNodes:
            return self.rangeMatch(nodes,firstNodes[0],lastNodes[-1])
        else:
            return None

    def getAtomRelRange(self,item,nodes):

        levItList=string.split(item, '-')
        #now the hard part: need to call pLI w/ each set of parent nodes
        selNodes = None
        parentNodes = ResidueSet(nodes.parent.uniq())
        for par in parentNodes:
            nds = AtomSet(filter(lambda x, par=par: x.parent==par, nodes))
            firstNodes = self.processListItem(nds, levItList[0], self.atomFD)
            lastNodes = self.processListItem(nds, levItList[-1], self.atomFD)
            if firstNodes and lastNodes:
                newNodes= self.rangeMatch(nds,firstNodes[0],lastNodes[-1])
            if newNodes:
                if selNodes: selNodes=selNodes + newNodes
                else: selNodes = newNodes
        return selNodes


    def getNamedAtomSet(self,item,nodes):
        #here get all atoms w/ name in  atomList_[item]
        #AND 8/2004: whose parents are std residues...
        alist = atomList_[item]
        #only get atoms in standard residues
        reslist = residueList_['std']
        res_atoms = filter(lambda x, nodes=nodes: x.parent.type in reslist, nodes)
        ans = filter(lambda x, alist=alist, res_atoms=res_atoms: x.name in alist, res_atoms)
        #previously:
        #ans = filter(lambda x, alist=alist, nodes=nodes: x.name in alist, nodes)
        return AtomSet(ans)


    def getAtomIndex(self, item, nodes):
        ans= filter(lambda x, item = item: x.number ==item,nodes)
        return AtomSet(ans)
            

    def getAtomRelIndex(self,item,nodes):
        try:
            number = int(item)
        except:
            msgStr = item + " is invalid relative index"
            #self.vf.warningMsg(msgStr)   
            return None
        #indices start w/ 1:
        parentNodes = nodes[0].parent.setClass(nodes.parent.uniq())
        l=[]
        for item in parentNodes:
           if len(item.children)>=number:
                l.append(item.children[number-1]) 
        return AtomSet(l)


    def getAtomMatch(self,item, nodes):
        #use self.procFunction to build a regexp
        reItem = self.procFunction(item)
        #FOR THE MOMENT: match to name field (should it be element?)
        return nodes.ReturnType(nodes.objectsFromString(reItem))


    #Functions for converting input for regular expressions:
    def processStringcIWEC(self,someString):
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
                    continue 
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
                        if c =='^'or c =='[': newExp = newExp + c
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


    def processStringcI(self,someString):
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
        someString = string.replace(someString, '-', ':')
        someString = string.replace(someString, '*', '.*')
        return someString


class Mv102StringSelector(MVStringSelector):


    def __init__(self, moleculeSet, selString, userPref = 'cS'):
        self.form = None
        if len(selString)<4:
            self.selString =self.getMvStrings(selString)
        else:
            self.selString = string.split(selString, ':')
        MVStringSelector.__init__(self, moleculeSet, self.selString,'cIWEC')
        

    def getMvStrings(self,selString):
        #nb  selString looks like ('::ALA:')
        selList = string.split(selString[0], ':')
        z = ''
        if len(selList) == 1:
           selList.append(z)
        if len(selList) == 2:
           selList.append(z)
        if len(selList)<4:
            pt1 = selList[:1]
            pt1.append(z)
            selList = pt1+ selList[1:]
        self.selList = selList
        return selList


    #THIS is one of two method overwritten in Mv102StringSelector
    def getResidueMatch(self,item, nodes):
        #use self.procFunction to build a regexp
        #any item starting w/ a digit is matched to number field
        #Residue strings match to 'type', so use 'objectsFromStringField'
        reItem = self.procFunction(item)
        try:
            t = int(item[0])
            return nodes.ReturnType(nodes.objectsFromString(reItem, 'number'))
            #return nodes.ReturnType(nodes.objectsFromStringField(reItem, 'number'))
        except ValueError:
            return nodes.ReturnType(nodes.objectsFromString(reItem, 'type'))
            #return nodes.ReturnType(nodes.objectsFromStringField(reItem, 'type'))


    # mv102 specific processListItem:
    # differences: callable items and residue set Names  not included
    #also, userpref has to be cIWEC
    def processListItem(self, nodes, item, FD):
        #classifies each item in selList and calls appropriate function

        #first detect 'lastitem' special char
        if item=='$':
            newNodes = nodes[-1]
            return nodes.ReturnType(newNodes)

        #detect ranges:
        if string.find(item, '-')!=-1 and string.find(item, '[')==-1:        
            #call range w/ nodes here:
            newNodes = FD['range'](item, nodes )
            return newNodes

        #next detect relative numbers
        if item[0]=='#':
            item = item[1:]
            newNodes = FD['relative'](item, nodes)
            return newNodes

        #next detect numbers
        try:
            item = int(item)
            newNodes = FD['index'](item, nodes)
            return newNodes        
        except:
            newNodes = FD['regexp'](item, nodes)
            return newNodes        

