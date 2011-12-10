#
#Authors: Michel F Sanner,Sowjanya Karnati
#

#NOTES:This class finds and returns  dependencies of all packages or one
#patricular package.When run with -V option finds dependencies at file level
#
#
import os, sys, re,string
from os.path import walk
import getopt,warnings
os.system('rm -rf depresult')
vinfo = sys.version_info
pythonv = "python%d.%d"%(vinfo[0], vinfo[1])
######################################################################
#   COMMAND LINE OPTIONS
######################################################################

class DEPENDENCYCHECKER:

########################################################################
#   Findind packages
########################################################################


 def rundependencychecker(self,pack,v=True,V=False,loc=None):    
  cwd = os.getcwd()
  import string
  packages = {}
  pn =[]
  for i in sys.path:
        s =i.split('/')    
        if s[-1]=="MGLToolsPckgs" or s[-1]== 'site-packages':
        #if s[-1]=='FDepPackages' or s[-1]=='RDepPackages' or s[-1]=='FSharedPackages' or s[-1]=='RSharedPackages' or s[-1]== 'site-packages':
            pn.append(i)
  
  for p in pn:
    os.chdir(p)
    files = os.listdir(p)
    for f in files:
        if not os.path.isdir(f):
            continue
        pwdir = os.path.join(p, f)
        if os.path.exists( os.path.join( pwdir, '__init__.py')):
            if not packages.has_key(f):
                packages[f] = pwdir
        elif os.path.exists( os.path.join( p, '.pth') ):
            if not packages.has_key(f):
                packages[f] = pwdir
  os.chdir(cwd)

################################################################
#     packages in site-packages
###################################################################


  pack_site=packages.keys()


##########################################################################
#   Finding list of import statements used in all packages
###########################################################################

  Total_packs = []
  if V == True:
   if pack!=None:
    packn = string.split(pack,'.')
    pack =packn[0]
    exec('packages={"%s":"%s"}'%(pack,packages[pack]))
    if packn[-1]!='py':
        pack_file =packn[-1] 
    if packn[-1]=='py':
        pack_file =packn[-2]
   else:
    pack_file=None
   for pack in packages:
    files = []
    pat = re.compile('import')
    print "please wait ...."
    for root, dirs, files in os.walk(packages[pack]):
        # remove directories not to visit
        for rem in ['CVS', 'regression', 'Tutorial', 'test','Doc','doc']:
            if rem in dirs:
                dirs.remove(rem)
        # look for files that contain the string 'import'
                
        for fi in files:
            
            if fi[-3:]!='.py':
                continue
            if fi[-3:]=='.py':
                #finds pattern "import" match in that particular file
                if pack_file!=pack:
                   if pack_file!=None: 
                    if fi !=pack_file+'.py': 
                        continue
                    else:    
                        candidates = []
                        f = open( os.path.join(root, fi) )
                        data = f.readlines()
                        f.close()
                        found = 0
                        for line in data:
                            match = pat.search(line)
                            if match:
                                candidates.append( (root, fi, line) )
                               
                #finds pattern "import" for packages given with option p at file level
                if pack_file==pack:
                         
                        candidates = []
                        f = open( os.path.join(root, fi) )
                        data = f.readlines()
                        f.close()
                        found = 0
                        for line in data:
                            match = pat.search(line)
                            if match:
                                candidates.append( (root, fi, line) )    
                #finds pattern "import" match for all packages at file level
                else:
                         
                        candidates = []
                        f = open( os.path.join(root, fi) )
                        data = f.readlines()
                        f.close()
                        found = 0
                        for line in data:
                            match = pat.search(line)
                            if match:
                                candidates.append( (root, fi, line) )    
                #######################################
                #finding dependencies
                #######################################
                 
                result= []
                import string
                if len(candidates)>0:
                    
                    for candidate_num in candidates:
                        p, f, imp = candidate_num
                        path =string.split(p,'site-packages')[-1]
                        implist =[]
                        fromlist=[] 
                        y =string.split(imp)
                        #omitting commemted imports
                        if '.' not in imp and y[0]=="import":
                            
                            len_space = len(imp.split(' '))
         
                            len_comma=len(imp.split(','))
                            if (len_space -1) > len_comma:
                                continue
                        # as im import statement
                        if "as" in y:
                            
                            for a in y:
                                   if a=='as':
                                    aind = y.index(a)
                                     
                                    if '.' not in y[aind-1] :
                                        implist.append(y[aind-1])
                                        continue 
                                    else:
                                        newa = y[aind-1].split('.')
                                        implist.append(newa[0])
                                        continue
                        if '#' in y:
                            continue
                        
                        #if first word is import in the list
                        if y[0]=='import':
                            
                            for i in range(1,len(y)):
                                if y[i][-1]==";":
                                    y[i]=y[i][:-1]
                                    if y[i] not in implist:
                                        implist.append(y[i])
                                        break
                                    
                                if 'as' in y:
                                    break 
                                if y[i][-1]==',':
                                    y[i]=y[i][:-1]
                                    if  ',' in y[i]:
                                        srg = string.split(y[i],',')    
                                        for j in srg:
                                            if j not in implist:
                                                implist.append(j)
                                                continue
                                
                                    elif len(y[i])<=2:
                                        continue 
                                    #if import statement is like a.q
                                    elif len(string.split(y[i],'.'))!=1:
                                        sr = string.split(y[i],'.')
                                        if sr[0] not in implist:
                                            #if module doesn't starts with __
                                            #append to list
                                            if sr[0][0]!='__':
                                                implist.append(sr[0])
                                    #if import statement with out '.'
                                    else:
                                        if y[i] not in implist:
                                            #print y[i]
                                            if y[i][0]!='__':
                                                implist.append(y[i]) 
                                #import statement with out ',' in the end
                                else:
                                    if len(y[i])==1:
                                        continue
                                    elif  ',' in y[i]:
                                        srg = string.split(y[i],',')    
                                        for j in srg:
                                            if j not in implist:
                                                implist.append(j)
                                                continue
                                    #import statement as a.b.c.d
                                    elif len(string.split(y[i],'.'))>1:
                                        sr = string.split(y[i],'.')
                                        if sr[0] not in implist:
                                            if sr[0][0]!='__': 
                                                implist.append(sr[0])
                                                continue
                                    #import statement without '.'
                                    elif y[i] not in implist:
                                        if y[i][0]!='__':
                                            implist.append(y[i])
                                            continue                                            
                            for im in implist:                                
                                #try importing module in implist
                                try:        
                                    exec('import %s'%im)
                                    if im == 'Pmw':
                                        if im not in result:
                                            if im!=pack:
                                                result.append(im)
                                                continue
                                        else:
                                            continue
                                    #if module.__file__ exists check in 
                                    #site-packages and append to result
                                    exec('fi = %s.__file__'%im)
                                    fil = os.path.abspath('%s'%fi)
                                    if os.path.exists(fil):
                                        file = string.split(str(fil),'/')
                                        if file[-2] in pack_site:
                                            if file[-2] not in result:     
                                                if file[-2] !=pack:
                                                    result.append(file[-2])    
                                        elif file[-2]=='Numeric':
                                            if 'Numeric' not in result:     
                                                if 'Numeric'!=pack:
                                                    result.append('Numeric')
                                        elif file[-2] not in ['lib-dynload', pythonv,'lib-tk']:                 
                                            if im not in result:     
                                                if im!=pack:
                                                    result.append(im)                
                                except:
                                    if im in ['sys','gc','thread','exceptions']:
                                        continue
                                    else:
                                        if im not in result:
                                            result.append(im)
                        #if first word is from in list
                        if y[0]=='from':
                            #if from statement is like a.b.c
                            if len(string.split(y[1],'.'))!=1:
                                sr = string.split(y[1],'.')
                                if sr[0] not in fromlist:
                                    fromlist.append(sr[0])    
                            else:
                                if y[1]!=pack:
                                        if y[1] not in fromlist:
                                            if y[1][0]!='__':
                                                fromlist.append(y[1]) 
                            for i in fromlist:
                                #checks importing module
                                try:        
                                    exec('import %s'%i)
                                    if i == 'Pmw':
                                        if i not in result:
                                            if i !=pack:
                                                result.append(i)
                                                continue
                                        else:
                                            continue
                                    #if __file exixts check in site-pacakges
                                    #and append to result
                                    exec('fi = %s.__file__'%i)
                                    fil = os.path.abspath('%s'%fi)
                                    if os.path.exists(fil):
                                        file = string.split(str(fil),'/')
                                        if file[-2] in pack_site:
                                                if file[-2] not in result:     
                                                    if file[-2] !=pack:
                                                        result.append(file[-2])
                                        elif file[-2]=='Numeric':
                                            if 'Numeric' not in result:     
                                                if 'Numeric'!=pack:
                                                    result.append('Numeric')
                                        elif file[-2] not in ['lib-dynload', pythonv,'lib-tk']:
                                            if i not in result:
                                                if i!=pack :
                                                    result.append(i)
                                except:
                                    if i in ['sys','gc','thread','exceptions']:
                                        continue
                                    else:
                                        if i not in result:
                                            result.append(i)
                    listdf=[]
                    for r,d,lf in os.walk(packages[pack]):
                            for rem in ['CVS', 'regression', 'Tutorial', 'test','Doc','doc']:
                                if rem in d:
                                    d.remove(rem)
                            #when files in pack are imported
                            listd = os.listdir(r)
                            for ld in listd:
                                if ld.endswith('.py')==True or ld.endswith('.so')==True:
                                    for res in result:
                                        if res == ld[:-3]:
                                            if res in result:
                                                result.remove(res)
                            for files in lf:
                                for res in result: 
                                    pat1 = re.compile('"%s"' %res)
                                    pat2 = re.compile("%s.py" %res)
                                    fptr=open("%s/%s" %(r,files))
                                    lines = fptr.readlines()
                                    for line in lines:
                                        match1 = pat1.search(line)
                                        match2 = pat2.search(line)
                                        if match1 or match2:
                                            if res in result:
                                                ind = result.index(res)
                                                if result[ind] not in pack_site:
                                                    del result[ind] 
                                                    continue
                                    #In case of Pmv multires,pdb2qr etc
                                    if res in files:
                                      if res in result:  
                                        ind = result.index(res)
                                        if result[ind] not in pack_site:
                                            del result[ind]
                            for res in result:
                                if res[:3] in ['tmp','TMP','win']:
                                    result.remove(res)
                                    continue
                    exec('l = len("%s")'%f)
                    if l>60:
                    
                        exec('md = string.split("%s",",")[0]'%f)
                        exec('f = string.split("%s","/")[-1][:-1]'%md)
                        Total_packs.append('result_%s %s %s %s' %(pack,path,f,result))
                        #return Total_packs
                    else:
                        Total_packs.append('result_%s %s %s %s' %(pack,path,f,result))    
                    print "result_%s %s %s %s" %(pack,path,f,result)
    if Total_packs:
        return Total_packs


  else:
   if pack != None:
     #print pack
     pack_list = pack.split(',')
     if len(pack_list)>=1:
      packs = {}
      for p in pack_list:
         packs[p]=packages[p]
      packages = packs
   print "please wait ......."
   for pack in packages:
    files = []
    pat = re.compile('import')
    candidates = [] 
    for root, dirs, files in os.walk(packages[pack]):
        # remove directories not to visit        
        for rem in ['CVS', 'regression', 'Tutorial', 'test','Doc','doc']:
            if rem in dirs:
                dirs.remove(rem)
        # look for files that contain the string 'import'
            for fi in files:
                if fi[-3:]!='.py':
                    continue
                #finding pattern "import" match 
                if fi[-3:]=='.py':
                    
                    f = open( os.path.join(root, fi) )
                    data = f.readlines()
                    f.close()
                    found = 0
                    for line in data:
                        match = pat.search(line)
                        if match:
                            candidates.append( (root, fi, line) )

##############################
#finding dependencies
##############################
    result= []
    import string
    for candidate_num in candidates:
        #print candidate_num
        p, f, imp = candidate_num
        implist =[]
        fromlist=[] 
        y =string.split(imp)
        #omitting commemted imports 
        if '.' not in imp and y[0]=="import":
            len_space = len(imp.split(' '))
            len_comma=len(imp.split(','))
            if (len_space -1) > len_comma:
                continue
        if "as" in y:
            for a in y:
                if a=='as':
                    aind = y.index(a)
                    if '.' not in y[aind-1] :
                        if y[aind-1] not in implist:
                            implist.append(y[aind-1])
                            continue
                    else:
                        newa = y[aind-1].split('.')
                        if newa[0] not in implist:
                            implist.append(newa[0])
                            continue
        if '#' in y:
            continue
        #if first word is import in the list
        if y[0]=='import':
            for i in range(1,len(y)):
                if y[i][-1]==";":
                    y[i]=y[i][:-1]
                    if y[i] not in implist:
                       implist.append(y[i])
                       break
                if "as" in y:
					break 
                
                if y[i][-1]==',':
                    y[i]=y[i][:-1]
                    if  ',' in y[i]:
                        srg = string.split(y[i],',')    
                        for j in srg:
                            if j not in implist:
                                implist.append(j)
                                continue
                    elif len(y[i])<=2:
                       continue 
                    elif len(string.split(y[i],'.'))!=1:
                        sr = string.split(y[i],'.')
                        if sr[0]!=pack:
                            if sr[0] not in implist:
                                if sr[0][0]!='__':
                                    implist.append(sr[0])
                    else:
                        if y[i]!=pack:
                            if y[i] not in implist:
                                if y[i][0]!='__':
                                    implist.append(y[i]) 
                else:
                    if len(y[i])==1:
                        continue
                    if  len(string.split(y[i],','))!=1:
                        srg = string.split(y[i],',')    
                        for j in range(0,len(srg)):
                            if srg[j] not in implist:
                                implist.append(srg[j])
                            else:
                                continue
                    elif len(string.split(y[i],'.'))!=1:
                        sr = string.split(y[i],'.')
                        if sr[0] not in implist:
                            if sr[0][0]!='__': 
                                implist.append(sr[0])
                                
                    elif y[i] not in implist:
                            if y[i][0]!='__':
                                implist.append(y[i])
                                                   
            
            for im in implist:
                try:        
                    exec('import %s'%im)
                    if im == 'Pmw':
                        if im not in result:
                            if im!=pack:
                                result.append(im)
                                continue
                        else:
                            continue
                    
                    exec('fi = %s.__file__'%im)
                    fil = os.path.abspath('%s'%fi)
                    if os.path.exists(fil):
                        file = string.split(str(fil),'/')
                        if file[-2] in pack_site:
                            if file[-2] not in result:     
                               if file[-2] !=pack:
                                    result.append(file[-2])    
                        elif file[-2]=='Numeric':
                            if 'Numeric' not in result:     
                               if 'Numeric'!=pack:
                                    result.append('Numeric')    
                        elif file[-2] not in ['lib-dynload', pythonv,'lib-tk']:
                            if im not in result:     
                               if im!=pack:
                                    result.append(im)                
                except:
                    if im in ['sys','gc','thread','exceptions']:
                        continue
                    if im not in result:
                        if im!=pack:
                            result.append(im)
        #if first word is from in list
        if y[0]=='from':
           if len(string.split(y[1],'.'))!=1:
                sr = string.split(y[1],'.')
                if sr[0] != pack:
                    if sr[0] not in fromlist:
                        fromlist.append(sr[0])  
           else:             
                if y[1]!=pack:
                    if y[1] not in fromlist:
                        fromlist.append(y[1]) 
           for i in fromlist:
                try:        
                    exec('import %s'%i)
                    if i == 'Pmw':
                        if i not in result:
                            if i !=pack:
                                result.append(i)
                                continue
                        else:
                            continue
                    exec('fi = %s.__file__'%i)
                    fil = os.path.abspath('%s'%fi)
                    if os.path.exists(fil):
                        file = string.split(str(fil),'/')
                        if file[-2] in pack_site:
                            if file[-2] not in result:     
                               if file[-2] !=pack:
                                    result.append(file[-2])
                        elif file[-2]=='Numeric':
                            if 'Numeric' not in result:     
                               if 'Numeric'!=pack:
                                    result.append('Numeric')
                        elif file[-2] not in ['lib-dynload', pythonv,'lib-tk']:
                            if i not in result:
                               if i!=pack :
                                    result.append(i)                
                except:
                    if i in ['sys','gc','thread','exceptions']:
                        continue
                    if i not in result:
                        if i!=pack:
                            result.append(i)    
                    
    listdf =[]
    for r,d,lf in os.walk(packages[pack]): 
        for rem in ['CVS', 'regression', 'Tutorial', 'test','Doc','doc']:
            if rem in d:
                d.remove(rem)
        #when directory is imported eg: import extent(opengltk/extent) 
        for res in result:
            if res in d:
                if res in result:
                    result.remove(res)
                    continue
        
        #when files in pack are imported
        listd = os.listdir(r)
        for ld in listd:
            if ld.endswith('.py')==True or ld.endswith('.so')==True:
                 for res in result:
                    if res == ld[:-3]:
                        if res in result:
                            result.remove(res)
                                    
        #This is for cases in which (mainly tests) some file is created and
        #for tests purpose and removed .but imported in testfile
        #like tmpSourceDial.py in mglutil
        for files in lf:
          if files[-3:]!=".py":
            lf.remove(files)   
            continue
          for res in result: 
            pat1 = re.compile('%s' %res)
            pat2 = re.compile('%s.py' %res)
            fptr=open("%s/%s" %(r,files))
            lines = fptr.readlines()
            for line in lines:
                match1 = pat1.search(line)
                match2 = pat2.search(line)
                if match1 or match2:
                        if res in result:
                            ind = result.index(res)
                            if result[ind] not in pack_site:
                                del result[ind] 
                                continue
            
            #In case of Pmv multires,pdb2qr etc
            if res in files:
               if res in result:
                    ind = result.index(res)
                    if result[ind] not in pack_site:
                        del result[ind]
                        continue
        for res in result:
            if res[:3] in ['tmp','TMP','win']:
                result.remove(res)
                continue 
    print "result_%s %s" %(pack,result)
    Total_packs.append('result_%s %s' %(pack,result))
   return Total_packs
############################################################
#   PRINTING DEPENDENCIES AND COPYING THEM TO FILE 
############################################################
    
  
