#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autodockHosts.py,v 1.1.1.1 2001/04/03 19:47:50 gillet Exp $
#
# $Id: autodockHosts.py,v 1.1.1.1 2001/04/03 19:47:50 gillet Exp $
#
#
#

import UserDict

class AutoDockHosts(UserDict.UserDict):

	def __init__(self, localHostDict):
		UserDict.UserDict.__init__(self)
		self.update(localHostDict)

	def buildEntry(self, host=None,agPath=None, adPath=None, qType='int', userSpecific=0):
		d={}
		d['host']=host
		d['autogrid']=agPath
		d['autodock']=adPath
		d['queuetype']=qType
		d['userSpecific']=userSpecific
		return d

	def addHost(self,macroName, hostdict=None, **kw):
		host = kw['host']
		agPath = kw['autogrid']
		adPath = kw['autodock']
		qType = kw['queuetype']
		userSpecific = kw['userSpecific']
		if not hostdict:
			hostdict=self.buildEntry(host=host, agPath=agPath, adPath=adPath,
				qType=qType,userSpecific=userSpecific)
		self[macroName]=hostdict
		#nb: preexisting macroName entry is overwritten

	def saveHostFile(self, filename, whichOnes='all'):
		#will be in file called .adthosts.py
		#and consist of python code for dictionary called 'newhosts'
		fptr = open(filename, 'w')
		outstr = 'hostMacros={'
		#outstr = 'hosts={'
		fptr.write(outstr)
		#always write a localHost line
		#get the correct macroList here
		if whichOnes=='all':
			macroList=self.keys()
		elif whichOnes=='userSpecific':
			#get the one with userSpecific=1, only
			macroList=[]
			for item in self.items():
				if item[1]['userSpecific']: 
					macroList.append(item)
		else:
			macroList=[]
			for item in self.items():
				if not item[1]['userSpecific']: 
					macroList.append(item)
			#get the other ones...
		for i in range(len(macroList)):
			h=macroList[i][0]
			self.writeEntry(h,fptr)
			if i<len(macroList)-1:
				fptr.write('\t\t},\n')
			else:
				fptr.write('\t\t}\n')
		#close the whole thing
		outstr = '\t}\n'
		fptr.write(outstr)
		fptr.close()

	def writeEntry(self, macroName,fptr):
		outstr='\t\''+macroName+'\': {\n'
		#outstr='\t\''+hostName+'\': {\n'
		fptr.write(outstr)
		d = self[macroName]
		#d = self[hostName]
		klist = d.keys()
		for i in range(len(klist)):
			k=klist[i]
			if k=='userSpecific':
				outstr=	'\t\t\''+k + '\': '+ str(d[k])
			else:
				outstr=	'\t\t\''+k + '\': \''+ str(d[k])+'\''
			fptr.write(outstr)
			if i<len(klist)-1:
				outstr= ',\n'
			else:
				outstr= '\n'
			fptr.write(outstr)
		

	def loadHostFile(self, filename):
		newStuff=__import__(filename)
		self.update(filename.adhosts)
		#self.update(adhosts)

	
