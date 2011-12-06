#
# $Id: GromacsEM.py 290 2011-06-03 09:51:17Z nielsen $
#
# Python class for running Gromacs
#
import os

class GROMACS:

    def __init__(self,neededfiles,createdir=True):
        #
        # Find Gromacs
        #
        GMXBIN=None
        import os
        path=os.environ['PATH']
        for dirname in path.split(':'):
            mdrun=os.path.join(dirname,'mdrun')
            if os.path.isfile(mdrun):
                GMXBIN=dirname
        if not GMXBIN:
            GMXBIN='/usr/bin'
        environment={'GMXBIN':GMXBIN} #,'GMXLIB':'/usr/share/gromacs/top/'}
        #
        # Set the environment
        #
        for key in environment.keys():
            os.putenv(key,environment[key])
        if createdir:
            #
            # Create the tempdir
            #
            import tempfile
            self.topdir=os.getcwd()
            tempfile.tempdir=self.topdir
            self.tmpdir=tempfile.mktemp()
            os.mkdir(self.tmpdir)
            print 'Running in :',self.tmpdir
            os.chdir(self.tmpdir)
            #
            # Copy all the input files to the tmpdir
            #
            for file in neededfiles:
                os.system('cp '+file+' ./')
        #
        # Filenames
        #
        self.files={}
        return

    #
    # ------
    #

    def cleanup(self):
        """Remove the tempdir"""
        import shutil
        os.chdir(self.topdir)
        shutil.rmtree(self.tmpdir)
        self.tmpdir=None
        return

    #
    # ---
    #

    def pdb2gmx(self,pdbfile,forcefield=0,ignore_Hs=True,auto_select_his=True):
        #
        # forcefield: 0 for Gromacs FF, 1 for Gromacs FF with all hyds.
        #
        print auto_select_his
        self.files['laststructure']='protein.gro'
        fd=open('pdb2gmx.scr','w')
        if auto_select_his:
            if ignore_Hs:
                fd.write('pdb2gmx -f '+pdbfile+' -ignh -o '+self.files['laststructure']+' << EOF\n')
            else:
                fd.write('pdb2gmx -f '+pdbfile+' -o '+self.files['laststructure']+' << EOF\n')
        else:
            if ignore_Hs:
                fd.write('pdb2gmx -f '+pdbfile+' -ignh -his -o '+self.files['laststructure']+' << EOF\n')
            else:
                fd.write('pdb2gmx -f '+pdbfile+' -his -o '+self.files['laststructure']+' << EOF\n')
        #
        # Select the force field
        #
        fd.write('\n'+str(forcefield)+'\n\n')
        #
        # Select histidine protonation states
        if not auto_select_his:
            #
            # Histidine
            # 0: on ND1 only (4lyt)
            # 1: on NE2 only  (7lyz, 2lzm)
            # 2: Double protonated (2lzt & 4lzt)
            fd.write('1\n')
        #
        # For protonating E35 and D52 in HEWL
        #
        #fd.write('0\n0\n1\n0\n0\n0\n0\n')
        
        #fd.write('0\n1\n')
        #
        # Water model
        #
        fd.write('1\n\n')
        fd.write('\n\n')
        fd.close()
        status=os.system('/bin/tcsh -i -c "source ./pdb2gmx.scr"')
        if status!=0:
            raise Exception('arrg')
        return

    def center(self,dist):
        #
        # Center molecule in box and add dist space between box and protein
        #
        if not self.files.has_key('laststructure'):
            raise "No valid input file for editconf"
        status=os.system('editconf -f '+self.files['laststructure']+' -o out.gro -bt cubic -d '+str(dist))
        self.files['laststructure']='out.gro'
        print 'Editconf returned',status
        return

    def solvate(self):
        #
        # Solvate the system
        #
        print 'Solvating system'
        status=os.system('genbox -cp '+self.files['laststructure']+' -cs -p topol.top -o solvated.gro')
        self.files['laststructure']='solvated.gro'
        print 'genbox returned:',status
        return

    def EM(self,emtol=2000,nsteps=1000,emstep=0.01,nstenergy=1):
        #
        # Define the standard em.mdp file
        #
        params={'cpp':'/lib/cpp','define':'-DFLEX_SPC','constraints':'none','integrator':'steep','nsteps':str(nsteps),'emtol':str(emtol),'emstep':str(emstep),'nstcomm':'1','ns_type':'grid','rlist':'2','rlist':'1','rcoulomb':'1.0','epsilon_r':'1000.0','rvdw':'1.0','Tcoupl':'no','Pcoupl':'no','gen_vel':'no','nstenergy':str(nstenergy)}
        fd=open('em.mdp','w')
        for key in params.keys():
            fd.write('%s\t\t=  %s\n' %(key,params[key]))
        fd.close()
        #
        # Preprocessing
        #
        status=os.system('grompp_d -v -f em.mdp -c '+self.files['laststructure']+' -o em.tpr -p topol.top')
        print 'grompp returned:',status
        #
        # Do minimisation
        #
        self.files['energytraj']='energy.ene'
        self.files['laststructure']='after_em.pdb'
        status=os.system('mdrun_d -v -s em.tpr -o em.trr -c '+self.files['laststructure']+' -g emlog -e '+self.files['energytraj'])
        print 'mdrun returned:',status
        return

    def PR_MD(self,nsteps=5000,genseed='173529',uparams={}):
        #
        # Standard pr.mdp file
        #
        params={'title':'Position_restrained_MD',
                'cpp':'/lib/cpp','define':'-DPOSRES',
                'constraints':'all-bonds',
                'integrator':'md','tinit':'0.0',
                'dt':'0.002','nsteps':'5000',
                'nstcomm':'1','nstxout':'50','nstvout':'1000','nstfout':'0','nstlog':'10','nstenergy':'10','nstlist':'10','ns_type':'grid',
                'rlist':1.0,
                'coulombtype':'PME','fourierspacing':'0.1','pme_order':'4','rcoulomb':'1.0',  
                'rvdw':1.0,
                'Tcoupl':'v-rescale','tc-grps':'System','tau_t':'0.1','ref_t':'300',
                'energygrps':'Protein SOL',
                'Pcoupl':'Berendsen','tau_p':'0.5','compressibility':'4.5e-5','ref_p':'1.0',
                'gen_vel':'yes','gen_temp':'300.0','gen_seed':genseed}

        params['nsteps']=str(nsteps)
        #
        # Did the user pass anything?
        #
        for p in uparams.keys():
            params[p]=uparams[p]
            print p,params[p]
        #
        # Write the file
        fd=open('pr.mdp','w')
        for key in params.keys():
            fd.write('%s\t\t=  %s\n' %(key,params[key]))
        fd.close()
        #
        # Preprocessing
        #
        status=os.system('grompp -v -f pr.mdp -c '+self.files['laststructure']+' -o pr.tpr -p topol.top')
        print 'grompp returned:',status

        #
        # Do minimisation
        #
        self.files['energytraj']='pr.ene'
        self.files['laststructure']='after_pr.pdb'
        status=os.system('mdrun -v -s pr.tpr -o pr.trr -c '+self.files['laststructure']+' -g prlog -e '+self.files['energytraj'])
        print 'mdrun returned:',status
        return

    #
    # ------
    #

    def MD(self,steps=100000,uparams={}):
        #
        # Standard full.mdp file from http://rugmd4.chem.rug.nl/~gmx/online2.0/getting_started.html#full
        #
        # Particle-Particle Particle-Mesh for long range electrostatics
        #
        params={'title':'Full_MD','cpp':'/lib/cpp','constraints':'all-bonds',
                'integrator':'md','tinit':'0.0','dt':'0.002','nsteps':str(steps),
                'nstcomm':'1','nstxout':'250','nstvout':'1000','nstfout':'0','nstlog':'100','nstenergy':'100','nstlist':'10','ns_type':'grid',
                'Tcoupl':'v-rescale','tc-grps':'System','tau_t':'0.1','ref_t':'300','energygrps':'Protein SOL',
                'Pcoupl':'Berendsen','tau_p':'0.5','compressibility':'4.5e-5','ref_p':'1.0','gen_vel':'no','gen_temp':'300.0','gen_seed':'173529',  

                # Neighbour searching
                'pbc':'xyz',  
                # VDW
                'rlist':'1.0','rvdw':'1.0',  
                # Electrostatics
                'coulombtype':'PME','fourierspacing':'0.1','pme_order':'4','rcoulomb':'1.0'  
                }
        #
        # Did the user pass anything?
        #
        for p in uparams.keys():
            params[p]=uparams[p]
        #
        # Write the control file
        #
        fd=open('full.mdp','w')
        for key in params.keys():
            fd.write('%s\t\t=  %s\n' %(key,params[key]))
        fd.close()
        #
        # Preprocessing
        #
        self.files['mdpfile']='full.mdp'
        self.files['runinput']='full.tpr'
        command='grompp -v -f '+self.files['mdpfile']+' -c '+self.files['laststructure']+' -o '+self.files['runinput']+' -p topol.top'
        print 'Using command'
        print command
        status=os.system(command)
        print 'grompp returned:',status
        if status!=0:
            raise Exception('Exiting with error after grompp')
        #
        # Do Molecular Dynamics
        #
        backup=self.files.copy()
        self.files['energytraj']='MD.ene'
        self.files['laststructure']='after_MD.pdb'
        self.files['trajectory']='MD.trr'
        self.MDtime=steps/500.0
        command='mdrun -v -s '+self.files['runinput']+' -o '+self.files['trajectory']+' -c '+self.files['laststructure']+' -g mdlog -e '+self.files['energytraj']
        status=os.system(command)
        print 'mdrun returned:',status
        if status!=0:
            self.files=backup.copy()
        return
        
    #
    # ---
    # 
    
    def get_snapshots(self,numsnapshots):
        """Get the Gromacs snapshots"""
        import os
        fd=open('trjconv.scr','w')
        self.files['snapshots']='MD.snapshot.pdb'
        command='trjconv -s %s -f %s -o %s.pdb -fit rot+trans -sep -dt %3d << EOF\n' %(self.files['runinput'],self.files['trajectory'],self.files['snapshots'],max(1,self.MDtime/numsnapshots))
        fd.write(command)
        fd.write('1\n') # Selecting protein for superpositioning
        fd.write('1\n') # Selecting protein for output
        fd.close()
        #
        # Run the script
        #
        status=os.system('/bin/tcsh -i -c "source ./trjconv.scr"')
        if status!=0:
            raise Exception('trjconv did not run')
        #
        # Find the snapshots
        #

        files=os.listdir(os.getcwd())
        snapshots=[]
        for fn in files:
            if self.files['snapshots'] in fn:
                snapshots.append(os.path.join(os.getcwd(),fn))
        if len(snapshots)<numsnapshots:
            print 'Got only %d snapshots...' %(len(snapshots))
        return snapshots
