from Tkinter import *

class charge_mon(Frame):

    def __init__(self):
        Frame.__init__(self)
        self.master.title('Charge monitor')
        height=800
        width=2000
        self.width=width
        self.cv=Canvas(self.master,bd=5,bg='white',
                       width=width,
                       height=height,
                       scrollregion=(0,0,width,height))
        self.cv.grid(row=0,column=0)
        self.calc=0
        self.seqstart=200
        self.text=''
        return

    def init_protein(self,pkaroutines):
        # Initialize the protein
        # count the number of residues
        self.numres=0
        for chain in pkaroutines.protein.chains:
            for residue in chain.residues:
                self.numres=self.numres+1
        self.numres=float(self.numres)
        #
        self.cv.create_text(0,self.calc,text='Setup',anchor='nw')
        x_count=self.seqstart
        self.res_pos={}
        for chain in pkaroutines.protein.chains:
            for residue in chain.residues:
                if residue.name=='ASP' or residue.name=='GLU':
                    fill='red'
                elif residue.name=='LYS' or residue.name=='ARG':
                    fill='blue'
                else:
                    fill='black'
                self.cv.create_text(x_count,self.calc,text='%3d' %residue.resSeq,anchor='nw',fill=fill)
                self.res_pos[residue.resSeq]=x_count
                x_count=x_count+int((self.width-100)/self.numres)
        self.calc=self.calc+15
        self.master.update()
        return

    def set_calc(self,text):
        self.text=text
        return

    def display_charges(self,charge_list):
        #
        # Print the calc
        #
        self.cv.create_text(0,self.calc,text=self.text,anchor='nw')
        charges={}
        for resnum,atomname,charge in charge_list:
            if not charges.has_key(resnum):
                charges[resnum]=[]
            charges[resnum].append(charge)
        #
        # Sum all charges
        #
        for res in charges.keys():
            non_zero=None
            sum=0.0
            for crg in charges[res]:
                sum=sum+crg
                if crg!=0.0:
                    non_zero=1
            if non_zero:
                charges[res]=sum
            else:
                charges[res]=None
        #
        #
        #
        later=[]
        for resid in charges.keys():
            x_count=self.res_pos[resid]
            if charges[resid] is None:
                fill='white'
            elif abs(charges[resid])<0.001:
                fill='grey'
            elif abs(charges[resid]-1.0)<0.001:
                fill='blue'
            elif abs(charges[resid]+1.0)<0.001:
                fill='red'
            else:
                fill='yellow'

            self.cv.create_rectangle(x_count,self.calc,x_count+int((self.width-100)/self.numres),self.calc+10,fill=fill)
            if fill=='yellow':
                later.append([x_count,'%4.2f' %charges[resid],resid])
        #
        # Print all the wrong charges
        #
        for x_count,text,resid in later:
            self.cv.create_text(x_count,self.calc,text=text,anchor='nw',fill='black')
            print '!!Wrong charge: %s %s' %(text,str(resid))
        #
        # Update and increment row
        #
        self.master.update()
        self.calc=self.calc+15
        return
