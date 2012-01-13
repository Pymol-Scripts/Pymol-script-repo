import sys
sys.path.append('/Users/nielsen/bin/APBS/tools/python/vgrid')
from vgrid import *
NMAX=5

class phimap_class:

    def __init__(self,xmap,ymap,zmap):
        #print
        #print 'Reading initial Epsmaps'
        #print
        self.maps={}
        for name,emap in [['x',xmap],['y',ymap],['z',zmap]]:
            print 'Reading',emap
            data=[]
            value=0.0
            
            startVio()
            import sys
            from sys import stdout, stderr
            grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data)
            Vgrid_readDX(grid, "FILE", "ASC", "", emap)

            nx = grid.nx
            ny = grid.ny
            nz = grid.nz
            hx = grid.hx
            hy = grid.hy
            hzed = grid.hzed
            xmin = grid.xmin
            ymin = grid.ymin
            zmin = grid.zmin
            # Get the grid data
            grid_data=[]
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):                    
                        inval=0.0
                        pt = [xmin + i*hx, ymin + j*hy, zmin + k*hzed]
                        ret, value = Vgrid_value(grid, pt, inval)
                        grid_data.append(value)
                        if not ret:
                            raise Exception('Grid point not found')
        return grid_data


    def writemap(self,map,filename='tester.dx'):
        """Write this map in dx format"""
        mydata=[]
        #count=0
        for k in range(map.nz):
            for j in range(map.ny):
                for i in range(map.nx):
                    mydata.append(map.get_value(i,j,k))
        #
        # Create the grid
        #
        #print mydata
        mygrid = Vgrid_ctor(map.nx, map.ny, map.nz,
                            map.hx, map.hy, map.hz,
                            map.xmin, map.ymin, map.zmin, mydata)
        Vgrid_writeDX(mygrid, "FILE", "ASC", "", filename,'What is this', null_array())
        delete_vgrid(mygrid)
        return
