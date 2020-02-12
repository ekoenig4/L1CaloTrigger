from ROOT import *
from objects import *

def getBox(iphi,ieta,geometry=(7,7),color=kBlack):
    dphi,deta=[d/2 for d in geometry]
    box = TBox(iphi-dphi-1,ieta-deta-1,iphi+dphi,ieta+deta)
    box.SetFillStyle(0)
    box.SetLineColor(color)
    box.SetLineWidth(2)
    return box


class Grid:
    def __init__(self,name,gridmap,center=None):
        self.cphi,self.ceta=0,0
        if center is not None: self.cphi,self.ceta = getPhiEta(center)
        if name=='t22':xmax,ymax=22,22
        elif name=='t19':xmax,ymax=19,19
        elif name=='t7':xmax,ymax=7,7
        else:
            xmax = max( iphi for iphi,ieta in keyiter(gridmap) )
            ymax = max( ieta for iphi,ieta in keyiter(gridmap) )
        grid = TH2F("grid_%s"%name,";iPhi;iEta",xmax,0,xmax,ymax,0,ymax)
        for i in irange(1,xmax): grid.GetXaxis().SetBinLabel(i,str(i))
        for i in irange(1,ymax): grid.GetYaxis().SetBinLabel(i,str(i))
        
        for y in irange(1,ymax):
            for x in irange(1,xmax):
                if getKey(x,y) in gridmap:
                    grid.SetBinContent(x,y,gridmap[getKey(x,y)].et)
        self.grid = grid
        self.boxes = []
   
heirarchy = ['t7','t19','t22','gct','calo']
def t19_to_t22(start,center):
    # (10,10) -> center
    sx,sy = start
    cx,cy = center
    return sx - 10 + cx,sy - 10 + cy
def t22_to_gct(start,center):
    # (11,12) -> center
    sx,sy = start
    cx,cy = center
    return sx - 11 + cx,sy - 12 + cy
def gct_to_calo(start,center):
    sx,sy = start
    cx,cy = center
    
    # (5,1) -> center (for 34x32)
    # return sx - 5 + cx,sy - 1 + cy
    # (10,1) -> center (for 34x42)
    return sx - 10 + cx,sy - 1 + cy
spacemap = {
    't19':t19_to_t22,
    't22':t22_to_gct,
    'gct':gct_to_calo
}

class Debug:
    gridmap = {}
    canvasmap = {}
    def add(self,name,gridmap,center=None):
        if not config.debug: return
        if name in self.gridmap:
            index = heirarchy.index(name)+1
            if index < len(heirarchy):
                for box in self.gridmap[name].boxes:
                    for parent in heirarchy[index:]:
                        if any(self.gridmap[parent].boxes):
                            self.gridmap[parent].boxes.pop()
        self.gridmap[name] = Grid(name,gridmap,center)
        if name == 't22': self.addBox(name,(11.5,11.5,3.0,3.0))
        if name == 't19': self.addBox(name,(10,10,7,7))
        
        self.canvasmap[name] = None
    def addBox(self,name,box):
        if not config.debug: return
        (cx,cy,wx,wy) = box
        grid = self.gridmap[name]
        grid.boxes.append( getBox(cx,cy,geometry=(wx,wy)) )
        # print name,box,(grid.cphi,grid.ceta)
        index = heirarchy.index(name)+1
        if index < len(heirarchy):
            parent=heirarchy[index]
            nx,ny = spacemap[name]( (cx,cy),(grid.cphi,grid.ceta) )
            self.addBox(parent,box=(nx,ny,wx,wy))
    def draw(self,name,gridmap=None,pause=False,box=None):
        if not config.debug: return
        if name not in self.gridmap and gridmap is not None: self.add(name,gridmap)
        
        grid = self.gridmap[name]
        
        center = None

        # if type(box) is not tuple:
        
        # if type(box) is tuple: self.addBox(name,box)

        if any(grid.boxes) or name == 't7':
            index = heirarchy.index(name)+1
            if index < len(heirarchy):
                parent=heirarchy[index]
                self.draw(parent)
        hs = grid.grid
        c = TCanvas("canvas_%s"%name,"")
        c.cd()
        c.SetGrid()
        gStyle.SetOptStat(0)
        hs.Draw("COLZ")
        if any(grid.boxes):
            for box in grid.boxes: box.Draw('same')
        c.grid=grid
        self.canvasmap[name]=c
        
        if pause:
            if config.play: time.sleep(1)
            else: raw_input("%s Paused."%name)
        
debugger = Debug()
