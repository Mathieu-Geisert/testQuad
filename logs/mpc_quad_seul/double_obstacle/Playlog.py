# -*- coding: utf-8 -*-
import numpy as np
from pinocchio.robot_wrapper import RobotWrapper
from hpp import Quaternion
import time
from gepetto.corbaserver import Client

class PlayFromFile2:
    def __init__(self):
        self.client=Client()
        self.client.gui.setRate(100)
        self.client.gui.createWindow("w")
        self.client.gui.createScene("world")
        self.client.gui.addSceneToWindow("world",0L)
        self.client.gui.addMesh("world/firefly","/local/mgeisert/models/firefly_dae/firefly.dae")
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = []
        list1 = []
        tmp3 = 0.
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[7]-1.57,-tmp[8],tmp[9])
            tmp2.extend(x.array[0:4])
            tmp2.extend([tmp[0]-tmp3])
            tmp3 = tmp[0]  
            list1.append(tmp2)
        self.list.append(list1)
    def play(self, i=0):
        for config in self.list[i]:
            self.client.gui.applyConfiguration("world/firefly",config[0:7])
            self.client.gui.refresh()
            time.sleep(config[7])
    def display(self,k ,i=0):
        self.client.gui.applyConfiguration("world/firefly",self.list[i][k][0:7])
        self.client.gui.refresh()
    def reload(self):
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        list2 = []
        tmp3 = 0.
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[7]-1.57,-tmp[8],tmp[9])
            tmp2.extend(x.array[0:4])
            tmp2.extend([tmp[0]-tmp3])
            tmp3 = tmp[0]
            list2.append(tmp2)
        self.list.insert(0,list2)

p = PlayFromFile2()
p.client.gui.addLandmark("world",0.5)
p.client.gui.addCylinder("world/cylinder11",1.2,10.,[1.,0.,0.,1.])
p.client.gui.addCylinder("world/cylinder12",.5,10.,[1.,0.,0.,1.])
p.client.gui.addCylinder("world/cylinder13",.5,10.,[1.,0.,0.,1.])
p.client.gui.addCylinder("world/cylinder21",1.2,10.,[0.,1.,0.,1.])
p.client.gui.addCylinder("world/cylinder22",.5,10.,[0.,1.,0.,1.])
p.client.gui.addCylinder("world/cylinder23",.5,10.,[0.,1.,0.,1.])

p.client.gui.applyConfiguration("world/cylinder11",[0.,-5.,3.,0.7,0.7,0.,0.])
p.client.gui.applyConfiguration("world/cylinder12",[1.3,-5.,3.,0.7,0.7,0.,0.])
p.client.gui.applyConfiguration("world/cylinder13",[-1.3,-5.,3.,0.7,0.7,0.,0.])
p.client.gui.applyConfiguration("world/cylinder21",[2.,-5.,6.,0.7,0.7,0.,0.])
p.client.gui.applyConfiguration("world/cylinder22",[3.3,-5.,6.,0.7,0.7,0.,0.])
p.client.gui.applyConfiguration("world/cylinder23",[0.7,-5.,6.,0.7,0.7,0.,0.])

p.client.gui.refresh()

p.play()
