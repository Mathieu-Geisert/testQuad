#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pinocchio.robot_wrapper import RobotWrapper
from hpp import Quaternion
import time
from gepetto.corbaserver import Client

class PlayFromFile2:
    def __init__(self,Lpend,mq,mp,t1,t2):
        self.GQ = Lpend*mp/(mp+mq)
        self.GP = Lpend*mq/(mp+mq)
        self.t1 = t1
        self.t2 = t2
        self.client=Client()
        self.client.gui.setRate(100)
        self.client.gui.createWindow("w")
        self.client.gui.createScene("world")
        self.client.gui.addSceneToWindow("world",0L)
        self.client.gui.addSphere("world/object", 0.2,[1,0,0,1])
        self.client.gui.applyConfiguration("world/object",[5,0,0,1,0,0,0])
        self.client.gui.createGroup("world/pend")
        self.client.gui.addMesh("world/pend/firefly","/local/mgeisert/models/firefly_dae/firefly.dae")
        self.client.gui.addLine("world/pend/line",[0.,0.,self.GP],[0.,0.,-self.GQ],[1.,0.,0.,1.])
        self.client.gui.addSphere("world/pend/sphere",0.05,[0.,0.,0.,1.])
        self.client.gui.addSphere("world/pend/object",0.2,[1.,0.,0.,1.])
        self.client.gui.applyConfiguration("world/pend/sphere",[0.,0.,self.GP,1.,0.,0.,0.])
        self.client.gui.applyConfiguration("world/pend/object",[0.,0.,self.GP,1.,0.,0.,0.])
        self.client.gui.setVisibility("world/pend/object","OFF")
        self.qi_pend = Quaternion().fromRPY(3.14,0,1.57)
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = []
        list1 = []
        tmp3 = 0.
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend([tmp[0]-tmp3])
            tmp2.extend(tmp[1:4])
            qr = Quaternion().fromRPY(tmp[9],0,0)
            qp = Quaternion().fromRPY(0,tmp[8],0)
            qy = Quaternion().fromRPY(0,0,tmp[7])
            qmesh = Quaternion().fromRPY(-1.57,0,0)
            q2 = Quaternion().fromRPY(tmp[17],0,0)
            q3 = Quaternion().fromRPY(0,tmp[18],0)
            x2 = self.qi_pend*q2*q3
            x = (self.qi_pend*q2*q3).inv()*qr*qp*qy*qmesh
            tmp2.extend(x.array[0:4])
            tmp2.extend(x2.array[0:4])
            tmp4 = 0
            if self.t1 <= tmp[0] <= self.t2 :
                tmp4 = 1 
            if tmp[0] >= self.t2 :
                tmp4 = 2
            tmp3 = tmp[0]
            tmp2.append(tmp4)
            list1.append(tmp2)
        self.list.append(list1)
    def play(self, i=0):
        for config in self.list[i]:
            a = config[1:4]
            a.extend(config[8:12])
            b = [0,0,-self.GQ]
            b.extend(config[4:8])
            self.client.gui.applyConfiguration("world/pend",a)
            self.client.gui.applyConfiguration("world/pend/firefly",b)
            if config[12] == 0:
                self.client.gui.setVisibility("world/pend/object","OFF")
                self.client.gui.applyConfiguration("world/object",[5,0,0,1,0,0,0])
                self.client.gui.setVisibility("world/object","ON")
            if config[12] == 1:
                self.client.gui.setVisibility("world/object","OFF")
                self.client.gui.setVisibility("world/pend/object","ON")
            if config[12] == 2:
                self.client.gui.setVisibility("world/pend/object","OFF")
                self.client.gui.applyConfiguration("world/object",[10,0,0,1,0,0,0])
                self.client.gui.setVisibility("world/object","ON")
            self.client.gui.refresh()
            time.sleep(config[0])
    def display(self,k ,i=0):
        config = self.list[i][k]
        a = config[1:4]
        a.extend(config[8:12])
        b = [0,0,-self.GQ]
        b.extend(config[4:8])
        self.client.gui.applyConfiguration("world/pend",a)
        self.client.gui.applyConfiguration("world/pend/firefly",b)
        self.client.gui.refresh()
    def reload(self):
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        list2 = []
        tmp3 = 0.
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend([tmp[0]-tmp3])
            tmp2.extend(tmp[1:4])
            qr = Quaternion().fromRPY(tmp[9],0,0)
            qp = Quaternion().fromRPY(0,tmp[8],0)
            qy = Quaternion().fromRPY(0,0,tmp[7])
            qmesh = Quaternion().fromRPY(-1.57,0,0)
            q2 = Quaternion().fromRPY(tmp[17],0,0)
            q3 = Quaternion().fromRPY(0,tmp[18],0)
            x2 = self.qi_pend*q2*q3
            x = (self.qi_pend*q2*q3).inv()*qr*qp*qy*qmesh
            tmp2.extend(x.array[0:4])
            tmp2.extend(x2.array[0:4])
            tmp3 = tmp[0]
            list2.append(tmp2)
        self.list.insert(0,list2)


p = PlayFromFile2(4.,0.9,0.45,3.0217,5.999)
p.client.gui.addLandmark("world",0.5)
p.client.gui.createGroup("world/pick")
p.client.gui.createGroup("world/place")
p.client.gui.addLandmark("world/pick",0.5)
p.client.gui.addLandmark("world/place",0.5)
p.client.gui.applyConfiguration("world/pick",[5,0,0,1,0,0,0])
p.client.gui.applyConfiguration("world/place",[10,0,0,1,0,0,0])

p.play()

