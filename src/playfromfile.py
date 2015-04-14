# -*- coding: utf-8 -*-
import numpy as np
from pinocchio.robot_wrapper import RobotWrapper
from hpp import Quaternion
import time
from gepetto.corbaserver import Client

class PlayFromFile:
    def __init__(self):
        self.wrapper = RobotWrapper("/local/mgeisert/models/ur5_quad.urdf")
        self.client=Client()
        self.client.gui.createWindow("w")
        self.client.gui.createScene("world")
        self.client.gui.addSceneToWindow("world",0L)
        self.client.gui.addURDF("world/ur5","/local/mgeisert/models/ur5_quad.urdf","/local/mgeisert/models/universal_robot/")
        self.wrapper.initDisplay("world/ur5/")
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = range(80)
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[4],tmp[5],tmp[6])
            tmp2.extend(x.array[1:4])
            tmp2.extend([x.array[0]])
            tmp2.extend(tmp[7:13])
            self.list[i] = tmp2
    def play(self):
        for config in self.list:
            q = np.matrix([[config[0]],[config[1]],[config[2]],[config[3]],[config[4]],[config[5]],[config[6]],[config[7]],[config[8]],[config[9]],[config[10]],[config[11]],[config[12]]])
            self.wrapper.display(q)
            time.sleep(0.2)
    def display(self,k):
        q = np.matrix([[self.list[k][0]],[self.list[k][1]],[self.list[k][2]],[self.list[k][3]],[self.list[k][4]],[self.list[k][5]],[self.list[k][6]],[self.list[k][7]],[self.list[k][8]],[self.list[k][9]],[self.list[k][10]],[self.list[k][11]],[self.list[k][12]]])
        self.wrapper.display(q)
    def reload(self):
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = range(80)
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[4],tmp[5],tmp[6])
            tmp2.extend(x.array[1:4])
            tmp2.extend([x.array[0]])
            tmp2.extend(tmp[7:13])
            self.list[i] = tmp2

class PlayFromFile2:
    def __init__(self):
        self.client=Client()
        self.client.gui.createWindow("w")
        self.client.gui.createScene("world")
        self.client.gui.addSceneToWindow("world",0L)
        self.client.gui.addMesh("world/firefly","/local/mgeisert/models/firefly_dae/firefly.dae")
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = range(1188)
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[7]-1.57,-tmp[8],tmp[9])
            tmp2.extend(x.array[0:4])
            self.list[i] = tmp2
    def play(self):
        for config in self.list:
            self.client.gui.applyConfiguration("world/firefly",config)
            self.client.gui.refresh()
            time.sleep(0.001)
    def display(self,k):
        self.client.gui.applyConfiguration("world/firefly",self.list[k])
        self.client.gui.refresh()
    def reload(self):
        file = open("/tmp/log.txt", "r")
        configs = file.readlines()
        self.list = range(1188)
        for i in range(len(configs)):
            tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
            tmp2 = []
            tmp2.extend(tmp[1:4])
            x = Quaternion().fromRPY(tmp[7]-1.57,-tmp[8],tmp[9])
            tmp2.extend(x.array[0:4])
            self.list[i] = tmp2




p = PlayFromFile2()
p.client.gui.addLandmark("world",0.5)
p.client.gui.addCylinder("world/cylinder",1.8,10.,[1.,1.,1.,1.])
p.client.gui.applyConfiguration("world/cylinder",[0.,0.,3.,0.7,0.7,0.,0.])
p.play()


q = np.matrix([[p.list[0][0]],[p.list[0][1]],[p.list[0][2]],[p.list[0][3]],[p.list[0][4]],[p.list[0][5]],[p.list[0][6]],[p.list[0][7]],[p.list[0][8]],[p.list[0][9]],[p.list[0][10]],[p.list[0][11]],[p.list[0][12]]])
