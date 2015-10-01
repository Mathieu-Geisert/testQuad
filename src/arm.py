import numpy as np
from pinocchio.robot_wrapper import RobotWrapper
from hpp import Quaternion
import time
from gepetto.corbaserver import Client
 
class PlayFromFile:
     def __init__(self,log="/tmp/log_pickup.txt"):
         self.wrapper = RobotWrapper("/local/mgeisert/models/ur5_quad.urdf")
         self.client=Client()
         self.client.gui.createWindow("w")
         self.client.gui.createScene("world")
         self.client.gui.addSceneToWindow("world",0L)
         self.client.gui.addURDF("world/ur5","/local/mgeisert/models/ur5_quad.urdf","/local/mgeisert/models/universal_robot/")
         self.client.gui.addSphere("world/sphere",0.07,[0.7,0.,0.,1.])
         self.client.gui.createGroup("world/gripper")
         self.client.gui.addLandmark("world/gripper",0.1)
         self.wrapper.initDisplay("world/ur5")
         self.list = []
         file = open(log, "r")
         configs = file.readlines()
         list1 = []
         tmp3 = 0.
         for i in range(len(configs)):
             tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
             tmp2 = []
             if tmp[0]-tmp3>0.:
               tmp2.extend([tmp[0]-tmp3])
               tmp2.extend(tmp[1:4])
               x = Quaternion().fromRPY(tmp[4],tmp[5],tmp[6])
               tmp2.extend(x.array[1:4])
               tmp2.extend([x.array[0]])
               tmp2.extend(tmp[7:13])
               tmp2.append(tmp[0])
               tmp3 = tmp[0]
               list1.append(tmp2)
         self.list.append(list1)
     def play(self, i=0, speed=1):
         for config in self.list[i]:
             q = np.matrix([[config[1]],[config[2]],[config[3]],[config[4]],[config[5]],[config[6]],[config[7]],[config[8]],[config[9]],[config[10]],[config[11]],[config[12]],[config[13]]])
             self.wrapper.display(q)
             x = np.matrix([[0],[0.15],[0]])
             x = self.wrapper.data.oMi[7]*x
             x = x.transpose()
             x = x.tolist()[0]
             print x
             quat = Quaternion(self.wrapper.data.oMi[7].rotation)
             x.extend(quat.array[0:4])
             self.client.gui.applyConfiguration("world/gripper",x)
             if config[14]<3:
                 self.client.gui.applyConfiguration("world/sphere",[2,-1,-1,1,0,0,0])
             elif config[14]<7:
                 self.client.gui.applyConfiguration("world/sphere",x)
             self.client.gui.refresh()
             time.sleep(config[0]*speed)
     def display(self,k, i=0):
         config = self.list[i][k]
         q = np.matrix([[config[1]],[config[2]],[config[3]],[config[4]],[config[5]],[config[6]],[config[7]],[config[8]],[config[9]],[config[10]],[config[11]],[config[12]],[config[13]]])
         self.wrapper.display(q)
     def reload(self,file="/tmp/log_pickup.txt"):
         file = open(file, "r")
         configs = file.readlines()
         list1 = []
         tmp3 = 0.
         for i in range(len(configs)):
             tmp = np.fromstring(configs[i][2:], dtype=float, sep="\t")
             tmp2 = []
             if tmp[0]-tmp3>0.:
               tmp2.extend([tmp[0]-tmp3])
               tmp2.extend(tmp[1:4])
               x = Quaternion().fromRPY(tmp[4],tmp[5],tmp[6])
               tmp2.extend(x.array[1:4])
               tmp2.extend([x.array[0]])
               tmp2.extend(tmp[7:13])
               tmp2.append(tmp[0])
               tmp3 = tmp[0]
               list1.append(tmp2)
         self.list.insert(0,list1)
 
 #p = PlayFromFile("/local/mgeisert/control/tests/logs/quad_arm/log_pickup_1000it.txt")
p = PlayFromFile("/tmp/log_pickup_ocp.txt")
p.client.gui.createGroup("world/final")
p.client.gui.createGroup("world/pickup")
p.client.gui.applyConfiguration("world/final",[7,0,0,1,0,0,0])
p.client.gui.applyConfiguration("world/pickup",[2,-1,-1,1,0,0,0])
p.client.gui.createGroup("world/place")
p.client.gui.applyConfiguration("world/place",[5,1,-1,1,0,0,0])
p.client.gui.addLandmark("world",0.5)
p.client.gui.addLandmark("world/final",0.5)
p.client.gui.addLandmark("world/pickup",0.25)
p.client.gui.addLandmark("world/ur5/wrist_3_link",0.1)
p.client.gui.addLandmark("world/place",0.25)
p.play()

