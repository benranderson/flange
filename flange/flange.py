from design_cond import DesignCond
from material import Material

class Scenario(object):
    def __init__(self, props):
        self.ID = props[0]
        self.kind = props[1] # "ASME"
        self.flange = props[2] #Flange()
        self.des_cond = props[3] #DesignCond()
        self.bolt_des_fac = props[4]

class Flange(object):
    def __init__(self, props):
        self.ID = props[0] 
        self.body = props[1] # Body(self.kind)
        self.body_mat = props[2]
        self.bolt = props[3] # Bolt(self.kind)
        self.bolt_mat = props[4]
        self.ringjoint = props[5] # RingJoint()
        
class Body(object):
    def __init__(self, props):
        '''
        API 6A and API 17D
        '''
        
        self.ID = props[0]
        self.A = props[1] #685.8 # flange OD
        self.B = props[2] #385.8 # flange bore (pipe bore)
        self.C = props[3] #603.3 # bolt circle dia
        self.h = props[4] #76 # hub length
        self.g_1 = props[5] #57.1 # hub thickness at back of flange
        self.g_0 = props[6] #12.7 # hub thickness at small end
        self.t = props[7] #76.2 # flange thickness
        #self.grade = props[2] #"F55"
        #self.material = props[8] #1 #Material(self.grade, kind)
        #self.des_fac = props[8]

        
class Bolt(object):
    def __init__(self, props):
        '''
        API 6A Annex D
        '''
        self.ID = props[0]
        self.num = props[1] #20
        self.A_r = props[2] #907.9
        #self.material = props[3] #2
        #self.des_fac = props[3] #0.5

        

class RingJoint(object):
    def __init__(self, props):
        '''
        Ring joint size from Flange dimension tables
        BX (SBX): API 6A Table 65 and API 17D Table 6
        - Gasket factors from ASME VII Div 2 Table 4.16.1
        '''
        self.ID = props[0]
        self.w = props[1] # 11.1 # gasket width
        self.y = props[2] # 179 # gasket seating pressure load
        self.m = props[3] # 6.5 # gasket factor
        self.G = props[4] # 469.9 # mean gasket diameter
        
        



        