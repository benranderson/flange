class DesignCond(object):
    def __init__(self, props):
        
        self.ID = props[0]
        self.P_d = props[1] # 4.7          #Design Pressure (bar)
        #self.T_d = props[2] # 50           #Design Temperature
        #self.T_a = 4
        self.F = props[2] # 0 # kN
        self.M = props[3] # 97 # kNm
        