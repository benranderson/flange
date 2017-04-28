class Material(object):
    def __init__(self, props):
        
        self.ID = props[0]
        self.grade = props[1]
        self.s_y_api = props[2]
        self.s_y_asme = props[3]
        self.E = props[4]

        yield_strengths = {'X52': {'API': 450, 'ASME': 152}, 
        'X60': {'API': 465, 'ASME': 152},
        'L7': {'API': 725, 'ASME': 172},
        'F51': {'API': 450, 'ASME': 152},
        'F55': {'API': 480, 'ASME': 152}}
        
        #self.s_y = yield_strengths[grade][kind]          #Minimum Yield Stress
