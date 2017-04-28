import math

class Asme8_Div1(object):
    def __init__(self, scenario, flange, body, body_mat, bolt, bolt_mat, ringjoint, des_cond, bolt_des_fac):
        '''
        Class containing calculations to determine internal flange stresses
        resulting from the internal pressure, external tensile force and
        external bending moment for weld neck flanges in accordance with
        ASME Boiler and Pressure Vessel Code, Section VIII, Division 1,
        Appendix 2.
        '''
        # Analysis description
        self.ID = scenario.ID
        self.flange_ID = scenario.flange
        self.des_cond_ID = scenario.des_cond
        self.bolt_des_fac = scenario.bolt_des_fac
        
        # Run methods
        self.gasket_dim(ringjoint)
        self.total_effective_pressure(des_cond, ringjoint)
        self.flange_design_strength(body_mat, bolt_mat)
        self.bolt_load(bolt, ringjoint, des_cond, bolt_mat, bolt_des_fac)
        self.force_moment(body, ringjoint)
        self.flange_factors(body)
        self.stresses_calc(body, bolt)
        self.stresses_allow()
        #self.bolt_ur()
        self.unity_check()
        self.rigid_i(body, body_mat)
    
    def gasket_dim(self, ringjoint):
        '''
        Function to determine effective gasket seating width
        '''
        self.b_o = ringjoint.w / 8 # (Table 2-5.2)
        if self.b_o > 6:
            self.b = 2.52 * math.sqrt(self.b_o) # (Table 2-5.2)
        else:
            self.b = self.b_o
    
    def total_effective_pressure(self, des_cond, ringjoint):
        P_d = des_cond.P_d
        F = des_cond.F
        M = des_cond.M
        G = ringjoint.G

        self.P_eff = P_d + ((4*0.001*F)/(math.pi*(0.001*G)**2)) + ((16*0.001*M)/(math.pi*(0.001*G)**3))
        
        print "P_eff = %s " % self.P_eff
    
    def flange_design_strength(self, body_mat, bolt_mat):
        '''
        Place holder for flange kind specific functions to determine flange
        design strength.
        - Asme_Div1_Api()
        - Asme_Div1_Asme()
        '''
        pass
    
    def bolt_load(self, bolt, ringjoint, des_cond, bolt_mat, bolt_des_fac):
        '''
        Function to determine the required bolt load for the operational and 
        gasket seating conditions.
        ASME VIII - 2-5
        '''
        # Extract inputs
        self.A_b = bolt.num*bolt.A_r
        G = ringjoint.G
        m = ringjoint.m
        y = ringjoint.y
        P_d = des_cond.P_d
        
        # Operational condition
        self.H = (math.pi/4)*G**2*self.P_eff # Total hydrostatic end force
        H_p = 2*self.b*math.pi*G*m*self.P_eff # Total joint compressive force
        self.W_m1 = self.H + H_p # Required bolt load (1)
        self.S_a = bolt_des_fac*self.s_y_a_b # Allowable bolt stress (note for ASME flange, bolt_des_fac = 1)
        A_m1 = self.W_m1/self.S_a # Required bolt area
        
        # Gasket seating condition
        W_m2 = math.pi*self.b*G*y # Required bolt load (2)
        A_m2 = W_m2/self.S_a # Required bolt area
        
        # Total required bolt area
        W_m = max(self.W_m1, W_m2) # Total required bolt load
        self.A_m = max(A_m1, A_m2) # Total required bolt area
        #self.ur_bolt_area = self.A_m/self.A_b # Bolt area unity check
        self.bolt_ur = self.A_m/self.A_b
        
        print "A_b = %s " % self.A_b
        print "H = %s " % self.H
        print "H_p = %s " % H_p
        print "W_m1 = %s " % self.W_m1
        print "S_a = %s " % self.S_a
        print "A_m1 = %s " % A_m1
        print "A_m2 = %s " % A_m2

    
    def force_moment(self, body, ringjoint):
        '''
        Function to calculate the forces and moments acting on the flange
        due to the loading.
        ASME VIII - 2-5, 2-6
        '''
        # Operational condition-------------------------------------------------
        # Forces
        W_o = self.W_m1 # Flange design bolt load (3)
        B = body.B
        H_d = (math.pi/4)*B**2*self.P_eff # Hydrostatic end load inside flange
        H_g = W_o - self.H # Gasket load
        H_t = self.H - H_d # Difference between hydrostatic end force and gasket load
        
        # Lever arms
        C = body.C
        g_1 = body.g_1
        G = ringjoint.G
        
        h_d = (C - B - g_1) / 2.0 # Radial distance from bolt circle
        h_g = (C - G) / 2.0 # Radial distance from gasket force to bolt circle
        h_t = (C - B) / 4.0 + h_g/2.0 # Radial distance from bolt circle to H_t circle
        
        # Moments
        M_d = H_d*h_d
        M_g = H_g*h_g
        M_t = H_t*h_t
        self.M_tot_o = M_d + M_g + M_t # Total flange moment
        
        # Gasket seating condition----------------------------------------------
        W_g_min = 0.5*(self.A_m + self.A_b)*self.S_a # Flange design load (4)
        W_g_act = self.A_b * self.S_a # Actual applied flange load
        W_g = max(W_g_min, W_g_act)
        
        self.M_tot_g = W_g_act*h_g # Total flange moment (5)
        # Alternative conservative method would be to use the actual applied 
        # bolt load (W_g_act). Not necessarily representative however as gasket 
        # deforms.
        
        # Determine governing load condition
        if self.M_tot_o > self.M_tot_g:
            self.govern = "Operational"
        else:
            self.govern = "Gasket Seating"
            
        print "H_d = %s " % H_d
        print "H_g = %s " % H_g
        print "H_t = %s " % H_t
        print "h_d = %s " % h_d
        print "h_g = %s " % h_g
        print "h_t = %s " % h_t
        print "M_d = %s " % M_d
        print "M_g = %s " % M_g
        print "M_t = %s " % M_t
        print "M_tot_o = %s " % self.M_tot_o
        print "W_g_min = %s " % W_g_min
        print "W_g_act = %s " % W_g_act
        print "W_g = %s " % W_g
        print "M_tot_g = %s " % self.M_tot_g
        
    def flange_factors(self, body):
        '''
        Function to calculate flange factors.
        ASME VIII - Fig. 2-7.1, Table 2-7.1
        '''
        # Shape constants (Fig. 2-7.1)
        K = body.A / float(body.B)
        T = (K**2 * (1 + 8.55246 * math.log10(K)) - 1) / ((1.04720 + 1.9448 * K**2) * (K - 1))
        U = (K**2 * (1 + 8.55246 * math.log10(K)) - 1) / (1.36136 * (K**2 - 1) * (K - 1))
        self.Y = (1 / (K - 1)) * (0.66845 + 5.7169 * (K**2 * math.log10(K)) / (K**2 - 1))
        self.Z = (K**2 + 1) / (K**2 - 1)
        
        # Flange factors (Table 2-7.1)
        A = body.g_1/float(body.g_0) - 1
        self.h_0 =math.sqrt(body.B*body.g_0)
        C = 43.68*(body.h/self.h_0)**4

        C_1 = (1/3.0) + (A/12.0)
        C_2 = (5/42.0) + 17*(A/336.0)
        C_3 = (1/210.0) + (A/360.0)
        C_4 = (11/360.0) + 59*(A/5040.0) + (1+3*A)/C
        C_5 = (1/90.0) + 5*(A/1008.0) - (1+A)**3/C
        C_6 = (1/120.0) + 17*(A/5040.0) + (1/C)
        C_7 = (215/2772.0) + 51*(A/1232.0) + ((60/7.0) + 225*(A/14.0) + 75*(A**2/7.0) + 5*(A**3/2.0))/C
        C_8 = (31/6930.0) + 128*(A/45045.0) + ((6/7.0) + 15*(A/7.0) + 12*(A**2/7) + 5*(A**3/11.0))/C
        C_9 = (533/30240.0) + 653*(A/73920.0) + (0.5 + 33*(A/14.0) + 39*(A**2/28.0) + 25*(A**3/84.0))/C
        C_10 = (29/3780.0) + 3*(A/704.0) - (0.5 + 33*(A/14.0) + 81*(A**2/28.0) + 13*(A**3/12.0))/C
        C_11 = (31/6048.0) + 1763*(A/665280) + (0.5 + 6*(A/7.0) + 15*(A**2/28.0) + 5*(A**3/42.0))/C
        C_12 = (1/2925.0) + 71*(A/300300.0) + ((8/35.0) + 18*(A/35.0) + 156*(A**2/385.0) + 6*(A**3/55))/C
        C_13 = (761/831600.0) + 937*(A/1663200) + ((1/35.0) + 6*(A/35.0) + 11*(A**2/70.0) + 3*(A**3/70.0))/C
        C_14 = (197/415800.0) + 103*(A/332640.0) - ((1/35.0) + 6*(A/35.0) + 17*(A**2/70.0) + (A**3/10.0))/C
        C_15 = (233/831600.0) + 97*(A/554400.0) + ((1/35.0) + 3*(A/35.0) + (A**2/14.0) + 2*(A**3/105.0))/C
        C_16 = C_1*C_7*C_12 + C_2*C_8*C_3 + C_3*C_8*C_2 - (C_3**2*C_7 + C_8**2*C_1 + C_2**2*C_12)
        C_17 = (C_4*C_7*C_12 + C_2*C_8*C_13 + C_3*C_8*C_9 - (C_13*C_7*C_3 + C_8**2*C_4 + C_12*C_2*C_9))/C_16
        C_18 = (C_5*C_7*C_12 + C_2*C_8*C_14 + C_3*C_8*C_10 - (C_14*C_7*C_3 + C_8**2*C_5 + C_12*C_2*C_10))/C_16
        C_19 = (C_6*C_7*C_12 + C_2*C_8*C_15 + C_3*C_8*C_11 - (C_15*C_7*C_3 + C_8**2*C_6 + C_12*C_2*C_11))/C_16
        C_20 = (C_1*C_9*C_12 + C_4*C_8*C_3 + C_3*C_13*C_2 - (C_3**2*C_9 + C_13*C_8*C_1 + C_12*C_4*C_2))/C_16
        C_21 = (C_1*C_10*C_12 + C_5*C_8*C_3 + C_3*C_14*C_2 - (C_3**2*C_10 + C_14*C_8*C_1 + C_12*C_5*C_2))/C_16
        C_22 = (C_1*C_11*C_12 + C_6*C_8*C_3 + C_3*C_15*C_2 - (C_3**2*C_11 + C_15*C_8*C_1 + C_12*C_6*C_2))/C_16
        C_23 = (C_1*C_7*C_13 + C_2*C_9*C_3 + C_4*C_8*C_2 - (C_3*C_7*C_4 + C_8*C_9*C_1 + C_2**2*C_13))/C_16
        C_24 = (C_1*C_7*C_14 + C_2*C_10*C_3 + C_5*C_8*C_2 - (C_3*C_7*C_5 + C_8*C_10*C_1 + C_2**2*C_14))/C_16
        C_25 = (C_1*C_7*C_15 + C_2*C_11*C_3 + C_6*C_8*C_2 - (C_3*C_7*C_6 + C_8*C_11*C_1 + C_2**2*C_15))/C_16
        C_26 = -(C/4.0)**(1/4.0)
        C_27 = C_20 - C_17 - (5/12.0) - (C_17*(C/4.0)**(1/4.0))
        C_28 = C_22 - C_19 - (1/12.0) - (C_19*(C/4.0)**(1/4.0))
        C_29 = -(C/4.0)**(1/2.0)
        C_30 = -(C/4.0)**(3/4.0)
        C_31 = 3*(A/2.0) + C_17*(C/4.0)**(3/4.0)
        C_32 = 0.5 + C_19*(C/4.0)**(3/4.0)
        C_33 = 0.5*C_26*C_32 + C_28*C_31*C_29 - (0.5*C_30*C_28 + C_32*C_27*C_29)
        C_34 = (1/12.0) + C_18 - C_21 + C_18*(C/4.0)**(1/4.0)
        C_35 = -C_18*(C/4.0)**(3/4.0)
        C_36 = (C_28*C_35*C_29 - C_32*C_34*C_29)/C_33
        C_37 = (0.5*C_26*C_35 + C_34*C_31*C_29 - (0.5*C_30*C_34 + C_35*C_27*C_29))/C_33

        E_1 = C_17*C_36 + C_18 + C_19*C_37
        E_2 = C_20*C_36 + C_21 + C_22*C_37
        E_3 = C_23*C_36 + C_24 + C_25*C_37
        E_4 = (1/4.0) + (C_37/12.0) + (C_36/4.0) - (E_3/5.0) - 3*(E_2/2) - E_1
        E_5 = E_1*(0.5 + (A/6.0)) + E_2*((1/4.0) + 11*(A/84.0)) + E_3*((1/70.0) + (A/105.0))
        E_6 = E_5 - C_36*((7/120.0) + (A/36.0) + 3*(A/C)) - (1/40.0) - (A/72.0) - C_37*((1/60.0) + (A/120.0) + (1/C))
        
        F = -E_6 / ((C/2.73)**(1/4.0)*((1+A)**3/C))
        self.V = E_4 / ((2.73/C)**(1/4.0)*(1+A)**3)
        f_1 = C_36 / (1 + A)
        
        if f_1 > 1:
            self.f = f_1
        else:
            self.f = 1
            
        d = (U/self.V)*self.h_0*body.g_0**2
        self.e = F/self.h_0
        self.L = (body.t*self.e + 1)/T + (body.t**3/d)

        print "K = %s " % K
        print "T = %s " % T
        print "U = %s " % U
        print "Y = %s " % self.Y
        print "Z = %s " % self.Z
        print "F = %s " % F
        print "V = %s " % self.V
        print "f = %s " % self.f
        print "d = %s " % d
        print "e = %s " % self.e
        print "L = %s " % self.L
        
    
    def stresses_calc(self, body, bolt):
        '''
        Function to calculate the stresses in the flange for the operational
        gasket seating conditions.
        ASME VIII - 2-7
        '''
        # Dictionary of conditions to allow looping
        M_0s = {'Operational': self.M_tot_o, 'Gasket Seating': self.M_tot_g}
        self.calculated_stresses = {}
        
        for case in M_0s:
            
            M_0 = M_0s[case]
        
            S_H = (self.f*M_0) / (self.L*body.g_1**2*body.B) # Longitudinal hub stress (6)
            S_R = ((1.33*body.t*self.e + 1)*M_0) / (self.L*body.t**2*body.B) # Radial flange stress (7)
            S_T = (self.Y*M_0)/(body.t**2*body.B) - self.Z*S_R # Tangential flange stress (8)
            S_add1 = 0.5*(S_H + S_R) # Additional stress check 1
            S_add2 = 0.5*(S_H + S_T) # Additional stress check 2

            # Populate dictionary with results
            self.calculated_stresses[case] = {}
            self.calculated_stresses[case]['S_H'] = S_H
            self.calculated_stresses[case]['S_R'] = S_R
            self.calculated_stresses[case]['S_T'] = S_T
            self.calculated_stresses[case]['S_add1'] = S_add1
            self.calculated_stresses[case]['S_add2'] = S_add2
            
            self.stress_bolt = self.S_a
            
        
    def stresses_allow(self):
        '''
        Function to determine allowable flange stresses.
        ASME VIII - 2-8
        '''
        self.allow_bolt = self.s_y_d_b
        self.allowable_stresses = {}
        self.allowable_stresses['S_H'] = 1.5*self.s_y_d_f
        self.allowable_stresses['S_R'] = self.s_y_d_f
        self.allowable_stresses['S_T'] = self.s_y_d_f
        self.allowable_stresses['S_add1'] = self.s_y_d_f
        self.allowable_stresses['S_add2'] = self.s_y_d_f


    def bolt_ur(self):
        pass
    
    
    def unity_check(self):
        '''
        Function to determine unity ratios based on calculated flange stresses
        for the operational and gasket seating conditions and the allowable
        flange stresses.
        '''

        self.urs = {}

        self.ur_max = self.bolt_ur

        for case in self.calculated_stresses:
            self.urs[case] = {}
            for stress in self.calculated_stresses[case]:
                
                self.urs[case][stress] = self.calculated_stresses[case][stress]/self.allowable_stresses[stress]
                if self.urs[case][stress] > self.ur_max:
                    self.ur_max = self.urs[case][stress]
                    
    
    def rigid_i(self, body, body_mat):
        '''
        Rigidity check to ensure the flange will not rotate and leak for the
        operational or gasket seating conditions.
        ASME VIII - 2-14
        '''
        self.K_I = 0.3 # Rigity factor for integral or optional flange
        M = max(self.M_tot_o, self.M_tot_g) # Maximum total moment acting on flange
        self.J = 52.14 * self.V * M / (self.L * body_mat.E * body.g_0**2 * self.K_I * self.h_0)
        
        print "K_I = " + str(self.K_I)
        print "J = " + str(self.J)
        
    
class Asme_Div1_Api(Asme8_Div1):
    def flange_design_strength(self, body_mat, bolt_mat):
        '''
        Function to determine allowable stresses based on API criteria
        API flanges use full yield strength values with design factor equal to
        0.5 along with the following deratings
        '''
        self.kind = "API"
        
        self.s_y_d_f = (2.0/3.0)*body_mat.s_y_api
        self.s_y_d_b = 0.83*bolt_mat.s_y_api
        self.s_y_a_b = bolt_mat.s_y_api
        
        print "s_y_d_f = %s MPa" % self.s_y_d_f
        print "s_y_d_b = %s MPa" % self.s_y_d_b
        
    def bolt_ur(self):
        #self.bolt_ur = self.stress_bolt/self.allow_bolt
        self.bolt_ur = self.A_m/self.A_b

        
class Asme_Div1_Asme(Asme8_Div1):
    def flange_design_strength(self, body_mat, bolt_mat):
        '''
        Function to determine allowable stresses based on API criteria
        ASME flanges use lower yield strength values from code with design 
        factors equal to 1
        '''
        self.kind = "ASME"
        
        self.s_y_d_f = body_mat.s_y_asme
        self.s_y_d_b = bolt_mat.s_y_asme
        self.s_y_a_b = bolt_mat.s_y_asme
        
        print "s_y_d_f = %s MPa" % self.s_y_d_f
        print "s_y_d_b = %s MPa" % self.s_y_d_b
        
    def bolt_ur(self):
        self.bolt_ur = self.A_m/self.A_b