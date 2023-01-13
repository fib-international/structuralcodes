"""Creep formulas MC2010 for concrete"""




class fib_MC2010_model:
    def __init__(self, concrete: dict, geometry: dict, environment: dict, load: dict):
        """ fib Model Code 2010 model according to [1].

        Args:
            concrete: dict() with concrete related parameters.
            geometry: dict() with specimen geometry related parameters.
            environment: dict() with environment related parameters.
            load: dict() with load related parameters.
        """
        self.concrete = concrete
        self.geometry = geometry
        self.environment = environment
        self.load = load

    def check_input(self):
        """
            Check if input can be handled by the model, and if certain parameter values are within the range of
            application (i.e. the range the model has been calibrated).
        """
        # Check if model considers the specified aggregate type (only used for calculating E at 28 days!):
        agg_types = ['Basalt', 'Quartzite', 'Limestone', 'Sandstone', '']
        if self.concrete['agg_type'] not in agg_types:
            raise ValueError(f"The specified cement type {self.concrete['agg_type']} is not implemented. It can only "
                             f"be one of these: {list(agg_types)}.")

        # Check if the mean concrete compressive strength is in the range of applicability:
        if self.concrete['fcm'] < 20*MPa or self.concrete['fcm'] > 130*MPa:
            print('OUT OF APPLICABILITY RANGE! The model is only applicable for a fcm-value between 20-130 MPa.')

        # Check if the compressive stress (load at t0) is in the range of applicability:
        if self.load['sigma'] < -0.4*self.concrete['fcm']:
            print('OUT OF APPLICABILITY RANGE! The concrete cannot be considered anymore as a concrete is considered '
                  'as an aging linear visco-elastic material. Use subclause 5.1.9.4.3 (d) of [1] for larger compressive'
                  ' stress values.')

        # Check if the age of loading is in the range of applicability:
        if self.load['tp'] < 1:
            print('OUT OF APPLICABILITY RANGE! The age of loading must at least be one day.')

        # Check if the humidity is in the range of applicability:
        if self.environment['RH'] < 0.4 or self.environment['RH'] > 1:
            raise ValueError('OUT OF APPLICABILITY RANGE! The model is only applicable for a RH-value between 40-100%.')

        # Check if the environmental temperature is in the range of applicability:
        if self.environment['T'] < 5 or self.environment['T'] > 30:
            print('OUT OF APPLICABILITY RANGE! The model is only applicable for a temperature between 5 and 30 degr.. '
                  'Use subclause 5.1.10 of [1] for a temperature range between 0 and 80 degr..')

    def get_tabulated_parmvalues(self):
        """
            Get calibrated (autogenous) shrinkage and creep model parameters belonging to the specified cement type.
        """
        # Get cement type dependent parameters:
        if self.concrete['cem_type'] == 'R':  # Assumption: strength classes 32.5 R, 42.5 N correspond with Regular

            # time correction factor (see Eq. 5.1-73 of [1]):
            alpha = 0

            # E-modulus related (see table 5.1-9 of [1]):
            if self.concrete['fcm'] / MPa > 60:
                s = 0.20
            else:
                s = 0.25

            # shrinkage related (see table 5.1-12 of [1]):
            alpha_as = 700
            alpha_ds1 = 4
            alpha_ds2 = 0.012

        elif self.concrete['cem_type'] == 'RS':  # Assumption: strength classes 42.5 R, 52.5 N, 52.5 R correspond with Rapid Hardening

            # time correction factor (see Eq. 5.1-73 of [1]):
            alpha = 1

            # E-modulus related (see table 5.1-9 of [1]):
            s = 0.20

            # shrinkage related (see table 5.1-12 of [1]):
            alpha_as = 600
            alpha_ds1 = 6
            alpha_ds2 = 0.012

        elif self.concrete['cem_type'] == 'SL':  # Assumption: strength class 32.5 N corresponds with Slow Hardening

            # time correction factor (see Eq. 5.1-73 of [1]):
            alpha = -1

            # E-modulus related (see table 5.1-9 of [1]):
            if self.concrete['fcm'] / MPa > 60:
                s = 0.20
            else:
                s = 0.38

            # shrinkage related (see table 5.1-12 of [1]):
            alpha_as = 800
            alpha_ds1 = 3
            alpha_ds2 = 0.013

        # Get aggregate dependent parameter; scaling factors for shrinkage (see Table 6 of [1]):
        if self.concrete['agg_type'] == 'Basalt':
            alpha_e = 1.2
        elif self.concrete['agg_type'] == 'Quartzite':
            alpha_e = 1.0
        elif self.concrete['agg_type'] == 'Limestone':
            alpha_e = 0.9
        elif self.concrete['agg_type'] == 'Sandstone':
            alpha_e = 0.7
        else:
            alpha_e = 1.0

        self.parm_cemtype = {'alpha': alpha, 's': s, 'alpha_as': alpha_as, 'alpha_ds1': alpha_ds1, 'alpha_ds2': alpha_ds2}
        self.parm_aggrtype = {'alpha_e': alpha_e}

    def calc_modified_tp(self, tp, T_cur):
        """
            Calculate the modified age at loading (tp) to account for the effect of the type of cement and curing
            temperature on the degree of hydration and - in turn - on creep.
        """
        # Temperature corrected concrete age in days at tp (see Eq. 5.1-85 of [1]):
        dt = np.ones(tp)
        T_dt = np.ones(tp) * T_cur
        tpT = np.sum(dt * np.exp(13.65 - (4000 / (273 + T_dt))))

        # Cement type and temperature corrected concrete age in days at tp (see Eq. 5.1-73 of [1]):
        # For slow hardening concrete, the creep coefficient is increased due to the lower modified age at loading.
        tp_adj = tpT * ((9/(2 + tpT**1.2)) + 1)**self.parm_cemtype['alpha']

        return tp_adj

    def calc_E28(self):
        """
            Calculate the modulus of elasticity for normal weight concrete at 28 days using Eq. 5.1-21 of [1].
        """
        E_28 = 21500 * self.parm_aggrtype['alpha_e'] * (self.concrete['fcm']/(MPa*10))**(1/3)

        return E_28

    def calc_E(self, time):
        """
            Calculate the modulus of elasticity for normal weight concrete at time 'time' (not 28 days) using Eqs.
            5.1-51, 5.1-56 and 5.1-57 of [1].
        """
        beta_cc = np.exp(self.parm_cemtype['s'] * (1 - np.sqrt(28/time)))
        beta_e = np.sqrt(beta_cc)
        E_time = beta_e*self.calc_E28()

        return E_time

    def calc_drying_shrinkage(self, time):

        # Calculate the drying shrinkage (eps_sh) (see Eqs. 5.1-77, 5.1-80 - 5.1-83):
        eps_sh0 = 10**(-6) * ((220 + 110*self.parm_cemtype['alpha_ds1']) *
                               np.exp(-self.parm_cemtype['alpha_ds2'] * self.concrete['fcm']/MPa))
        beta_sh = np.sqrt((time - self.environment['t0']) / (0.035*(self.geometry['notional_size']/mm)**2 + (time - self.environment['t0'])))
        beta_s1 = min((35 / (self.concrete['fcm']/MPa))**0.1, 1.0)
        if self.environment['RH'] >= 0.99*beta_s1:
            beta_RH = 0.25
        elif 0.4*beta_s1 <= self.environment['RH'] < 0.99*beta_s1:
            beta_RH = -1.55 * (1 - self.environment['RH']**3)
        else:
            raise ValueError("The specified relative humidity * beta_s1 is not in the range of application.")
        eps_sh = eps_sh0 * beta_RH * beta_sh

        return eps_sh

    def calc_autogenous_shrinkage(self, time):

        # Calculate the autogenous shrinkage (eps_au) (see Eqs. 5.1-76, 5.1-78, 5.1-79):
        eps_au0 = -self.parm_cemtype['alpha_as'] * 10**(-6) * \
                   ((0.1*self.concrete['fcm']/MPa) / (6 + 0.1*self.concrete['fcm']/MPa))**2.5
        beta_au = 1 - np.exp(-0.2 * np.sqrt(time))
        eps_au = eps_au0 * beta_au

        return eps_au

    def calc_creep(self, time):

        # Get the adjusted time for creep:
        tp_adj = self.calc_modified_tp(self.load['tp'], self.environment['T_cur'])

        # Calculate the basic creep coefficient (phi_bc) (see Eqs. 5.1-64 - 5.1-66 of [1]):
        beta_bc_fcm = 1.8 / (self.concrete['fcm']/MPa)**0.7
        beta_bc_t = np.log(((30/tp_adj + 0.035)**2) * (time - self.load['tp']) + 1)
        phi_bc = beta_bc_fcm * beta_bc_t

        # Calculate the drying creep coefficient (phi_dc) (see Eqs. 5.1-67 - 5.1-71 of [1]):
        beta_dc_fcm = 412 / (self.concrete['fcm']/MPa)**1.4
        beta_RH = (1 - self.environment['RH']) / ((0.1*self.geometry['notional_size']/(mm*100))**(1/3))
        beta_dc_tp = 1 / (0.1 + tp_adj**0.2)

        alfa_fcm = np.sqrt(35 / (self.concrete['fcm']/MPa))
        beta_h = min(1.5*self.geometry['notional_size']/mm + 250*alfa_fcm, 1500*alfa_fcm)
        gamma_tp = 1 / (2.3 + 3.5/np.sqrt(tp_adj))
        beta_dc_t = ((time - self.load['tp']) / (beta_h + (time - self.load['tp'])))**gamma_tp

        phi_dc = beta_dc_fcm * beta_RH * beta_dc_tp * beta_dc_t

        # Calculate the total creep coefficient (see Eq. 5.1-63 of [1]):
        phi = phi_bc + phi_dc

        # Include the effect of high stress if needed (see Eq. 5.1-74 of [1]):
        k_sigma = self.load['sigma'] / self.concrete['fcm']
        if -0.6 <= k_sigma <= -0.4:
            phi = np.exp(1.5*(abs(k_sigma)-0.4)) * phi
        elif k_sigma < -0.6:
            raise ValueError(f"The stress level exceeds the range of application.")

        # Calculate the creep compliance (see Eq. 5.1-61 of [1]):
        Eci_tp = self.calc_E(self.load['tp'])
        J = (1 / Eci_tp) + (phi / self.calc_E28())  # in 1/MPa

        return J, phi