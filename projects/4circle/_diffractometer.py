
class Diffractometer():

    def __init__(self):
        
        self.chi = 0
        self.omega = 0
        self.phi = 0
        self.theta = 0
        
    def _get_chi(self):
    
        return self.chi
        
    def _set_chi(self, chi):
        
        self.chi = chi
        
    def _get_omega(self):
    
        return self.omega
        
    def _set_omega(self, omega):
        
        self.omega = omega
        
    def _get_phi(self):
    
        return self.phi
        
    def _set_phi(self, phi):
        
        self.phi = phi
        
    def _get_theta(self):
    
        return self.theta
        
    def _set_theta(self, theta):
        
        self.theta = theta
        
    def __repr__(self):
    
        return "Instance of the class Diffractometer"
        