# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:38:40 2018

@author: shafigh.bahamin
"""
import sys

class kinematic_object:
    
    def __init__(self, distance = None, initial_velocity = None, final_velocity = None, acceleration = None, time = None):
        self.distance = distance
        self.initial_velocity = initial_velocity
        self.final_velocity = final_velocity
        self.acceleration = acceleration
        self.time = time
        
        
class math_lib:
    
    def raise_to_power(self, x, y):
        return x**y
    
    def log(self, base, x):
        exp = 0;
        pow_ten = 0
        ans = ""
        for i in range(0, 20):
            while(pow_ten < x):
                exp+=1
                pow_ten = self.raise_to_power(base, exp)
            exp-=1
            pow_ten = self.raise_to_power(base, exp)
            ans +=str(exp)
            if len(ans) == 1:
                ans+='.'
            res = x/float(pow_ten)
            x = res**float(base)
            exp = 0
            pow_ten = 0
        return ans
               
    def sqrt(self, number):
        ans = (number/2.0) / 2.0
        maximum_val = number/2.0
        lookup = ans * ans
        while(lookup!=number):
            lookup = ans * ans
            if((lookup < number + .00000000001) and (lookup > number - .0000000001)):
                for i in range(0, (number/2)/2):
                    if((ans + .0000001 < i + .000001) and (ans + .0000001 > i - .000001)):
                        return i
                return ans
            if(lookup > number):
                maximum_val = ans
                ans = ans / 2.0
            if(lookup < number):
                ans = (ans+maximum_val)/2.0 
        
class physics_and_astronomy:
    
    def __init__(self):
        self.math = math_lib()
        
    def luminosity(self, absolute_magnitude):
        if(absolute_magnitude is None):
            return None
        return (3.0128*self.math.raise_to_power(10, 28)*self.math.raise_to_power(10, -.4*absolute_magnitude))
    
    
    def apparent_magnitude(self, absolute_magnitude, distance):
        distance_pc = 0
        if(distance is not None):
            distance_pc = .306601 * distance
        if(absolute_magnitude is None or distance is None):
            return None
        return 5*float(self.math.log(10,distance_pc/10))+absolute_magnitude
          
    def absolute_magnitude(self, apparent_magnitude, distance):
        distance_pc = 0
        if(distance is not None):
            distance_pc = .306601 * distance
        if(apparent_magnitude is None or distance is None):
            return None
        return apparent_magnitude - (5*float(self.math.log(10,distance_pc/10)))
    
    
    def distance(self, apparent_magnitude, absolute_magnitude):
        if(apparent_magnitude is None or absolute_magnitude is None):
            return None
        return (self.math.raise_to_power(10,(apparent_magnitude - absolute_magnitude)/5)*10)/.306601
    
    #e=hv       e = hc/wavelength
    def energy(self, frequency, wavelength):
        planck_constant = 6.6260693 * self.math.raise_to_power(10, -34)
        speed_of_light = 299792458.0
        if(frequency is None and wavelength is None):
            return None
        if(frequency is None):
            return (planck_constant * speed_of_light) / wavelength
        if(wavelength is None):
            return (frequency * planck_constant)
        else:
            return (frequency * planck_constant)
    
    def frequency(self, energy, wavelength):
        planck_constant = 6.6260693 * self.math.raise_to_power(10, -34)
        speed_of_light = 299792458.0
        if(energy is None and wavelength is None):
            return None
        if(wavelength is None):
            return energy / planck_constant
        else:
            return speed_of_light / wavelength
        
    def kinetic_energy(self, mass, velocity):
        if(mass is None and velocity is None):
            return None
        return 1.0/2*mass*velocity*velocity
    
    def kinetic_energy_atom(self, temprature):
        k = 1.28 * self.math.raise_to_power(10, -23)
        if(temprature is None):
            return None
        return 3.0/2*temprature*k
    
    def velocity_from_energy(self, energy, mass):
        if(mass is None or energy is None):
            return None
        return self.math.sqrt(2.0 * energy / mass)
    
    def velocity_of_atom(self, temprature, mass):
        k = 1.28 * self.math.raise_to_power(10, -23)
        if(temprature is None or mass is None):
            return None
        return self.math.sqrt((3 * k * temprature) / mass)
    
    #working on the kinematics equations
    # distance, time, intitial_velocity, final_velocity, acceleration
    def kinematic_equations(self, kinematic_variables):
        if(isinstance(kinematic_variables, kinematic_object) == False):
            return None
        if(kinematic_variables.distance is None and ):
            
            
            
    
def main():
    m = physics_and_astronomy()
    
    
if __name__ == "__main__":
    main()






