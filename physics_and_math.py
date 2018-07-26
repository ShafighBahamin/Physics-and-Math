# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:38:40 2018

@author: shafigh.bahamin
"""
import sys
from random import randint


"""
remember to put a try catch function in every method in this script
"""
class kinematic_object:
    
    def __init__(self, distance = None, initial_velocity = None, final_velocity = None, acceleration = None, time = None):
        self.nmbrs_vars_set = 0;
        self.distance = distance
        if(distance != None):
            self.nmbrs_vars_set+=1
        self.initial_velocity = initial_velocity
        if(initial_velocity != None):
            self.nmbrs_vars_set+=1
        self.final_velocity = final_velocity
        if(final_velocity != None):
            self.nmbrs_vars_set+=1
        self.acceleration = acceleration
        if(acceleration != None):
            self.nmbrs_vars_set+=1
        self.time = time
        if(time != None):
            self.nmbrs_vars_set+=1
        
        
class math_lib:
    
    #the reason there is no tangent function is becuase tan(r) = sin(r) / cos(r)
    
    def sin(self, radians):
        ans = radians
        # 0 means add and 1 means subtract
        add_sub = 1
        for i in range(3, 41, 2):
            r = self.raise_to_power(radians, i)/self.factorial(i)
            if add_sub == 1:
                ans-=r
                add_sub = 0
            else:
                ans+=r
                add_sub = 1
        return ans
    
    def cos(self, radians):
        ans = 1
        # 0 means add and 1 means subtract
        add_sub = 1
        for i in range(2, 41, 2):
            r = self.raise_to_power(radians, i)/self.factorial(i)
            if add_sub == 1:
                ans -= r
                add_sub = 0
            else:
                ans += r
                add_sub = 1
        return ans
    
    #maybe implement a monte carlo multithreaded algorithm to calc pi
    def pi(self):
        return 3.1415926535897932384626433
    
    def raise_to_power(self, x, y):
        return x**y
    
    def factorial(self, number):
        if(isinstance(number, int) == False):
            return None
        else:
            ret = 1
            for i in range(number, 0, -1):
                ret = ret * i
        return ret
    
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
                for i in range(0, int((number/2)/2)):
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
    
    #ke = p^2/2m
    def kinetic_energy_power(self, mass, power):
        if(mass is None or power is None):
            return None
        return (power * power)/(2.0*mass)
    
    def kinetic_energy(self, mass, velocity):
        if(mass is None or velocity is None):
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
    
    #unfinished you need to implement a cos function in the math library
    def power_velocity(self, force, velocity, radians=None, degrees=None):
        if degrees is None and radians is None:
            return None
        if degrees is not None:
            radians = (degrees/180.0)*self.math.pi
        return force * velocity * self.math.cos(radians)
    
    def energy_from_mass_conversion(self, mass):
        speed_of_light = 299792458.0
        return mass * speed_of_light * speed_of_light
    
    def mass_from_energy_conversion(self, energy):
        speed_of_light = 299792458.0
        return energy / (speed_of_light * speed_of_light)
    
    def create_combinations(self, list_of_items, length_of_each_combination):
        number_of_combinations = self.math.factorial(len(list_of_items)) / self.math.factorial(len(list_of_items)-length_of_each_combination)
        list_of_combinations = []
        while(len(list_of_combinations) != number_of_combinations):
            first_idx = randint(0, len(list_of_items)-1)
            second_idx = randint(0, len(list_of_items)-1)
            third_idx = randint(0, len(list_of_items)-1)
            if(second_idx != first_idx and third_idx != first_idx and second_idx != third_idx):
                combination = list_of_items[first_idx] + " " + list_of_items[second_idx] + " " + list_of_items[third_idx]
                if(combination not in list_of_combinations):
                    list_of_combinations.append(combination)
        return list_of_combinations           
    
    def kinematic_equations(self, kinematic_variables):
        if(isinstance(kinematic_variables, kinematic_object) == False):
            return None
        if(kinematic_variables.nmbrs_vars_set<3):
            return None
        else:
            d = kinematic_variables.distance
            t = kinematic_variables.time
            v_i = kinematic_variables.initial_velocity
            v_f = kinematic_variables.final_velocity
            a = kinematic_variables.acceleration
            calculations_done = False
            
        while(calculations_done == False):
            
            if(a is None):
                if(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(d, (float, int)) != False):
                    a = (v_f * v_f - v_i * v_i) / (2.0 * d)
                elif(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(t, (float, int)) != False):
                    a = (v_f - v_i) / float(t)
                elif(isinstance(d, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(d, (float, int)) != False):
                    a = (d - v_i * t) * 2.0 / (t * t)
                    
            if(d is None):
                if(isinstance(a, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(t, (float, int)) != False):
                    d = v_i * t + 1 / 2.0 * a * t * t
                elif(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(t, (float, int)) != False):
                    d = (v_i + v_f) / 2.0 * t
                elif(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(a, (float, int)) != False):
                    d = (v_f * v_f - v_i * v_i) / (2.0 * a)
            
            if(t is None):
                if(isinstance(d, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(a, (float, int)) != False):
                    if((v_i * v_i) - 4 * a * d >= 0):
                        t = ((v_i * -1.0) - self.math.sqrt((v_i * v_i) - 4 * a * d)) / (2 * a)
                        if(t < 0):
                            t = ((v_i * -1.0) + self.math.sqrt((v_i * v_i) - 4 * a * d)) / (2 * a)
                if(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(a, (float, int)) != False):
                    if((v_f - v_i) / float(a) >= 0):
                        t = (v_f - v_i) / float(a)
                if(isinstance(v_f, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(d, (float, int)) != False):
                    if((2.0 * d) / (v_i + v_f) >= 0):
                        t = (2.0 * d) / (v_i + v_f) 
                    
            if(v_i is None):
                if(isinstance(v_f, (float, int)) != False and isinstance(a, (float, int)) != False and isinstance(d, (float, int)) != False):
                    if((v_f * v_f) - 2.0 * a * d >= 0):
                        v_i = self.math.sqrt((v_f * v_f) - 2.0 * a * d)
                    else:
                        v_i = self.math.sqrt(((v_f * v_f) - 2.0 * a * d) * -1) * -1
                if(isinstance(v_f, (float, int)) != False and isinstance(a, (float, int)) != False and isinstance(t, (float, int)) != False):
                    v_i = v_f - a * float(t)
                if(isinstance(a, (float, int)) != False and isinstance(t, (float, int)) != False and isinstance(d, (float, int)) != False):
                    v_i = (d - 1 / 2.0 * a * t * t) / t
                
            if(v_f is None):
                if(isinstance(a, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(d, (float, int)) != False):
                    if((v_i * v_i) + 2.0 * a * d >= 0):
                        v_f = self.math.sqrt((v_i * v_i) + 2.0 * a * d)
                    else:
                        v_f = self.math.sqrt(((v_i * v_i) + 2.0 * a * d) * -1) * -1
                if(isinstance(a, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(t, (float, int)) != False):
                    v_f = v_i + a * float(t)
                if(isinstance(t, (float, int)) != False and isinstance(v_i, (float, int)) != False and isinstance(d, (float, int)) != False):
                    v_f = (2.0 * d / t) - v_i 
                            
            if(d is not None and t is not None and v_i is not None and v_f is not None and a is not None):
                kinematic_variables.distance = d
                kinematic_variables.time = t
                kinematic_variables.initial_velocity = v_i
                kinematic_variables.final_velocity = v_f
                kinematic_variables.acceleration = a
                calculations_done = True
                return kinematic_variables
    
    def static_fric(self, coeficient, normal):
        return coeficient * float(normal)
    
    def kinetic_friction(self, coeficient, normal):
        return coeficient * float(normal)
    
    def momentum(self, mass, velocity):
        return mass * float(velocity)
    
    def impulse(self, force, initial_time = None, final_time = None, change_in_time = None):
        if((final_time is None and change_in_time is None) or (initial_time is None and change_in_time is None)):
            return None
        if(final_time < 0 or initial_time < 0):
            return None
        if(initial_time is not None and final_time is not None and change_in_time is None):
            change_in_time = float(final_time) - initial_time
        return force * change_in_time
    
    def work(self, time = None, power = None):
        if(power is not None and time is not None):
            return power * float(time)
        else:
            return None
    
    def universal_gravitation(self, mass_one, mass_two, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return (gravitational_constant * mass_one * mass_two)/(float(radius)*radius)
    
    def gravitational_field(self, mass, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return (gravitational_constant * mass)/(radius*float(radius))
        
    def gravitational_potential_energy(self, mass_one, mass_two, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return (gravitational_constant * mass_one * mass_two)/(float(radius))
    
    def gravitational_potential(self, mass, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return (gravitational_constant * mass)/(float(radius))
    
    def orbital_speed(self, mass, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return self.math.sqrt((gravitational_constant * mass)/float(radius))
    
    def escape_speed(self, mass, radius):
        gravitational_constant = 6.67408 * self.math.raise_to_power(10.0, -11)
        return self.math.sqrt((2.0 * gravitational_constant * mass)/float(radius))
    
    def elastic_potential_energy(self, k, distance):
        return 1/2.0*k*distance*distance
    
    def density(self, mass, volume):
        return float(mass)/volume
    
    def pressure(self, force, area_covered):
        return float(force)/area_covered
    
    def pressure_in_fluid(self, initial_pressure, density_of_fluid, gravity, depth):
        return initial_pressure + density_of_fluid*gravity*float(depth)
    
    def buoyancy(self, density, gravity, volume):
        return density*gravity*float(volume)
    
    
def main():
    m = physics_and_astronomy()
    k  = kinematic_object(distance = 3, final_velocity = 1, acceleration = 2)
    k = m.kinematic_equations(k)
if __name__ == "__main__":
    main()






