from scramjet import scramjet

class parametrize_scramjet():
    def __init__(self, alpha_array, mach_array, scale_factors_nominal, input_vals_nominal):
        self.alpha_array = alpha_array
        self.mach_array = mach_array
        self.scale_factors = scale_factors_nominal
        self.input_vals = input_vals_nominal
    
    def parametrize_alpha(self):
        data_array = []

        for i in range(len(self.alpha_array)):
            self.input_vals['alpha'] = self.alpha_array[i]
            data = scramjet(self.scale_factors, self.input_vals)
            data.run()
            data_array.append(data)

        return data_array

    def parametrize_mach(self):
        data_array = []

        for i in range(len(self.mach_array)):
            self.input_vals['M_freestream'] = self.mach_array[i]
            data = scramjet(self.scale_factors, self.input_vals)
            data.run()
            data_array.append(data)

        return data_array