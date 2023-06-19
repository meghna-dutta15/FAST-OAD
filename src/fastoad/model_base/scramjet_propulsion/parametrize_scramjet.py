from scramjet import scramjet

class parametrize_scramjet():
    def __init__(self, scale_factors_nominal, input_vals_nominal, alpha_array=None, mach_array=None, x_scale_array = None, y_scale_array=None):
        self.alpha_array = alpha_array
        self.mach_array = mach_array
        self.scale_factors = scale_factors_nominal
        self.input_vals = input_vals_nominal
        self.x_scale_array = x_scale_array
        self.y_scale_array = y_scale_array
    
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
    

    def parametrize_scale_factors(self):
        data_array = []

        for i in range(len(self.x_scale_array)):
            self.scale_factors['x_scale'] = self.x_scale_array[i]
            temp = []
            for i in range(len(self.y_scale_array)):
                self.scale_factors['y_scale'] = self.y_scale_array[i]
            
                data = scramjet(self.scale_factors, self.input_vals)
                data.run()
                temp.append(data)


            data_array.append(temp)

        return data_array


    def parametrize_x(self):
        data_array = []

        for i in range(len(self.x_scale_array)):
            self.scale_factors['x_scale'] = self.x_scale_array[i]
            data = scramjet(self.scale_factors, self.input_vals)
            data.run()
            data_array.append(data)

        return data_array
    
    def parametrize_y(self):
        data_array = []

        for i in range(len(self.y_scale_array)):
            self.scale_factors['y_scale'] = self.y_scale_array[i]
            data = scramjet(self.scale_factors, self.input_vals)
            data.run()
            data_array.append(data)

        return data_array
    
