"""Script to parse output files"""
import pandas as pd


def parse_single_assignment(line):
    """Parses a line of form x = y, and returns y"""
    return float(line.split('=')[1].strip())


def parse_table(lines, header=True):
    """Parses a table into a list of lists

    Each list contains a row of the table
    """
    if header:
        these_lines = lines[1:]
    else:
        these_lines = lines
    data = [line.split() for line in these_lines]
    return data

    
class DeltaEdOutput():
    """Class to hold output from Delta Eddington Model"""

    def __init__(self, output_file):

        with open(filename, 'r') as f:
            lines = f.readlines()

        self.day_of_year = parse_single_assignment(lines[5])
        self.latitude = parse_single_assignment(lines[6])
        self.surface_pressure = parse_single_assignment(lines[26])
        self.atmospheric_co2_vmr = parse_single_assignment(lines[27])
        self.surface_air_temperature = parse_single_assignment(lines[28])
        self.surface_skin_temperature = parse_single_assignment(lines[29])
        self.snow_depth = parse_single_assignment(lines[30])
        self.snow_density = parse_single_assignment(lines[31])
        self.snow_grain_radius = parse_single_assignment(lines[32])
        self.pond_depth = parse_single_assignment(lines[33])
        self.pond_tuning_parameter = parse_single_assignment(lines[34])
        self.sea_ice_thickness = parse_single_assignment(lines[35])
        self.sea_ice_tuning_parameter = parse_single_assignment(lines[36])
        self.atmosphere_profile = parse_table(lines[7:25])


    def print_inputs(self):
        print(f'Day of Year: {self.day_of_year}')
        print(f'Latitude: {self.latitude}')
        print(f'Surface Pressure: {self.surface_pressure}')
        print('')
        print(self.atmosphere_profile)
        

def main(filename):

    results = DeltaEdOutput(filename)

    results.print_inputs()
    
    
if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
