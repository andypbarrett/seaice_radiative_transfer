"""Script to parse output files"""
import pandas as pd


def parse_single_assignment(line):
    """Parses a line of form x = y, and returns y"""
    return float(line.split('=')[1].strip())


def parse_dble_assignment(line):
    """Parses a line of form x = y0 y1"""
    pass


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

    
def parse_atmospheric_profiles(lines):
    """Parses atmospheric profile input to pandas DataFrame"""
    data = parse_table(lines, header=True)
    df = pd.DataFrame(data)
    df.columns = ['level',  'pressure(mb)',  'temperature(k)',
                  'h2ommr(g/g)', 'o3mmr(g/g)', 'cld_cover',  'cld_lwp(g/m2)']
    df.set_index('level', inplace=True)
    return df


class DeltaEdOutput():
    """Class to hold output from Delta Eddington Model"""

    def __init__(self, output_file):

        with open(filename, 'r') as f:
            lines = f.readlines()

        # Input parameters
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
        self.atmosphere_profile = parse_atmospheric_profiles(lines[7:25])

        # Results
        self.cosine_solar_zenith_angle = parse_single_assignment(lines[40])
        
        
    def print_inputs(self):
        print(f'Day of Year: {self.day_of_year}')
        print(f'Latitude: {self.latitude}')
        print(f'Snow depth (m): {self.snow_depth}')
        print(f'Snow density (kg/m3): {self.snow_density}')
        print(f'Pond depth (m): {self.pond_depth}')
        print(f'Pond tuning parameter: {self.pond_tuning_parameter}')
        print(f'Sea ice thickness (m): {self.sea_ice_thickness}')
        print(f'Sea ice tuning parameter: {self.sea_ice_tuning_parameter}')
        print('')
        print('Atmospheric Profile')
        print(self.atmosphere_profile)
        

def main(filename):

    results = DeltaEdOutput(filename)

    results.print_inputs()
    
    
if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
