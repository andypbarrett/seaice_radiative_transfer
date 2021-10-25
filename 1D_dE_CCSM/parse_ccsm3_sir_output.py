"""Script to parse output files"""
import pandas as pd


def parse_single_assignment(line):
    """Parses a line of form x = y, and returns y"""
    return float(line.split('=')[1].strip())


def parse_dble_assignment(line):
    """Parses a line of form x = y0 y1"""
    return [float(v.strip()) for v in line.split('=')[-1].split()]


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


def parse_spectral_albedos(lines):
    """Parses direct and diffuse spectral albedos"""
    spectral_albedos = []
    for line in lines:
        spectral_albedos.append(parse_dble_assignment(line))
    df = pd.DataFrame(spectral_albedos, index=[1, 2, 3], columns=['Direct', 'Diffuse'])
    df.index.name = 'Interval'
    return df


def parse_normalized_fluxes(lines):
    """Parses the normalized radiative flux table
    lines[0:16] contain visible fluxes and extinction_coef
    lines[17:] contain near-infrared
    """
    vs_lines = lines[:16]
    header = vs_lines[0].split()
    flux = parse_table(vs_lines[1::2], header=False)
    extinction_coef = [parse_single_assignment(line) for line in vs_lines[2::2]]
    df_flux = pd.DataFrame(flux, columns=header)
    df_flux.set_index('level', inplace=True)
    df_extinction = pd.DataFrame(extinction_coef, columns=['Extinction_coef'])
    df_extinction.index = df_extinction.index + 0.5
    df_extinction.index.name = 'level'

    ni_lines = lines[16:]
    header = ni_lines[0].split()
    flux = parse_table(ni_lines[1:], header=False)
    df_tmp = pd.DataFrame(flux, columns=header)
    df_tmp.set_index('level', inplace=True)

    df_flux = df_flux.join(df_tmp)
    return df_flux, df_extinction


class DeltaEdOutput():
    """Class to hold output from Delta Eddington Model

    
    Spectral albedo intervals are (I think) 0.2-0.7, 0.7-1.19, and
    1.19-5.0 micro-meters
    """

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
        self.spectral_albedos = parse_spectral_albedos(lines[41:44])
        self.albedo_visible_direct = parse_single_assignment(lines[45])
        self.albedo_visible_diffuse = parse_single_assignment(lines[46])
        self.albedo_nir_direct = parse_single_assignment(lines[47])
        self.albedo_nir_diffuse = parse_single_assignment(lines[48])
        flux_table, extinction_coefs = parse_normalized_fluxes(lines[51:76])
        self.normalized_flux_table = flux_table
        self.extinction_coeficients = extinction_coefs
        
        
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
        

    def print_results(self):
        print(f'Cosine Solar Zenith Angle: {self.cosine_solar_zenith_angle}')
        print(f'Albedo visible direct: {self.albedo_visible_direct}')
        print(f'Albedo visible diffuse: {self.albedo_visible_diffuse}')
        print(f'Albedo near-infrared direct: {self.albedo_nir_direct}')
        print(f'Albedo near-infrared diffuse: {self.albedo_nir_diffuse}')
        print('')
        print('Spectral Albedo for intervals')
        print(self.spectral_albedos)
        print('')
        print('Normalized Fluxes')
        print(self.normalized_flux_table)
        print(self.extinction_coeficients)
        
        
def main(filename):

    results = DeltaEdOutput(filename)

    results.print_results()
    
    
if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
