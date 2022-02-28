"""Script to parse output files"""
import pandas as pd


def parse_single_assignment(line, with_units=False):
    """Parses a line of form x = y, and returns y"""
    result = line.split('=')[1].strip()
    if with_units:
        result = result.split()[0]
    return float(result)


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


def parse_cloud_particle_table(lines):
    """Parses table of cloud particle radii and ice fraction"""
    header = lines[0].split()
    data = parse_table(lines[1:], header=False)
    df = pd.DataFrame(data, columns=header)
    df.set_index('level', inplace=True)
    return df


def parse_heating_rates(lines):
    """Parses heating rate table"""
    header = lines[0].split()
    data = parse_table(lines[1:], header=False)
    df = pd.DataFrame(data, columns=header)
    df.set_index('level', inplace=True)
    return df


def parse_absorption(lines):
    #header = lines[0].split()
    header = ['level_type', 'level_id', 'Q_vs', 'Q_ni', 'Q_total']
    absorption = parse_table(lines[1::2], header=False)
    absorption[0].insert(1, None)
    absorption[-1].insert(1, None)
    return pd.DataFrame(absorption, columns=header)


def parse_transmittance(lines):
    header = ['depth', 'Tr_vs', 'Tr_ni']
    transmittance = parse_table(lines[2::2], header=False)
    return pd.DataFrame(transmittance, columns=header)


def parse_transmittance_absorption(lines):
    absorption_df = parse_absorption(lines)
    transmittance_df = parse_transmittance(lines)
    return absorption_df, transmittance_df


class DeltaEdOutput():
    """Class to hold output from Delta Eddington Model

    
    Spectral albedo intervals are (I think) 0.2-0.7, 0.7-1.19, and
    1.19-5.0 micro-meters
    """

    def __init__(self, output_file):

        with open(output_file, 'r') as f:
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
        
        # Atmosphere/Surface Radiation
        self.earth_sun_distance_factor = parse_single_assignment(lines[81])
        self.local_solar_time = parse_single_assignment(lines[85])
        self.solar_toa_insolation = parse_single_assignment(lines[86], with_units=True)
        self.solar_toa_albedo = parse_single_assignment(lines[87])
        self.solar_toa_absorbed = parse_single_assignment(lines[88], with_units=True)
        self.solar_absorbed_atmosphere = parse_single_assignment(lines[89], with_units=True)
        self.solar_absorbed_surface = parse_single_assignment(lines[90], with_units=True)
        self.solar_clear_toa_absorbed = parse_single_assignment(lines[91], with_units=True)
        self.solar_clear_atm_absorbed = parse_single_assignment(lines[92], with_units=True)
        self.solar_clear_srf_absorbed = parse_single_assignment(lines[93], with_units=True)
        self.solar_cloud_forcing = parse_single_assignment(lines[94], with_units=True)
        self.solar_clear_sky_albedo = parse_single_assignment(lines[95], with_units=True)
        self.longwave_net_toa_up = parse_single_assignment(lines[97], with_units=True)
        self.longwave_net_surface_up = parse_single_assignment(lines[98], with_units=True)
        self.longwave_surface_down = parse_single_assignment(lines[99], with_units=True)
        self.longwave_clear_outgoing = parse_single_assignment(lines[100], with_units=True)
        self.longwave_clear_net_srf = parse_single_assignment(lines[101], with_units=True)
        self.longwave_cloud_forcing = parse_single_assignment(lines[102], with_units=True)
        self.net_cloud_forcing = parse_single_assignment(lines[104], with_units=True)

        # Atmospheric tables
        self.cloud_particle = parse_cloud_particle_table(lines[107:126])
        self.heating_rates = parse_heating_rates(lines[128:147])

        self.solar_vs_direct_surface_irrad = parse_single_assignment(lines[149], with_units=True)
        self.solar_vs_diff_surface_irrad = parse_single_assignment(lines[150], with_units=True)
        self.vs_fraction_of_direct_irrad = parse_single_assignment(lines[151])
        self.vs_fraction_of_diff_irrad = parse_single_assignment(lines[152])
        self.solar_ni_direct_surface_irrad = parse_single_assignment(lines[153], with_units=True)
        self.solar_ni_diff_surface_irrad = parse_single_assignment(lines[154], with_units=True)
        self.ni_fraction_of_direct_irrad = parse_single_assignment(lines[155])
        self.ni_fraction_of_diff_irrad = parse_single_assignment(lines[156])
        self.total_solar_surface_irrad = parse_single_assignment(lines[157], with_units=True)
        self.vs_fraction_total_irrad = parse_single_assignment(lines[158])
        self.ni_fraction_total_irrad = parse_single_assignment(lines[159])
        self.solar_absorbed_at_surface = parse_single_assignment(lines[160], with_units=True) 
        self.broad_band_albedo = parse_single_assignment(lines[161])
        self.solar_vs_absorbed_seaice = parse_single_assignment(lines[162], with_units=True)
        self.solar_ni_absorbed_seaice = parse_single_assignment(lines[163], with_units=True)
        self.solar_total_absorbed_surface = parse_single_assignment(lines[164], with_units=True)
        self.frac_vs_abs_pentrt_srf = parse_single_assignment(lines[165])
        self.frac_ni_abs_pentrt_srf = parse_single_assignment(lines[166])
        
        absorption_df, transmit_df = parse_transmittance_absorption(lines[169:187])
        self.absorption = absorption_df
        self.transmittance = transmit_df

        # Estimate all sky spectral albedos
        self.estimate_total_spectral_albedo()


    def estimate_total_spectral_albedo(self):
        fraction = pd.DataFrame(
            {
                'Direct': [self.vs_fraction_of_direct_irrad,
                           self.ni_fraction_of_direct_irrad,
                           self.ni_fraction_of_direct_irrad],
                'Diffuse': [self.vs_fraction_of_diff_irrad,
                            self.ni_fraction_of_diff_irrad,
                            self.ni_fraction_of_diff_irrad],
            },
            index=self.spectral_albedos.index
            )
        self.spectral_albedos['Total'] = (self.spectral_albedos * fraction).sum(axis=1)


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
        print('Atmosphere/Surface Radiation')
        print(self.earth_sun_distance_factor)
        print(self.local_solar_time)
        print(self.solar_toa_insolation)
        print('')
        print('Cloud particle radii and ice fraction')
        print(self.cloud_particle)
        print('')
        print('Snow/Sea Ice Absorption')
        print(self.absorption)
        print('')
        print('Snow/Sea Ice Transmittance')
        print(self.transmittance)
        
        
def main(filename):

    results = DeltaEdOutput(filename)

    results.print_results()
    
    
if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
