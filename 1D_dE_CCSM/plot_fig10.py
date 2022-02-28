"""Generates Figure 10 from Briegleib and Light page 90"""
from pathlib import Path

import matplotlib.pyplot as plt

from parse_ccsm3_sir_output import DeltaEdOutput


spectral_bands = [(300., 700.), (700., 1200.), (1200., 2100.)]  # Need to fix

OUTPUT_DIRPATH = Path.home() / 'src' / 'seaice_radiative_transfer' / '1D_dE_CCSM' / '1D_DE_test'
FILE_PREFIX = 'ccsm3_sir_de_output'
FIGURE = '10'

EXPERIMENTS = ['bare1sd',
               'bare2sd',
               'bareminus1sd',
               'bareminus2sd',
               'pond1sd',
               'pondminus1sd',
               ]


linecolor = ['k', 'k', 'k', 'k', 'b', 'b']
linestyle = ['--', '-.', '--', '-.', '--', '--']

def main():

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(300., 2100.)
    ax.set_ylim(0., 1.)
    ax.set_ylabel('Albedo')
    ax.set_xlabel('Wavelength (nm)')

    bare = DeltaEdOutput(OUTPUT_DIRPATH /
                               f'{FILE_PREFIX}_fig{FIGURE}bare.dat')
    salbedo = bare.spectral_albedos
    salbedo.index = [300., 700., 1190.]
    salbedo['Total'].plot(color='k', drawstyle='steps-post', linewidth=2,
                           label='bare')

    pond = DeltaEdOutput(OUTPUT_DIRPATH /
                         f'{FILE_PREFIX}_fig{FIGURE}pond.dat')
    salbedo = pond.spectral_albedos
    salbedo.index = [300., 700., 1200.]
    salbedo['Total'].plot(color='b', drawstyle='steps-post', linewidth=2,
                           label='pond')
    
    for experiment, c, lst in zip(EXPERIMENTS, linecolor, linestyle):
        output = DeltaEdOutput(OUTPUT_DIRPATH /
                               f'{FILE_PREFIX}_fig{FIGURE}{experiment}.dat')
        salbedo = output.spectral_albedos
        #salbedo.index = [300., 700., 1200.]
        #salbedo['Direct'].plot(drawstyle='steps-post', linewidth=2,
        #                       label=experiment)
        for index, row in salbedo.iterrows():
            ax.plot(spectral_bands[index-1], [row['Total'], row['Total']],
                    color=c, linestyle=lst)

    ax.legend()
    plt.show()


if __name__ == "__main__":
    main()
