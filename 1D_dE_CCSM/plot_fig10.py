"""Generates Figure 10 from Briegleib and Light page 90"""
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt

from parse_ccsm3_sir_output import DeltaEdOutput
from filepath import TEST_OUTPUT_DIRPATH

spectral_bands = [(300., 700.), (700., 1190.), (1190., 2100.)]  # Need to fix

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

def plot_albedo_as_step(df, ax=None, label=None, color='k'):
    """Helper function to make nice step plot"""
    if not ax:
        ax = plt.gca()

    x = [300., 700., 1190., 5000.]
    y = np.append(df['Total'].values, df['Total'][3])
    ax.step(x, y, where='post', color=color,
            linewidth=2, label=label)


def main():

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(300., 2100.)
    ax.set_ylim(0., 1.)

    ax.set_yticks(np.arange(0., 1.05, .1))
    ax.set_yticks(np.arange(0.05, 1.05, .1), minor=True)

    ax.set_xticks([500., 1000., 1500., 2000.])
    ax.set_xticks(np.arange(300., 2150., 100.), minor=True)

    ax.set_ylabel('Albedo', fontsize=20)
    ax.set_xlabel('Wavelength (nm)', fontsize=20)

    ax.tick_params(axis='both', which='major', labelsize=15, length=10)
    ax.tick_params(which='minor', length=7)

    bare = DeltaEdOutput(TEST_OUTPUT_DIRPATH /
                         f'{FILE_PREFIX}_fig{FIGURE}bare.dat')
    plot_albedo_as_step(bare.spectral_albedos, ax=ax, color='k', label='Bare')

    pond = DeltaEdOutput(TEST_OUTPUT_DIRPATH /
                         f'{FILE_PREFIX}_fig{FIGURE}pond.dat')
    plot_albedo_as_step(pond.spectral_albedos, ax=ax, color='b', label='Pond')

    for experiment, c, lst in zip(EXPERIMENTS, linecolor, linestyle):
        output = DeltaEdOutput(TEST_OUTPUT_DIRPATH /
                               f'{FILE_PREFIX}_fig{FIGURE}{experiment}.dat')
        salbedo = output.spectral_albedos
        for index, row in salbedo.iterrows():
            ax.plot(spectral_bands[index-1], [row['Total'], row['Total']],
                    color=c, linestyle=lst)

    ax.legend(fontsize=20)
    fig.savefig(TEST_OUTPUT_DIRPATH / f"{FILE_PREFIX}_fig{FIGURE}.png")
    plt.show()


if __name__ == "__main__":
    main()
