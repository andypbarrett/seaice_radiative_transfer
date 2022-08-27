"""Makes snow depth and density with time and geolocation for Polarstern MOSAiC drift track"""

from pathlib import Path

import pandas as pd

DATAPATH = Path(__file__).resolve().parent / ".." / "data"
SNOWPATH = DATAPATH / "polarstern_mosaic_magnaprobe.csv"
DRIFTPATH = DATAPATH / "PSTrack20200817.csv"

# Load snow depth and density data
snow_df = pd.read_csv(SNOWPATH, index_col=0, header=0, parse_dates=True)
snow_df.index = snow_df.index + pd.Timedelta(12, 'h')  # make snow data times midday
snow_df[snow_df["Depth"] < 0.005] = 0.0  # Set depth and density < 0.5 cm to 0.0 cm

# Load track data and extract just midday position - might need to be changed
drift_df = pd.read_csv(DRIFTPATH, index_col=2, header=0, parse_dates=True)
drift_df = drift_df[(drift_df.index.hour == 12) & (drift_df.index.minute == 0)]

joined_df = snow_df.join(drift_df)
joined_df = joined_df.dropna()
joined_df.to_csv(DATAPATH / "rt_snow_input.csv")