"""
Example script to plot the before and after k-correction data for a given slice and frame.
"""

import glob
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

data_dir = Path("EQSANS_125707_elastic_correction") / "slice_0" / "frame_0"
before_files = sorted(glob.glob(os.path.join(data_dir, "IQ_*_before_k_correction.dat")))
after_files = sorted(glob.glob(os.path.join(data_dir, "IQ_*_after_k_correction.dat")))

output_dir = Path("plots")
output_dir.mkdir(exist_ok=True)

for before_file, after_file in zip(before_files, after_files):
    wavelength = Path(before_file).name.split("_")[1].split("_")[0]
    before_data = np.loadtxt(before_file)
    after_data = np.loadtxt(after_file)
    plt.plot(before_data[:, 0], before_data[:, 1], label="Before")
    plt.plot(after_data[:, 0], after_data[:, 1], label="After")
    plt.legend()
    plt.savefig(f"plots/plot_{wavelength}.png")
    plt.close()
