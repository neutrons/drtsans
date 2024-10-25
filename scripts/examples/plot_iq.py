"""
Example script to plot I(q) before and after corrections.
This script is intended to be run in the same directory as the I(q) files,
for example, after the `test_elastic_and_inelastic_corrections.py` integration test.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

before_file = "EQSANS_125707_no_correction_Iq.dat"
after_k_file = "EQSANS_125707_elastic_correction_Iq.dat"
after_b_file = "EQSANS_125707_inelastic_correction_Iq.dat"
after_both_file = "EQSANS_125707_elastic_inelastic_correction_Iq.dat"

output_dir = Path("plots")
output_dir.mkdir(exist_ok=True)

before_data = np.loadtxt(before_file)
after_k_data = np.loadtxt(after_k_file)
after_b_data = np.loadtxt(after_b_file)
after_both_data = np.loadtxt(after_both_file)
plt.plot(before_data[:, 0], before_data[:, 1], label="Before")
plt.plot(after_k_data[:, 0], after_k_data[:, 1], label="Elastic")
plt.plot(after_b_data[:, 0], after_b_data[:, 1], label="Inelastic")
plt.plot(after_both_data[:, 0], after_both_data[:, 1], label="Both")
plt.legend()
plt.savefig("plots/plot_Iq.png")
plt.close()
