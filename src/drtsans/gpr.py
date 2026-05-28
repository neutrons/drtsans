"""
Gaussian Process Regression (GPR) analysis for EQSANS I(Q) profiles.

This module provides GPR fitting and uncertainty quantification for I(Q) data,
generating both interactive Plotly plots and PNG images for autoreduction reports.

Original algorithm developed by Changwoo Do (ORNL).
Adapted for drtsans integration: May 2026.

Command-Line Interface
----------------------
The module can be run as a standalone CLI tool:

    python -m drtsans.gpr input_file.dat [options]

Shell Tab Completion (argcomplete):
    To enable bash/zsh tab completion, add to ~/.bashrc:
        eval "$(register-python-argcomplete drtsans-gpr)"

Interactive TUI (argparse-tui):
    Launch an interactive terminal UI with:
        python -m drtsans.gpr --tui

References
----------
- RBF kernel-based GPR with log-space transformations
- Iterative background estimation via kernel smoothing
- Noise scale optimization via maximum likelihood
"""

# standard imports
import argparse
import os
import sys
from typing import Callable, Dict, List, Optional, Tuple, Union

# third party imports
import matplotlib

matplotlib.use("Agg")  # Use non-GUI backend for headless environments
import matplotlib.pyplot as plt
import numpy as np
from mantid.simpleapi import logger

# local imports
from drtsans.dataobjects import IQmod

__all__ = [
    "generate_gpr_analysis",
    "run_gpr_from_file",
    "plotly_gpr_plot",
]


# ============================================================================
# Core GPR Algorithm Functions (from original eqsans_gpr.py)
# ============================================================================


def f_loglin(x: np.ndarray, x_c: float = 1) -> np.ndarray:
    """Log-linear transformation function."""
    return np.where(x < x_c, (x - x_c) / x_c + np.log(x_c), np.log(x))


def f_loglin_inv(y: np.ndarray, y_c: float = 0) -> np.ndarray:
    """Inverse log-linear transformation."""
    return np.where(y < y_c, (1 + y - y_c) * np.exp(y_c), np.exp(y))


def f_loglin_deriv(x: np.ndarray, x_c: float = 1) -> np.ndarray:
    """Derivative of log-linear transformation."""
    return np.where(x < x_c, 1 / x_c, 1 / x)


def f_log_likelihood(
    y: np.ndarray,
    y_err: np.ndarray,
    x: np.ndarray,
    x_err: np.ndarray,
    z_list: np.ndarray = np.linspace(0, 0, 1),
    importance: Optional[np.ndarray] = None,
    index: Optional[np.ndarray] = None,
) -> float:
    """
    Calculate log-likelihood of y in the distribution of x.

    Parameters
    ----------
    y : np.ndarray
        Observed data values
    y_err : np.ndarray
        Uncertainties in y
    x : np.ndarray
        Model predictions
    x_err : np.ndarray
        Uncertainties in model predictions
    z_list : np.ndarray, optional
        Integration points for marginalizing over observation uncertainty
    importance : np.ndarray, optional
        Importance weights for each data point
    index : np.ndarray, optional
        Boolean mask for which data points to include

    Returns
    -------
    float
        Log-likelihood value
    """
    log_likelihood = 0
    sum_weight = 0

    if importance is None:
        importance = np.ones(len(y))
    if index is None:
        index = np.arange(len(y))

    for z in z_list:
        weight_z = np.exp(-(z**2) / 2) / np.sqrt(2 * np.pi)
        y_z = y + y_err * z
        log_likelihood_z = -0.5 * np.sum(
            (((y_z - x) ** 2 / x_err**2 + np.log(2 * np.pi * x_err**2)) * importance)[index]
        )
        log_likelihood += log_likelihood_z * weight_z
        sum_weight += weight_z

    log_likelihood = log_likelihood / sum_weight
    return log_likelihood


def gpr_core(
    q: np.ndarray,
    I: np.ndarray,
    I_err: np.ndarray,
    I_bg: np.ndarray,
    sigma_scale: float,
    f_I: Callable[[np.ndarray], np.ndarray],
    f_I_deriv: Callable[[np.ndarray], np.ndarray],
    f_Q: Callable[[np.ndarray], np.ndarray],
    lmbda: float,
    index_eval: Optional[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Core GPR computation: posterior mean and standard deviation.

    Parameters
    ----------
    q : np.ndarray
        Q values (momentum transfer)
    I : np.ndarray
        Intensity values
    I_err : np.ndarray
        Intensity uncertainties
    I_bg : np.ndarray
        Background intensity estimate
    sigma_scale : float
        Noise scale parameter
    f_I : callable
        Intensity transformation function
    f_I_deriv : callable
        Derivative of intensity transformation
    f_Q : callable
        Q transformation function
    lmbda : float
        RBF kernel length scale
    index_eval : np.ndarray
        Boolean mask for evaluation points

    Returns
    -------
    tuple
        (mu_s, std_s, y, y_err, sig_var) - posterior mean, std, transformed data, error, variance
    """
    y_data = f_I(I)
    y_bg = f_I(I_bg)
    y = y_data - y_bg
    y_err = I_err * f_I_deriv(I)

    if index_eval is None:
        index_eval = np.ones(len(y), dtype=bool)

    sig_var = np.mean(y[index_eval] ** 2)
    prior_sigma = np.sqrt(sig_var)

    def f_K_rbf(q1, q2, lmbda, std_1=None, std_2=None):
        if std_1 is None:
            std_1 = np.ones_like(q1)
        if std_2 is None:
            std_2 = np.ones_like(q2)
        return np.exp(-((q1[:, None] - q2[None, :]) ** 2) / (2 * lmbda**2)) * std_1[:, None] * std_2[None, :]

    def f_y_err(err):
        return sigma_scale**2 * (err**2 / sig_var)

    K_sigma = np.diag(f_y_err(y_err))

    q_tr = f_Q(q)
    K_xx = f_K_rbf(q_tr, q_tr, lmbda) + K_sigma
    K_xs = f_K_rbf(q_tr, q_tr, lmbda)
    k_ss = f_K_rbf(q_tr, q_tr, lmbda)

    L = np.linalg.cholesky(K_xx)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y))
    v = np.linalg.solve(L, K_xs)

    mu_s = K_xs.T @ alpha
    cov_s = k_ss - v.T @ v
    var_s = np.diag(cov_s) * prior_sigma**2
    std_s = np.sqrt(var_s)

    return mu_s, std_s, y, y_err, sig_var


def optimize_m_factor_gpr(
    q: np.ndarray,
    I: np.ndarray,
    I_err: np.ndarray,
    I_bg: np.ndarray,
    f_I: Callable[[np.ndarray], np.ndarray],
    f_inv_I: Callable[[np.ndarray], np.ndarray],
    f_I_deriv: Callable[[np.ndarray], np.ndarray],
    f_Q: Callable[[np.ndarray], np.ndarray],
    lmbda: float,
    index_eval: Optional[np.ndarray],
) -> float:
    """
    Optimize GPR noise scale parameter via maximum likelihood.

    Parameters
    ----------
    q : np.ndarray
        Q values
    I : np.ndarray
        Intensity values
    I_err : np.ndarray
        Intensity uncertainties
    I_bg : np.ndarray
        Background estimate
    f_I : callable
        Intensity transformation
    f_inv_I : callable
        Inverse intensity transformation
    f_I_deriv : callable
        Derivative of intensity transformation
    f_Q : callable
        Q transformation
    lmbda : float
        RBF kernel length scale
    index_eval : np.ndarray
        Boolean mask for evaluation

    Returns
    -------
    float
        Optimal noise scale parameter
    """
    sigma_list = np.logspace(-1, 3, 25)
    log_likelihood_list = []

    for sigma in sigma_list:
        mu_s, std_s, y, y_err, _ = gpr_core(
            q=q,
            I=I,
            I_err=I_err,
            I_bg=I_bg,
            sigma_scale=sigma,
            f_I=f_I,
            f_I_deriv=f_I_deriv,
            f_Q=f_Q,
            lmbda=lmbda,
            index_eval=index_eval,
        )
        log_likelihood = f_log_likelihood(y, y_err, mu_s, std_s, index=index_eval)
        log_likelihood_list.append(log_likelihood)

    log_likelihood_list = np.array(log_likelihood_list)
    log_likelihood_list -= log_likelihood_list[-1]
    log_likelihood_list /= len(q)

    sigma_best = sigma_list[np.argmax(log_likelihood_list)]
    return sigma_best


def gpr_posterior_predictive(
    q: np.ndarray,
    I: np.ndarray,
    I_err: np.ndarray,
    I_bg: np.ndarray,
    sigma_best: float,
    f_I: Callable[[np.ndarray], np.ndarray],
    f_inv_I: Callable[[np.ndarray], np.ndarray],
    f_I_deriv: Callable[[np.ndarray], np.ndarray],
    f_Q: Callable[[np.ndarray], np.ndarray],
    lmbda: float,
    index_eval: Optional[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute GPR posterior predictive distribution.

    Parameters
    ----------
    q : np.ndarray
        Q values
    I : np.ndarray
        Intensity values
    I_err : np.ndarray
        Intensity uncertainties
    I_bg : np.ndarray
        Background estimate
    sigma_best : float
        Optimized noise scale parameter
    f_I : callable
        Intensity transformation
    f_inv_I : callable
        Inverse intensity transformation
    f_I_deriv : callable
        Derivative of intensity transformation
    f_Q : callable
        Q transformation
    lmbda : float
        RBF kernel length scale
    index_eval : np.ndarray
        Boolean mask for evaluation

    Returns
    -------
    tuple
        (I_mean_GP, I_std_GP) - posterior mean and uncertainty in original scale
    """
    mu_s, std_s, _, _, _ = gpr_core(
        q=q,
        I=I,
        I_err=I_err,
        I_bg=I_bg,
        sigma_scale=sigma_best,
        f_I=f_I,
        f_I_deriv=f_I_deriv,
        f_Q=f_Q,
        lmbda=lmbda,
        index_eval=index_eval,
    )

    y_bg = f_I(I_bg)
    I_mean_GP = f_inv_I(mu_s + y_bg)

    # Convert uncertainty from transformed space back to intensity space
    # For identity transform, f_I_deriv(...) == 1 so this is a no-op
    I_std_GP = std_s / np.abs(f_I_deriv(I_mean_GP))

    return I_mean_GP, I_std_GP


# ============================================================================
# Main GPR Analysis Function
# ============================================================================


def run_gpr(
    q: np.ndarray,
    I: np.ndarray,
    I_err: np.ndarray,
    dq: Optional[np.ndarray] = None,
    use_log_I: bool = False,
    use_log_Q: bool = True,
    lmbda: float = 0.25,
    background_filter_width: float = 0.2,
    background_intensity_offset: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    Run GPR analysis on I(Q) data.

    Parameters
    ----------
    q : np.ndarray
        Q values (momentum transfer in 1/Å)
    I : np.ndarray
        Intensity values (1/cm)
    I_err : np.ndarray
        Intensity uncertainties (1/cm)
    dq : np.ndarray, optional
        Q uncertainties (1/Å)
    use_log_I : bool, optional
        Apply log-linear transformation to intensity (default: False)
    use_log_Q : bool, optional
        Apply log transformation to Q (default: True)
    lmbda : float, optional
        RBF kernel length scale (default: 0.25)
    background_filter_width : float, optional
        Kernel width for background estimation (default: 0.2)
    background_intensity_offset : float, optional
        Background level for log transform (default: 1.0)

    Returns
    -------
    dict
        Dictionary containing:
        - 'q': Q values
        - 'I_gpr': GPR mean prediction
        - 'I_gpr_err': GPR uncertainty
        - 'I_exp': Experimental data
        - 'I_err_exp': Experimental uncertainties
        - 'dq': Q uncertainties (if provided)
        - 'sigma_best': Optimized noise scale parameter

    Raises
    ------
    ValueError
        If input arrays have incompatible shapes or contain invalid values
    """
    # Input validation
    if len(q) != len(I) or len(q) != len(I_err):
        raise ValueError("Input arrays q, I, and I_err must have the same length")

    if dq is not None and len(dq) != len(q):
        raise ValueError("dq array must have the same length as q")

    if np.any(q <= 0):
        raise ValueError("Q values must be positive")

    if lmbda <= 0:
        raise ValueError("lmbda must be positive")

    if background_filter_width <= 0:
        raise ValueError("background_filter_width must be positive")

    if np.any(I_err <= 0):
        logger.warning("Found zero or negative uncertainties, replacing with minimum positive value")
        I_err = I_err.copy()
        min_err = np.min(I_err[I_err > 0]) if np.any(I_err > 0) else 1.0
        I_err[I_err <= 0] = min_err

    # Define Q transformation functions
    if use_log_Q:

        def f_Q(q):
            return np.log(q)

    else:

        def f_Q(q):
            return q

    # Define I(Q) transformation functions
    if use_log_I:
        x_c = background_intensity_offset * 2
        y_c = np.log(x_c)

        def f_I(I):
            return f_loglin(I, x_c=x_c)

        def f_inv_I(fI):
            return f_loglin_inv(fI, y_c=y_c)

        def f_I_deriv(I):
            return f_loglin_deriv(I, x_c=x_c)

    else:

        def f_I(I):
            return I

        def f_inv_I(fI):
            return fI

        def f_I_deriv(I):
            return np.ones_like(I)

    # Define outlier mask (exclude first and last 2 points)
    index_outlier_all = np.zeros(len(q), dtype=bool)
    if len(q) > 4:
        index_outlier_all[:2] = True
        index_outlier_all[-2:] = True
    index_outlier_sm = index_outlier_all
    index_outlier = ~index_outlier_all

    # Estimate background using iterative kernel smoothing
    I_q_gf = np.zeros_like(q)
    for _ in range(5):
        I_q_gf_i = np.zeros_like(q)
        I_iter = I - I_q_gf
        for i in range(len(q)):
            dq_kern = (f_Q(q[i]) - f_Q(q)) / background_filter_width
            weights = np.exp(-0.5 * dq_kern**2)
            weights[index_outlier_sm] = 0
            if np.sum(weights) > 0:
                weights /= np.sum(weights)
                I_q_gf_i[i] = np.sum(I_iter * weights)
        I_q_gf += I_q_gf_i

    # Optimize GPR noise scale parameter
    sigma_best = optimize_m_factor_gpr(
        q=q,
        I=I,
        I_err=I_err,
        I_bg=I_q_gf,
        f_I=f_I,
        f_inv_I=f_inv_I,
        f_I_deriv=f_I_deriv,
        f_Q=f_Q,
        lmbda=lmbda,
        index_eval=index_outlier,
    )

    # GPR prediction: posterior mean and uncertainty
    I_mean_GP, I_std_GP = gpr_posterior_predictive(
        q=q,
        I=I,
        I_err=I_err,
        I_bg=I_q_gf,
        sigma_best=sigma_best,
        f_I=f_I,
        f_inv_I=f_inv_I,
        f_I_deriv=f_I_deriv,
        f_Q=f_Q,
        lmbda=lmbda,
        index_eval=index_outlier,
    )

    result = {
        "q": q,
        "I_gpr": I_mean_GP,
        "I_gpr_err": I_std_GP,
        "I_exp": I,
        "I_err_exp": I_err,
        "sigma_best": sigma_best,
    }

    if dq is not None:
        result["dq"] = dq

    return result


# ============================================================================
# Plotting Functions
# ============================================================================


def create_png_plot(gpr_result: Dict[str, np.ndarray], output_path: str, title: str = "") -> str:
    """
    Create PNG plot of GPR results using matplotlib.

    Parameters
    ----------
    gpr_result : dict
        Dictionary from run_gpr() containing q, I_gpr, I_gpr_err, I_exp, I_err_exp
    output_path : str
        Path to save PNG file
    title : str, optional
        Plot title

    Returns
    -------
    str
        Path to saved PNG file
    """
    q = gpr_result["q"]
    I_gpr = gpr_result["I_gpr"]
    I_gpr_err = gpr_result["I_gpr_err"]
    I_exp = gpr_result["I_exp"]
    I_err_exp = gpr_result["I_err_exp"]

    fig, ax = plt.subplots(figsize=(5, 5))

    # Experimental data
    ax.errorbar(q, I_exp, yerr=I_err_exp, fmt="o", color="k", ms=6, label=r"$I_\mathrm{Exp}$", alpha=0.3)

    # GPR prediction + uncertainty band
    ax.plot(q, I_gpr, "-r", label=r"$I_\mathrm{GPR}$")
    ax.fill_between(q, I_gpr - I_gpr_err, I_gpr + I_gpr_err, color="r", alpha=0.3)

    # Axis and formatting
    ax.set_xlabel(r"$Q$ (1/Å)", fontsize=20)
    ax.set_ylabel(r"$I(Q)$ (1/cm)", fontsize=20)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(direction="in", labelsize=16)

    # Square aspect in log-log scale
    xlim_log = np.log10(ax.get_xlim())
    ylim_log = np.log10(ax.get_ylim())
    ax.set_aspect((xlim_log[1] - xlim_log[0]) / (ylim_log[1] - ylim_log[0]))

    # Legend and layout
    ax.legend(frameon=False, fontsize=16, loc="lower center", bbox_to_anchor=(0.75, 0.6))

    if title:
        ax.set_title(title, fontsize=18)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    return output_path


def plotly_gpr_plot(
    gpr_result: Dict[str, np.ndarray], title: str = "", labels: Optional[Dict[str, str]] = None
) -> str:
    """
    Generate interactive Plotly HTML plot of GPR results.

    Parameters
    ----------
    gpr_result : dict
        Dictionary from run_gpr() containing q, I_gpr, I_gpr_err, I_exp, I_err_exp
    title : str, optional
        Plot title
    labels : dict, optional
        Dictionary with keys 'exp' and 'gpr' for custom labels

    Returns
    -------
    str
        HTML div string containing the Plotly figure

    Notes
    -----
    Uses plot_publisher.plot1d for consistency with existing drtsans plots.
    """
    try:
        from plot_publisher import plot1d
    except ImportError:
        logger.error("plot_publisher not available, cannot create Plotly plot")
        return ""

    q = gpr_result["q"]
    I_gpr = gpr_result["I_gpr"]
    I_gpr_err = gpr_result["I_gpr_err"]
    I_exp = gpr_result["I_exp"]
    I_err_exp = gpr_result["I_err_exp"]

    if labels is None:
        labels = {"exp": "I_Exp", "gpr": "I_GPR"}

    # Prepare data for plot_publisher
    # Format: [[x, y, yerr], [x, y, yerr], ...]
    data_list = [
        [q, I_exp, I_err_exp],  # Experimental data
        [q, I_gpr, I_gpr_err],  # GPR prediction
    ]

    data_names = [labels["exp"], labels["gpr"]]

    # Create plot with log-log scale
    html = plot1d(
        run_number=None,
        data_list=data_list,
        data_names=data_names,
        x_title="Q (1/Å)",
        y_title="I(Q) (1/cm)",
        title=title,
        x_log=True,
        y_log=True,
        publish=False,
    )

    return html


# ============================================================================
# Entry Point Functions for drtsans Integration
# ============================================================================


def generate_gpr_analysis(
    iqmod_objects: Union[IQmod, List[IQmod]],
    output_dir: str,
    base_filename: str,
    **kwargs,
) -> Tuple[str, List[str], List[str]]:
    """
    Generate GPR analysis for one or more IQmod objects.

    This is the main entry point for autoreduction/livereduction integration.
    Generates both PNG files and Plotly HTML plots for each I(Q) profile.

    Parameters
    ----------
    iqmod_objects : IQmod or list of IQmod
        Single IQmod or list of IQmod objects (e.g., from slicing or wedges)
    output_dir : str
        Directory to save output files
    base_filename : str
        Base filename for outputs (e.g., "EQSANS_172839")
    **kwargs : dict
        Optional GPR parameters:
        - use_log_I : bool (default: False)
        - use_log_Q : bool (default: True)
        - lmbda : float (default: 0.25)
        - background_filter_width : float (default: 0.2)
        - background_intensity_offset : float (default: 1.0)

    Returns
    -------
    tuple
        (html_report, png_files, dat_files) where:
        - html_report : str - HTML string with all Plotly plots
        - png_files : list of str - Paths to generated PNG files
        - dat_files : list of str - Paths to generated .dat files

    Examples
    --------
    >>> from drtsans.dataobjects import IQmod
    >>> iq = IQmod.read_csv("EQSANS_172839_Iq.dat")
    >>> html, pngs, dats = generate_gpr_analysis(iq, "/output", "EQSANS_172839")
    >>> # html can be embedded in autoreduction report
    >>> # PNG files saved for archival purposes

    Notes
    -----
    - Handles multiple IQmod objects (wedges, time slices)
    - Gracefully handles errors - returns empty strings/lists on failure
    - Logs all operations for debugging
    """
    # Ensure iqmod_objects is a list
    if not isinstance(iqmod_objects, list):
        iqmod_objects = [iqmod_objects]

    # Validate inputs
    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            logger.error(f"Failed to create output directory {output_dir}: {e}")
            return "", [], []

    html_report = ""
    png_files = []
    dat_files = []

    # Process each IQmod object
    for idx, iqmod in enumerate(iqmod_objects):
        try:
            # Determine suffix for multiple profiles
            suffix = f"_wedge_{idx}" if len(iqmod_objects) > 1 else ""

            # Extract data from IQmod
            q = iqmod.mod_q
            I = iqmod.intensity
            I_err = iqmod.error
            dq = iqmod.delta_mod_q if iqmod.delta_mod_q is not None else None

            logger.information(f"Running GPR analysis for profile {idx} ({len(q)} data points)")

            # Run GPR analysis
            gpr_result = run_gpr(q, I, I_err, dq, **kwargs)

            # Generate output filenames
            png_path = os.path.join(output_dir, f"{base_filename}{suffix}_Iq_gpr.png")
            dat_path = os.path.join(output_dir, f"{base_filename}{suffix}_Iq_gpr.dat")

            # Create PNG plot
            create_png_plot(gpr_result, png_path, title=f"{base_filename}{suffix}")
            png_files.append(png_path)
            logger.information(f"Saved GPR PNG plot: {png_path}")

            # Save .dat file
            q_out = gpr_result["q"]
            I_gpr = gpr_result["I_gpr"]
            I_gpr_err = gpr_result["I_gpr_err"]
            dq_out = gpr_result.get("dq", np.zeros_like(q_out))

            data_to_save = np.column_stack((q_out, I_gpr, I_gpr_err, dq_out))
            header = "Q (1/A)\tI_GPR (1/cm)\tI_GPR_err (1/cm)\tdQ (1/A)"
            np.savetxt(dat_path, data_to_save, fmt="%.6e", delimiter="\t", header=header)
            dat_files.append(dat_path)
            logger.information(f"Saved GPR data file: {dat_path}")

            # Create Plotly HTML plot
            plot_title = f"GPR Analysis: {base_filename}{suffix}"
            html_plot = plotly_gpr_plot(gpr_result, title=plot_title)
            html_report += html_plot + "\n"

        except Exception as e:
            logger.error(f"GPR analysis failed for profile {idx}: {e}")
            # Continue with other profiles even if one fails
            continue

    if not html_report:
        logger.warning("No GPR plots generated")

    return html_report, png_files, dat_files


def run_gpr_from_file(input_filename: str, output_dir: Optional[str] = None, **kwargs) -> Tuple[str, str]:
    """
    Run GPR analysis on a .dat file (for CLI usage).

    Parameters
    ----------
    input_filename : str
        Path to input .dat file with columns: Q, I, I_err, dQ
    output_dir : str, optional
        Output directory (default: same as input file)
    **kwargs : dict
        Optional GPR parameters (see run_gpr)

    Returns
    -------
    tuple
        (png_path, dat_path) - Paths to generated PNG and .dat files

    Raises
    ------
    FileNotFoundError
        If input file does not exist
    ValueError
        If input file format is invalid
    """
    if not os.path.isfile(input_filename):
        raise FileNotFoundError(f"Input file not found: {input_filename}")

    # Load data
    try:
        data = np.loadtxt(input_filename, skiprows=2)
        # Check if data is 1D (single column) or 2D
        if data.ndim == 1:
            raise ValueError("Input file must have at least 3 columns: Q, I, I_err")
        if data.shape[1] < 3:
            raise ValueError("Input file must have at least 3 columns: Q, I, I_err")
        q = data[:, 0]
        I = data[:, 1]
        I_err = data[:, 2]
        dq = data[:, 3] if data.shape[1] > 3 else None
    except ValueError:
        # Re-raise ValueError with our message
        raise
    except Exception as e:
        raise ValueError(f"Failed to read input file {input_filename}: {e}")

    # Run GPR
    gpr_result = run_gpr(q, I, I_err, dq, **kwargs)

    # Determine output paths
    if output_dir is None:
        output_dir = os.path.dirname(input_filename)
    base_name, _ = os.path.splitext(os.path.basename(input_filename))

    png_path = os.path.join(output_dir, f"{base_name}_gpr.png")
    dat_path = os.path.join(output_dir, f"{base_name}_gpr.dat")

    # Save outputs
    create_png_plot(gpr_result, png_path)

    q_out = gpr_result["q"]
    I_gpr = gpr_result["I_gpr"]
    I_gpr_err = gpr_result["I_gpr_err"]
    dq_out = gpr_result.get("dq", np.zeros_like(q_out))

    data_to_save = np.column_stack((q_out, I_gpr, I_gpr_err, dq_out))
    header = "Q (1/A)\tI_GPR (1/cm)\tI_GPR_err (1/cm)\tdQ (1/A)"
    np.savetxt(dat_path, data_to_save, fmt="%.6e", delimiter="\t", header=header)

    logger.information(f"GPR analysis complete: {png_path}, {dat_path}")
    return png_path, dat_path


# ============================================================================
# Command-Line Interface
# ============================================================================


def main() -> int:
    """Command-line interface for GPR analysis."""
    parser = argparse.ArgumentParser(
        description="Gaussian Process Regression analysis for EQSANS I(Q) profiles",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help="Input .dat file with columns: Q, I, I_err, dQ",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory (default: same as input file)",
        default=None,
    )

    parser.add_argument(
        "--use-log-I",
        action="store_true",
        help="Apply log-linear transformation to intensity",
    )

    parser.add_argument(
        "--no-log-Q",
        action="store_true",
        help="Disable log transformation of Q (default: enabled)",
    )

    parser.add_argument(
        "--lmbda",
        type=float,
        default=0.25,
        help="RBF kernel length scale",
    )

    parser.add_argument(
        "--background-filter-width",
        type=float,
        default=0.2,
        help="Kernel width for background estimation",
    )

    parser.add_argument(
        "--background-intensity-offset",
        type=float,
        default=1.0,
        help="Background level for log transform",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    # Set up bash completion (optional)
    try:
        import argcomplete

        argcomplete.autocomplete(parser)
    except ImportError:
        pass  # argcomplete is optional

    # Set up TUI (optional, activated with --tui)
    try:
        from argparse_tui import add_tui_argument

        add_tui_argument(parser)
    except ImportError:
        pass  # argparse_tui is optional

    args = parser.parse_args()

    # Check if input_file is provided (required unless --tui was used)
    if args.input_file is None:
        if hasattr(args, "tui") and getattr(args, "tui", False):
            # TUI was run, exit gracefully
            return 0
        parser.error("input_file is required (or run with --tui)")

    # Run GPR analysis
    try:
        png_path, dat_path = run_gpr_from_file(
            args.input_file,
            output_dir=args.output_dir,
            use_log_I=args.use_log_I,
            use_log_Q=not args.no_log_Q,
            lmbda=args.lmbda,
            background_filter_width=args.background_filter_width,
            background_intensity_offset=args.background_intensity_offset,
        )
        print("Success! Generated files:")
        print(f"  PNG: {png_path}")
        print(f"  DAT: {dat_path}")
        return 0
    except Exception as e:
        logger.error(f"GPR analysis failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
