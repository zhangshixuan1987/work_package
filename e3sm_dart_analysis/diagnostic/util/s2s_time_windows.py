# utils/time_windows.py
import numpy as np

def build_windows_from_start(start_date, lead_bins):
    """
    Convert lead-day bins into calendar-date windows.

    Parameters
    ----------
    start_date : str
        Forecast start date, e.g. "2012-01-01"
    lead_bins : dict
        {"week1-2": (1,14), "week3-4": (15,28), ...}

    Returns
    -------
    dict
        {"week1-2": (np.datetime64, np.datetime64), ...}
    """
    start = np.datetime64(start_date)
    windows = {}
    for label, (d0, d1) in lead_bins.items():
        windows[label] = (
            start + np.timedelta64(d0 - 1, "D"),
            start + np.timedelta64(d1 - 1, "D"),
        )
    return windows
