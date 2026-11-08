%matplotlib inline

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from IPython.display import Image, display
from matplotlib import animation
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from esp_lab.utils.filename_utils import figure_filename

from util import rtenso_case as rc
from util import strat_case as sc
from util import vortex_plot as vp

warnings.filterwarnings("ignore", message="invalid value encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")