#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
