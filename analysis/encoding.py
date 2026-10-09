"""Lossless JSON preparation for exact deterministic examples."""
from fractions import Fraction as F
import numpy as np

def clean(x):
    """Represent rational values as exact fraction strings, never rounded floats."""
    if isinstance(x, F):
        return str(x)
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    return x
