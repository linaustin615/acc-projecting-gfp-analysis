"""Repository paths and the frozen analysis settings."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
MATCH_RADIUS_UM = 9.
DENSE_RADIUS_UM = 250.
DENSE_MINIMUM = 20
DENSE_GRID_PX = 32
DENSE_SEPARATION_UM = 500.
