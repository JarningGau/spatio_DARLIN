import os

_legacy_matlab_dir = os.path.dirname(os.path.dirname(__file__))
_repo_root = os.path.dirname(os.path.dirname(_legacy_matlab_dir))

script_dir = os.path.join(_legacy_matlab_dir, "bin")
QC_dir = os.path.join(_repo_root, "QC")
CARLIN_dir = os.path.join(os.path.dirname(_repo_root), "CARLIN_pipeline")
