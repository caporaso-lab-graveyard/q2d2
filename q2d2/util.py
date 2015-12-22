import os
from appdirs import AppDirs
import q2d2

def get_install_path():
    return os.path.abspath(os.path.split(__file__)[0])

data_dirs = AppDirs("q2d2", "biocore", version=q2d2.__version__)
