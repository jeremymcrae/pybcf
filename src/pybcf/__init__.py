from importlib.metadata import version

__name__ = 'pybcf'
__version__ = version(__name__)

from pybcf.reader import BcfReader
