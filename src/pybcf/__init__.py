from importlib.metadata import version

__name__ = 'bcf'
__version__ = version(__name__)

from bcf.reader import BcfReader
