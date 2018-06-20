
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from prlworkflows.structure_builders.sqs import SQS, enumerate_sqs
