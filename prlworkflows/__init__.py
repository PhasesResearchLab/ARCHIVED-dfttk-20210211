
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from prlworkflows.sqs import SQS, enumerate_sqs
from prlworkflows.sqs_db import SQSDatabase, get_structures_from_database
