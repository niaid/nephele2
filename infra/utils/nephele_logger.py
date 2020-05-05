import logging
from logging.config import dictConfig


def create_logger():
    logging_config = dict(
        version=1,
        formatters={
            'f': {'format':
                  '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
        },
        handlers={
            'h': {'class': 'logging.StreamHandler',
                  'formatter': 'f',
                  'level': logging.ERROR}
        },
        root={
            'handlers': ['h'],
            'level': logging.ERROR,
        },
    )

    dictConfig(logging_config)
    #logger = logging.getLogger()
    return logging.getLogger()
