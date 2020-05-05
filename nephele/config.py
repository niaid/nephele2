"""web server specific config"""
import os
from itsdangerous import URLSafeTimedSerializer
SECRET_KEY = os.environ['SECRET_KEY']
SALT = os.environ['SALT']
RECAPTCHA_PUBLIC_KEY = os.environ['RC_SITE_KEY']
RECAPTCHA_PRIVATE_KEY = os.environ['RC_SECRET_KEY']

# Timed Serializer: for securely transferring user email between client/server
TS = URLSafeTimedSerializer(SECRET_KEY)


def salt_string(string):
    """salt_string

    :param string:
    """
    return TS.dumps(string, salt=SALT)


def unsalt_string(string):
    """unsalt_string

    :param string:
    """
    return TS.loads(string,
                    salt=SALT,
                    max_age=86400)
