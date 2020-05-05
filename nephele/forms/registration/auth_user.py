"""
An implementation of the user class that is required by Flask-Login.
See https://flask-login.readthedocs.io/en/latest/ for details.
"""
from flask_login import UserMixin
from nephele2.rds.db_utils import DBUtils


class AuthUser(UserMixin):
    """
    Based on the Mixin provided by flask-login, this class
    defines the authenticated user object required by flask-login.
    """
    def __init__(self, email_addr):
        self._email_addr = email_addr
        self._is_registered = False
        self._is_bad = False
        self._is_confirmed = False
        self._is_active = False
        self._has_compute = False

        (user, user_address) = DBUtils.get_user_and_email_by_email(email_addr)
        if user and user_address:
            self._is_registered = True
            if user_address.is_bad:
                self._is_bad = True
            if user_address.is_confirmed:
                self._is_confirmed = True
            if user.compute_remaining > 0:
                self._has_compute = True
            if (self._is_registered and
                    self._is_confirmed and
                    not self._is_bad and
                    self._has_compute):
                self._is_active = True

    def __str__(self):
        return ("AuthUser :\n"
                "id:{}, "
                "is_registered:{}, "
                "is_bad:{}, "
                "is_confirmed:{}, "
                "has_compute:{}, "
                "is_active:{}, ".format(
                    self._email_addr,
                    self._is_registered,
                    self._is_bad,
                    self._is_confirmed,
                    self._has_compute,
                    self._is_active))

    @property
    def is_registered(self):
        """is_registered"""
        return self._is_registered

    @property
    def is_bad(self):
        """is_bad"""
        return self._is_bad

    @property
    def is_confirmed(self):
        """is_confirmed"""
        return self._is_confirmed

    @property
    def has_compute(self):
        """
        Checks remaining compute hours allotted to the user.
        This is not required by Flask-Login, but is used by
        us to determine whether or not to allow the user to
        log in to the system and submit jobs.
        """
        return self._has_compute

    @property
    def is_active(self):
        """Needed by Flask Login:
        must return boolean! (exceptions = false)"""
        return self._is_active

    @property
    def is_authenticated(self):
        """Needed by Flask Login:
        must return boolean! (exceptions = false)"""
        return self._is_active

    @property
    def is_anonymous(self):
        """Needed by Flask Login:
        must return boolean! (exceptions = false)"""
        return False

    def get_id(self):
        """Method is needed by Flask Login:
        https://flask-login.readthedocs.io/en/latest/#your-user-class"""
        return self._email_addr
