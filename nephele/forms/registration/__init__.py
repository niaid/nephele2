#!/usr/bin/env python
# -*- coding: utf-8 -*-
from flask_login import LoginManager

from nephele2.nephele.forms.registration.auth_user import AuthUser

login_manager = LoginManager()
login_manager.login_view = "login"
login_manager.session_protection = "strong"


@login_manager.user_loader
def load_user(user_id):
    """
    The user_loader required by flask_login.
    https://flask-login.readthedocs.io/en/latest/

    According to the docs, this must return None if the user
    isn't found.  It cannot raise an exception.

    Args:
        user_id (str): an ID that uniquely identifies a user

    Returns:
        AuthUser: an AuthUser or None if the user is not found
    """
    user = AuthUser(user_id)
    if user.is_registered:
        return user
    return None
