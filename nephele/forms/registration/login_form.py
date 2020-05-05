#!/usr/bin/env python
# -*- coding: utf-8 -*-
import flask
import flask_login
import flask_wtf
import wtforms
from wtforms.validators import InputRequired, Email, ValidationError
from nephele2.nephele.forms.registration.auth_user import AuthUser
from nephele2.nephele import config as server_conf


def _check_if_auth_user(_, field):
    """
    Validator: Checks to see if the user is authorized to submit jobs.

    Raises:
        ValidationError
    """
    email_addr = field.data
    salted_email = server_conf.salt_string(email_addr)
    resend_url = flask.url_for('resend_registration_email', email=salted_email)
    reg_url = flask.url_for('register')
    resend_link = "<a href="'{}'">re-send</a>".format(resend_url)
    reg_link = "<a href="'{}'">register</a>".format(reg_url)
    horrible_link = """<a href=
    "mailto:nephelesupport?Subject=Nephele Help - Request Compute"
     onmouseover="this.href=this.href.replace(\'@@\',\'.\')"
     onclick="this.href=this.href.replace(\'@@\',\'.\')">contact us</a>"""
    no_compute_err = ('You have used your current allocation of compute time, '
                      ' please {} to request more.'.format(horrible_link))
    bad_email_err = ("We have had issues sending email to this address. "
                     "Please {} a new address.".format(reg_link))
    no_confirm_err = ('Please reply to the confirmation email we sent you. '
                      'If you can\'t find that mail we can {} one.'.
                      format(resend_link))
    no_reg_err = ('You must be a registered user to use this system. '
                  'Please {} your email address.'.format(reg_link))
    user = AuthUser(email_addr)
    # these Exceptions seem to have to be thrown here
    # It doesn't like them being raised from AuthUser constructor,
    if user.is_active:
        flask_login.login_user(user)
    elif not user.is_registered:
        raise ValidationError(no_reg_err)
    elif not user.is_confirmed:
        raise ValidationError(no_confirm_err)
    elif user.is_bad:
        raise ValidationError(bad_email_err)
    elif not user.has_compute:
        raise ValidationError(no_compute_err)


class LoginForm(flask_wtf.FlaskForm):
    """
    Defines the form rendered with the login.html template.
    """
    email_addr = wtforms.StringField(
        'Please enter your registered email address:',
        [InputRequired("Please enter your registered email address"),
         Email("This field requires a valid email address"),
         _check_if_auth_user])
