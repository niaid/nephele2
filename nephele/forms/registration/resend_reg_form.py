#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wtforms import SubmitField
from flask_wtf import FlaskForm, RecaptchaField


class ResendEmailForm(FlaskForm):
    """ResendRegistrationForm
    handles resending of confirmation emails"""
    recaptcha = RecaptchaField()
    submit = SubmitField('Resend email')
