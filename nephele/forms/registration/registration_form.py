"""
WTForms definition of the registration for that is rendered with the
register.html template.
"""
from flask_wtf import FlaskForm, RecaptchaField
from wtforms import BooleanField, StringField, SelectField, TextAreaField
from wtforms.validators import InputRequired, DataRequired, Email, EqualTo,\
    ValidationError
from sqlalchemy.orm.exc import MultipleResultsFound
from nephele2.infra.utils.neph_utils import N2Manager
from nephele2 import NepheleError


def check_if_registered(_, field):
    """
    Custom validator that raises errors if the email is already
    in the database, or if no email address was provided.
    """
    try:
        user = N2Manager.get_user(field.data)
        if user:
            raise ValidationError('This email address is already registered')
    except NepheleError.NepheleMissingArgument:
        raise ValidationError('This field requires a valid email address')
    except MultipleResultsFound:
        raise ValidationError('This email address is already registered')
    except NepheleError.NepheleRowNotFound:
        pass


class RegistrationForm(FlaskForm):
    """
    Defines the user registration form
    """
    fname = StringField('First Name',
                        [InputRequired("Please enter your first name")])
    lname = StringField('Last Name',
                        [InputRequired("Please enter your last name")])
    email = StringField('Email Address',
                        [InputRequired("Please enter your email address"),
                         Email("This field requires a valid email address"),
                         check_if_registered]
                        )
    email_chk = StringField(
        'Re-enter Email Address',
        [InputRequired("Please re-enter your email address"),
         Email("This field requires a valid email address"),
         EqualTo('email', "This field must match the original email address")])
    affil = StringField(
        'Affiliation',
        [InputRequired("We want to know who's using our service. \
                              Please let us know where you're coming from.")])
    affil_cat = SelectField(
        u'How would you categorize yourself?',
        choices=[('NIH', 'NIH Employee or Contractor'),
                 ('Government (non-NIH)', 'Government (non-NIH)'),
                 ('University', 'University'),
                 ('Research Institute (Private)',
                  'Research Institute (Private)'),
                 ('Citizen Scientist (no-affiliation)',
                  'Citizen Scientist (no-affiliation)'),
                 ('Other', 'Other')
                 ])
    ref = SelectField(
        u'How did you hear about us?',
        choices=[
            ('NA', '--'),
            ('NIH Bioinformatics listserv', 'NIH Bioinformatics listserv'),
            ('NIH Microbiome listserv', 'NIH Microbiome listserv'),
            ('Internet', 'Internet search'),
            ('Colleague', 'My colleague'),
            ('Other', 'Other')
        ])
    analysis = TextAreaField('Please provide a brief explanation of the \
                              analysis you plan to conduct using Nephele')
    subscribe = BooleanField('I would like to subscribe to Nephele newsletter',
                             default="True")
    priv_label = 'I agree to the \
        <a href="/privacy_policy">Nephele Privacy Policy</a>'
    priv_pol = BooleanField(
        priv_label,
        [DataRequired("You must accept the privacy policy to continue")])
    recaptcha = RecaptchaField()
