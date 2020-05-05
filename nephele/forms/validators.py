#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
from urllib.parse import urlparse

from wtforms.validators import InputRequired, ValidationError


def url_file_allowed(allowed_prefixes, message=None):
    """
    Validator that checks if a filename has an allowed prefix.

    Args:
        allowed_prefixes (array): list of allowed protocols
        message (str, optional): The error message to display if the value has
        an invalid protocol.
    """
    message = 'Invalid URL type.' if not message else message

    def _url_file_allowed(_, field):
        input_url = field.data
        parsed_url = urlparse(input_url)
        if parsed_url[0] not in allowed_prefixes:
            raise ValidationError(message)
        if not re.match(r'.*\/$', input_url):
            raise ValidationError(
                "Please provide the path to the folder where your \
                files reside (your path should end in a forward slash '/'). \
                We'll download the files requested by your mapping file.")
    return _url_file_allowed


def url_valid(_, field):
    """
    Checks that the provided value looks like a valid URL.
    """
    putative_url = urlparse(field.data)
    if putative_url.scheme is None or putative_url.netloc is None:
        raise ValidationError('The URL you provided appears invalid.\
                              Please provide a well formatted URL.')


def greater_than(minimum, message=None):
    """
    If a value is provided, checks to see if it is
    strictly greater than the minimum value.

    Args:
        minimum (int): a number or None
        message (str, optional): the error message to display

    Raises:
        ValidationError
    """
    message = "Value must be greater than " + \
        str(minimum) if not message else message

    def _greater_than(_, field):
        # we can't use if field.data here because it won't pass the if when
        # field.data = 0
        if (field.data is not None and
                field.data != '' and
                field.data <= minimum):
            raise ValidationError(message)
    return _greater_than


class RequiredIf(InputRequired):
    """
    A validator which makes a field required if
    another field is set and has a truthy value

    from Mehdi Sadeghi:
    https://stackoverflow.com/questions/8463209/how-to-make-a-field-conditionally-optional-in-wtforms
    """

    def __init__(self, other_field_name, *args, **kwargs):
        self.other_field_name = other_field_name
        super(RequiredIf, self).__init__(*args, **kwargs)

    def __call__(self, form, field):
        other_field = form._fields.get(self.other_field_name)
        if other_field is None:
            raise Exception('no field named "%s" in form' %
                            self.other_field_name)
        if bool(other_field.data):
            super(RequiredIf, self).__call__(form, field)
