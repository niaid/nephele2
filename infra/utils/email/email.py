#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import collections
import html2text
import jinja2
import boto3
from nephele2 import config, tfvars
from nephele2.rds.db_utils import DBUtils

FAQ_URL = tfvars.SERVER_ADDR+"/user_guide_pipes/#output"


def notify_user_job_queued(job_id):
    """notify_user_job_queued

    :param job_id:
    Notifies users if job is unable to be run immedeatly due to high
    usage.
    """
    user = DBUtils.get_user_by_job_id(job_id)
    mail_subject = DBUtils.make_email_subject(job_id)
    msg = _apply_template('job_queued',
                          base=tfvars.SERVER_ADDR,
                          job_id=job_id)
    _send(user.email_address, mail_subject, msg)


def notify_user_job_completion(job_id):
    """
    Sends an email to the user with a completion status. If the job
    succeeds, a success email is sent. If the job fails, any error message
    associated with the failure is retreived from the database and a
    failure email is sent.

    Args:
        job_id (str): the unique idenifier for the completed job

    Raises:
        Exception: any exception raised by either DBUtils or aws_utils
    """
    user = DBUtils.get_user_by_job_id(job_id)
    job_details = DBUtils.get_job_type_details(job_id)
    job_type = job_details.get('data_type')
    mail_subject = DBUtils.make_email_subject(job_id)
    final_stat = DBUtils.get_job_status(job_id)

    if final_stat and final_stat == 'Succeeded':
        msg = _apply_template('pipeline_completed',
                              base=tfvars.SERVER_ADDR,
                              job_id=job_id,
                              job_type=job_type,
                              faq_url=FAQ_URL)
    else:
        err_msg = DBUtils.get_error_msg(job_id)
        msg = _apply_template('pipeline_error',
                              base=tfvars.SERVER_ADDR,
                              job_id=job_id,
                              job_type=job_type,
                              err_msg=err_msg)
    return _send(user.email_address, mail_subject, msg)


def _apply_template(template_name, **kwargs):
    """
    Generates a properly formatted SES email body dictionary from a jinja2
    template.

    Args:
        template_name (str): the name of the template to use (no .html
        extension)
        **kwargs (dict): the key, value arguments that are needed by the
        template being requested.

    Returns:
        dict: a properly formatted SES email body dictionary

    Raises:
        TemplateNotFound: the requested template wasn't found in the
        resources folder
    """
    # the FileSystemLoader demands a full path
    # relative paths only work if caller never changes loc
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    template_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=template_dir))
    template = template_env.get_template(template_name+'.html')
    html = template.render(kwargs)

    msg = {'Text':  {'Data': html2text.html2text(html)},
           'Html':  {'Data': html}}
    return msg


def _send(email_address, subject, msg):
    destination = {'ToAddresses': [email_address]}
    message = {'Subject': {'Data': subject},
               'Body': msg}
    session = boto3.Session(region_name='us-east-1')
    ses = session.client('ses')
    return ses.send_email(Source=config.get_sender_address(),
                          Destination=destination,
                          Message=message)


def send_admin_info_email(subject, msg):
    msg = _apply_template('infra_info', msg=str(msg))
    return _send(config.INFRA_MAIL, subject, msg)


def send_infra_failure_email(msg, job_id=None):
    subject = "Infrastructure failure -- {env}".format(
        env=tfvars.environment)
    if job_id:
        stack_trace = DBUtils.get_stack_trace(job_id)
        if stack_trace:
            msg = msg + "\n" + stack_trace
    msg = _apply_template('infra_failure', jobid=job_id, stack=str(msg))
    return _send(config.INFRA_MAIL, subject, msg)


def send_started_email(job_id):
    """
    Send an email to the user notifying them that their job started.
    Uses the pipeline_started.html template to render the message body.

    Args:
        job_id (str): a valid job ID

    Raises:
        NepheleError.NepheleEmailError: on failure
    """
    user = DBUtils.get_user_by_job_id(job_id)
    job_args = DBUtils.get_job_arguments(job_id)
    job_desc = DBUtils.get_job_description(job_id)
    email_address = user.email_address
    mail_subject = DBUtils.make_email_subject(job_id)
    sorted_args = collections.OrderedDict(job_args)
    job_details = DBUtils.get_job_type_details(job_id)
    msg = _apply_template('pipeline_started',
                          job_args=sorted_args,
                          pipe_name=job_details['name'],
                          job_id=job_id,
                          job_desc=job_desc,
                          neph_version=config.NEPHELE_VERSION,
                          logUrl=tfvars.SERVER_ADDR+"/view_log/"+job_id)
    return _send(email_address, mail_subject, msg)


def send_registration_email(email_addr, confirm_url):
    subject = "Please complete your Nephele registration"
    body = _apply_template('registration_confirmation',
                           confirm_url=confirm_url)
    return _send(email_addr, subject, body)
