#!/usr/bin/env python
"""
deals with the control of the pending queue
"""
# -*- coding: utf-8 -*-
from collections import namedtuple
from datetime import datetime
import boto3
from nephele2 import tfvars

Job = namedtuple('Job', ['job_id', 'backoff', 'user_notified'])

_SQS = boto3.resource('sqs', region_name='us-east-1')
_PENDING_Q = _SQS.Queue(tfvars.PENDING_Q_URL)


def _create_q(qname):
    """
    Creates a queue with the given name.

    Args:
        qname (str): the name to give the queue

    Raises:
        NepheleError.NepheleQueueCreationException: queue creation failed
    """
    return _SQS.create_queue(QueueName=qname, Attributes={'DelaySeconds': '5'})


def add_to_pending(job_id, user_notified=None, backoff=None):
    """add_to_pending

    :param job_id:
    :param user_notified:
    :param backoff:
    """
    if not user_notified:
        user_notified = 'false'
    if not backoff:
        backoff = datetime.now()
    _PENDING_Q.send_message(
        MessageBody=job_id,
        MessageAttributes={
            'job_id': {
                'DataType': 'String',
                'StringValue': job_id
            },
            'backoff': {
                'DataType': 'Number',
                'StringValue': str(backoff.timestamp())
            },
            'user_notified': {
                'DataType': 'String',
                'StringValue': user_notified
            }
        })


def receive_pending():
    """receive_pending
    deletes anything in the queue that doesn't have a job attribute
    Looks for messages that are due to be run, that is messages
    whose timestamp has passed.
    If message backoff timestamp has passed, job will attempt to be started
    """
    jobs = list()
    msgs = _PENDING_Q.receive_messages(
        MessageAttributeNames=['All'],
        VisibilityTimeout=0,
        WaitTimeSeconds=20
    )
    for msg in msgs:
        job_id = msg.message_attributes.get('job_id')
        if not job_id:
            print(msg + 'has no job_id - deleting')
            msg.delete()
        backoff = msg.message_attributes.get('backoff')
        backoff_dt = datetime.fromtimestamp(float(backoff['StringValue']))
        if backoff_dt > datetime.now():
            continue
        user_notified = msg.message_attributes.get('user_notified')
        u_note = bool(user_notified['StringValue'] == 'true')
        job = Job(job_id=job_id['StringValue'],
                  backoff=backoff_dt,
                  user_notified=u_note)
        jobs.append(job)
        msg.delete()
    return jobs
