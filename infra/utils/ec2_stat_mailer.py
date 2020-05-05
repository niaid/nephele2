#!/usr/bin/env python3

import datetime
import boto3
from nephele2.infra.utils.email import email


MAIL_SUBJCT='Amazon instances running for a long time'

def list_long_running():
    session = boto3.Session( region_name='us-east-1')
    ec2 = session.resource('ec2')
    instance_iterator = ec2.instances.all()
    long_running = list()
    for i in instance_iterator:
        if i.state['Name'] == 'running':
            lt = i.launch_time
            td = datetime.datetime.now(datetime.timezone.utc) - lt
            line = ['Days alive : '+str(td.days),
                    'Tag : '+i.tags[0]['Value'],
                    'ID : '+i.id]
            long_running.append(' ; '.join(line))
    return long_running


if __name__ == "__main__":
    print('Trying to run EC2 status mailer')
    long_runs = list_long_running()
    list_elts = ''
    for long_run in long_runs:
        list_elts += '<li>'+long_run+'</li>'
    email.send_admin_info_email(MAIL_SUBJCT, '<ul>'+list_elts+'</ul>')
