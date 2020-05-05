"""AWS lambda function for adding USERS to the Bastion server, to allow them
ssh access. The reason there's two names is that the internal BCBB account 
uses diff unames to NIH. This lambda will run periodically to ensure bastion 
has keys if machine was nuked."""

# -*- coding: utf-8 -*-

import boto3
from nephele2 import config


#  BASTION_ID = 'i-092154924f6d72258'
BASTION_ID = ''
WEB_SERVER = ''


def _run_command(user, uname):
    key = config.get_pub_keys(user)
    if key:
        cmds = config.gen_user_authkeys_cmd(uname, key)
        client = boto3.client('ssm')
        return client.send_command(
            InstanceIds=[BASTION_ID, WEB_SERVER],
            DocumentName='AWS-RunShellScript',
            Parameters={'commands': cmds})

def lambda_handler(event, context):
    users = config.get_ec2_users()
    for user, uname in users.items():
        _run_command(user, uname)

