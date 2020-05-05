import boto3
"""system wide vars"""

NEPHELE_LOC_ON_SRVR = '/usr/local/src/nephele2'
PIPELINES_LOC_ON_WRKR = NEPHELE_LOC_ON_SRVR + '/pipelines/'
UPLOAD_PATH = '/mnt/EFS/user_uploads/'
INPUTS_DNAME = 'inputs'
OUTPUTS_DNAME = 'outputs'
ALLOWED_EXTS = '*.sff, *.fa, *.fasta, *.fas, *.faa, *.fna, *.qual, *.fastq, '\
    '*.fq, *.oligos, *.oligo, *.txt, *.sff.gz, *.fa.gz, *.fasta.gz, *.fas.gz, '\
    '*.faa.gz, *.fna.gz, *.qual.gz, *.fastq.gz, *.fq.gz, *.oligos.gz, '\
    '*.oligo.gz, *.txt'
GZ_COL_NAMES = ['ForwardFastqFile', 'ReverseFastqFile']
NEPHELE_VERSION = '2.4.0'
NEPHELE_TAG = 'Nephele_2020_Jan_29'
DEV_ACC_ID = ''
PROD_ACC_ID = ''
DEV_USERS = {}
PROD_USERS = {}
SHA_LONG = "892b1dd31a841aac59fe8b9a5defa3beb6995c65"
SHA_SHORT = "892b1dd3"

def get_sender_address():
    if get_account() == DEV_ACC_ID:
        return "NIAID BCBB Nephele Team"
    return "NIAID BCBB Nephele Team"


def gen_user_authkeys_cmd(uname, key):
    """gen_user_authkeys_cmd

    Creates a string with when executed creates an authorized_keys
    file in user home dir.
    :param uname:
    :param key:
    """
    cmds = [
        'id {0} || /usr/sbin/useradd -d /home/{0} -m -s /bin/bash -U {0}'.
        format(uname),
        'chmod 0700 /home/{}'.format(uname),
        'mkdir -p /home/{}/.ssh'.format(uname),
        'printf "{key}\n" > /home/{uname}/.ssh/authorized_keys '.
        format(uname=uname, key=key)
    ]
    return cmds


def get_pub_keys(user):
    """
    looks up users' keys
    Returns:
        list of public keys
    """
    keys = []
    iam = boto3.client('iam')
    try:
        user_response = iam.list_ssh_public_keys(UserName=user)
    except iam.exceptions.NoSuchEntityException:
        return None
    aws_keys = user_response.get('SSHPublicKeys')
    for aws_key in aws_keys:
        sshpublickeyid = aws_key.get('SSHPublicKeyId')
        if sshpublickeyid:
            key_response = iam.get_ssh_public_key(
                UserName=user,
                SSHPublicKeyId=sshpublickeyid,
                Encoding='SSH')
            keys.append(key_response.get(
                'SSHPublicKey').get('SSHPublicKeyBody'))
    return "\n".join(keys)


def get_account():
    """get_account
    asks the name of the AWS account - eg dev or prod."""
    sts = boto3.client('sts')
    ident = sts.get_caller_identity()
    return ident.get('Account')


def get_ec2_users():
    """ asks if we are running in a dev account or not,
    returns appropriate user names
    Returns:
        list
    """
    account_id = get_account()
    if account_id == DEV_ACC_ID:
        return DEV_USERS
    if account_id == PROD_ACC_ID:
        return PROD_USERS
