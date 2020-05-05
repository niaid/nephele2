import os
import boto3
from botocore.exceptions import ClientError
from nephele2 import NepheleError, tfvars, config

EFS_LOC = '/mnt/EFS'

# from AWS
MNT_OPTS = 'nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport'


def _gen_auth_keys():
    """
    calls gen_user_authkeys_cmd() for all users
    Returns:
        list of commands that create the authorized_keys for all users
    """
    cmds = []
    users = config.get_ec2_users()
    for user, uname in users.items():
        key = config.get_pub_keys(user)
        if key:
            cmds.extend(config.gen_user_authkeys_cmd(uname, key))
    return "\n".join(cmds)


def _gen_user_data(job_id):
    TF_VAR_neph_db_read_pw = os.environ.get('TF_VAR_neph_db_read_pw')
    TF_VAR_neph_db_write_pw = os.environ.get('TF_VAR_neph_db_write_pw')
    if not TF_VAR_neph_db_read_pw or not TF_VAR_neph_db_write_pw:
        msg = 'Unable to find read_pw:{r} or write:{w}'.format(
            r=TF_VAR_neph_db_read_pw,
            w=TF_VAR_neph_db_write_pw)
        raise NepheleError.UnableToStartEC2Exception(msg=msg)

    tmp_d = config.UPLOAD_PATH + job_id + '/tmp/'
    mnt_str = 'mount -t nfs -o {opts} {mnt_trgt}:/ {mnt_pt}/'.format(
        opts=MNT_OPTS,
        mnt_trgt=tfvars.EFS_IP,
        mnt_pt=EFS_LOC)
    env_str = 'export PYTHONPATH=/usr/local/src;'\
        'export TF_VAR_neph_db_read_pw={r};'\
        'export TF_VAR_neph_db_write_pw={w};'\
        'export TMPDIR={tmp_d};'\
        'export TMP={tmp_d};'\
        'export TEMP={tmp_d};'\
        'export JOB_ID={job_id};'\
        'export SECRET_KEY=;'\
        'export SALT=;'.format(job_id=job_id,
                               r=TF_VAR_neph_db_read_pw,
                               w=TF_VAR_neph_db_write_pw,
                               tmp_d=tmp_d)

    gen_auth_keys = _gen_auth_keys()
    pub_keys_str = "\n" + gen_auth_keys + "\n"
    # putting the job_id environment for the check EFS service to run
    u_data = '#!/bin/bash\n'\
             'printf "\nnet.ipv6.conf.all.disable_ipv6 = 1\n" >> /etc/sysctl.conf;'\
             'timedatectl set-timezone America/New_York;'\
             'printf "JOB_ID={job_id}" >> /etc/environment;'\
             '{pub_keys_str}'\
             'sysctl -p;'\
             'cd /usr/local/src/ && aws s3 cp s3://{pipe_src}/current.tgz .\n'\
             'tar xzf current.tgz \n'\
             'rm -rf current.tgz \n'\
             '{mnt_str}\n'\
             '{env_str}'\
             'mkdir -p {tmp_d};'\
             'chmod -R g+s {base_d};'\
             'chmod -R 777 {base_d};'\
             'sudo -Eu www-data bash -c "touch {base_d}/outputs/logfile.txt";'\
             'chown -R www-data:www-data {base_d};'\
             'chmod -R g+w {base_d};'\
             'mv /usr/local/src/nephele2/resources/misc_files_for_worker/chk_efs.* /etc/systemd/system/;'\
             'chmod a+x /etc/systemd/system/chk_efs.service;'\
             'systemctl start chk_efs.timer;'\
             'python3 /usr/local/src/nephele2/infra/utils/worker.py --process_job;'\
             '/sbin/shutdown -h "now"'.format(
                 pipe_src=tfvars.SRC_BUCKET_ID,
                 mnt_str=mnt_str,
                 env_str=env_str,
                 base_d=config.UPLOAD_PATH + job_id,
                 tmp_d=tmp_d,
                 job_id=job_id,
                 pub_keys_str=pub_keys_str)
    return u_data


def _find_available_subnet():
    """_find_available_subnet
    get all subnets, returns subnet with the highest number of
    free addresse.
    (will catch failure to launch)"""
    client = boto3.client('ec2', region_name='us-east-1')
    subnet_ids = tfvars.internal_subnets.split(',')
    subnets = client.describe_subnets(SubnetIds=subnet_ids)
    id_to_num_addr = dict()
    for subnet in subnets['Subnets']:
        id_to_num_addr[subnet['SubnetId']] = \
            int(subnet['AvailableIpAddressCount'])
    # list of subnet names, sorted by num addresses available.
    return sorted(id_to_num_addr, key=id_to_num_addr.get, reverse=True)


def start_worker_EC2(job_id, ami_id, instance_type):
    """
    returns a boto3 Instance object or raises UnableToStartEC2Exception
    TagSpecifications=[{'Tags' : [{'Key':'Name', 'Value': job_id}]}],
    # States can be:
    #     0 : pending
    #     16 : running
    #     32 : shutting-down
    #     48 : terminated
    #     64 : stopping
    #     80 : stopped
    Attempts to start EC2s until a suitable subnet is found.
    Will raise exception if no instance can be started
    """
    ec2 = boto3.resource('ec2', region_name='us-east-1')
    u_data = _gen_user_data(job_id)
    subnet_ids = _find_available_subnet()
    for i, subnet_id in enumerate(subnet_ids):
        try:
            instances = ec2.create_instances(
                DryRun=False,
                SecurityGroupIds=[tfvars.INTERNAL_SECURITY_GROUP,
                                  tfvars.ecs_cluster_security_group_id],
                IamInstanceProfile={'Arn': tfvars.N2_WORKER_INSTANCE_PROFILE},
                InstanceType=instance_type,
                KeyName=tfvars.KEY_NAME,
                ImageId=ami_id,
                MinCount=1,
                MaxCount=1,
                InstanceInitiatedShutdownBehavior='terminate',
                SubnetId=subnet_id,
                UserData=u_data)
            break
        except ClientError as ex:
            if i == len(subnet_ids) - 1:
                raise ex
    if len(instances) != 1:
        raise RuntimeError('Instances launched: %s' % str(instances))
    instance = instances[0]
    instance.wait_until_running()
    instance.create_tags(Tags=[{'Key': 'Name', 'Value': job_id}])
    return instance
