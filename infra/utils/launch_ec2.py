#!/usr/bin/env python3
import os
import boto3
import botocore.exceptions
import argparse
import yaml

from nephele2 import NepheleError
mand_vars = ['AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY']
perm_error = """\n\nIt seems you have not set up your AWS correctly.
Should you be running this with Awssume? Or have profile with appropriate role?
Exiting now.\n"""



def main(args):
    """Launch ec2 instance"""
    if args.profile is None:
            ec2_resource = boto3.Session(region_name='us-east-1').resource('ec2')
    else:
        ec2_resource = boto3.Session(region_name='us-east-1', profile_name=args.profile).resource('ec2')
    test_sanity(ec2_resource, args)
    envs = load_stack_vars(args.yaml_env.name)
    start_EC2(ec2_resource, args.ami_id, args.instance_type,
              args.key_path, args.label, envs, args.dry_run)


def load_stack_vars(fname):
    try:
        with open(fname) as f:
            data_map = yaml.safe_load(f)
        return data_map
    except FileNotFoundError as fnf:
        print(fnf)
        print('Unable to find yaml file, exiting.')
        exit(1)
    except:
        raise


def gen_mnt_str(efs_ip):
    mnt_opts = 'nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport'  # from AWS
    return 'mount -t nfs -o {opts} {trgt}:/ {mnt}/'.format(opts=mnt_opts,
                                                           trgt=efs_ip,
                                                           mnt='/mnt/EFS')

def read_key(key_path):
    try:
        with open(key_path, 'r') as f:
            key = f.read()
            return key
    except:
        raise


def test_sanity(ec2_resource, args):
    """Test if env vars are set, key exists, and can access ec2"""
    if args.profile is None:
        for var in mand_vars:
            if os.environ.get(var) is None:
                print(var + ' must be set as an evironment variable. \nExiting.')
                exit(1)
    if not os.path.exists(args.key_path):
        print('Unable to see your key: {}, exiting now :-('.format(args.key_path))
        exit(1)
    try:
        ec2_resource.instances.all().__iter__().__next__()
    except botocore.exceptions.ClientError as expn:
        print(expn)
        print(perm_error)
        exit(1)

def create_EC2(ec2_resource, ami_id, i_type, envs, u_data='', dry_run=True):
    """create ec2 instance. by default DryRun is T, and only checks perms."""
    inst = ec2_resource.create_instances(
            DryRun=dry_run,
            SecurityGroupIds=[envs['INTERNAL_SECURITY_GROUP'],
                              envs['ecs_cluster_security_group_id']],
            IamInstanceProfile={'Arn': envs['N2_WORKER_INSTANCE_PROFILE']},
            InstanceType=i_type,
            ImageId=ami_id,
            MinCount=1,
            MaxCount=1,
            InstanceInitiatedShutdownBehavior='terminate',
            SubnetId=envs['VPC_SUBNET'],
            UserData=u_data
        )
    return inst

def start_EC2(ec2_resource, ami_id, i_type, key_path, label, envs, dry_run):
    """check if have perms to create instance.
       https://boto3.amazonaws.com/v1/documentation/api/latest/guide/ec2-example-managing-instances.html#start-and-stop-instances
       if so, the create instance and tag with label.
    """
    try:
        create_EC2(ec2_resource, ami_id, i_type, envs)
    except botocore.exceptions.ClientError as e:
        if 'DryRunOperation' not in str(e):
            print(e.response['Error']['Message'])
            print(perm_error)
            exit(1)
        elif dry_run:
            print(e.response['Error']['Message'])
            exit(0)
        else:
            pass

    mnt_str = gen_mnt_str(envs['EFS_IP'])
    key_str = read_key(key_path)
    auth_key_str = 'printf "{}" >> /home/admin/.ssh/authorized_keys;'.format(
        key_str)
    u_data = '#!/bin/bash\n{mnt_str}\n{auth_key_str}\n'.format(mnt_str=mnt_str,
                                                               auth_key_str=auth_key_str)
    print('Creating EC2...')
    try:
        instances = create_EC2(ec2_resource, ami_id, i_type, envs, u_data, False)

    except botocore.exceptions.ClientError as bce:
        print(bce)
        print('\nUnable to launch EC2. \nExiting.')
        exit(1)
    if len(instances) is not 1:
        msg = 'Instances launched: %s' % str(instances)
        raise NepheleError.UnableToStartEC2Exception(msg=msg)
    instance = instances[0]
    instance.wait_until_running()
    instance.create_tags(Tags=[{'Key': 'Name', 'Value': label}])
    print(str(instance) + ' has been created.')
    print('To connect type:\nssh {ip_addr}'.format(
        ip_addr=instance.instance_id))
    print('To terminate instance type:')
    print('awssume aws ec2 terminate-instances --instance-ids ' + instance.instance_id)


if __name__ == "__main__":
    usage = 'Eg:\nsource ~/code/neph2-envs/dev/environment_vars\n'\
            'awssume launch_ec2.py -e ../../neph2-envs/dev/dev_outputs.yaml -a ami-0ae1b7201f4a236f9 -t m5.4xlarge -k ~/.ssh/id_rsa.pub --label instance_name_tag\n\n'\
            'Alternately, pass profile which has correct role/permissions:\n'\
            'launch_ec2.py -e dev_outputs.yaml -a ami-003eed27e5bf2ef91 -t t2.micro -k ~/.ssh/id_rsa.pub -l name_tag --profile aws_profile_name'
    parser = argparse.ArgumentParser(
        description='CLI Interface to N2.', usage=usage)
    req = parser.add_argument_group('required args')
    req.add_argument("-e", "--yaml_env",
                        type=argparse.FileType('r'), required=True)
    req.add_argument("-t", "--instance_type", type=str, required=True)
    req.add_argument("-a", "--ami_id", type=str, required=True)
    req.add_argument("-k", "--key_path", type=str, required=True)
    req.add_argument("-l", "--label", type=str, required=True)
    parser.add_argument("-p", "--profile", type=str)
    parser.add_argument("-d", "--dry_run", action='store_true')
    args = parser.parse_args()
    main(args)
