# Connection EG:
# mysql -u fancy_uname -h
# terraform-0089b44840ace5118ca2481414.c4dtf09rvl4e.us-east-1.rds.amazonaws.com
# -P 3306 -p

import os
from nephele2 import tfvars

_DB_NAME = 'nephele_db'
# if the the env is not set (eg for tests), give it a fake val
if os.environ.get('NEPH_DB_LOCAL'):
    _ROOT_UNAME = 'root'
    _ROOT_PW = 'root_pw'
    _READ_PW = 'read_pw'
    _WRITE_PW = 'write_pw'
    _INIT_PW = 'init_pw'
    _HOST = 'localhost'
else:
    _ROOT_UNAME = os.environ.get('TF_VAR_neph_db_root_uname')
    _ROOT_PW = os.environ.get('TF_VAR_neph_db_root_pw')
    _READ_PW = os.environ.get('TF_VAR_neph_db_read_pw')
    _WRITE_PW = os.environ.get('TF_VAR_neph_db_write_pw')
    _INIT_PW = os.environ.get('TF_VAR_neph_db_init_pw')
    _HOST = tfvars.DB_ENDPOINT

_RO_URI = ('mysql+mysqlconnector://{uname}:{pw}@{host}/{db}'
           .format(uname='nephele_ro',
                   pw=_READ_PW,
                   host=_HOST,
                   db='nephele_db'))

_RW_URI = ('mysql+mysqlconnector://{uname}:{pw}@{host}/{db}'
           .format(uname='nephele_rw',
                   pw=_WRITE_PW,
                   host=_HOST,
                   db='nephele_db'))
# used in __init__
READ = {'SQL_DATABASE_URI': _RO_URI,
        'SQL_POOL_SIZE': 2,
        'SQL_ECHO': False}

# used in __init__
WRITE = {'SQL_DATABASE_URI': _RW_URI,
         'SQL_POOL_SIZE': 2,
         'SQL_ECHO': False}
#
DB_CONFIG = {
    'user': _ROOT_UNAME,
    'password': _ROOT_PW,
    'host': _HOST,
    'port': '3306',
    'database': _DB_NAME,
    'raise_on_warnings': True
}

INIT_URI = ('mysql+mysqlconnector://{uname}:{pw}@{host}/{db}'
            .format(uname='inituser',
                    pw=_INIT_PW,
                    host=_HOST,
                    db='nephele_db'))
