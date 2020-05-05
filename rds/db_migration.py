"""
Database migration script. To update the schema::
    cd rds
    python3 db_migration.py db migrate
    python3 db_migration.py db upgrade

The migration script, which is added to the `migrations` directory, should be reviewed
and can be manually edited before running upgrade.

Note that the migration script requires more privileges than the default nephele\_user has,
so the script needs to be run with the privileged `inituser` instead, which it should do
by default. Both of these users should exist within the database if the
`sql/create_users.sql` script was run after the database was created."""

from flask import Flask
from flask_script import Manager
from flask_migrate import Migrate, MigrateCommand
from sqlservice import SQLClient
from nephele2.rds import rds_config
from nephele2.rds.db_models import Model, UserEmail, UserInfo, User, MachineInfo, Job, Checkpoints, JobType

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = rds_config.INIT_URI
db = SQLClient({'SQL_DATABASE_URI':rds_config.INIT_URI}, model_class=Model)

migrate = Migrate(app, db, compare_type=True)
manager = Manager(app)

manager.add_command('db', MigrateCommand)


if __name__ == '__main__':
    manager.run()
