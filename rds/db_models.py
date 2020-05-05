from datetime import datetime
from sqlalchemy import ForeignKey, Column, Integer, String, Boolean, JSON, TIMESTAMP, Enum, FLOAT, Text
from sqlalchemy.orm import relationship
from sqlalchemy import event
from nephele2.rds.orm_base import Model

# NOTE if we want JSON representation of these objects for some reason see:
# http://blog.sampingchuang.com/setup-user-authentication-in-flask/
# SEE http://docs.sqlalchemy.org/en/rel_1_0/orm/tutorial.html

class UserEmail(Model):
    __tablename__ = 'user_email'

    address = Column(String(256), primary_key=True)
    is_confirmed = Column(Boolean, nullable=False, default=0)
    is_bad = Column(Boolean, nullable=False, default=0)
    subscribed = Column(Boolean, nullable=False, default=0)
    subscr_date = Column(TIMESTAMP, nullable=False, default=datetime.now)
    user = relationship("User", uselist=False, back_populates="user_address")

    def __repr__(self):
        return "{UserEmail:{address:'%s', confirmed:'%s', bad:'%s', subscribed:'%s', subscr_date:'%s'}}"\
            % (self.address, self.is_confirmed, self.is_bad, self.subscribed, self.subscr_date)

class UserInfo(Model):
    __tablename__ = 'user_info'

    info_id = Column(Integer, primary_key=True)
    affiliation = Column(String(150))
    referrer = Column(String(30))
    analysis_type = Column(String(150))
    contacted_support = Column(Boolean, nullable=False, default=0)
    affiliation_category = Column(String(30))

    user = relationship("User", uselist=False, back_populates="user_info")

    def __repr__(self):
        return "{UserInfo:{info_id:'%s', affiliation: '%s', affiliation_category: '%s', referrer: '%s', analysis_type: '%s', contacted_support: '%s'}}"\
            % (self.info_id, self.affiliation, self.affiliation_category, self.referrer, self.analysis_type, self.contacted_support)

class User(Model):
    __tablename__ = 'users'

    user_id = Column(Integer, primary_key=True)
    email_address = Column(String(256), ForeignKey('user_email.address'), nullable=False)
    first_name = Column(String(25))
    last_name = Column(String(50))
    info_id = Column(Integer, ForeignKey('user_info.info_id'), nullable=False)
    compute_remaining = Column(FLOAT, default=1800000, nullable=False)

    user_address = relationship("UserEmail", back_populates="user")
    user_info = relationship("UserInfo", back_populates="user")
    user_jobs = relationship("Job", back_populates="user")

    def __repr__(self):
        return "{User:{user_id:'%s', email_address: '%s', first_name: '%s', last_name: '%s', info_id: '%s', compute_remaining: '%s'}}"\
            % (self.user_id, self.email_address, self.first_name, self.last_name, self.info_id, self.compute_remaining)


class MachineInfo(Model):
    __tablename__ = 'machine_info'

    instance_id = Column(String(20), primary_key=True)
    job_id = Column(String(16), ForeignKey('jobs.job_id'), nullable=False)
    instance_type = Column(String(20))
    ami = Column(String(20))

    job = relationship("Job", back_populates="machine_info")

    def __repr__(self):
        return "{MachineInfo:{instance_id: '%s', job_id: '%s', instance_type: '%s', AMI: '%s'}}"\
            % (self.instance_id, self.job_id, self.instance_type, self.ami)

class JobType(Model):
    """
    The job_type table contains a fixed list of all job types that the system knows
    how to run.  This table should be accessed in a read only manner, write
    operations should not be performed on this table.
    """
    __tablename__ = 'job_type'

    type_id = Column(Integer, primary_key=True)
    name = Column(String(56), nullable=False)
    data_type = Column(String(25), nullable=True)
    ami_id = Column(String(25), nullable=True)
    default_instance_type = Column(String(25), nullable=True)
    package = Column(String(25), nullable=True)
    script_name = Column(String(56), nullable=True)

    def __repr__(self):
        return "{JobType: "\
            "{type_id: '%s', name: '%s', data_type: '%s', ami_id: '%s', default_instance_type: '%s', package: '%s', script_name: '%s'}}"\
            % (self.type_id, self.name, self.data_type, self.ami_id, self.default_instance_type, self.package, self.script_name)

class Job(Model):
    __tablename__ = 'jobs'
    job_id = Column(String(16), primary_key=True)
    user_id = Column(Integer, ForeignKey('users.user_id'), nullable=False)
    type_id = Column(Integer, ForeignKey('job_type.type_id'), nullable=True)
    status = Column(Enum('Initializing', 'Pending', 'Pre-Processing', 'Running', 'Failed', 'Succeeded'), nullable=False)
    transferred = Column(Boolean, nullable=False, default=0)
    args = Column(JSON)
    description = Column(String(100))
    error_msg = Column(String(512))
    stack_trace = Column(Text())
    created = Column(TIMESTAMP, nullable=False, default=datetime.now)
    submitted = Column(TIMESTAMP)
    started = Column(TIMESTAMP)
    completed = Column(TIMESTAMP)

    machine_info = relationship("MachineInfo", uselist=False, back_populates="job")
    user = relationship("User", back_populates="user_jobs")
    checkpoints = relationship("Checkpoints")
    job_type = relationship("JobType")

    def __repr__(self):
        return "{Job: {job_id: '%s', args: '%s', description: '%s', error_msg: '%s', created: '%s', submitted: '%s', started: '%s', completed: '%s'}}"\
            %(self.job_id, self.args, self.description, self.error_msg, self.created, self.submitted, self.started, self.completed)

class Checkpoints(Model):
    __tablename__ = 'checkpoints'

    job_id = Column(String(16), ForeignKey('jobs.job_id'), nullable=False, primary_key=True)
    checkpoint = Column(String(100), nullable=False, primary_key=True)
    transition_in = Column(TIMESTAMP, nullable=False, primary_key=True)

    def __repr__(self):
        return "{Checkpoints: {job_id:'%s', checkpoint: '%s', transition_in: '%s'}}"\
            % (self.job_id, self.checkpoint, self.transition_in)

@event.listens_for(Job.completed, 'set')
def receive_set(target, value, old, initiator):  # old, initiator not being used but need to be here
    """
    Updates the compute_remaining for a user when the jobs.completed
    field is set for a job.
    """
    ### Added because started can be not set (only .created is defaulted)
    ### and you get a problem with you try to take None from a timestamp.
    if target.started is None:
        return
    delta = value - target.started
    target.user.compute_remaining = target.user.compute_remaining - delta.total_seconds()
