"""Used session variables, get cleared on clear_session

job_id - The ID for the current job, gets set in upload_data or
transfer_s3_to_efs
previous_jobid - The ID for the job the user wants to re-run
input_from - local or url, tells the system how to let the user add data,
    set on the input_selection page (doesn't get set for job re-submission)
data_type - type of data files the user is uploading, set on one of the data
selection or init pages
data_type_humanized - human readable version of the data_type, set at the
    same time as data_type
map_file - full path to the mapping file, set in upload mapping, re-set in
    checkValidation after file is corrected
job_type - the pipeline that's going to be run, set on one of the pipeline
    selection pages

session variable set by us but not cleared on clear_session
user_id - the identifier for the current user, set on login, removed on logout
"""

import json
import traceback
import os.path
from os import listdir
from functools import wraps
from urllib.parse import urlparse, urljoin
from werkzeug.datastructures import MultiDict
from pygtail import Pygtail
from flask import render_template, redirect, request, url_for, Markup,\
    session, abort, send_from_directory, Response, stream_with_context
from flask_login import logout_user, login_required, current_user
from flask_wtf import FlaskForm
from wtforms import RadioField
from wtforms.validators import DataRequired
from nephele2.nephele import app, models
from nephele2.nephele.forms import forms
from nephele2.nephele.forms.pipeline_form_map import PIPELINE_NAME_TO_WEB_FORM
from nephele2.nephele.forms.registration import login_form as neph_login_form
from nephele2.nephele.forms.registration.registration_form \
    import RegistrationForm
from nephele2.nephele.forms.registration import resend_reg_form
from nephele2.infra.utils.nephele_logger import create_logger
from nephele2.infra.utils.map_validator import map_validator
from nephele2.infra.utils.neph_utils import N2Manager
from nephele2.infra.job import Job
from nephele2 import NepheleError
from nephele2.nephele import config as server_conf
from nephele2.rds.db_utils import DBUtils
from nephele2.infra.utils import fs_utils

LOGGER = create_logger()
CONTACT_US = '<a href="mailto:nephelesupport?Subject=Nephele Help" '\
    'onmouseover="this.href=this.href.replace(\'@@\',\'.\')" '\
    'onclick="this.href=this.href.replace(\'@@\',\'.\')">contact us</a>'


ALLOWED_EXTS = ['fastq', 'fq', 'fastq.gz', 'fq.gz']


@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404


@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html'), 500


@app.errorhandler(410)
def gone_error(e):
    return render_template('410.html'), 410

# TODO: figure out how to call the errorhandler instead of this from
# transfer_wait.html


@app.route('/410/')
@app.route('/410')
def gone_error_route():
    return render_template('410.html')


@app.route('/504/')
@app.route('/504')
def timeout_error_route():
    return render_template('504.html')


@app.route('/')
@app.route('/index')
def index():
    """
    Displays the home page.
    """
    # in case we're trying to start a new job
    clear_session()
    return render_template('index.html')


# This is apparently ununsed right now
@app.route('/warning/')
@app.route('/warning')
def display_warning():
    """
    A call to a simple HTML page that displays a title and a message.

    I set this up with a route so that we can redirect to it from the html
    templates.
    """
    return render_template('warning.html',
                           title="Compute hours exceeded",
                           msg="You've reached your maximum allowed compute \
                           hours (50 hrs.) for this year, and will be unable \
                           to submit any additional jobs.  If you would like \
                           to request additional compute time, please "
                           + CONTACT_US+".")


@app.route('/login/<email_addr>')
@app.route('/login', methods=['GET', 'POST'])
def login(email_addr=None):
    """
    Handles user login using the flask_login library.  The method will
    initially display the login page (we're only requiring a valid email
    for authorization, we don't have anything that requires authentication
    right now), and will attempt to verify the user upon submission of the
    form.
    """
    login_form = neph_login_form.LoginForm()
    if email_addr:
        login_form.email_addr.data = email_addr
    if login_form.validate_on_submit():
        # validate_on_submit is called on submit by WTF
        # checks valid email, good email, compute, exists, etc.
        target = request.args.get('next')
        if not is_safe_url(target):
            return abort(400)
        return redirect(target or url_for('index'))
    return render_template('login.html',
                           title='Sign In',
                           form=login_form)


def login_required_save_post(f):
    """
    Wrapper for the login method that saves any POST data to the session before
    running login. Use this redirect method for any page that requires POST
    data to operate.

    Methods that refer to this login as required, should pop session variables
    'next_form_data' and 'next_form_path' from the session on GET (which is the
    request method login calls).

    from https://stackoverflow.com/questions/32215548/store-post-data-for-use-after-authenticating-with-flask-login
    """
    @wraps(f)
    def decorated(*args, **kwargs):
        if app.login_manager._login_disabled or current_user.is_authenticated:
            # auth disabled or already logged in
            return f(*args, **kwargs)

        # store data before handling login
        session['next_form_data'] = request.form.to_dict(flat=False)
        session['next_form_path'] = request.path
        return app.login_manager.unauthorized()

    return decorated


@app.route('/logout/')
@app.route('/logout')
# @login_required  # if not logged in, this route cannot run.
def logout():
    """
    Clears out the user's session and runs the logout function from the
    flask_login library.
    """
    clear_session()
    logout_user()
    return redirect(url_for('index'))


@app.route('/view_log/<jobid>')
def view_log(jobid=None):
    """
    Displays the contents of the job's log file, if it exists.  The log file
    will not exist for jobs that have not yet started, or that have already
    completed.

    To display addition lines in the file as the job progresses the user will
    have to refresh the page manually, this does not do auto-update, although
    we may want to implement that feature in the future.
    """
    if jobid is None:
        return render_template(
            'warning.html', title="Log not found",
            msg="You must provide a job ID to view a log file.")
    job = Job(jobid)
    logfile = job.log_file
    if not os.path.exists(logfile):
        return render_template(
            'warning.html',
            title="Log not found",
            msg="We couldn't find the log file you requested. \
            If your job is just starting, it can take several minutes \
            for the log file to be created. If your job has been running, \
            it's possible that it has already completed. \
            Please check your email for an email containing the link to \
            your job results or "+CONTACT_US+" if you did not receive \
            the completion email.")
    try:
        offset_file = logfile+".offset"
        if os.path.exists(offset_file):
            os.remove(offset_file)
        return render_template('job_progress.html', log_file=logfile)
    except OSError:
        LOGGER.error("Couldn't remove offset file for log")
        return render_template('job_progress.html', log_file=logfile)
    except Exception:
        LOGGER.error(str(traceback.format_exc()))
        abort(500)


@app.route('/stream')
def streamed_response():
    def generate(log_file):
        for line in Pygtail(log_file):
            yield line.replace('\n', '<br>')
    return Response(stream_with_context(
        generate(request.args.get("log_file"))))


@app.route('/input_selection/')
@app.route('/input_selection', methods=['GET', 'POST'])
@login_required
def input_selection():
    """
    Sets session variable 'input_from' and redirects to the appropriate data
    upload page.
    """
    if request.method == 'POST':
        session['input_from'] = request.form['input_from']
        if request.form['input_from'] == 'local':
            return redirect(url_for('upload_data'))
        return redirect(url_for('upload_url_data'))
    return render_template('upload_selection.html')


@app.route('/data_selection_16s/')
@app.route('/data_selection_16s', methods=['GET', 'POST'])
@login_required
def data_selection_16s():
    """
    Sets session variables 'data_type', 'data_type_humanized'
    """
    select_form = forms.DataTypeSelectionForm()
    if select_form.validate_on_submit():
        session['data_type'] = select_form.data_type.data
        session['data_type_humanized'] = dict(
            select_form.data_type.choices).get(select_form.data_type.data)
        return redirect(url_for('input_selection'))
    return render_template('datatype_selection.html',
                           type="16s", form=select_form)


@app.route('/data_selection_qc/')
@app.route('/data_selection_qc', methods=['GET', 'POST'])
@login_required
def data_selection_qc():
    """
    Sets session variables 'data_type', 'data_type_humanized'
    """
    select_form = forms.DataTypeQCSelectionForm()
    if select_form.validate_on_submit():
        session['data_type'] = select_form.data_type.data
        session['data_type_humanized'] = dict(
            select_form.data_type.choices).get(select_form.data_type.data)
        return redirect(url_for('input_selection'))
    return render_template('datatype_selection.html',
                           type="16s", form=select_form)


@app.route('/init_its/')
@app.route('/init_its', methods=['GET', 'POST'])
@login_required
def init_ITS():
    """
    Sets session variables 'data_type', 'data_type_humanized'
    """
    session['data_type'] = "ITS_PE"
    session['data_type_humanized'] = "ITS Paired End FASTQ"
    return redirect(url_for('input_selection'))


@app.route('/data_selection_wgs/')
@app.route('/data_selection_wgs', methods=['GET', 'POST'])
@login_required
def data_selection_wgs():
    """
    Sets session variables 'data_type', 'data_type_humanized'
    """
    select_form = forms.DataTypeWGSSelectionForm()
    if select_form.validate_on_submit():
        session['data_type'] = select_form.data_type.data
        session['data_type_humanized'] = dict(
            select_form.data_type.choices).get(select_form.data_type.data)
        return redirect(url_for('input_selection'))
    return render_template('datatype_selection.html',
                           type="WGS", form=select_form)


@app.route('/retrieve_data/')
@app.route('/retrieve_data', methods=['GET', 'POST'])
@login_required_save_post
def get_resubmit_data():
    """
    Landing point from index for re-run previous job.
    We'll hit the GET method if we're redirected from the login page.
    """
    # TODO -  get all error messages in one conf file.
    job_not_found_err = "We couldn't find any data for job {}. This happens "\
        "when we are unable to find the data associated with a previous "\
        "run. Maybe it has expired?"
    if request.method == "GET":
        try:
            if session.pop('next_form_path', None) == request.path:
                prev_id = MultiDict(session.pop(
                    'next_form_data')).get('previous_jobid')
                # if we do this outside of the GET/POST methods it causes
                # exceptions on resubmit...but is this safe?
                session['previous_jobid'] = prev_id
        except Exception:
            return redirect(url_for('index'))
    if request.method == "POST":
        prev_id = request.form['previous_jobid']
        job = Job(prev_id)
        session['previous_jobid'] = prev_id
        # check if job data expired
        if not job.exists():
            return render_template('warning.html',
                                   title="Job not found",
                                   msg=job_not_found_err.format(prev_id))
    return render_template('transfer_wait.html')


@app.route('/transfer_data/')
@app.route('/transfer_data', methods=['GET', 'POST'])
@login_required
def transfer_data():
    """
    The POST is automatically triggered from the html template and transfers
    the old job data from S3 into the new job directory in EFS.

    Sets session variable 'job_id', 'data_type', 'data_type_humanized',
    'job_type' Requires session variable 'previous_jobid'
    """
    # TODO: implement error handling
    def transfer():
        job_id = N2Manager.retrieve_prev_job(session['_user_id'],
                                             session['previous_jobid'])
        job_details = DBUtils.get_job_type_details(session['previous_jobid'])
        session['job_id'] = job_id
        session['data_type'] = job_details['data_type']
        session['data_type_humanized'] = job_details['name']
        session['job_type'] = job_details['name']
    return Response(transfer(), mimetype='text/event-stream')


@app.route('/upload_data/')
@app.route('/upload_data', methods=['GET', 'POST'])
@login_required
def upload_data():
    """
    Renders the blueimp multifile upload page.  Redirects to login
    if the user isn't logged-in properly, and creates a job ID if
    one doesn't already exist.

    requires session variables 'job_id', 'data_type', 'data_type_humanized'
    """
    # Create the job ID if we don't have one
    if 'job_id' not in session:
        try:
            session['job_id'] = N2Manager.init_job(session['_user_id'])
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            abort(500)
    next_pg = redirect_from_file_upload()
    if 'previous_jobid' in session:
        title = "Do you want to update any of your files?"
    else:
        title = "Please upload your " + \
            session['data_type_humanized']+" sequence files"
    return render_template('multifile_upload.html',
                           next=next_pg,
                           title=title,
                           allowed_exts=", ".join(ALLOWED_EXTS))

# TODO: implement error handling
# TODO: should this be a login_required route? Will login required break blueimp?
# This is blueimp code, DO NOT change method signature


@app.route('/upload/')
@app.route('/upload', methods=['GET', 'POST'])
def upload():
    """
    Requires session variables 'job_id', 'data_type'
    """
    # TODO: figure out what to do here if we don't have a job ID, this is
    # called for every file

    if request.method == 'POST':
        files = request.files['file']
        # These three lines will print the attributes of the file object
        # (useful for debugging)
        # for attr in dir(files):
        #    if hasattr( files, attr ):
        #        print( "obj.%s = %s" % (attr, getattr(files, attr)))
        # NOTE: This will return false if the file has improper permissions,
        # but we don't get any exceptions being raised.
        # The front end displays an "Unknown error" msg.
        if files:
            inputs_dir = N2Manager.get_job_inputs_d(session['job_id'])
            return models.save_file(inputs_dir, files, ALLOWED_EXTS)

    if request.method == 'GET':
        # get all file in ./data directory
        inputs_dir = N2Manager.get_job_inputs_d(session['job_id'])
        return models.get_file_list(inputs_dir)

    return redirect(url_for('upload_data'))

# To get this to work with CSRF protection, I had to set the X-CSRF-TOKEN
# in the .ajaxSetup headers in the JS.  CSRF token is stored in a meta tag
# in the HTML.  There are some notes on-line about this being unsafe as
# it may pass the token to a 3rd party, but it's also the recommended
# way to do things in the documentation for many libraries that use AJAX.
# This is blueimp, DO NOT change method signature


@app.route("/delete/<string:filename>", methods=['DELETE'])
def delete(filename):
    """
    Implementation of blueimp delete function.

    Requires session variable 'job_id'
    """
    job_base_dir = N2Manager.get_job_inputs_d(session['job_id'])
    return models.delete_file(job_base_dir, filename)

# This is blueimp, DO NOT change method signature
# TODO: look into security here


@app.route("/data/<string:filename>", methods=['GET'])
def get_file(filename):
    """
    Implementation of blueimp get_file function.

    Requires session variable 'job_id'
    """
    inputs_dir = N2Manager.get_job_inputs_d(session['job_id'])
    return send_from_directory(inputs_dir, filename=filename)


@app.route('/upload_url_data/')
@app.route('/upload_url_data', methods=['GET', 'POST'])
@login_required
def upload_url_data():
    """
    Check the remote host to ensure it's accessible.
    Create a file listing

    Requires the session variables job_id, data_type
    Adds the session variables remote_files_list, remote_url
    """
    try:
        form = forms.UploadSeqUrlForm()         # request.POST
    except Exception:
        abort(404)

    # I think if there's no job_id here it's ok to assign one at this point
    if 'job_id' not in session:
        session['job_id'] = N2Manager.init_job(session['_user_id'])

    if form.validate_on_submit():
        try:
            input_url = request.form['data_file']
            inputs_dir = N2Manager.get_job_inputs_d(session['job_id'])
            listing_file = fs_utils.curl_check_url(input_url, inputs_dir)
            if not listing_file:
                LOGGER.error(str(traceback.format_exc()))
                form.data_file.errors.append(
                    "There was an error reaching your remote server. " +
                    str(traceback.format_exc()))
            else:
                session['remote_files_list'] = listing_file
                session['remote_url'] = input_url
            return redirect(redirect_from_file_upload())
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            error = "Unable to reach server. " + str(traceback.format_exc())
            form.data_file.errors.append(error)
    help_1 = "Provide the path to the remote directory where your sequence \
        files are located. Only these extensions <strong>{}</strong> are \
        allowed".format(", ".join(ALLOWED_EXTS))
    help_2 = "We'll check to make sure we can reach the folder for now, \
        and we'll fetch every file in the top level of \
        that folder (we don't do a recursive get) that have the \
        appropriate file extensions when we're ready to run \
        your job."
    help_3 = "Do not add escape characters to the ftp path you're providing. \
        If there's a space in your folder name, just put the space in the \
        path."
    return render_template('upload.html',
                           title="Upload Data",
                           help_text=[help_1, help_2, help_3],
                           form=form)


@app.route("/verify_file_reuse")
@app.route("/verify_file_reuse/", methods=['GET', 'POST'])
@login_required
def check_upload_or_reuse_map():
    """
    Displayed only for job re-submission, after the files are displayed in the
    upload page.  Let's the user indicate if they want to upload new mapping
    (and special) file(s) or if they want to run with the same one(s) from the
    original job.

    Requires the session 'data_type' variable
    """
    if request.method == 'POST':
        try:
            base_dir = N2Manager.get_job_inputs_d(session['job_id'])
            old_args = DBUtils.get_job_arguments(session['previous_jobid'])
            if request.form['input_from'] == 'old':
                # set the map file and other variables here so we have them at
                # the end
                if old_args and 'map_file' in old_args:
                    map_file = models.get_filename_from_path(old_args
                                                             ['map_file'])
                session['map_file'] = base_dir+map_file
                if session['data_type'] == 'DS_Analysis':
                    session['biom_fp'] = base_dir +\
                        models.get_filename_from_path(old_args['biom_fp'])
                return redirect(url_for('start_job'))
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            abort(500)
        # propmt the user for new files if they said they want to change them,
        # or if we can't find the old ones
        # if the user wants to upload new files, continue as usual
        next_pg = url_for('upload_mapping')
        return redirect(next_pg)
    # Display the page
    title = "Do you want to re-use your mapping file, or upload a new one?"
    return render_template('verify_file_reuse.html', title=title)


@app.route('/select_data_file/')
@app.route('/select_data_file', methods=['GET', 'POST'])
@login_required
def select_data_file():
    """
    Page that prompts the user to identify essential, non-sequence files
    needed to run a pipeline, such as an index or qual file.

    The form that is displayed is dynamically rendered because we don't
    know the contents of the radio selection until after the user uploads
    their data, so we can't generate the form at application start. It
    must be defined in this method.
    """
    if 'job_id' not in session:
        # TODO: we may want a more informative message here
        abort(500)
    try:
        job_base_dir = N2Manager.get_job_inputs_d(session['job_id'])
        # define the form first - we have to do this here
        # so it isn't rendered at app start
        # the form needs to be dynamic because we don't
        # know the list of files until we fetch them

        class DynamciForm(FlaskForm):
            """
            This is a dynamically defined WTForm class, so we must create it
            here or it will be rendered once at application start, before we
            know the contents.
            """
            pass

        # define the structure of the form before we instantiate it
        # we have to do this after we create the file class and before
        # we instantiate it
        file_type = "supplementary data"
        session_arg = 'data_file'
        # collect the file names from the input dir
        choices_list = []
        if session['input_from'] == 'url':
            # get list from remote .listing file (generated by wget)
            file_list = models.get_remote_file_list(job_base_dir, ALLOWED_EXTS)
        else:
            file_list = listdir(job_base_dir)
        for filename in file_list:
            choices_list.append((filename, filename))
        # create the RadioField and add it to the form
        # TODO: consider moving field creating into forms.py so we don't have
        # to import
        # the field and validation libs
        setattr(
            DynamciForm,
            'data_file',
            RadioField(
                "Which file is the " +
                file_type +
                " file?",
                choices=choices_list,
                validators=[
                    DataRequired(
                        "Please select the " +
                        file_type +
                        " file.")]))
        # instantiate the form so we can use it
        # we have to do this after we define the form
        form = DynamciForm()

        # If we have a valid form submission, upload the file and add it's name
        # to the session
        if form.validate_on_submit():
            file = form.data_file.data
            job_base_dir = N2Manager.get_job_inputs_d(session['job_id'])
            filepath = job_base_dir+""+file
            session[session_arg] = filepath
            return redirect(url_for('upload_mapping'))
        # This is the default behavior for this page, it displays the generated
        # form
        return render_template(
            'supplementary_file_selection.html', title="Please select the " +
            file_type + " file that you uploaded:",
            help_text=[
                "The quality file is required for FASTA and Quality jobs."],
            form=form)
    except OSError:
        LOGGER.error(str(traceback.format_exc()))
        render_template(
            'warning.html', title="Input file error",
            msg="We either couldn't find the input directory for your job, \
            or we encountered an error while trying to read the data.")
    except Exception:
        LOGGER.error(str(traceback.format_exc()))
        abort(500)


@app.route('/upload_biom/')
@app.route('/upload_biom', methods=['GET', 'POST'])
@login_required
def upload_biom():
    """
    Uploads a biom file from local.

    requires session variables 'job_id', 'input_from', 'data_type_humanized'
    sets session variable 'map_file'
    """

    form = forms.UploadBiomForm()
    # If there's no job id, make one
    if 'job_id' not in session:
        session['job_id'] = N2Manager.init_job(session['_user_id'])

    # If we have a valid form submission, upload the map file and add it's
    # name to the session
    if form.validate_on_submit():
        try:
            file = request.files[form.data_file.name]
            job_base_dir = N2Manager.get_job_inputs_d(session['job_id'])
            filepath = N2Manager.save_user_file(job_base_dir, file)
            session['biom_fp'] = filepath
            session['data_type_humanized'] = 'Downstream Analysis'
            session['data_type'] = 'DS_Analysis'
            session['job_type'] = 'Downstream Analysis'
            return redirect(url_for('upload_mapping'))
        except OSError:
            LOGGER.error(str(traceback.format_exc()))
            form.data_file.errors.append(
                "There was an error uploading your file. " +
                str(traceback.format_exc()))
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            form.data_file.errors.append(
                "There was an error uploading your file. " +
                str(traceback.format_exc()))
    help_text = ['This pipeline accepts <strong>BIOM V1 or QIIME\'s '
                 'BIOM V2 file format only</strong>. All biom files '
                 'from Nephele\'s amplicon pipelines are of these '
                 'formats and should work.']
    return render_template('upload.html',
                           title="Upload Biom File for Downstream Analysis.",
                           form=form,
                           help_text=help_text)


@app.route('/upload_mapping/')
@app.route('/upload_mapping', methods=['GET', 'POST'])
@login_required
def upload_mapping():
    """
    Uploads a mapping file from local.

    requires session variables 'job_id', 'data_type', 'input_from',
        'data_type_humanized'
    sets session variable 'map_file'
    """
    # FIXME this function needs to be refactored.
    form = forms.UploadMapForm()

    map_warn = ('The mapping file should contain all information about the '
                'samples essential to perform the data analysis. In general, '
                'it should contain the name of each sample and the '
                'corresponding sequence data filenames exactly as they have '
                'been uploaded.')
    red_warn = ('<span style="color: red;">If you are uploading an excel '
                'file, please make sure that it contains only one sheet. '
                'Multiple sheet workbooks will not be validated properly.'
                '<span>')
    # FIXME: If we're missing a job_id we should let the user know what's
    # going on
    if 'job_id' not in session:
        if session['input_from'] == 'local':
            return redirect(url_for('upload_data'))
        return redirect(url_for('upload_url_data'))
    # If we have a valid form submission, upload the map file and add it's
    # name to the session
    if form.validate_on_submit():
        try:
            file = request.files[form.data_file.name]
            job_base_dir = N2Manager.get_job_inputs_d(session['job_id'])
            filepath = N2Manager.save_user_file(job_base_dir, file)
            session['map_file'] = filepath
            return redirect(url_for('validate_mapping'))
        except OSError:
            LOGGER.error(str(traceback.format_exc()))
            form.data_file.errors.append(
                "There was an error uploading your file. " +
                str(traceback.format_exc()))
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            form.data_file.errors.append(
                "There was an error uploading your file. " +
                str(traceback.format_exc()))
    map_template = N2Manager.get_map_file_template(session['data_type'])
    if session['data_type'] == 'DS_Analysis':
        se_url = N2Manager.get_map_file_template('SE')
        pe_url = N2Manager.get_map_file_template('PE')
        pe_link = ("<a href='{}' target='_blank' rel='noopener noreferrer'>"
                   "PE</a>".format(pe_url))
        se_link = ("<a href='{}' target='_blank' rel='noopener noreferrer'>"
                   "SE</a>".format(se_url))
        download_link = ('The mapping file format is the same as that used by '
                         'the amplicon pipelines: either {} or '
                         '{} templates. The FASTQ file columns will be '
                         'ignored.'.format(se_link, pe_link))
    else:
        download_link = ("<a href='{}' target='_blank' "
                         "rel='noopener noreferrer'>Download</a> template file"
                         " for instructions.".format(map_template))
    return render_template(
        'upload.html', title="Upload Mapping File for " +
        session['data_type_humanized'],
        info_link=map_template,
        info_link_title="Download mapping file template",
        help_text=[map_warn, download_link, red_warn],
        form=form)


@app.route('/validate_mapping/')
@app.route('/validate_mapping', methods=['GET', 'POST'])
@login_required
def validate_mapping():
    """
    Initializes the map validation unit with either the uploaded map file,
    or the changes made in the validation UI.

    Requires session variables 'data_type', 'map_file', 'job_id'
    Optionally checks session variables 'remote_files_list'
    """
    analysis = models.get_analysis_type(session['data_type'])

    if request.method == 'POST':
        table_data = json.loads(request.form['tableData'])
        headers_list = json.loads(request.form['tableHeaders'])
        col_map = models.get_column_map(table_data)

        try:
            if 'remote_files_list' in session:
                # we're doing a URL upload
                validator = map_validator(
                    session['map_file'],
                    analysis, ALLOWED_EXTS, column_map=col_map,
                    headers_list=headers_list,
                    listing_file=session['remote_files_list'])
            else:
                validator = map_validator(session['map_file'],
                                          analysis,
                                          ALLOWED_EXTS,
                                          column_map=col_map,
                                          headers_list=headers_list)
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            return render_template("warning.html",
                                   title="Parsing Error",
                                   msg="There was a problem constructing your mapping file. \
                                   Please "+CONTACT_US+" and provide job ID " +
                                   session['job_id']+" for reference.")
        return check_validation(validator)
    else:

        if 'map_file' in session and os.path.isfile(session['map_file']):
            try:
                if 'remote_files_list' in session:
                    validator = map_validator(
                        session['map_file'],
                        analysis, ALLOWED_EXTS,
                        listing_file=session['remote_files_list'])
                else:
                    validator = map_validator(
                        session['map_file'], analysis, ALLOWED_EXTS)
            except Exception:
                LOGGER.error(str(traceback.format_exc()))
                # TODO this sucks, make it a template.
                return render_template(
                    "warning.html", title="Parsing Error", msg="There was a problem reading your mapping file. \
                    This is most often caused by an invalid Excel format or a corrupted file. \
                    Here are some suggestions for fixing your file:\
                    <br><br>\
                    <ol>\
                    <li>If you are uploading an Excel file, please ensure that it contains <strong>only one sheet</strong>. \
                    Excel workbooks containing multiple sheets will not parse.</li>\
                    <li>If you are uploading an Excel file, please ensure that it is <strong>saved as \"Excel Workbook (.xlsx)\" \
                    format</strong>. Other .xlsx formats such as \"Strict Open XML Spreadsheet (.xlsx)\" will not parse.</li>\
                    <li>If you are having problems with an Excel file, try saving it as a .csv file. Ensure that the Excel file was not saved \
                    in \"Strict Open XML Spreadsheet (.xlsx)\" format first, as the format may be held over on the text file.</li>\
                    <li>If you still have an issue after trying the above steps, please "+CONTACT_US+" and include your mapping file.</li>\
                    </ol>")
            return check_validation(validator)
        # This is here so that if there is no map_file in the session
        # the user will be requested to upload the map file again.
        # Used to skip the error page if the browser back button is
        # hit on the start_job page.
        return redirect(url_for('upload_mapping'))


@app.route('/job_select/')
@app.route('/job_select', methods=['GET', 'POST'])
@login_required
def select_job():
    """
    Sets session variable 'job_type'
    Requires session variable 'data_type'
    """
    try:
        if request.method == 'POST':
            session['job_type'] = request.form['job_type']
            return redirect(url_for('start_job'))
        data_type = session['data_type']
        return render_template(data_type+"_jobs.html")
    except Exception:
        LOGGER.error(str(traceback.format_exc()))
        if 'data_type' in session:
            msg = "We couldn't find any pipelines for data of type " + \
                session['data_type']
        else:
            msg = "Missing required value 'data_type'.  Cannot display pipeline options \
            if we don't know what we're looking for."
        return render_template(
            'warning.html', title="Error loading pipeline options", msg=msg)


@app.route('/start_job/')
@app.route('/start_job', methods=['GET', 'POST'])
@login_required
def start_job():
    """
    All job inputs have been OK'd and we're ready to submit a job to start.
    Calls into N2Manager.submit_job() which handles this process.

    Requires session variables 'job_type', 'job_id', 'map_file'
    Checks for optional session variables 'previous_job','oligo_file'
    """
    # we've somehow lost something from the session here
    # and can't continue (or we've already submitted the job)
    if 'job_id' not in session:
        return redirect(url_for('index'))
    # now define how to build the form
    job_name = session['job_type']

    form_class = PIPELINE_NAME_TO_WEB_FORM[job_name]
    if 'previous_jobid' in session:
        data = DBUtils.get_job_arguments(session['previous_jobid'])
        form = form_class(data=data)
    else:
        form = form_class()

    if request.method == 'POST' and form.validate_on_submit():
        # get all of the user's argument values and submit the job
        try:
            job_id = session['job_id']
            details_list = models.get_non_job_argslist(forms.JobDetailsForm())
            job_args = {}
            for fieldname, value in form.data.items():
                # skip all of the non-pipeline argument form fields and any
                # non-boolean empty values
                if (fieldname not in details_list and
                        value is not None and
                        value != ''):
                    job_args[fieldname] = value
            # add the map file, any extra files, and the data_type and job_id
            job_args['map_file'] = session['map_file']
            job_args['data_type'] = session['data_type']
            job_args['job_id'] = session['job_id']
            # add ftp flag if files are from remote
            if 'remote_url' in session:
                job_args['ftp'] = session['remote_url']
            if session.get('biom_fp'):
                job_args['biom_fp'] = session.get('biom_fp')
            # if we enable selection in the front end
            # print(job_args)
            N2Manager.submit_job(job_id,
                                 job_name,
                                 form.job_desc.data,
                                 job_args)
            clear_session()

            return redirect(url_for(
                'show_job_submission_success', job_id=job_id))
        except Exception:
            LOGGER.error(str(traceback.format_exc()))
            abort(500)
    page_name = job_name.lower()
    page_name = page_name.replace(" ", "_")
    return render_template(
        page_name+"_opts.html", form=form, job_name=job_name)


@app.route('/submissionSuccess/<string:job_id>')
@login_required
def show_job_submission_success(job_id):
    """
    Renders the job submission success page, which displays the Job ID
    and the email address that communications about the job will be sent to.
    """
    try:
        return render_template("job_submission_success.html",
                               job_id=job_id,
                               email=current_user.get_id())
    except NepheleError.NepheleError:
        LOGGER.error(str(traceback.format_exc()))
        return render_template("job_submission_success.html",
                               job_id=job_id,
                               email=current_user.get_id(),
                               resources=False)


@app.route('/register', methods=['GET', 'POST'])
def register():
    """
    Handles user registration.
    """
    form = RegistrationForm()
    if request.method == 'POST' and form.validate_on_submit():
        user_email = form.email.data
        salted_email = server_conf.salt_string(user_email)
        confirm_url = url_for('confirm_registration',
                              token=salted_email,
                              _external=True)
        try:
            N2Manager.create_user(form.email.data,
                                  form.fname.data,
                                  form.lname.data,
                                  form.affil.data,
                                  form.affil_cat.data,
                                  form.ref.data,
                                  form.subscribe.data,
                                  analysis=form.analysis.data)
            N2Manager.send_registration_email(user_email, confirm_url)
            return redirect(url_for('show_registration_success',
                                    email=salted_email))
        except NepheleError.NepheleDatabaseWriteError as rnf:
            LOGGER.error(rnf)
            return render_template(
                "registration_confirm.html",
                error="Oops!  Looks like something went wrong.")
        except NepheleError.NepheleMissingArgument:
            # We really should never see this, since we're checking that we
            # have the data in the form validation, but it's here just in case
            LOGGER.error(traceback.format_exc())
            return render_template(
                "registration_confirm.html",
                error="Oops!  Looks like something went wrong.")
        except Exception:
            try:
                resend_url = url_for('resend_registration_email',
                                     email=salted_email)
                resend_link = "<a href='{}'>Retry</a>".format(resend_url)
                return render_template(
                    "registration_confirm.html",
                    error="We were unable to send your confirmation email. {}".
                    format(resend_link))
            except Exception:
                LOGGER.error(str(traceback.format_exc()))
                abort(500)
    return render_template('register.html', form=form)


@app.route('/registrationSuccess/<string:email>')
def show_registration_success(email):
    """
    Displays the success template
    """
    try:
        unsalted_email = server_conf.unsalt_string(email)
        return render_template("registration_success.html",
                               email=unsalted_email)
    except Exception:
        LOGGER.error("Couldn't find template for registration_success:")
        LOGGER.error(traceback.format_exc())
        abort(404)

# based on http://exploreflask.com/en/latest/users.html


@app.route("/registrationConfirm/<token>")
def confirm_registration(token):
    """
    Confirms user's registration.
    """
    try:
        user_email = server_conf.unsalt_string(token)
    except Exception:
        LOGGER.error(str(traceback.format_exc()))
        abort(500)
    try:
        N2Manager.confirm_user(user_email, True)
        return render_template("registration_confirm.html", email=user_email)
    except Exception:
        return render_template(
            "registration_confirm.html",
            error="We encountered an unexpected error "
            "while confirming your email. Confirmation failed.")
    return render_template("registration_confirm.html",
                           error="No email address provided as input. "
                           "We can't process your registration.")


@app.route("/registrationResend/<string:email>", methods=['GET', 'POST'])
def resend_registration_email(email):
    """resend_registration_email

    Decrypts the user_email received from the client and resends the
    registration email
    :param email:
    """
    try:
        resend_form = resend_reg_form.ResendEmailForm()
        unsalted_email = server_conf.unsalt_string(email)
        if resend_form.validate_on_submit():
            confirm_url = url_for('confirm_registration',
                                  token=email,
                                  _external=True)
            N2Manager.send_registration_email(unsalted_email, confirm_url)
            # redirect to a Thank you page on success
            return render_template("registration_success.html",
                                   email=unsalted_email)
        return render_template('resend_registration.html',
                               form=resend_form,
                               email=unsalted_email)
    except Exception:
        LOGGER.error("Error re-sending confirmation email")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route("/results")
@app.route("/results/")
@app.route("/results/<jobid>")
def show_results(jobid=None):
    """show_results

    Displays the results page.

    :param jobid:
    """
    if jobid is None:
        err_msg = "Please provide a valid Job ID to retrieve results."
        return render_template('results.html',
                               jobid=jobid,
                               errmsg=err_msg)
    job = Job(jobid)
    if not job.exists():
        err_msg = ("We couldn't find any data for job: {}. Either you "
                   "provided an incorrect job identifier, or the job has "
                   " expired.".format(jobid))
        return render_template('results.html',
                               jobid=jobid,
                               errmsg=err_msg)

    results_url = job.get_results_url()
    log_url = job.get_logfile_url()

    if not results_url or not log_url:
        err_msg = ("You job seems incomplete. Please {} for assistance."
                   .format(Markup(CONTACT_US)))
        return render_template('results.html',
                               jobid=jobid,
                               errmsg=err_msg)

    return render_template('results.html',
                           jobid=jobid,
                           expiration=job.get_job_expiration(),
                           download=results_url,
                           download_size=job.get_results_size(),
                           logfile=log_url)


@app.route('/about')
@app.route('/about/')
def show_about():
    try:
        return render_template("nephele-about-main-page.html")
    except Exception:
        LOGGER.error("Couldn't find template for about page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/about_microbes')
@app.route('/about_microbes/')
def show_about_mirobes():
    try:
        return render_template("nephele-about-microbiome.html")
    except Exception:
        LOGGER.error("Couldn't find template for about microbiome page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/references')
@app.route('/references/')
def show_references():
    try:
        return render_template("references.html")
    except Exception:
        LOGGER.error("Couldn't find template for references page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/about_us')
@app.route('/about_us/')
def show_team():
    try:
        return render_template("about-team.html")
    except Exception:
        LOGGER.error("Couldn't find template for about team page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/about_release')
@app.route('/about_release/')
def show_release_notes():
    try:
        return render_template("release_notes.html")
    except Exception:
        LOGGER.error("Couldn't find template release_notes.html")
        LOGGER.error(traceback.format_exc())
        abort(400)


@app.route('/user_guide')
@app.route('/user_guide/')
def show_user_guide():
    try:
        try:
            seqs_link = N2Manager.get_test_files(
                "pe_user_guide/N2_16S_example_data.zip")
        except Exception:
            LOGGER.error("Couldn't find test data:")
            LOGGER.error(traceback.format_exc())
            seqs_link = None
        try:
            map_link = N2Manager.get_test_files(
                "pe_user_guide/N2_16S_example_mapping_file.xlsx")
        except Exception:
            LOGGER.error("Couldn't find test map file:")
            LOGGER.error(traceback.format_exc())
            map_link = None
        return render_template("user_guide_main.html",
                               test_seqs_url=seqs_link,
                               test_map_url=map_link)
    except Exception:
        LOGGER.error("Couldn't find template for main user guide page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_pipes')
@app.route('/user_guide_pipes/')
def show_pipes_guide():
    try:
        wgs_seqs_link = N2Manager.get_test_files(
            "wgs_user_guide/N2_wgs_example_data.zip")
        wgs_map_link = N2Manager.get_test_files(
            "wgs_user_guide/N2_wgs_example_mapping_file.xlsx")
        pe_map_link = N2Manager.get_map_file_template("PE")
        pe_wgs_map_link = N2Manager.get_map_file_template("WGS_PE")
        qc_pe_map_link = N2Manager.get_map_file_template("QC_PE")
        se_map_link = N2Manager.get_map_file_template("SE")
        se_wgs_map_link = N2Manager.get_map_file_template("WGS_SE")
        qc_se_map_link = N2Manager.get_map_file_template("QC_SE")
        qiime_output_url = N2Manager.get_test_files("qiime_example_output.zip")
        mothur_output_url = N2Manager.get_test_files(
            "mothur_example_output.zip")
        dada_output_url = N2Manager.get_test_files("dada_example_output.zip")
        wgs_output_url = N2Manager.get_test_files("wgs_example_output.zip")
        qc_output_url = N2Manager.get_test_files("qc_example_output.zip")
        da_output_url = N2Manager.get_test_files("da_example_output.zip")
        return render_template("user_guide_pipes.html",
                               wgs_test_seqs_url=wgs_seqs_link,
                               wgs_test_map_url=wgs_map_link,
                               pe_map_url=pe_map_link,
                               se_map_url=se_map_link,
                               wgs_pe_map_url=pe_wgs_map_link,
                               wgs_se_map_url=se_wgs_map_link,
                               qc_pe_map_url=qc_pe_map_link,
                               qc_se_map_url=qc_se_map_link,
                               qiime_output_url=qiime_output_url,
                               mothur_output_url=mothur_output_url,
                               dada_output_url=dada_output_url,
                               wgs_output_url=wgs_output_url,
                               qc_output_url=qc_output_url,
                               da_output_url=da_output_url)
    except Exception:
        LOGGER.error("Couldn't find template for pipelines user guide page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_tutorials')
@app.route('/user_guide_tutorials/')
def show_tutorials():
    try:
        download_url = N2Manager.get_resource_files("datavis_tutorial.zip")
        return render_template("datavis_tutorial.html",
                               download_url=download_url)
    except Exception:
        LOGGER.error("Couldn't find template or data files for visualization \
                     tutorials page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_phyloseq_tutorial')
@app.route('/user_guide_phyloseq_tutorial/')
def show_phyloseq_tutorial():
    try:
        return render_template("phyloseq_tutorial.html")
    except Exception:
        LOGGER.error("Couldn't find template for phyloseq tutorials page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_qiime2_tutorial')
@app.route('/user_guide_qiime2_tutorial/')
def show_qiime2_tutorial():
    try:
        return render_template("qiime2_tutorial.html")
    except Exception:
        LOGGER.error("Couldn't find template for QIIME2 tutorials page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/webinars')
@app.route('/webinars/')
def show_webinars():
    try:
        slides_url = N2Manager.get_resource_files(
            "Nephele_webinar_Nov2018_compress.pdf")
        return render_template("webinars_info.html", slides_url=slides_url)
    except Exception:
        LOGGER.error("Couldn't find template for webinars page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/external_resources')
@app.route('/external_resources/')
def show_external_resources():
    try:
        return render_template("external_resources.html")
    except Exception:
        LOGGER.error("Couldn't find template for external resources page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_visualization_details')
@app.route('/user_guide_visualization_details/')
def show_vis_details():
    try:
        return render_template("datavis16s_pipeline.html")
    except Exception:
        LOGGER.error("Couldn't find template for visualization dtails page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/user_guide_data')
@app.route('/user_guide_data/')
def show_data_comp_guide():
    try:
        return render_template("user_guide_datasets.html")
    except Exception:
        LOGGER.error(
            "Couldn't find template for user guide public data sets page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/privacy_policy')
@app.route('/privacy_policy/')
def show_privacy_policy():
    try:
        return render_template("privacy_policy.html")
    except Exception:
        LOGGER.error("Couldn't find template for privacy policy page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/faq')
@app.route('/faq/')
def show_FAQ():
    try:
        return render_template("faq.html")
    except Exception:
        LOGGER.error("Couldn't find template for FAQ page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/details_mothur')
@app.route('/details_mothur/')
def show_mothur_details():
    try:
        return render_template("mothur_pipeline.html")
    except Exception:
        LOGGER.error(
            "Couldn't find template for mothur pipeline details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/details_qiime')
@app.route('/details_qiime/')
def show_qiime_details():
    try:
        return render_template("qiime_pipeline.html")
    except Exception:
        LOGGER.error("Couldn't find template for QIIME pipeline details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/details_dada2')
@app.route('/details_dada2/')
def show_dada2_details():
    try:
        return render_template("dada2_pipeline.html")
    except Exception:
        LOGGER.error("Couldn't find template for DADA2 pipeline details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/details_biobakery')
@app.route('/details_biobakery/')
def show_biobakery_details():
    try:
        return render_template("biobakerywgs_pipeline.html")
    except Exception:
        LOGGER.error(
            "Couldn't find template for bioBakery pipeline details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/details_qc')
@app.route('/details_qc/')
def show_qc_details():
    try:
        return render_template("qc_details.html")
    except Exception:
        LOGGER.error("Couldn't find template for QC pipeline details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/result_details')
@app.route('/result_details/')
def show_result_details():
    try:
        return render_template("analysis_details.html")
    except Exception:
        LOGGER.error("Couldn't find template for analysis details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/using_ftp')
@app.route('/using_ftp/')
def using_ftp():
    try:
        return render_template("using_ftp.html")
    except Exception:
        LOGGER.error("Couldn't find template for using_ftp page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


@app.route('/da_details')
@app.route('/da_details/')
def show_da_details():
    """route for Downstream analysis details page"""
    try:
        pe_map_url = N2Manager.get_map_file_template("PE")
        se_map_url = N2Manager.get_map_file_template("SE")
        return render_template("amplicon_da_pipeline.html",
                               se_map_url=se_map_url,
                               pe_map_url=pe_map_url)
    except Exception:
        LOGGER.error("Couldn't find template for downstream analysis '\
                     'details page:")
        LOGGER.error(traceback.format_exc())
        abort(404)


def is_safe_url(target):
    """
    A function that checks if the url is safe for redirectes by
    ensuring that a redirect target will lead to the same server.
    See http://flask.pocoo.org/snippets/62/ for an example.
    Must remain in views because it requires the request object.
    """
    ref_url = urlparse(request.host_url)
    test_url = urlparse(urljoin(request.host_url, target))
    return test_url.scheme in ('http', 'https') and \
        ref_url.netloc == test_url.netloc


def clear_session():
    """
    Clears all non-login, non-csrf security keys from the session.
    NOTE: this must be in the views.py even though it doesn't have
    a route, because it needs access to the session object.
    """
    required_keys = ['_user_id', '_fresh', '_flashes', '_id', 'csrf_token']
    session_keys = [x for x in session.keys() if x not in required_keys]
    for key in session_keys:
        session.pop(key, None)


def redirect_from_file_upload():
    """
    Determines the correct redirect from a data file upload page.  Since we now
    have several ways to upload data files and they all have to use the same
    logical redirects, they should all call this function.

    NOTE: Because this function accesses the session, it must
    be in this file, or we have to pass the session into the function.
    """
    # use this to skip map file upload if adding files to correct an error,
    # we already have the map file
    # redirect to pg asking to reuse mapping file on job resubmit
    # go to extra file selection page for data types that use them
    # else, ask for the map file
    if 'map_errors' in session:
        next_pg = url_for('validate_mapping')
    elif 'previous_jobid' in session:
        next_pg = url_for('check_upload_or_reuse_map')
    else:
        next_pg = url_for('upload_mapping')
    return next_pg


def check_validation(validator):
    """
    Runs the map file validator and takes the appropriate action
    based on the validation results and the data type.

    Updates session variable 'map_file'
    Sets session variable 'map_errors' if mapping file has errors else clear it
    Checks for session variable 'previous_jobid'
    remote_files_list
    """
    try:
        response = validator.validate_mapping()
        validation_result = json.loads(response)
        session['map_file'] = validation_result['corrected_file']
    except Exception:
        LOGGER.error(str(traceback.format_exc()))
        abort(500)
    if validation_result['is_valid']:
        # remove flag, use this to skip the map upload page when adding files
        # to correct validation errors
        session.pop('map_errors', None)
        # clean up the remote file listings if remote upload
        # remove remote_files_list from session
        remote_list_file = session.pop('remote_files_list', None)
        # remove .listing from dir
        if remote_list_file:
            models.delete_file(N2Manager.get_job_inputs_d(session['job_id']),
                               os.path.basename(remote_list_file))
        if not validation_result['Warnings']:
            # NOTE: if this if/else is modified, modify the Continue button in
            # map_file_errors.html as well
            if 'previous_jobid' in session:
                return redirect(url_for('start_job'))  # TODO: test this
            return redirect(url_for('select_job'))
    else:
        # set flag, we use this to skip the map upload page when adding files
        # to correct validation errors
        session['map_errors'] = True
    map_template = N2Manager.get_map_file_template(session['data_type'])
    return render_template('map_file_errors.html',
                           template=map_template,
                           is_valid=validation_result['is_valid'],
                           errors=validation_result['Errors'],
                           warnings=validation_result['Warnings'],
                           table=validator.render_table())
