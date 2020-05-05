from flask import Flask
from flask_bootstrap import Bootstrap
from flask_wtf.csrf import CSRFProtect
from nephele2.nephele.forms.registration import login_manager
from nephele2.config import SHA_LONG, SHA_SHORT, NEPHELE_VERSION
# these imports have to be done after we define app
# because they import app
from nephele2.infra.utils import timer

# Set default logging handler to avoid "No handler found" warnings.
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

app = Flask(__name__)
@app.context_processor
def inject_git_sha():
    return dict(git_sha=SHA_LONG, git_sha_short=SHA_SHORT, tag=NEPHELE_VERSION)

app.config.from_object('nephele2.nephele.config')
login_manager.init_app(app)
csrf = CSRFProtect(app)
bootstrap = Bootstrap(app)
app.email_timer = timer.Timer()

# this import needs to be down here to prevent circ imports
from nephele2.nephele import views

