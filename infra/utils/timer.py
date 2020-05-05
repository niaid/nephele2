from threading import Lock
from datetime import datetime, timedelta
lock = Lock()

MINS_GAP = 60

class Timer:
    def __init__(self):
        self.next_tick = datetime.now()

    def ready_to_send(self):
        with lock:
            if datetime.now() > self.next_tick:
                self.next_tick = datetime.now() + timedelta(minutes=MINS_GAP)
                return True
            else:
                return False

