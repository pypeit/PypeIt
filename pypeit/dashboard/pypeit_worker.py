
# this is the worker process and it will process the logs



import logging
from logging.handlers import QueueHandler
from pypeit_wrapper import PypeItWrapper

def PypeItWorker(setup_file_path,log_queue):

    # Configure logging
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    # Remove any handlers inherited from parent
    # this is apparently needed for macos or windows 
    for h in root.handlers[:]:
        root.removeHandler(h)

    root.addHandler(QueueHandler(log_queue))

    # Run the worker
    return PypeItWrapper(setup_file_path).run()

