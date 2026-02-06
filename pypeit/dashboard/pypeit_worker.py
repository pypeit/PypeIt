
# this is the worker process and it will process the logs



import logging
from logging.handlers import QueueHandler
from pypeit_wrapper import PypeItWrapper
from pypeit import log

def PypeItWorker(setup_file_path,log_queue,msg_queue):

    log.addHandler(QueueHandler(log_queue))

    # Run the worker
    return PypeItWrapper(setup_file_path,msg_queue).run()

