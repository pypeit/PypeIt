import sys
import zmq
from PyQt6.QtWidgets import QApplication
from pypeit.setup_gui.controller import SetupGUIController

# start a zmq context and bind it to a socket
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:5555")

print("socket binded") # in the future I will have it so that the dashboard continuously send signals until it gets a reply which will then tell the dashboard that the socket is binded

#define the start controller function
app = QApplication(sys.argv)
verbosity = 1 # in the future I will store the verbosity in its own file and then read the verbosity from that file maybe
controller = SetupGUIController(app,verbosity)

def open_setup():
    controller.start()

"""
TODO: add functions depending on what the dashboard sends via the socket

need:
    - a function that will take in the pypeit setup file and edit the controllers data to set it up as such. (this is sent to it via the controller)
    - a function that will take in the science files and do the same as above
    - a function that will take in the run all instruction and "click" the run button which will automatically run all
        = this function should needs to publish the state of the running if at all possible
        = this is especially true for the "run" button on the dashboard which needs to see which step it is on

NOTE: edit setup will change the pypeit setup file from the dashboard (so open up a prompt for the user) and also change the science files
    - this will be done within dashboard.py and NOT within controller_runner.py

BUGS:
    - this isn't really pertinent to the controller but the dashboard will stop responding when you click close after you have opened the setup
        = this could be maybe because the socket is still open and it needs to be closed?
"""

def obtain_setup_file(file_path):
    # open file and obtain contents from file path
    setup_file = file_path
    return setup_file

def obtain_science_or_standard_file(file_path): # need to change the name of the function because it is kind of icky
    # obtains the science or standard file by reading the contents from the file path
    pass

def change_spectrograph(setup_file):
    # changes the meta data of the controller so it uses the setup_file to start
    pass

def change_raw_data_file(file_path):
    # changes the metadata of the controller to have the raw data file
    pass

def run_setup(): # should change this name to be more descriptive
    # this will call the function withing the controller that will run
    pass

valid_requests = [b'open setup',b'edit setup',b'run all']

def close_socket():
    socket.close()
    context.term()

while True:
    # continually recieve messages
    message = socket.recv()
    print(message)
    # message = message.decode('utf-8')
    if message == b'open setup':  #TODO change this so that it reflects the valid requests list
        open_setup()
        socket.send(b'setup opened')

    elif message not in valid_requests:
        socket.send(b'f"{message} not recognized"')
    print(message)

