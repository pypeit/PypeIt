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

