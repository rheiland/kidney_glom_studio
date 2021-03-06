"""
Authors:
Randy Heiland (heiland@iu.edu)
Adam Morrow, Michael Siler, Grant Waldrow, Drew Willis, Kim Crevecoeur
Dr. Paul Macklin (macklinp@iu.edu)

"""

import sys
import os
import time
from pathlib import Path
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from PyQt5 import QtCore, QtGui
# from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QFormLayout,QLineEdit, QHBoxLayout,QVBoxLayout,QRadioButton,QLabel,QCheckBox,QComboBox,QScrollArea, QPushButton,QPlainTextEdit

from PyQt5.QtCore import QProcess

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class RunModel(QWidget):
    def __init__(self):
        super().__init__()

        #-------------------------------------------
        self.vis_tab = None

        # following set in studio.py
        self.homedir = ''   
        self.config_file = None
        self.tree = None

        self.config_tab = None
        self.microenv_tab = None
        self.celldef_tab = None
        self.user_params_tab = None

        self.sim_output = QWidget()
        self.main_layout = QVBoxLayout()

        self.scroll = QScrollArea()

        self.p = None
        # self.xmin = 0.0
        # self.xmax = 1.0
        # self.ymin = 0.0
        # self.ymax = 1.0

        self.control_w = QWidget()

        self.vbox = QVBoxLayout()

        #------------------
        hbox = QHBoxLayout()

        self.run_button = QPushButton("Run Simulation")
        self.run_button.setStyleSheet("background-color: lightgreen")
        hbox.addWidget(self.run_button)
        self.run_button.clicked.connect(self.run_model_cb)

        self.cancel_button = QPushButton("Cancel")
        # self.cancel_button.setStyleSheet("background-color: red")
        self.cancel_button.setStyleSheet("background-color: rgb(250,100,100)")
        hbox.addWidget(self.cancel_button)
        self.cancel_button.clicked.connect(self.cancel_model_cb)

        # self.cancel_button = QPushButton("Cancel")
        # hbox.addWidget(self.cancel_button)
        # self.new_button.clicked.connect(self.append_more_cb)

        hbox.addWidget(QLabel("Exec:"))
        self.exec_name = QLineEdit()
        self.exec_name.setText('mymodel')
        # self.exec_name.setText('biorobots')
        hbox.addWidget(self.exec_name)

        hbox.addWidget(QLabel("Config:"))
        self.config_xml_name = QLineEdit()
        self.config_xml_name.setText('mymodel.xml')
        # self.config_xml_name.setText('copy_PhysiCell_settings.xml')
        hbox.addWidget(self.config_xml_name)

        # self.vbox.addStretch()

        self.vbox.addLayout(hbox)

        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.resize(400,900)  # nope

        self.vbox.addWidget(self.text)

        #==================================================================
        self.control_w.setLayout(self.vbox)

        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setWidgetResizable(True)

        self.scroll.setWidget(self.control_w) 
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.scroll)

#------------------------------
    def update_xml_from_gui(self):
        self.xml_root = self.tree.getroot()
        self.config_tab.xml_root = self.xml_root
        self.microenv_tab.xml_root = self.xml_root
        self.celldef_tab.xml_root = self.xml_root
        self.user_params_tab.xml_root = self.xml_root

        self.config_tab.fill_xml()
        self.microenv_tab.fill_xml()
        self.celldef_tab.fill_xml()
        self.user_params_tab.fill_xml()

        
    def message(self, s):
        self.text.appendPlainText(s)

    def run_model_cb(self):
        print("===========  run_model_cb():  ============")

        auto_load_params = True
        # if auto_load_params:

        if self.vis_tab:
            # self.vis_tab.reset_axes()
            self.vis_tab.reset_model_flag = True

        # for f in Path('./output').glob('*.*'):
        #     try:
        #         f.unlink()
        #     except OSError as e:
        #         print("Error: %s : %s" % (f, e.strerror))
        print("  deleting /output")
        os.system('rm -rf output/*')
        time.sleep(1)

        # if os.path.isdir('tmpdir'):
        #     # something on NFS causing issues...
        #     tname = tempfile.mkdtemp(suffix='.bak', prefix='output_', dir='.')
        #     shutil.move('output', tname)
        # os.makedirs('output')

        self.update_xml_from_gui()

        print("\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # print("run_tab.py: ----> writing modified model to ",new_config_file)
        # self.tree.write(new_config_file)  # saves modified XML to tmpdir/config.xml 
        # self.tree.write("model.xml")  # saves modified XML to tmpdir/config.xml 

        self.tree.write(self.config_file)  # saves modified XML to tmpdir/config.xml 

        # update axes ranges on Plots


        if self.p is None:  # No process running.
            # self.tree = ET.parse(self.config_file)
            self.tree = ET.parse("mymodel.xml")
            # tree = ET.parse(read_file)
            # self.tree = ET.parse(read_file)
            self.xml_root = self.tree.getroot()

            self.message("Executing process")
            self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.stateChanged.connect(self.handle_state)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.
            # self.p.start("mymodel", ['biobots.xml'])
            exec_str = self.exec_name.text()
            xml_str = self.config_xml_name.text()
            print("running: ",exec_str,xml_str)
            self.p.start(exec_str, [xml_str])
            # self.p = None  # No, don't do this
        else:
            print("self.p is not None???")

    def cancel_model_cb(self):
        print("===========  cancel_model_cb():  ============")
        if self.p:  # process running.
            self.p.kill()

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def handle_state(self, state):
        states = {
            QProcess.NotRunning: 'Not running',
            QProcess.Starting: 'Starting',
            QProcess.Running: 'Running',
        }
        state_name = states[state]
        self.message(f"State changed: {state_name}")

    def process_finished(self):
        self.message("Process finished.")
        print("-- process finished.")
        self.p = None