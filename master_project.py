from PyQt4 import QtGui
import sys
import gui
import json
import main as mn


class MasterProject(QtGui.QMainWindow, gui.Ui_projectgui):
    def __init__(self, parent=None):
        super(MasterProject, self).__init__(parent)
        self.setupUi(self)


def main():
    app = QtGui.QApplication(sys.argv)
    form = MasterProject()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
