from PyQt4 import QtGui
import sys
import gui
import json
import main as mn
import numpy as np


class MasterProject(QtGui.QMainWindow, gui.Ui_projectgui):
    def __init__(self, parent=None):
        super(MasterProject, self).__init__(parent)
        self.setupUi(self)
        #Configuration parameters:
        self.configs['graph_path'] = ''
        self.terminations['old'] = False
        self.terminations['new'] = False
        self.terminations['dmark'] = False
        self.algorithms['mk1'] = False
        self.algorithms['mk2'] = False
        self.algorithms['dmark'] = False
        self.dmark['ini'] = 0
        self.dmark['end'] = 0
        self.lrange['ini'] = 0
        self.lrange['end'] = 0
        self.alpharange['ini'] = 0.0
        self.alpharange['end'] = 0.0
        self.configs['tag'] = ''
        self.config_file_path = ''
        self.configs['lmax'] = 0
        #Connection functions:
        self.saveconfig.clicked.connect(self.save)
        self.create_term.clicked.connect(self.call_create_term)
        self.create_dmark.clicked.connect(self.call_create_dmark)
        self.apply_algo.clicked.connect(self.call_apply_algo)
        self.gen_seq.clicked.connect(self.call_gen_seq)
        self.analyze.connect(self.call_analyze_seq)
        self.plot.connect(self.call_plot)

    def save(self):
        self.configs['graph_path'] = self.graph_path.text()
        self.configs['lmax'] = int(self.lmax.text())
        self.terminations['old'] = self.old_term.isChecked()
        self.terminations['new'] = self.new_term.isChecked()
        self.terminations['dmark'] = self.dmark_term.isChecked()
        self.algorithms['dmark'] = self.dmark.isChecked()
        self.algorithms['mk1'] = self.mk1.isChecked()
        self.algorithms['mk2'] = self.mk2.isChecked()
        self.dmark['ini'] = int(self.d_ini.text())
        self.dmark['end'] = int(self.d_end.text())
        self.lrange['ini'] = int(self.l_ini.text())
        self.lrange['end'] = int(self.l_end.text())
        self.alpharange['ini'] = float(self.alpha_ini.text())
        self.alpharange['end'] = float(self.alpha_end.text())
        self.configs['tag'] = self.tag.text()
        self.config_file_path = self.params['graph_path'] + '/configs/config_file_' + self.params['tag'] + '.json'
        lrange = range(self.lrange['ini'], self.lrange['end'] + 2, 2)
        drange = range(self.dmark['ini'], self.dmark['end'] + 1)
        if self.alpharange['end'] == 0.99:
            alpharange = list(arange(self.alpharange['ini'], 0.95, 0.05))
            alpharange.append(0.99)
        else:
            alpharange = list(arange(self.alpharange['ini'], self.alpharange['end'], 0.05))
        self.configs['lrange'] = lrange
        self.configs['alpharange'] = alpharange
        self.configs['algorithms'] = [x for x in self.algoritms.keys() if self.alogorithms[x]]
        self.configs['terminations'] = [x for x in self.terminations.keys() if self.terminations[x]]
        self.configs['drange'] = drange
        self.configs['test'] = 'chi-squared'
        synch_path = self.configs['graph_path'] + '/synch_words/' + self.configs['tag'] + '.json'
        with open(synch_path, 'r') as f:
            self.configs['synch_words'] = json.load(f)
        with open(self.config_file_path, 'w') as f:
            json.dump(self.configs, f)

    def call_create_term(self):
        self.create_term_bar.value = 0
        mn.main(self.config_file_path, terminate=True)
        self.create_term_bar.value = 100

    def call_create_dmark(self):
        self.create_dmark_bar.value = 0
        mn.main(self.config_file_path, dmark=True)
        self.create_dmark_bar.value = 100

    def call_apply_algo(self):
        self.apply_algo_bar.value = 0
        mn.main(self.config_file_path, generate=True)
        self.apply_algo_bar.value = 100

    def call_gen_seq(self):
        seq_len = int(self.seq_len.value())
        self.gen_seq_bar.value = 0
        mn.main(self.config_file_path, gen_seq=True, seq_len=seq_len)
        self.gen_seq_bar.value = 100

    def call_analyze_seq(self):
        seq_len = int(self.an_seq_len.value())
        to_analyze = {}
        to_analyze['probabilities'] = self.probs.isChecked()
        to_analyze['cond_probabilities'] = self.cond_probs.isChecked()
        to_analyze['cond_entropy'] = self.cond_entropy.isChecked()
        to_analyze['autocorrelation'] = self.autocorr.isChecked()
        to_analyze['kld'] = self.kld.isChecked()
        to_analyze['l1metric'] = self.l1m.isChecked()

        other_params = {}
        other_params['L'] = int(self.probs_l.text())
        other_params['upto'] = int(self.autocorr_upto.text())
        other_params['K'] = int(self.kld_l.text())
        other_params['l1'] = int(self.l1m_upto.text())

        params = {}
        params['to_analyze'] = to_analyze
        params['other_params'] = other_params

        p = self.configs['graph_path'] + '/configs/params.json'
        with open(p, 'w') as f:
            json.dump(params, f)

        self.analyze_bar.value = 0
        mn.main(self.config_file_path, an_seq=True, seq_len=seq_len)
        self.analyze_bar.value = 100

    def call_plot(self):
        params = {}
        params['cond_entropy'] = self.plot_entropy.isChecked()
        params['autocorrelation'] = self.plot_autocorr.isChecked()
        params['kld'] = self.plot_kld.isChecked()
        params['l1metric'] = self.plot_l1m.isChecked()
        params['eval_l'] = int(self.eval_l.text())
        params['up_to'] = int(self.plot_upto.text())
        p = self.configs['graph_path'] + '/configs/plotconfigs.json'
        with open(p, 'w') as f:
            json.dump(params, f)
        mn.main(self.config_file_path, plot=True)


def main():
    app = QtGui.QApplication(sys.argv)
    form = MasterProject()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
