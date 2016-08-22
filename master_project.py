from PyQt4 import QtGui
import sys
import gui
import yaml
import main as mn
import numpy as np


class MasterProject(QtGui.QMainWindow, gui.Ui_projectgui):
    def __init__(self, parent=None):
        super(MasterProject, self).__init__(parent) 
        self.setupUi(self)
        #Configuration parameters:
        self.configs = {}
        self.terminations = {}
        self.drange = {}
        self.algorithms = {}
        self.lrange = {}
        self.l2range = {}
        self.alpharange = {}
        self.config_file_path = ''
        self.config_tag = ''
        #Connection functions:
        self.saveconfig.clicked.connect(self.save)
        self.loadconfig.clicked.connect(self.load)
        self.create_term.clicked.connect(self.call_create_term)
        self.create_dmark.clicked.connect(self.call_create_dmark)
        self.apply_algo.clicked.connect(self.call_apply_algo)
        self.gen_seq.clicked.connect(self.call_gen_seq)
        self.analyze.clicked.connect(self.call_analyze_seq)
        self.plot.clicked.connect(self.call_plot)
        self.find_synch_words.clicked.connect(self.call_fsw)

    def save(self):
        self.configs['graph_path'] = str(self.graph_path.text())
        self.configs['lmax'] = int(self.max_l.text())
        self.terminations['old'] = self.old_term.isChecked()
        self.terminations['new'] = self.new_term.isChecked()
        self.terminations['dmark'] = self.dmark_term.isChecked()
        self.terminations['omega'] = self.omega_term.isChecked()
        self.terminations['omega_inverted'] = self.omega_inverted_term.isChecked()
        self.algorithms['dmark'] = self.dmark.isChecked()
        self.algorithms['mk1'] = self.mk1.isChecked()
        self.algorithms['mk2'] = self.mk2.isChecked()
        self.algorithms['mk2_moore'] = self.mk2_moore.isChecked()
        self.algorithms['mk3'] = self.mk3.isChecked()
        self.algorithms['crissis'] = self.crissis.isChecked()
        self.drange['ini'] = int(self.d_ini.text())
        self.drange['end'] = int(self.d_end.text())
        self.lrange['ini'] = int(self.l_ini.text())
        self.lrange['end'] = int(self.l_end.text())
        self.l2range['ini'] = int(self.l2_ini.text())
        self.l2range['end'] = int(self.l2_end.text())
        self.alpharange['ini'] = float(self.alpha_ini.text())
        self.alpharange['end'] = float(self.alpha_end.text())
        self.config_tag = str(self.tag.text())
        self.config_file_path = 'configs/' + self.configs['graph_path'] + '/config_file_' + self.config_tag + '.yaml'
        lrange = range(self.lrange['ini'], self.lrange['end'] + 2, 2)
        drange = range(self.drange['ini'], self.drange['end'] + 1)
        l2range = range(self.l2range['ini'], self.l2range['end'] + 1)
        if self.alpharange['ini'] ==  self.alpharange['end']:
            alpharange = [self.alpharange['ini']]
        else:
            if self.alpharange['ini'] == self.alpharange['end']:
                alpharange = [self.alpharange['ini']]
            else:
                if self.alpharange['end'] == 0.99:
                    alpharange = list(np.arange(self.alpharange['ini'], 1, 0.05))
                    alpharange.append(0.99)
                else:
                    alpharange = list(np.arange(self.alpharange['ini'], self.alpharange['end'] + 0.05, 0.05))
        self.configs['lrange'] = lrange
        self.configs['l2range'] = l2range
        self.configs['alpharange'] = alpharange
        self.configs['algorithms'] = [x for x in self.algorithms.keys() if self.algorithms[x]]
        self.configs['terminations'] = [x for x in self.terminations.keys() if self.terminations[x]]
        self.configs['drange'] = drange
        self.configs['test'] = 'chi-squared'
        synch_path = 'synch_words/' + self.configs['graph_path'] + '/sw.yaml'
        with open(synch_path, 'r') as f:
            self.configs['synch_words'] = yaml.load(f)
        with open(self.config_file_path, 'w') as f:
            yaml.dump(self.configs, f)

    def load(self):
        self.configs['graph_path'] = str(self.graph_path.text())
        self.config_tag = str(self.tag.text())
        self.config_file_path = 'configs/' + self.configs['graph_path'] + '/config_file_' + self.config_tag + '.yaml'

    def call_fsw(self):
        self.fsw_bar.setValue(0)
        fsw_params = {}
        fsw_params['w'] = int(self.fsw_w.text())
        fsw_params['alpha'] = float(self.fsw_alpha.text())
        fsw_params['test'] = self.fsw_test.text()

        p = 'configs/' + self.configs['graph_path'] + '/fsw_params.yaml'
        with open(p, 'w') as f:
            yaml.dump(fsw_params, f)

        mn.main(self.config_file_path, fsw=True, tag=self.config_tag)
        self.fsw_bar.setValue(100)

    def call_create_term(self):
        self.create_term_bar.value = 0
        mn.main(self.config_file_path, terminate=True, tag=self.config_tag)
        self.create_term_bar.value = 100

    def call_create_dmark(self):
        self.create_dmark_bar.value = 0
        mn.main(self.config_file_path, dmark=True, tag=self.config_tag)
        self.create_dmark_bar.value = 100

    def call_apply_algo(self):
        self.apply_algo_bar.value = 0
        mn.main(self.config_file_path, generate=True, tag=self.config_tag)
        self.apply_algo_bar.value = 100

    def call_gen_seq(self):
        seq_len = int(self.seq_len.text())
        self.gen_seq_bar.value = 0
        mn.main(self.config_file_path, gen_seq=True, seq_len=seq_len, tag=self.config_tag)
        self.gen_seq_bar.value = 100

    def call_analyze_seq(self):
        seq_len = int(self.an_seq_len.text())
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
        other_params['K'] = int(self.kld_l.text()) - 1
        other_params['l1'] = int(self.l1m_upto.text()) - 1

        params = {}
        params['to_analyze'] = to_analyze
        params['other_params'] = other_params

        p = 'configs/' + self.configs['graph_path'] + '/params.yaml'
        with open(p, 'w') as f:
            yaml.dump(params, f)

        self.analyze_bar.value = 0
        mn.main(self.config_file_path, an_seq=True, seq_len=seq_len, tag=self.config_tag)
        self.analyze_bar.value = 100

    def call_plot(self):
        params = {}
        params['cond_entropy'] = self.plot_entropy.isChecked()
        params['autocorrelation'] = self.plot_autocorr.isChecked()
        params['kld'] = self.plot_kld.isChecked()
        params['l1metric'] = self.plot_l1m.isChecked()
        params['eval_l'] = int(self.eval_l.text())
        params['upto'] = int(self.plot_upto.text())
        p = 'configs/' + self.configs['graph_path'] + '/plotconfigs.yaml'
        with open(p, 'w') as f:
            yaml.dump(params, f)
        mn.main(self.config_file_path, plot=True, tag=self.config_tag)


def main():
    app = QtGui.QApplication(sys.argv)
    form = MasterProject()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
