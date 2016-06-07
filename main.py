#!/usr/bin
import probabilisticgraph as pg
import graphgenerator as gg
import dmarkov as dm
import sequenceanalyzer as sa
import json
import matplotlib.pyplot as plt


def main(config_file, terminate=False, dmark=False, generate=False, gen_seq=False, an_seq=False, plot=False,
         seq_len=10000000, tag='default'):
    with open(config_file, 'r') as f:
        configs = json.load(f)
    graph_path = configs['graph_path']
    terminations = configs['terminations']
    lmax = configs['lmax']
    algorithms = configs['algorithms']
    lrange = configs['lrange']
    alpharange = configs['alpharange']
    drange = configs['drange']
    test = configs['test']
    synch_words = configs['synch_words']
    if terminate:
        terminate_graphs(graph_path, terminations, lrange, lmax, alpharange, test)
    if dmark:
        generate_dmarkov(graph_path, drange, lmax)
    if generate:
        generate_graphs(algorithms, terminations, lrange, alpharange, graph_path, synch_words, test)
    if gen_seq:
        generate_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len)
    if an_seq:
        p = 'configs/' + graph_path + '/params.json'
        with open(p, 'r') as f:
            params = json.load(f)
        analyze_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len,
                          params['to_analyze'], params['other_params'])
    if plot:
        p = 'configs/' +  graph_path + '/plotconfigs.json'
        with open(p, 'r') as f:
            params = json.load(f)
        if params['cond_entropy']:
            plot_entropies(graph_path, algorithms, terminations, drange, lrange, alpharange, params['eval_l'], tag)
        if params['autocorrelation']:
            plot_entropies(graph_path, algorithms, terminations, drange, lrange, alpharange, params['upto'], tag)
        if params['kld']:
            plot_others('kld', graph_path, algorithms, terminations, drange, lrange, alpharange, tag)
        if params['l1metric']:
            plot_others('l1metric', graph_path, algorithms, terminations, drange, lrange, alpharange, tag)


def terminate_graphs(graph_path, terminations, lrange, lmax, alpharange, test):
    g = pg.ProbabilisticGraph([], [])
    for t in terminations:
        for l in lrange:
            for alpha in alpharange:
                p = 'graphs/' + graph_path + '/rtp_L' + str(lmax) + '.json'
                g.open_graph_file(p)
                h = g.expand_last_level(l, t, alpha, test)
                path = 'graphs/' + graph_path + '/rtp_L' + str(l) + '_alpha' + str(alpha) + '_' + t + '.json'
                h.save_graph_file(path)


def generate_graphs(algorithms, terminations, lrange, alpharange, save_path, synch_words, test):
    for t in terminations:
        for l in lrange:
            for alpha in alpharange:
                p1 = 'graphs/' + save_path + '/rtp_L' + str(l) + '_alpha' + str(alpha) + '_' + t + '.json'
                p2 = 'graphs/' + save_path + '/L' + str(l) + '_alpha' + str(alpha) + '_' + t
                g = gg.GraphGenerator(p1, synch_words, p2)
                for algo in algorithms:
                    if algo == 'mk1':
                        g.mk1(test, alpha)
                    elif algo == 'mk2':
                        g.mk2()


def generate_dmarkov(graph_path, drange, lmax):
    g = pg.ProbabilisticGraph([], [])
    for d in drange:
        p = 'graphs/' + graph_path + '/rtp_L' + str(lmax) + '.json'
        g.open_graph_file(p)
        h = dm.DMarkov(g, d)
        path = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.json'
        h.save_graph_file(path)


def generate_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len):
    g = pg.ProbabilisticGraph([], [])
    for algo in algorithms:
        if algo == 'dmark':
            for d in drange:
                p = 'dmarkov_d' + str(d) + '.json'
                path = 'graphs/' + graph_path + '/' + p
                generate_sequences_core(g, graph_path, path, p, seq_len)
        else:
            for t in terminations:
                for l in lrange:
                    for alpha in alpharange:
                        p = 'L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + algo + '.json'
                        path = 'graphs/' + graph_path + '/' + p
                        generate_sequences_core(g, graph_path, path, p, seq_len)


def generate_sequences_core(g, graph_path, path, p, seq_len):
    g.open_graph_file(path)
    seq = g.generate_sequence(seq_len, g.states[0])
    p = 'sequences/' + graph_path + '/len_' + str(seq_len) + '_' + p
    with open(p, 'w') as f:
        json.dump(seq, f)


def analyze_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len, to_analyze, params):
    for algo in algorithms:
        if algo == 'dmark':
            kld = []
            l1 = []
            for d in drange:
                p = 'dmarkov_d' + str(d) + '.json'
                path = 'sequences/' + graph_path + '/len_' +str(seq_len) + '_' + p
                seq_an = sa.SequenceAnalyzer(path)
                kld_step, l1_step = analyze_sequences_core_1(graph_path, p, to_analyze, params, seq_an)
                kld.append(kld_step)
                l1.append(l1_step)
            if to_analyze['kld']:
                k_path = 'results/' +  graph_path + '/kld/' + p
                with open(k_path, 'w') as f:
                    json.dump(kld, f)
            if to_analyze['l1metric']:
                l_path = 'results/' + graph_path + '/l1/' + p
                with open(l_path, 'w') as f:
                    json.dump(l1, f)
        else:
            for t in terminations:
                for l in lrange:
                    for alpha in alpharange:
                        p = 'L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + algo + '.json'
                        path = 'sequences/' + graph_path + '/len_' +str(seq_len) + '_' + p
                        seq_an = sa.SequenceAnalyzer(path)
                        kld_step, l1_step = analyze_sequences_core_1(graph_path, p, to_analyze, params, seq_an)
                        kld.append(kld_step)
                        l1.append(l1_step)
                    if to_analyze['kld']:
                        k_path = 'results/' + graph_path + '/kld/' + p
                        with open(k_path, 'w') as f:
                            json.dump(kld, f)
                    if to_analyze['l1metric']:
                        l_path = 'results/' + graph_path + '/l1/' + p
                        with open(l_path, 'w') as f:
                            json.dump(l1, f)


def analyze_sequences_core_1(graph_path, path, to_analyze, params, seq_an):
    kld = 0
    l1 = 0
    if to_analyze['probabilities']:
        p, alph = seq_an.calc_probs(params['L'])
        p_path = 'results/'+ graph_path + '/probabilities/' + path
        with open(p_path, 'w') as f:
            json.dump([p, alph], f)
    if to_analyze['cond_probabilities']:
        check_probs(seq_an, graph_path, path)
        p_cond = seq_an.calc_cond_probs(params['L']-1)
        p_cond_path = 'results/'+ graph_path + '/probabilities/cond_' + path
        with open(p_cond_path, 'w') as f:
            json.dump(p_cond, f)
    if to_analyze['cond_entropy']:
        check_probs(seq_an, graph_path, path)
        check_cond_probs(seq_an, graph_path, path)
        h = seq_an.calc_cond_entropy(params['L']-1)
        h_path = 'results/'+ graph_path + '/cond_entropies/' + path
        with open(h_path, 'w') as f:
            json.dump(h, f)
    if to_analyze['autocorrelation']:
        a = seq_an.calc_autocorrelation(params['upto'])
        a_path = 'results/' + graph_path + '/autocorrelations/' + path
        with open(a_path, 'w') as f:
            json.dump(a, f)
    if to_analyze['kld']:
        check_probs(seq_an, graph_path, path)
        p = load_reference_probs(graph_path)
        kld = seq_an.calc_kldivergence(p, params['K'])
    if to_analyze['l1metric']:
        check_probs(seq_an, graph_path, path)
        p = load_reference_probs(graph_path)
        l1 = seq_an.calc_l1metric(p, params['l1'])
    return [kld, l1]


def check_probs(seq_an, graph_path, path):
    if not seq_an.probabilities:
        p_path = 'results/'+ graph_path + '/probabilities/' + path
        with open(p_path, 'r') as f:
                p, alph = json.load(f)
                seq_an.probabilities = p
                seq_an.alphabet = alph


def check_cond_probs(seq_an, graph_path, path):
    if not seq_an.conditional_probabilities:
        p_path = 'results/'+ graph_path + '/probabilities/cond_' + path
        with open(p_path, 'r') as f:
                pcond = json.load(f)
                seq_an.conditional_probabilities = pcond


def load_reference_probs(graph_path):
    path = 'results/' + graph_path + '/probabilities/original.json'
    with open(path, 'r') as f:
        p = json.load(f)
    return p[0]


def plot_entropies(graph_path, algorithms, terminations, drange, lrange, alpharange, eval_l, tag):
    path_original = 'results/' + graph_path + '/cond_entropies/original.json'
    with open(path_original, 'r') as f:
        h_original = json.load(f)
    h_base = h_original[eval_l]
    h = []
    states = []
    labels = []
    g = pg.ProbabilisticGraph([], [])
    for algo in algorithms:
        if algo == 'dmark':
            h_dmark = []
            states_dmark = []
            for d in drange:
                h_path = 'results/' + graph_path + '/cond_entropies/dmarkov_d' + str(d) + '.json'
                with open(h_path, 'r') as f:
                    h_eval = json.load(f)
                    h_dmark.append(h_eval[eval_l])
                g_path = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.json'
                g.open_graph_file(g_path)
                states_dmark.append(len(g.states))
            h.append(h_dmark)
            states.append(states_dmark)
            lbl = 'D-Markov, D from ' + str(drange[0]) + ' to ' + str(drange[-1])
            labels.append(lbl)
        else:
            for t in terminations:
                h_term = []
                states_term = []
                for l in lrange:
                    for alpha in alpharange:
                        p = 'L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + algo + '.json'
                        h_path = 'results/' + graph_path + '/cond_entropies/' + p
                        with open(h_path, 'r') as f:
                            h_eval = json.load(f)
                            h_term.append(h_eval[eval_l])
                        g_path = 'graphs/' + graph_path + '/' + p
                        g.open_graph_file(g_path)
                        states_term.append(len(g.states))
                lbl = algo + ', ' + t
                labels.append(lbl)
                h.append(h_term)
                states.append(states_term)
    i = 0
    for entropy in h:
        plt.semilogx(states[i], entropy, marker = 'o', label = labels[i])
        i += 1

    plt.axhline(y=h_base, color='k', linewidth = 3, label='Original sequence baseline')
    plt.legend(loc='upper right', shadow=False, fontsize='medium')
    plt.xlabel('Number of states')
    plt.ylabel('Conditional Entropy')
    save_path = 'plots/' + graph_path + '/cond_entropies_' + tag + '.png'
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def plot_others(kind, graph_path, algorithms, terminations, drange, lrange, alpharange, tag):
    h = []
    states = []
    labels = []
    g = pg.ProbabilisticGraph([], [])
    for algo in algorithms:
        if algo == 'dmark':
            h_dmark = []
            states_dmark = []
            for d in drange:
                h_path = 'results/' + graph_path + '/' + kind + '/dmarkov_d' + str(d) + '.json'
                with open(h_path, 'r') as f:
                    h_dmark.append(json.load(f))
                g_path = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.json'
                g.open_graph_file(g_path)
                states_dmark.append(len(g.states))
            h.append(h_dmark)
            states.append(states_dmark)
            lbl = 'D-Markov, D from ' + str(drange[0]) + ' to ' + str(drange[-1])
            labels.append(lbl)
        else:
            for t in terminations:
                h_term = []
                states_term = []
                for l in lrange:
                    for alpha in alpharange:
                        p = 'L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + algo + '.json'
                        h_path = 'results/' + graph_path + '/' + kind + '/' + p
                        with open(h_path, 'r') as f:
                            h_term.append(json.load(f))
                        g_path = 'graphs/' + graph_path + '/' + p
                        g.open_graph_file(g_path)
                        states_term.append(len(g.states))
                lbl = algo + ', ' + t
                labels.append(lbl)
                h.append(h_term)
                states.append(states_term)
    i = 0
    for value in h:
        plt.semilogx(states[i], value, marker='o', label=labels[i])
        i += 1
    plt.legend(loc='upper right', shadow=False, fontsize='medium')
    plt.xlabel('Number of states')
    if kind == 'l1metric':
        plt.ylabel('L1-Metric')
    elif kind == 'kld':
        plt.ylabel('Kullback-Leibler Divergence')
    save_path = 'plots/' + graph_path + '/' + kind + '_' + tag + '.png'
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def plot_autocorr(graph_path, algorithms, terminations, drange, lrange, alpharange, up_to, tag):
    path_original = 'results/' + graph_path + '/autocorrelations/original.json'
    with open(path_original, 'r') as f:
        h_base = json.load(f)
    h = []
    labels = []
    g = pg.ProbabilisticGraph([], [])
    for algo in algorithms:
        if algo == 'dmark':
            for d in drange:
                h_path = 'results/' + graph_path + '/autocorrelations/dmarkov_d' + str(d) + '.json'
                with open(h_path, 'r') as f:
                    h_eval = json.load(f)
                    h.append(h_eval)
                g_path = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.json'
                g.open_graph_file(g_path)
                lbl = 'D-Markov, D = ' + str(d) + ', ' + str(len(g.states)) + ' states'
                labels.append(lbl)
        else:
            for t in terminations:
                for l in lrange:
                    for alpha in alpharange:
                        p = 'L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + algo + '.json'
                        h_path = 'results/' + graph_path + '/autocorrelations/' + p
                        with open(h_path, 'r') as f:
                            h_eval = json.load(f)
                            h.append(h_eval)
                        g_path = 'graphs/' + graph_path + '/' + p
                        g.open_graph_file(g_path)
                        lbl = algo + ', ' + t + ', ' + str(len(g.states)) + ' states'
                        labels.append(lbl)
    i = 0
    x = range(1, up_to + 1)
    for autocorr in h:
        plt.plot(x, autocorr[1:], marker='o', label=labels[i])
        i += 1

    plt.plot(x, h_base, color='k', linewidth=3, label='Original sequence')
    plt.legend(loc='upper right', shadow=False, fontsize='medium')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation')
    save_path = 'plots/' + graph_path + '/autocorrelations_' + tag + '.png'
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
