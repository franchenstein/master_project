#!/usr/bin
import probabilisticgraph as pg
import graphgenerator as gg
import dmarkov as dm
import sequenceanalyzer as sa
import json


def terminate_graphs(graph_path, terminations, lrange, alpharange, save_path, test):
    g = pg.ProbabilisticGraph([], [])
    for t in terminations:
        for l in lrange:
            for alpha in alpharange:
                g.open_graph_file(graph_path)
                h = g.expand_last_level(l, t, alpha, test)
                path = 'graphs/' + save_path + '/rtp_L' + str(l) + '_alpha' + str(alpha) + '_' + t + '.json'
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
                    else:
                        g.mk2()


def generate_dmarkov(graph_path, drange, save_path):
    g = pg.ProbabilisticGraph([], [])
    for d in drange:
        g.open_graph_file(graph_path)
        h = dm.DMarkov(g, d)
        path = 'graphs/' + save_path + '/dmarkov_d' + str(d) + '.json'
        h.save_graph_file(path)


def generate_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len):
    g = pg.ProbabilisticGraph([], [])
    for algo in algorithms:
        if algo == 'dmark':
            for d in drange:
                path = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.json'
                generate_sequences_core(g, path, seq_len)
        else:
            for t in terminations:
                for l in lrange:
                    for alpha in alpharange:
                        path = 'graphs/' + graph_path + '/L' + str(l) + '_alpha' + str(alpha) + \
                               '_' + t + '_' + algo + '.json'
                        generate_sequences_core(g, path, seq_len)


def generate_sequences_core(g, path, seq_len):
    g.open_graph_file(path)
    seq = g.generate_sequence(seq_len, g.states[0])
    p = 'sequences/len_' + str(seq_len) + '_' + path[7:]
    with open(p, 'w') as f:
        json.dump(seq, f)


def analyze_sequences(graph_path, algorithms, drange, terminations, lrange, alpharange, seq_len, to_analyze, params):
    for algo in algorithms:
        if algo == 'dmark':
            kld = []
            l1 = []
            for d in drange:
                path = 'sequences/len_' +str(seq_len) + '_' + graph_path + '/dmarkov_d' + str(d) + '.json'
                seq_an = sa.SequenceAnalyzer(path)
                kld, l1 = analyze_sequences_core_1(path, to_analyze, params, seq_an)
        else:
            for t in terminations:
                for l in lrange:
                    for alpha in alpharange:
                        path = 'sequences/len_' +str(seq_len) + '_' '/L' + str(l) + '_alpha' +\
                               str(alpha) + '_' + t + '_' + algo + '.json'


def analyze_sequences_core_1(path, to_analyze, params, seq_an):
    if to_analyze['probabilities']:
        p, alph = seq_an.calc_probs(params['L'])
        p_path = 'results/probabilities_' + path[10:]
        with open(p_path, 'w') as f:
            json.dump([p, alph], f)
    if to_analyze['cond_probabilities']:
        if not seq_an.probabilities:
            p_path = 'results/probabilities_' + path[10:]
            with open(p_path, 'r') as f:
                p, alph = json.load(f)
                seq_an.probabilities = p
                seq_an.alphabet = alph
        p_cond = seq_an.calc_cond_probs(params[L])
        p_cond_path = 'results/cond_probabilities_' + path[10:]
        with open(p_cond_path, 'w') as f:
            json.dump(p_cond, f)
    if to_analyze['cond_entropy']:
        if not seq_an.probabilities:
            p_path = 'results/probabilities_' + path[10:]
            with open(p_path, 'r') as f:
                p, alph = json.load(f)
                seq_an.probabilities = p
            if not seq_an.conditional_probabilities:
                p_cond_path = 'results/cond_probabilities_' + path[10:]
                with open(p_cond_path], 'r') as f:
                    seq_an.conditional_probabilities = json.load(f)
        h = seq_an.calc_cond_entropy(params[L])
    if to_analyze['autocorrelation']:
        a = seq_an.calc_autocorrelation()
    if to_analyze['kld']:
        if not seq_an.probabilities:
            kld.append(seq_an.calc_kldivergence())
    if to_analyze['l1metric']:
        if not seq_an.probabilities:
            l1.append(seq_an.calc_l1metric())
    return kld, l1