#!/usr/bin
import sequenceanalyzer as seqan
import yaml
import sys
import getopt


def main(graph_path, length, l):
    sequence_path = 'sequences/' + graph_path + '/original_length_' + str(length) + '.yaml'
    sa = seqan.SequenceAnalyzer(sequence_path)
    p, alph = sa.calc_probs(l)
    p_cond = sa.calc_cond_probs(l - 1)
    h = sa.calc_cond_entropy(l - 1)
    sa.save_cond_probs_as_graph('graphs/' + graph_path + '/rtp_L' + str(l) + '.yaml')
    with open('results/' + graph_path + '/probabilities/original.yaml', 'w') as f:
        yaml.dump([p, alph], f)
    with open('results/' + graph_path + '/probabilities/original_cond.yaml', 'w') as f:
        yaml.dump(p_cond, f)
    with open('results/' + graph_path + '/cond_entropies/original.yaml', 'w') as f:
        yaml.dump(h, f)


def read_input(argv):
    graph_path = ''
    length = 10000000
    l = 0
    try:
        opts, args = getopt.getopt(argv, "hg:x:l:", ["graphpath=", "length=", "l="])
    except getopt.GetoptError:
        print "setuporiginals.py -g <graphpath> -x <length> -l <l>"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print "setuporiginals.py -g <graphpath> -x <length> -l <l>"
            sys.exit()
        if opt in ("-g", "--graphpath"):
            graph_path = arg
        if opt in ("-x", "--length"):
            length = arg
        if opt in ("-l", "--l"):
            l = arg
    return [graph_path, length, l]

if __name__ == "__main__":
    graph_path, length, l = read_input(sys.argv[1:])
    main(graph_path, length, l)
