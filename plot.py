import probabilisticgraph as pg
import yaml
import matplotlib.pyplot as plt
graph_path = 'vanderpol2'
lrange = range(4,10,2)
l2range = range(1,5)
drange = range(4,10)
alpha = 0.95
states = []
terms = ['dmark', 'omega_inverted']
algos = ['mk1', 'mk2_moore']
labels = ['Mk1, dmark', 'Mk1, $\Omega$', 'Mk2, dmark', 'Mk2, $\Omega$', 'CRISSiS', 'D-Markov']
h = []
k = []
fi = []
tag = 'v_log'

with open('results/' + graph_path + '/cond_entropies/original.yaml', 'r') as f:
    h_base = yaml.load(f)[-1]

for a in algos:
    for t in terms:
        aux = []
        h_aux = []
        for l in lrange:
            ending = '/L' + str(l) + '_alpha' + str(alpha) + '_' + t + '_' + a + '.yaml'
            p = 'graphs/' + graph_path + ending
            g = pg.ProbabilisticGraph(path = p)
            aux.append(len(g.states))
            with open('results/' + graph_path + '/cond_entropies/' + ending, 'r') as f:
                h_aux.append(yaml.load(f)[-1])

        states.append(aux)
        h.append(h_aux)

        with open('results/' + graph_path + '/kld/' + t + '_' + a + '.yaml', 'r') as f:
            k.append(yaml.load(f))

        with open('results/' + graph_path + '/l1metric/' + t + '_' + a + '.yaml', 'r') as f:
            fi.append(yaml.load(f))

aux = []
h_aux = []
for l in l2range:
    ending = '/L_2_' + str(l) + '_alpha' + str(alpha) + '_crissis.yaml'
    p = 'graphs/' + graph_path + ending
    g = pg.ProbabilisticGraph(path=p)
    aux.append(len(g.states))
    with open('results/' + graph_path + '/cond_entropies'+ ending, 'r') as f:
        h_aux.append(yaml.load(f)[-1])

with open('results/' + graph_path + '/kld/crissis.yaml', 'r') as f:
    k.append(yaml.load(f))

with open('results/' + graph_path + '/l1metric/crissis.yaml', 'r') as f:
    fi.append(yaml.load(f))

h.append(h_aux)
states.append(aux)
aux = []
h_aux = []
for d in drange:
    p = 'graphs/' + graph_path + '/dmarkov_d' + str(d) + '.yaml'
    g = pg.ProbabilisticGraph(path=p)
    aux.append(len(g.states))
    with open('results/' + graph_path + '/cond_entropies'+ ending, 'r') as f:
        h_aux.append(yaml.load(f)[-2])

h.append(h_aux)
states.append(aux)

with open('results/' + graph_path + '/kld/dmarkov.yaml', 'r') as f:
    k.append(yaml.load(f))

with open('results/' + graph_path + '/l1metric/dmarkov.yaml', 'r') as f:
    fi.append(yaml.load(f))


plt.clf()
i = 0
for s in states:
    plt.semilogx(s, h[i], marker = 'o', label = labels[i])
    i += 1

plt.axhline(y=h_base, color='k', linewidth = 3, label='Original sequence baseline')
plt.legend(loc='upper right', shadow=False, fontsize='medium')
plt.xlabel('Number of states')
plt.ylabel('Conditional Entropy')
save_path = 'plots/' + graph_path + '/cond_entropies_' + tag + '.png'
plt.savefig(save_path, bbox_inches='tight')
plt.show()

type = ['kld', 'l1metric']
for t in type:
    plt.clf()
    i = 0
    x = k if t == 'kld' else fi
    for s in states:
        plt.semilogx(s, x[i], marker='o', label=labels[i])
        i += 1
    plt.legend(loc='upper right', shadow=False, fontsize='medium')
    plt.xlabel('Number of states')
    ylbl = 'Kullback-Leibler Divergence' if t == 'kld' else '$\Phi$'
    plt.ylabel(ylbl)
    save_path = 'plots/' + graph_path + '/' + t + '_' + tag + '.png'
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
