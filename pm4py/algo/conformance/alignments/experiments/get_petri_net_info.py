import os

from pm4py.objects import petri


def get_petri_net_statistics(filepath, filename):
    net, marking, fmarking = petri.importer.pnml.import_net(os.path.join(filepath, filename))
    print(filename)

    number_silent_transitions = 0
    for t in net.transitions:
        if not t.label:
            number_silent_transitions += 1

    print("Number places: ", len(net.places))
    print("Number total transitions: ", len(net.transitions))
    print("Number silent transitions: ", number_silent_transitions)
    print("Number arcs: ", len(net.arcs), "\n")


def bpi2019():
    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "data", "bpi_challenge_2019")
    petri_nets = ['petri_net_1.pnml', 'petri_net_2.pnml', 'petri_net_3.pnml', 'petri_net_4.pnml', 'petri_net_5.pnml']
    for net in petri_nets:
        get_petri_net_statistics(path_to_files, net)


if __name__ == '__main__':
    bpi2019()
