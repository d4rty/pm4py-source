from pm4py.util.create_artificial_event_log import create_xes_string
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.objects.petri.petrinet import PetriNet
from pm4py.objects.petri import utils as petri_net_utils
from pm4py.visualization.petrinet import factory as petri_net_visualization_factory
from pm4py.objects.petri.petrinet import Marking
from pm4py.algo.conformance.alignments import factory
from pm4py.algo.conformance.alignments.utils import print_alignments


def execute_script(show_petri_nets):
    petri_net_1 = PetriNet("running_example")
    places = {}
    for i in range(1, 4):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        petri_net_1.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),

    }

    for transition in transitions:
        # print(transition)
        petri_net_1.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net_1)
    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net_1)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_B'], petri_net_1)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_C'], petri_net_1)

    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], petri_net_1)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_3'], petri_net_1)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_3']] = 1

    # petri net visualization
    if show_petri_nets:
        parameters = {"format": "svg", 'debug': True}
        # parameters = None
        gviz = petri_net_visualization_factory.apply(petri_net_1, initial_marking, final_marking, parameters=parameters)
        petri_net_visualization_factory.view(gviz)

    traces = [
        {"frequency": 1, "events": ["A", "A", "B"]}
    ]

    event_log = xes_importer.import_log_from_string(create_xes_string(traces))

    alignments = factory.apply(event_log[0], petri_net_1, initial_marking, final_marking,
                               find_all_opt_alignments=True)
    print("########################################################################\nPetri Net 1")
    print(traces[0]["events"])
    print_alignments(alignments)

    petri_net_2 = PetriNet("2")
    places = {}
    for i in range(1, 7):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        petri_net_2.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),
        't_D': PetriNet.Transition('t_D', 'D'),
        't_E': PetriNet.Transition('t_E', 'E'),
        't_F': PetriNet.Transition('t_F', 'F'),
        't_G': PetriNet.Transition('t_G', 'G')
    }

    for transition in transitions:
        # print(transition)
        petri_net_2.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net_2)
    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net_2)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_B'], petri_net_2)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_D'], petri_net_2)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_F'], petri_net_2)

    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], petri_net_2)
    petri_net_utils.add_arc_from_to(transitions['t_D'], places['p_4'], petri_net_2)
    petri_net_utils.add_arc_from_to(transitions['t_F'], places['p_5'], petri_net_2)

    petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_C'], petri_net_2)
    petri_net_utils.add_arc_from_to(places['p_4'], transitions['t_E'], petri_net_2)
    petri_net_utils.add_arc_from_to(places['p_5'], transitions['t_G'], petri_net_2)

    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_6'], petri_net_2)
    petri_net_utils.add_arc_from_to(transitions['t_E'], places['p_6'], petri_net_2)
    petri_net_utils.add_arc_from_to(transitions['t_G'], places['p_6'], petri_net_2)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_6']] = 1

    # petri net visualization
    if show_petri_nets:
        parameters = {"format": "svg", 'debug': True}
        # parameters = None
        gviz = petri_net_visualization_factory.apply(petri_net_2, initial_marking, final_marking, parameters=parameters)
        petri_net_visualization_factory.view(gviz)

    traces = [
        {"frequency": 1, "events": ["A"]}
    ]
    event_log = xes_importer.import_log_from_string(create_xes_string(traces))
    alignments = factory.apply(event_log[0], petri_net_2, initial_marking, final_marking,
                               find_all_opt_alignments=True)
    print("########################################################################\nPetri Net 2")
    print(traces[0]["events"])
    print_alignments(alignments)

    petri_net_3 = PetriNet("3")
    places = {}
    for i in range(1, 6):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        petri_net_3.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),
        't_D': PetriNet.Transition('t_D', 'D'),
        't_E': PetriNet.Transition('t_E', 'E'),
        't_F': PetriNet.Transition('t_F', 'F'),
    }

    for transition in transitions:
        # print(transition)
        petri_net_3.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net_3)
    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_B'], petri_net_3)
    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_C'], petri_net_3)

    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net_3)
    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], petri_net_3)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_3'], petri_net_3)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_4'], petri_net_3)

    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_D'], petri_net_3)
    petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_E'], petri_net_3)
    petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_F'], petri_net_3)
    petri_net_utils.add_arc_from_to(places['p_4'], transitions['t_F'], petri_net_3)

    petri_net_utils.add_arc_from_to(transitions['t_D'], places['p_5'], petri_net_3)
    petri_net_utils.add_arc_from_to(transitions['t_E'], places['p_5'], petri_net_3)
    petri_net_utils.add_arc_from_to(transitions['t_F'], places['p_5'], petri_net_3)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_5']] = 1

    # petri net visualization
    if show_petri_nets:
        parameters = {"format": "svg", 'debug': True}
        # parameters = None
        gviz = petri_net_visualization_factory.apply(petri_net_3, initial_marking, final_marking, parameters=parameters)
        petri_net_visualization_factory.view(gviz)

    traces = [
        {"frequency": 1, "events": ["B", "D"]}
    ]
    event_log = xes_importer.import_log_from_string(create_xes_string(traces))
    alignments = factory.apply(event_log[0], petri_net_3, initial_marking, final_marking,
                               find_all_opt_alignments=True)
    print("########################################################################\nPetri Net 3")
    print(traces[0]["events"])
    print_alignments(alignments)

    petri_net_4 = PetriNet("4")
    places = {}
    for i in range(1, 8):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        petri_net_4.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),
        't_D': PetriNet.Transition('t_D', 'D'),
        't_E': PetriNet.Transition('t_E', 'E'),
        't_F': PetriNet.Transition('t_F', 'F'),
        't_G': PetriNet.Transition('t_G', 'G'),
        't_H': PetriNet.Transition('t_H', 'H'),
        't_I': PetriNet.Transition('t_I', 'I'),
        't_J': PetriNet.Transition('t_J', 'J'),

    }

    for transition in transitions:
        # print(transition)
        petri_net_4.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_B'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_C'], petri_net_4)

    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_4'], petri_net_4)

    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_D'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_E'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_F'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_4'], transitions['t_G'], petri_net_4)

    petri_net_utils.add_arc_from_to(transitions['t_D'], places['p_5'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_E'], places['p_5'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_F'], places['p_5'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_G'], places['p_5'], petri_net_4)

    petri_net_utils.add_arc_from_to(places['p_5'], transitions['t_H'], petri_net_4)

    petri_net_utils.add_arc_from_to(transitions['t_H'], places['p_6'], petri_net_4)

    petri_net_utils.add_arc_from_to(places['p_6'], transitions['t_I'], petri_net_4)
    petri_net_utils.add_arc_from_to(places['p_6'], transitions['t_J'], petri_net_4)

    petri_net_utils.add_arc_from_to(transitions['t_I'], places['p_7'], petri_net_4)
    petri_net_utils.add_arc_from_to(transitions['t_J'], places['p_7'], petri_net_4)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_7']] = 1

    # petri net visualization
    if show_petri_nets:
        parameters = {"format": "svg", 'debug': True}
        # parameters = None
        gviz = petri_net_visualization_factory.apply(petri_net_4, initial_marking, final_marking, parameters=parameters)
        petri_net_visualization_factory.view(gviz)

    traces = [
        {"frequency": 1, "events": ["H"]}
    ]
    event_log = xes_importer.import_log_from_string(create_xes_string(traces))
    alignments = factory.apply(event_log[0], petri_net_4, initial_marking, final_marking,
                               find_all_opt_alignments=True)
    print("########################################################################\nPetri Net 4")
    print(traces[0]["events"])
    print_alignments(alignments)

    petri_net_5 = PetriNet("5")
    places = {}
    for i in range(1, 7):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        petri_net_5.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),
        't_D': PetriNet.Transition('t_D', 'D'),
        't_1': PetriNet.Transition('t_1', None),

    }

    for transition in transitions:
        # print(transition)
        petri_net_5.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net_5)

    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net_5)
    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_3'], petri_net_5)

    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_B'], petri_net_5)
    petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_C'], petri_net_5)

    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_4'], petri_net_5)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_5'], petri_net_5)

    petri_net_utils.add_arc_from_to(places['p_4'], transitions['t_1'], petri_net_5)
    petri_net_utils.add_arc_from_to(places['p_4'], transitions['t_D'], petri_net_5)
    petri_net_utils.add_arc_from_to(places['p_5'], transitions['t_D'], petri_net_5)
    petri_net_utils.add_arc_from_to(transitions['t_1'], places['p_2'], petri_net_5)

    petri_net_utils.add_arc_from_to(transitions['t_D'], places['p_6'], petri_net_5)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_6']] = 1

    # petri net visualization
    if show_petri_nets:
        parameters = {"format": "svg", 'debug': True}
        # parameters = None
        gviz = petri_net_visualization_factory.apply(petri_net_5, initial_marking, final_marking, parameters=parameters)
        petri_net_visualization_factory.view(gviz)

    traces = [
        {"frequency": 1, "events": ["A", "B", "D"]}
    ]
    event_log = xes_importer.import_log_from_string(create_xes_string(traces))
    alignments = factory.apply(event_log[0], petri_net_5, initial_marking, final_marking,
                               find_all_opt_alignments=True)
    print("########################################################################\nPetri Net 5")
    print(traces[0]["events"])
    print_alignments(alignments)

    # print()
    # for alignment in alignments["alignments"]:
    #     for step in alignment:
    #         print(step)
    #     print()


if __name__ == "__main__":
    execute_script(show_petri_nets=True)
