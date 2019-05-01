import matplotlib
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import os
import numpy as np

figure(num=None, figsize=(11, 6))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ','))


def plot_time_per_algorithm(time_a_star_computation_without_heuristic, time_heuristic_computation, description,
                            number_solved_lps,
                            path_to_store="", svg=False):
    # time_a_star_computation_without_heuristic = (20, 35, 30, 35, 27, 34, 78, 78)
    # time_heuristic_computation = (25, 32, 34, 200, 25, 0, 0, 0)
    fig = figure(num=None, figsize=(len(description) * 1.5, 6))

    ind = range(len(time_a_star_computation_without_heuristic))  # the x locations for the groups
    width = 0.35  # the width of the bars: can also be len(x) sequence

    plt.grid(zorder=0, color=(.9, .9, .9))
    p1 = plt.bar(ind, time_heuristic_computation, width, color=(0.6, 0.6, 0.6), zorder=3)
    p2 = plt.bar(ind, time_a_star_computation_without_heuristic, width, color=(0.1, 0.1, 0.1),
                 bottom=time_heuristic_computation, zorder=3)
    plt.ylabel('average time per trace (seconds)', fontsize=12)
    plt.title('average number of solved LPs per trace', fontsize=12, color="blue")

    plt.xticks(ind, description)
    plt.legend((p1[0], p2[0]),
               ('heuristic computation time', 'A* computation time\n(excluding heuristic computation time)'),
               fontsize=12)
    if max(time_heuristic_computation) > 1000 or max(time_a_star_computation_without_heuristic) > 1000:
        axes = fig.get_axes()
        axes[-1].get_yaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    for i in range(len(time_a_star_computation_without_heuristic)):
        y = time_heuristic_computation[i] + time_a_star_computation_without_heuristic[i]
        plt.text(x=i, y=y, s=number_solved_lps[i], size=11, color="blue", horizontalalignment='center')

    if path_to_store:
        if svg:
            plt.savefig(path_to_store + ".svg")
        else:
            plt.savefig(path_to_store + ".png")
    else:
        plt.show()
    plt.clf()
    plt.close()


def generate_simple_bar_plot(y_values, description, path_to_store="", attribute="", svg=False):
    # time_a_star_computation_without_heuristic = (20, 35, 30, 35, 27, 34, 78, 78)
    # time_heuristic_computation = (25, 32, 34, 200, 25, 0, 0, 0)
    fig = figure(num=None, figsize=(len(description) * 1.5, 6))

    ind = range(len(y_values))  # the x locations for the groups
    width = .5  # the width of the bars: can also be len(x) sequence

    plt.bar(ind, y_values, width, color=(0.1, 0.1, 0.1), zorder=2)
    plt.grid(zorder=0, color=(.9, .9, .9))
    plt.ylabel('average ' + attribute.replace("_", " ") + ' per trace', fontsize=12)
    # plt.title('Time to compute prefix-alignments for 100 traces', fontsize=12)
    plt.xticks(ind, description)
    axes = fig.get_axes()
    axes[-1].get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    if path_to_store:
        if svg:
            plt.savefig(path_to_store + ".svg")
        else:
            plt.savefig(path_to_store + ".png")
    else:
        plt.show()
    plt.clf()
    plt.close()


def plot_length_distribution(lengths, path_to_store):
    figure(num=None, figsize=(5, 3))
    print(lengths)

    names = range(1, 11)
    values = []
    for i in range(1, 11):
        values.append(lengths.count(i))
    plt.bar(names, values, color=(0.1, 0.1, 0.1), width=0.4)
    # plt.title("distribution trace length")
    plt.xlabel("trace length", fontsize=12)
    plt.ylabel("frequency", fontsize=12)
    path = os.path.join(path_to_store, "length_distribution_plot.svg")
    plt.savefig(path)
    plt.close()
    for i in range(1, 11):
        print(i, " ", lengths.count(i))


def plot_search_space_size_depending_on_prefix_length(keys, keys_to_label, data, attribute, path_to_store="",
                                                      svg=False):
    # t = [1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # s = [1, 2, 5, 3.5, 4, 5, 6, 7, 8, 9, 10]
    # s2 = [1, 2, 2, 1, 5, 5, 2, 7, 8, 9, 10]
    # s3 = [1, 2, 5, 3, 5, 5, 2, 7, 8, 9, 19]
    #
    # plt.plot(t, s, marker='o',linestyle='dashed', label="$a_1$")
    # plt.plot(t, s2,  marker='o',linestyle='dashed', label="b")
    # plt.plot(t, s3, marker='o',linestyle='dashed', label="b")

    fig = figure(num=None, figsize=(9, 5))
    for key in keys:
        if key in data:
            plt.plot(range(1, len(data[key]) + 1), data[key], marker='o', linestyle='dashed', label=keys_to_label[key],
                     zorder=3)
    plt.xlabel("prefix length", fontsize=12)
    plt.ylabel("cumulated number of " + attribute.replace('_', ' '), fontsize=12)
    plt.xticks([i + 1 for i in range(len(data[key]))])
    plt.legend(loc='upper left')
    plt.grid(zorder=0, color=(.9, .9, .9))
    axes = fig.get_axes()
    axes[-1].get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    if path_to_store:
        if svg:
            plt.savefig(path_to_store + ".svg")
        else:
            plt.savefig(path_to_store + ".png")
    else:
        plt.show()
    plt.clf()
    plt.close()


if __name__ == '__main__':
    # plot_time_per_algorithm((20, 35, 30, 35, 27, 34, 78, 78), (25, 0, 34, 0, 25, 0, 0, 0))
    a = (
        58882.53734087944, 565.6675276756287, 6085.766488313675, 66.61048531532288, 633.4740769863129, 4.20329737663269,
        4.312514305114746, 6.859708309173584)
    b = (0, 62.672616720199585, 0, 137.3296763896942, 71.92269444465637, 2.5155222415924072, 3.9844913482666016,
         4.968816757202148)
    # plot_time_per_algorithm(a, b, None)
    plot_search_space_size_depending_on_prefix_length()
