from __future__ import annotations

import matplotlib.pyplot as plt
import pybmad
from pytao import Tao

tao = Tao(init_file="$ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init")

print(tao)

print("\n".join(tao.cmd("show uni")))

s = pybmad.get_super_universe()

regions = [r for r in s.plot_page.region if r.plot.description and r.plot.graph.is_valid()]

plt.ion()

for region in regions:
    for graph in region.plot.graph:
        curves = [curve for curve in graph.curve if curve.valid]

        if not curves:
            continue

        plt.figure()
        for curve in curves:
            plt.plot(curve.x_line, curve.y_line, label=curve.legend_text)
            plt.title(graph.title)
            plt.legend()


plt.ioff()
