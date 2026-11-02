print("LOADING DIAGRAM FILE:", __file__)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle

def draw_loading_diagram(L, point_loads, distributed_loads):
    """Draw simply supported beam and loads. Forces positive downward."""
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot([0, L], [0, 0], color="black", lw=4, solid_capstyle="round")
    h = 0.065 * L
    ax.add_patch(Polygon([[0, -0.02*L], [-h/2, -h], [h/2, -h]],
                         closed=True, fill=False, lw=1.5))
    ax.add_patch(Polygon([[L, -0.02*L], [L-h/2, -h], [L+h/2, -h]],
                         closed=True, fill=False, lw=1.5))
    for dx in [-h/3, 0, h/3]:
        ax.add_patch(Circle((L+dx, -h-0.018*L), 0.012*L, fill=False))
    maxq = max([max(abs(qa), abs(qb)) for _, _, qa, qb in distributed_loads] + [1])
    for k, (a, b, qa, qb) in enumerate(distributed_loads):
        xs = np.linspace(a, b, 13)
        qs = qa + (qb-qa)*(xs-a)/(b-a)
        base = 0.19*L + k*0.16*L
        heights = 0.12*L * np.abs(qs)/maxq
        for xx, q, hh in zip(xs, qs, heights):
            if abs(q) < 1e-12: continue
            yy = base + hh
            ax.annotate("", xy=(xx, base), xytext=(xx, yy),
                        arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=1.1))
        ax.plot(xs, base+heights, color="tab:blue", lw=1.5)
        ax.text((a+b)/2, base+max(heights)+0.025*L,
                f"{qa/1000:g} → {qb/1000:g} kN/m",
                ha="center", fontsize=9, color="tab:blue")
        ax.plot([a, b], [base, base], color="tab:blue", alpha=.4)
    for i, (P, a) in enumerate(point_loads):
        top = (0.68 + 0.12*(i%2))*L
        ax.annotate("", xy=(a, 0.015*L), xytext=(a, top),
                    arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2))
        ax.text(a, top+0.025*L, f"{P/1000:g} kN",
                color="tab:red", ha="center", fontsize=9)
    ax.annotate("", xy=(0, -0.19*L), xytext=(L, -0.19*L),
                arrowprops=dict(arrowstyle="<->", lw=1))
    ax.text(L/2, -0.17*L, f"L = {L:g} m", ha="center")
    ax.set_xlim(-0.12*L, 1.12*L)
    ax.set_ylim(-0.27*L, max(0.95*L, (0.19+0.16*len(distributed_loads)+0.2)*L))
    ax.axis("off")
    fig.tight_layout()
    return fig

