import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from Beam_loading_diagram import draw_loading_diagram

st.set_page_config(page_title="Beam Analysis", layout="wide")
st.title("Simply Supported Beam Analysis")
st.caption("Point loads and linear distributed loads • FEM deflection")

L = st.number_input("Beam length (m)", min_value=0.1, value=6.0)
E = st.number_input("Young's modulus (GPa)", min_value=0.001, value=200.0) * 1e9
I = st.number_input("Second moment of area (m⁴)", min_value=1e-12, value=8e-6, format="%.8f")

st.subheader("Point loads")
n_point = st.number_input("Number of point loads", 0, 20, 2)
point_loads = []
for i in range(n_point):
    c1, c2 = st.columns(2)
    P = c1.number_input(f"Load {i+1} (kN, downward)", min_value=0.0, value=10.0 if i == 0 else 15.0, key=f"p{i}")
    a = c2.number_input(f"Position {i+1} (m)", min_value=0.0, max_value=float(L), value=min(float(L), 2.0 if i == 0 else 4.0), key=f"a{i}")
    point_loads.append((P*1000, a))

st.subheader("Distributed loads")
st.caption("Each load varies linearly between its starting and ending intensity. Equal values give a UDL.")
n_dist = st.number_input("Number of distributed loads", 0, 20, 2)
distributed_loads = []
for i in range(n_dist):
    c1, c2, c3, c4 = st.columns(4)
    a = c1.number_input(f"Start {i+1} (m)", min_value=0.0, max_value=float(L), value=min(float(L), 0.0 if i == 0 else 2.0), key=f"da{i}")
    b = c2.number_input(f"End {i+1} (m)", min_value=0.0, max_value=float(L), value=float(L) if i == 0 else min(float(L), 5.0), key=f"db{i}")
    qa = c3.number_input(f"Start intensity {i+1} (kN/m)", min_value=0.0, value=2.0 if i == 0 else 0.0, key=f"qa{i}")
    qb = c4.number_input(f"End intensity {i+1} (kN/m)", min_value=0.0, value=2.0 if i == 0 else 3.0, key=f"qb{i}")
    distributed_loads.append((a, b, qa*1000, qb*1000))

def shapes(s, le):
    return np.array([1-3*s*s+2*s**3, le*(s-2*s*s+s**3),
                     3*s*s-2*s**3, le*(-s*s+s**3)])

st.subheader("Beam Loading Diagram")

fig = draw_loading_diagram(
    L,
    point_loads,
    distributed_loads
)

st.pyplot(fig)
plt.close(fig)

if st.button("Calculate", type="primary"):
    if any(b <= a for a,b,_,_ in distributed_loads):
        st.error("Each distributed load must end after it starts.")
        st.stop()
    breaks = sorted(set([0.0, L] + [a for _,a in point_loads]
                        + [v for a,b,_,_ in distributed_loads for v in (a,b)]))
    nodes = [breaks[0]]
    for a,b in zip(breaks[:-1], breaks[1:]):
        count = max(1, int(np.ceil((b-a)/(L/24))))
        nodes.extend(np.linspace(a,b,count+1)[1:])
    nodes = np.array(nodes)
    n = len(nodes)
    K = np.zeros((2*n,2*n))
    F = np.zeros(2*n)
    gp,gw = np.polynomial.legendre.leggauss(3)
    def intensity(xx):
        return sum(qa+(qb-qa)*(xx-a)/(b-a)
                   for a,b,qa,qb in distributed_loads if a <= xx <= b)
    for j in range(n-1):
        le = nodes[j+1]-nodes[j]
        ke = E*I/le**3*np.array([
            [12,6*le,-12,6*le],
            [6*le,4*le**2,-6*le,2*le**2],
            [-12,-6*le,12,-6*le],
            [6*le,2*le**2,-6*le,4*le**2]])
        fe = np.zeros(4)
        for z,w in zip(gp,gw):
            s=(z+1)/2
            xx=nodes[j]+s*le
            fe -= w*le/2*intensity(xx)*shapes(s,le)
        ids=[2*j,2*j+1,2*j+2,2*j+3]
        K[np.ix_(ids,ids)] += ke
        F[ids] += fe
    for P,a in point_loads:
        j=int(np.argmin(abs(nodes-a)))
        F[2*j] -= P
    fixed=[0,2*(n-1)]
    free=np.array([i for i in range(2*n) if i not in fixed])
    U=np.zeros(2*n)
    U[free]=np.linalg.solve(K[np.ix_(free,free)],F[free])
    R=K@U-F
    RA,RB=R[fixed]
        # Engineering validation: check force and moment equilibrium
    total_load = sum(P for P, a in point_loads)
    total_moment = sum(P * a for P, a in point_loads)

    for a, b, qa, qb in distributed_loads:
        length = b - a
        resultant = (qa + qb) * length / 2
        total_load += resultant

        # First moment of a linearly varying distributed load
        total_moment += a * resultant + length**2 * (qa + 2 * qb) / 6

    force_error = RA + RB - total_load
    moment_error = RB * L - total_moment

    x=np.linspace(0,L,1501)
    v=np.zeros_like(x)
    for k,xx in enumerate(x):
        j=min(max(np.searchsorted(nodes,xx,side="right")-1,0),n-2)
        le=nodes[j+1]-nodes[j]
        v[k]=shapes((xx-nodes[j])/le,le)@U[2*j:2*j+4]
    V=np.full_like(x,RA)
    M=RA*x
    for P,a in point_loads:
        V-=P*(x>=a)
        M-=P*np.maximum(x-a,0)
    for a,b,qa,qb in distributed_loads:
        t=np.clip(x-a,0,b-a)
        slope=(qb-qa)/(b-a)
        V-=qa*t+slope*t*t/2
        M-=qa*(t*(x-a)-t*t/2)+slope*(t*t*(x-a)/2-t**3/3)
    st.subheader("Results")
    with st.expander("Engineering validation checks"):
            st.write(f"Force equilibrium error: **{force_error:.6f} N**")
            st.write(f"Moment equilibrium error: **{moment_error:.6f} Nm**")
    c1,c2,c3=st.columns(3)
    c1.metric("Left reaction",f"{RA/1000:.3f} kN")
    c2.metric("Right reaction",f"{RB/1000:.3f} kN")
    c3.metric("Maximum moment",f"{np.max(np.abs(M))/1000:.3f} kNm")
    # Find the maximum absolute bending moment
    moment_idx = int(np.argmax(np.abs(M)))
    moment_max = M[moment_idx] / 1000
    moment_x = x[moment_idx]

       # Find the maximum absolute bending moment
    moment_idx = int(np.argmax(np.abs(M)))
    moment_max = M[moment_idx] / 1000
    moment_x = x[moment_idx]

    st.write(
        f"Maximum bending moment: **{moment_max:.3f} kNm** "
        f"at **{moment_x:.3f} m**"
    )



    idx=int(np.argmin(v))
    st.write(f"Maximum downward deflection: **{v[idx]*1000:.3f} mm** at **{x[idx]:.3f} m**")
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
    for ax,y,title,unit in zip(axes,[V/1000,M/1000,v*1000],
                               ["Shear force","Bending moment","FEM deflection"],
                               ["kN","kNm","mm"]):
        ax.plot(x,y)
        ax.axhline(0,color="black",lw=0.7)
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(True)

        # Highlight maximum bending moment
    axes[1].scatter(
        moment_x, moment_max,
        color="red", s=65, zorder=5
    )

    axes[1].annotate(
        f"{moment_max:.3f} kNm at {moment_x:.3f} m",
        xy=(moment_x, moment_max),
        xytext=(10, 15),
        textcoords="offset points",
        fontsize=10,
        color="red"
    )    

        # Highlight maximum downward deflection
    axes[2].scatter(
        x[idx], v[idx] * 1000,
        color="red", s=65, zorder=5
    )

    axes[2].annotate(
        f"{v[idx]*1000:.3f} mm at {x[idx]:.3f} m",
        xy=(x[idx], v[idx] * 1000),
        xytext=(10, 20),
        textcoords="offset points",
        fontsize=10,
        color="red"
    )
    for ax in axes:
        ax.set_xlabel("Position along beam (m)")

    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)