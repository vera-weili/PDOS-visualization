#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
from pymatgen.io.vasp.outputs import Vasprun    # read vasprun.xml output file of VASP
from pymatgen.electronic_structure.core import Spin


# In[3]:


import chart_studio.plotly as pltly      # plotting functions
import plotly.tools as tls  # plotly tools
from plotly.subplots import make_subplots
import plotly.graph_objs as go 
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotterProjected


# In[ ]:





# In[ ]:


# dosplot = DosPlotter(sigma=0.1)
# dosplot.add_dos("Total DOS", dos)
# dosplot.add_dos_dict(dos.get_element_dos())
# plt = dosplot.get_plot(xlim=[-10,10])
# plt
# plt.grid()

# plt.show()


# In[5]:


dosrun = Vasprun("./vasprun.xml")
spd_dos = dosrun.complete_dos.get_element_spd_dos("Sn")
print(spd_dos)


# In[6]:


pdos_orbt = {}
for obt, value in spd_dos.items():
    pdos_orbt[obt.name] = value


# In[ ]:





# In[7]:


trace_tdos = go.Scatter(
    x=dosrun.tdos.densities[Spin.up],
    y=dosrun.tdos.energies - dosrun.efermi - 0.151,
    mode="lines",
    name="total DOS",
    line=go.scatter.Line(color="#444444"),
#     fill="tozeroy"
)

# 3s contribution to the total DOS
trace_3s = go.Scatter(
    x=pdos_orbt["s"].densities[Spin.up],
    y=dosrun.tdos.energies - dosrun.efermi - 0.151,
    mode="lines",
    name="s",
    line=go.scatter.Line(color="red")
)

# 3p contribution to the total DOS
trace_3p = go.Scatter(
    x=pdos_orbt["p"].densities[Spin.up],
    y=dosrun.tdos.energies - dosrun.efermi - 0.151,
    mode="lines",
    name="p",
    line=go.scatter.Line(color="green")
)
dosdata = go.Data([trace_tdos, trace_3s, trace_3p])


# In[8]:


# Customize axes and general aspect of the plot
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
dosxaxis = go.layout.XAxis(
    showgrid=True,
    showline=True,
    range=[.01, 5],
    mirror="ticks",
    ticks="inside",
    ticklen=10,
    linewidth=2,
    tickwidth=2
)
dosyaxis = go.layout.YAxis(
    title="$E - E_f \quad / \quad \\text{eV}$",
    showgrid=True,
    showline=True,
    ticks="inside",
    mirror='ticks',
    linewidth=2,
    tickwidth=2,
    ticklen=10,
    range=[-13,10],
    zerolinewidth=2
)
doslayout = go.Layout(
    width=550,
    height=700,
    title="Density of states",
    xaxis=dosxaxis,
    yaxis=dosyaxis,
    font=dict(
    family="Arial",
    size=22
        )
)

dosfig = go.Figure(data=dosdata, layout=doslayout)
plot_url = iplot(dosfig, filename="DOS")


# In[9]:


run = Vasprun("./vasprun.xml", parse_projected_eigen = True)
bands = run.get_band_structure("./KPOINTS", line_mode=True, efermi=dosrun.efermi)


# In[10]:


# Look for the boundaries of the band diagram in order to set up y axes range.
emin = 1e100
emax = -1e100
for spin in bands.bands.keys():
    for band in range(bands.nb_bands):
        emin = min(emin, min(bands.bands[spin][band]))
        emax = max(emax, max(bands.bands[spin][band]))
emin = emin - bands.efermi - 1 
emax = emax - bands.efermi + 1


# In[11]:


# Each band is plotted using a scatter plot.
kptslist = [k for k in range(len(bands.kpoints))]
bandTraces = list()
for band in range(bands.nb_bands):
    bandTraces.append(
        go.Scatter(
            x=kptslist,
            y=[(e - bands.efermi + 0.09) for e in bands.bands[Spin.up][band]],
            mode="lines",
            line=go.scatter.Line(color="#666666"),
            showlegend=False
        )
    )


# In[12]:


labels = [r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$", r"$R$", r"$X|M$",r"$R$"]
step = len(bands.kpoints) / (len(labels) - 1)
# vertical lines
vlines = list()
for i, label in enumerate(labels):
    vlines.append(
        go.Scatter(
            x=[i * step, i * step],
            y=[emin, emax],
            mode="lines",
            line=go.Line(color="#111111", width=1),
            showlegend=False
        )
    )
    
# Labels of highsymetry k-points are added as Annotation object
annotations = list()
for i, label in enumerate(labels):
    annotations.append(
        go.layout.Annotation(
            x=i * step, y=emin,
            xref="x1", yref="y1",
            text=label,
            xanchor="center", yanchor="top",
            showarrow=False
        )
    )


# In[13]:


# Customize axes and general aspect of the plot
bandxaxis = go.layout.XAxis(
    title="k-points",
    range=[0, len(bands.kpoints)],
    showgrid=True,
    showline=True,
    ticks="",
    showticklabels=False,
    mirror=True,
    linewidth=2
)

bandyaxis = go.layout.YAxis(
    title="$E - E_f \quad / \quad \\text{eV}$",
    range=[emin, emax],

    showgrid=True,
    showline=True,
    zeroline=True,
    mirror="ticks",
    ticks="inside",
    linewidth=2,
    tickwidth=2,
    zerolinewidth=2,
)

bandlayout = go.Layout(
    width=700,
    height=500,
    title="Bands diagram",
    xaxis=bandxaxis,
    yaxis=bandyaxis,
    annotations=go.Annotations(annotations)
)

bandfig = go.Figure(data=bandTraces + vlines, layout=bandlayout)
plot_url = iplot(bandfig, filename="Bands")


# In[14]:


# plot DOS and bands together
dosbandfig = make_subplots(rows=1, cols=2, shared_yaxes=True)
# add the bands
for btrace in bandTraces:
    dosbandfig.append_trace(btrace, 1, 1)
# add vlines for specific k-points
for vline in vlines:
    dosbandfig.append_trace(vline, 1, 1)
# add the densities
dosbandfig.append_trace(trace_tdos, 1, 2)
dosbandfig.append_trace(trace_3s, 1, 2)
dosbandfig.append_trace(trace_3p, 1, 2)

dosbandfig["layout"].update(
    go.Layout(
        title="Bands diagram and density of states",
        xaxis1=bandxaxis,
        yaxis1=dosyaxis,
        xaxis2=dosxaxis,
        annotations=go.Annotations(annotations)
    )
)

# adjust size of subplots
dosbandfig["layout"]["xaxis1"]["domain"] = [0., 0.7]
dosbandfig["layout"]["xaxis2"]["domain"] = [0.702, 1.]
# add some specific options
dosbandfig["layout"]["yaxis1"]["mirror"] = "allticks"
dosbandfig["layout"]["xaxis2"]["mirror"] = "allticks"

plot_url = iplot(dosbandfig, filename="DOS_bands_Si")


# In[15]:


# Atomic orbital contributions in bands using a color scale
name = "O"
pbands = bands.get_projections_on_elements_and_orbitals({name: ["s", "p", "d"]})


# In[16]:


# compute contributions
contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
for band in range(bands.nb_bands):
    for k in range(len(bands.kpoints)):
        sc = pbands[Spin.up][band][k][name]["s"]**2
        pc = pbands[Spin.up][band][k][name]["p"]**2
        dc = pbands[Spin.up][band][k][name]["d"]**2
        tot = sc + pc + dc
        if tot != 0.0:
            contrib[band, k, 0] = sc / tot
            contrib[band, k, 1] = pc / tot
            contrib[band, k, 2] = dc / tot


# In[17]:


colorBands = list() # will contain the list of all lines
nkpts = len(bands.kpoints)
for band in range(bands.nb_bands):
    eband = [e - bands.efermi for e in bands.bands[Spin.up][band]]
    for k in range(nkpts - 1):
        red, green, blue = [int(255 * (contrib[band, k, i] + contrib[band, k+1, i])/2) for i in range(3)]
        colorBands.append(
            go.Scatter(
                x=[k, k+1],
                y=[eband[k], eband[k+1]],
                mode="lines",
                line=go.Line(color="rgb({}, {}, {})".format(red, green, blue)),
                showlegend=False
            )
        )


# In[18]:


# set up a new figure with two subplots
colorbandfig = make_subplots(rows=1, cols=2, shared_yaxes=True)
# add the bands in the first subplot
for btrace in colorBands:
    colorbandfig.append_trace(btrace, 1, 1)
# add vlines for specific k-points in the first subplot
for vline in vlines:
    colorbandfig.append_trace(vline, 1, 1)
# add the densities in the second subplot
colorbandfig.append_trace(trace_tdos, 1, 2)
colorbandfig.append_trace(trace_3s, 1, 2)
colorbandfig.append_trace(trace_3p, 1, 2)
# Layout configuration
colorbandfig["layout"].update(
    go.Layout(
        title="Bands diagram and density of states",
        xaxis1=bandxaxis,
        yaxis1=dosyaxis,
        xaxis2=dosxaxis,
        annotations=go.Annotations(annotations)
    )
)


# In[19]:


# adjust size of subplots
colorbandfig["layout"]["xaxis1"]["domain"] = [0., 0.7]
colorbandfig["layout"]["xaxis2"]["domain"] = [0.702, 1.]
# add some specific options
colorbandfig["layout"]["yaxis1"]["mirror"] = "allticks"
colorbandfig["layout"]["xaxis2"]["mirror"] = "allticks"
# add a custom legend
legend = go.layout.Legend(
    x=.98, y=.98,
    xanchor="right", yanchor="top",
    bordercolor='#333', borderwidth=1
)
colorbandfig["layout"]["legend"] = legend


# In[20]:


plot_url = iplot(colorbandfig, filename="DOS_bands_color")


# In[78]:





# In[ ]:





# In[ ]:




