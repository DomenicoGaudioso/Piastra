import streamlit as st
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go
import pandas as pd

# ==========================================
# 1. FUNZIONI GEOMETRICHE
# ==========================================

def parse_input(testo):
    try: return [float(x.strip()) for x in testo.split(',')]
    except: return []

def discretizza_asse(interassi, dim_mesh_max):
    coordinate = [0.0]
    for campata in interassi:
        num_elementi = int(np.ceil(campata / dim_mesh_max))
        passo_reale = campata / num_elementi
        for _ in range(num_elementi):
            coordinate.append(coordinate[-1] + passo_reale)
    return np.array(coordinate)

def genera_griglia(interassi_x, interassi_y, dim_mesh_max):
    return discretizza_asse(interassi_x, dim_mesh_max), discretizza_asse(interassi_y, dim_mesh_max)

def trova_nodi_vincolati(x_coords, y_coords, assi_x_continui, assi_y_continui, punti_vincolati):
    nodi_bloccati = []
    tol = 1e-5
    for i, x in enumerate(x_coords):
        for j, y in enumerate(y_coords):
            is_supporto = False
            if any(abs(x - ax) < tol for ax in assi_x_continui): is_supporto = True
            if any(abs(y - ay) < tol for ay in assi_y_continui): is_supporto = True
            if any((abs(x - px) < tol and abs(y - py) < tol) for px, py in punti_vincolati): is_supporto = True
            if is_supporto: nodi_bloccati.append((i, j))
    return nodi_bloccati

# ==========================================
# 2. MOTORE FEM E RISULTATI
# ==========================================

def esegui_analisi_fem(x_coords, y_coords, nodi_bloccati, carichi_conc_df, carichi_uni_df, E, nu, spessore):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    ops.section('ElasticMembranePlateSection', 1, E, nu, spessore, 0.0)
    
    node_tags = {}
    tag_nodo = 1
    for i, x in enumerate(x_coords):
        for j, y in enumerate(y_coords):
            ops.node(tag_nodo, float(x), float(y), 0.0)
            node_tags[(i, j)] = tag_nodo
            tag_nodo += 1
            
    element_tags = {}
    tag_ele = 1
    for i in range(len(x_coords) - 1):
        for j in range(len(y_coords) - 1):
            n1, n2 = node_tags[(i, j)], node_tags[(i+1, j)]
            n3, n4 = node_tags[(i+1, j+1)], node_tags[(i, j+1)]
            ops.element('ShellMITC4', tag_ele, n1, n2, n3, n4, 1)
            element_tags[(i, j)] = tag_ele
            tag_ele += 1

    for (i, j) in nodi_bloccati:
        ops.fix(node_tags[(i, j)], 1, 1, 1, 0, 0, 0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    tol = 1e-5
    
    for _, row in carichi_conc_df.iterrows():
        xc, yc, L, B, F_tot = row['X'], row['Y'], row['L'], row['B'], row['Carico_Totale_kN']
        xmin, xmax = xc - L/2, xc + L/2
        ymin, ymax = yc - B/2, yc + B/2
        
        nodi_caricati = []
        for (i, j), tag in node_tags.items():
            if (xmin - tol <= x_coords[i] <= xmax + tol) and (ymin - tol <= y_coords[j] <= ymax + tol):
                nodi_caricati.append(tag)
        
        if not nodi_caricati:
            distanze = {tag: (x_coords[i]-xc)**2 + (y_coords[j]-yc)**2 for (i, j), tag in node_tags.items()}
            nodi_caricati.append(min(distanze, key=distanze.get))

        f_nodo = -(F_tot * 1000) / len(nodi_caricati) 
        for tag in nodi_caricati:
            ops.load(tag, 0.0, 0.0, f_nodo, 0.0, 0.0, 0.0)

    for _, row in carichi_uni_df.iterrows():
        xmin, xmax, ymin, ymax = row['X_min'], row['X_max'], row['Y_min'], row['Y_max']
        q_N_m2 = -(row['q_kN_m2'] * 1000)
        
        for i in range(len(x_coords) - 1):
            for j in range(len(y_coords) - 1):
                xc_el, yc_el = (x_coords[i] + x_coords[i+1]) / 2.0, (y_coords[j] + y_coords[j+1]) / 2.0
                if (xmin - tol <= xc_el <= xmax + tol) and (ymin - tol <= yc_el <= ymax + tol):
                    area_el = (x_coords[i+1] - x_coords[i]) * (y_coords[j+1] - y_coords[j])
                    f_nodo_uni = (q_N_m2 * area_el) / 4.0
                    for nodo in [node_tags[(i, j)], node_tags[(i+1, j)], node_tags[(i+1, j+1)], node_tags[(i, j+1)]]:
                        ops.load(nodo, 0.0, 0.0, f_nodo_uni, 0.0, 0.0, 0.0)

    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ok = ops.analyze(1)
    
    Nx, Ny = len(x_coords), len(y_coords)
    U_z = np.zeros((Ny, Nx))
    U_x = np.zeros((Ny, Nx))
    U_y = np.zeros((Ny, Nx))
    
    for i in range(Nx):
        for j in range(Ny):
            tag = node_tags[(i, j)]
            U_x[j, i] = ops.nodeDisp(tag, 1)
            U_y[j, i] = ops.nodeDisp(tag, 2)
            U_z[j, i] = ops.nodeDisp(tag, 3)

    sforzi = {}
    
    if ok == 0:
        dw_dy, dw_dx = np.gradient(U_z, y_coords, x_coords)
        d2w_dy2, d2w_dydx = np.gradient(dw_dy, y_coords, x_coords)
        d2w_dxdy, d2w_dx2 = np.gradient(dw_dx, y_coords, x_coords)
        
        dux_dy, dux_dx = np.gradient(U_x, y_coords, x_coords)
        duy_dy, duy_dx = np.gradient(U_y, y_coords, x_coords)
        
        D = (E * (spessore**3)) / (12 * (1 - nu**2))
        C = (E * spessore) / (1 - nu**2)
        
        sforzi['Nxx'] = C * (dux_dx + nu * duy_dy) / 1000.0
        sforzi['Nyy'] = C * (duy_dy + nu * dux_dx) / 1000.0
        sforzi['Nxy'] = C * ((1 - nu) / 2.0) * (dux_dy + duy_dx) / 1000.0
        
        sforzi['Mxx'] = -D * (d2w_dx2 + nu * d2w_dy2) / 1000.0
        sforzi['Myy'] = -D * (d2w_dy2 + nu * d2w_dx2) / 1000.0
        sforzi['Mxy'] = -D * (1 - nu) * d2w_dydx / 1000.0
        
        dMxx_dy, dMxx_dx = np.gradient(sforzi['Mxx'], y_coords, x_coords)
        dMxy_dy, dMxy_dx = np.gradient(sforzi['Mxy'], y_coords, x_coords)
        dMyy_dy, dMyy_dx = np.gradient(sforzi['Myy'], y_coords, x_coords)
        
        sforzi['Vxz'] = dMxx_dx + dMxy_dy
        sforzi['Vyz'] = dMyy_dy + dMxy_dx
    else:
        sforzi = {key: np.zeros((Ny, Nx)) for key in ['Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy', 'Vxz', 'Vyz']}

    return U_z, sforzi, ok

# ==========================================
# 3. INTERFACCIA STREAMLIT E PLOT
# ==========================================

st.set_page_config(layout="wide", page_title="FEM Soletta")
st.title("🏗️ Analisi FEM Soletta in Calcestruzzo")

st.sidebar.header("1. Geometria (Interassi)")
int_x_str = st.sidebar.text_input("Longitudinali X (m)", "1.0, 1.5, 4.0, 4.0, 4.0, 4.0, 1.5, 1.0")
int_y_str = st.sidebar.text_input("Trasversali Y (m)", "3.0, 3.0")

st.sidebar.header("2. Proprietà")
dim_mesh = st.sidebar.number_input("Dim. Mesh (m)", value=0.5, step=0.1)
spessore = st.sidebar.number_input("Spessore (m)", value=0.3, step=0.05)
E_cls = st.sidebar.number_input("E (N/m2)", value=3e10, format="%e")
nu_cls = st.sidebar.number_input("Poisson", value=0.2)

# --- NUOVA OPZIONE PER IL GRAFICO ---
st.sidebar.header("Impostazioni Grafiche")
mostra_isolinee = st.sidebar.checkbox("Visualizza Isolinee (Contour)", value=True)

int_x, int_y = parse_input(int_x_str), parse_input(int_y_str)

if int_x and int_y:
    x_coords, y_coords = genera_griglia(int_x, int_y, dim_mesh)
    appoggi_x, appoggi_y = np.cumsum([0.0] + int_x), np.cumsum([0.0] + int_y)

    st.sidebar.header("3. Vincoli")
    dict_assi_x = {f"Asse X a {val:.2f} m": val for val in appoggi_x}
    dict_assi_y = {f"Asse Y a {val:.2f} m": val for val in appoggi_y}
    sel_x = st.sidebar.multiselect("Muri (Assi X continui)", list(dict_assi_x.keys()))
    sel_y = st.sidebar.multiselect("Muri (Assi Y continui)", list(dict_assi_y.keys()))
    
    intersezioni = [(ax, ay) for ax in appoggi_x for ay in appoggi_y]
    dict_punti = {f"Punto (X: {pt[0]:.2f}, Y: {pt[1]:.2f})": pt for pt in intersezioni}
    sel_punti = st.sidebar.multiselect("Pilastri (Puntuali)", list(dict_punti.keys()))
    
    nodi_bloccati = trova_nodi_vincolati(x_coords, y_coords, 
                                         [dict_assi_x[k] for k in sel_x], 
                                         [dict_assi_y[k] for k in sel_y], 
                                         [dict_punti[k] for k in sel_punti])

    col1, col2 = st.columns([1, 1.5])
    
    with col1:
        st.subheader("Carichi Uniformi per Aree (kN/m²)")
        default_uni_loads = pd.DataFrame([{"X_min": 0.0, "X_max": appoggi_x[-1], "Y_min": 0.0, "Y_max": appoggi_y[-1], "q_kN_m2": 5.0}])
        carichi_uni_df = st.data_editor(default_uni_loads, num_rows="dynamic", use_container_width=True)

        st.subheader("Carichi Concentrati (Impronte in kN)")
        default_conc_loads = pd.DataFrame([{"X": 2.0, "Y": 1.5, "L": 0.5, "B": 0.5, "Carico_Totale_kN": 50.0}])
        carichi_conc_df = st.data_editor(default_conc_loads, num_rows="dynamic", use_container_width=True)
        
        visualizza_3d = st.button("🔍 Visualizza Modello 3D Pre-Analisi", use_container_width=True)
        esegui = st.button("🚀 Esegui Analisi FEM", type="primary", use_container_width=True)

    with col2:
        if visualizza_3d:
            st.subheader("Modello 3D (Mesh, Vincoli e Carichi)")
            Z_zero = np.zeros((len(y_coords), len(x_coords)))
            fig_pre = go.Figure(data=[go.Surface(z=Z_zero, x=x_coords, y=y_coords, colorscale='Greys', opacity=0.7, showscale=False)])
            
            if nodi_bloccati:
                vx = [x_coords[i] for (i, j) in nodi_bloccati]
                vy = [y_coords[j] for (i, j) in nodi_bloccati]
                vz = [0] * len(vx)
                fig_pre.add_trace(go.Scatter3d(x=vx, y=vy, z=vz, mode='markers', marker=dict(size=4, color='red', symbol='diamond'), name='Vincoli'))
            
            if not carichi_conc_df.empty:
                cx, cy, cz = carichi_conc_df['X'].tolist(), carichi_conc_df['Y'].tolist(), [0.5]*len(carichi_conc_df)
                fig_pre.add_trace(go.Scatter3d(x=cx, y=cy, z=cz, mode='markers', marker=dict(size=6, color='blue', symbol='x'), name='Carichi (Centri)'))
                
            fig_pre.update_layout(scene=dict(aspectratio=dict(x=2, y=1, z=0.5)), margin=dict(l=0, r=0, b=0, t=0))
            st.plotly_chart(fig_pre, use_container_width=True)
            
        else:
            st.subheader("Pianta 2D")
            fig, ax = plt.subplots(figsize=(10, 5))
            for y in y_coords: ax.axhline(y, color='lightgray', lw=0.5)
            for x in x_coords: ax.axvline(x, color='lightgray', lw=0.5)
            vx = [x_coords[i] for (i, j) in nodi_bloccati]
            vy = [y_coords[j] for (i, j) in nodi_bloccati]
            if vx: ax.plot(vx, vy, 'go', markersize=5, label='Vincoli')
            
            show_uni = st.checkbox("Mostra Zone Carico Uniforme", value=False)
            if show_uni:
                for _, row in carichi_uni_df.iterrows():
                    ax.add_patch(patches.Rectangle((row['X_min'], row['Y_min']), row['X_max'] - row['X_min'], row['Y_max'] - row['Y_min'], 
                                                   linewidth=2, edgecolor='orange', facecolor='orange', alpha=0.3, hatch='//'))
            for _, row in carichi_conc_df.iterrows():
                ax.add_patch(patches.Rectangle((row['X'] - row['L']/2, row['Y'] - row['B']/2), row['L'], row['B'], color='blue', alpha=0.5))
            ax.set_aspect('equal'); st.pyplot(fig)

    # --- RISULTATI ANALISI ---
    if esegui:
        if not nodi_bloccati:
            st.error("⚠️ Inserisci almeno un vincolo!")
        else:
            with st.spinner("Calcolo in corso..."):
                U_z, sforzi, ok = esegui_analisi_fem(x_coords, y_coords, nodi_bloccati, carichi_conc_df, carichi_uni_df, E_cls, nu_cls, spessore)
                
                if ok == 0:
                    st.success("Analisi completata!")
                    
                    sforzi_keys = ['Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy', 'Vxz', 'Vyz']
                    nomi_tabs = ["📉 Deformata (Uz)"] + sforzi_keys
                    tabs = st.tabs(nomi_tabs)
                    
                    with tabs[0]:
                        fig_uz = go.Figure(data=[go.Surface(z=U_z * 1000, x=x_coords, y=y_coords, colorscale='Turbo', colorbar_title='Uz (mm)')])
                        fig_uz.update_layout(scene=dict(zaxis=dict(autorange="reversed"), aspectratio=dict(x=2, y=1, z=0.5)))
                        st.plotly_chart(fig_uz, use_container_width=True)
                        
                    for idx, key in enumerate(sforzi_keys):
                        with tabs[idx + 1]:
                            unita_misura = "kNm/m" if "M" in key else "kN/m"
                            
                            # Configurazione per il Contour: 'fill' crea le classiche bande a gradini
                            # showlines disegna la riga scura di divisione
                            tipo_colore = 'fill' if mostra_isolinee else 'heatmap'
                            spessore_linea = 1 if mostra_isolinee else 0

                            fig_sforzi = go.Figure(data=[go.Contour(
                                z=sforzi[key], 
                                x=x_coords, 
                                y=y_coords, 
                                colorscale='Turbo',
                                contours=dict(
                                    coloring=tipo_colore, 
                                    showlines=mostra_isolinee
                                ),
                                line_width=spessore_linea,
                                colorbar=dict(title=f"{key} ({unita_misura})")
                            )])
                            
                            fig_sforzi.update_layout(
                                title=f"Mappa Sollecitazione: {key}",
                                yaxis=dict(scaleanchor="x", scaleratio=1), 
                                xaxis=dict(title="X (m)"),
                                margin=dict(l=20, r=20, t=40, b=20)
                            )
                            st.plotly_chart(fig_sforzi, use_container_width=True)
                            
                else:
                    st.error("Errore FEM. Controlla la labilità del modello.")