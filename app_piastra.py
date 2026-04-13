import streamlit as st
import openseespy.opensees as ops
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go
import pandas as pd
import io
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, white, black
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, Image, PageBreak)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from datetime import datetime

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
# FUNZIONE GENERAZIONE PDF - PIASTRA
# ==========================================

def _draw_contour_matplotlib(x_coords, y_coords, Z, title, cmap='Turbo', unit='mm'):
    """Disegna una mappa di sollecitazione con matplotlib e ritorna BytesIO."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    X, Y = np.meshgrid(x_coords, y_coords)
    v_min = np.min(Z)
    v_max = np.max(Z)
    if abs(v_max - v_min) < 1e-12:
        v_max = v_min + 1
    cs = ax.contourf(X, Y, Z, levels=20, cmap=cmap)
    fig.colorbar(cs, ax=ax, label=f"{title} ({unit})")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title(f"{title} ({unit})", fontsize=11, fontweight='bold')
    ax.set_aspect('equal')
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return buf


def _header_footer_piastra(canvas, doc):
    canvas.saveState()
    pagina = doc.page
    larghezza, altezza = A4
    canvas.setStrokeColor(HexColor("#2C3E50"))
    canvas.setLineWidth(0.8)
    canvas.line(15*mm, altezza - 12*mm, larghezza - 15*mm, altezza - 12*mm)
    canvas.setFont("Helvetica-Bold", 9)
    canvas.setFillColor(HexColor("#2C3E50"))
    canvas.drawString(15*mm, altezza - 10*mm, "Analisi FEM Soletta")
    canvas.setFont("Helvetica", 8)
    canvas.drawRightString(larghezza - 15*mm, altezza - 10*mm, "Relazione di Calcolo")
    canvas.line(15*mm, 10*mm, larghezza - 15*mm, 10*mm)
    canvas.setFont("Helvetica", 7)
    canvas.drawString(15*mm, 6*mm, f"Generato il {datetime.now().strftime('%d/%m/%Y %H:%M')}")
    canvas.drawRightString(larghezza - 15*mm, 6*mm, f"Pagina {pagina}")
    canvas.restoreState()


def genera_pdf_piastra(int_x_str, int_y_str, dim_mesh, spessore, E_cls, nu_cls,
                       x_coords, y_coords, nodi_bloccati,
                       carichi_uni_df, carichi_conc_df,
                       U_z, sforzi, ok_status):
    """Genera la relazione PDF per l'analisi FEM della piastra."""
    buf_pdf = io.BytesIO()
    doc = SimpleDocTemplate(buf_pdf, pagesize=A4,
                            leftMargin=15*mm, rightMargin=15*mm,
                            topMargin=18*mm, bottomMargin=18*mm)

    styles = getSampleStyleSheet()
    style_title = ParagraphStyle('SezTitle', parent=styles['Heading2'],
                                  backColor=HexColor("#2C3E50"), textColor=white,
                                  fontName='Helvetica-Bold', fontSize=12,
                                  spaceAfter=2*mm, spaceBefore=6*mm,
                                  leftIndent=2*mm, borderPadding=(2*mm, 2*mm, 2*mm, 2*mm))
    style_normal = styles['Normal']

    story = []

    # Frontespizio
    story.append(Spacer(1, 25*mm))
    story.append(Paragraph("Analisi FEM Soletta in Calcestruzzo", ParagraphStyle('BigTitle', parent=styles['Title'], fontSize=22, fontName='Helvetica-Bold', textColor=HexColor("#2C3E50"))))
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph("Relazione di Calcolo – Piastra sottile su appoggi", ParagraphStyle('Sub', parent=styles['Normal'], fontSize=14, textColor=HexColor("#555555"))))
    story.append(Spacer(1, 5*mm))
    story.append(Paragraph(f"Data: {datetime.now().strftime('%d/%m/%Y %H:%M')}", style_normal))
    story.append(PageBreak())

    # Sezione 1: Parametri
    story.append(Paragraph("1. Dati di Input", style_title))
    data_rows = [
        [Paragraph("<b>Parametro</b>", style_normal), Paragraph("<b>Valore</b>", style_normal), Paragraph("<b>Unità</b>", style_normal)],
        ["Interasse X", str(int_x_str), "m"],
        ["Interasse Y", str(int_y_str), "m"],
        ["Dimensione mesh", f"{dim_mesh:.2f}", "m"],
        ["Spessore", f"{spessore:.3f}", "m"],
        ["Modulo E", f"{E_cls:.3e}", "N/m²"],
        ["Coefficiente Poisson", f"{nu_cls:.3f}", "—"],
        ["Nodi vincolati", str(len(nodi_bloccati)), "—"],
        ["Nodi totali", str(len(x_coords) * len(y_coords)), "—"],
    ]
    tab_style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), HexColor("#2C3E50")),
        ('TEXTCOLOR', (0, 0), (-1, 0), white),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 9),
        ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
        ('GRID', (0, 0), (-1, -1), 0.4, HexColor("#CCCCCC")),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [white, HexColor("#F2F2F2")]),
    ])
    story.append(Table(data_rows, colWidths=[60*mm, 55*mm, 30*mm], style=tab_style))
    story.append(Spacer(1, 3*mm))

    # Carichi
    story.append(Paragraph("2. Carichi Applicati", style_title))
    carichi_rows = [[Paragraph("<b>Tipo</b>", style_normal), Paragraph("<b>Dettagli</b>", style_normal)]]
    if not carichi_uni_df.empty:
        for _, row in carichi_uni_df.iterrows():
            carichi_rows.append(["Distribuito", f"q={row.get('q_kN_m2',0):.2f} kN/m², area [{row.get('X_min',0):.1f}, {row.get('X_max',0):.1f}] x [{row.get('Y_min',0):.1f}, {row.get('Y_max',0):.1f}]"])
    if not carichi_conc_df.empty:
        for _, row in carichi_conc_df.iterrows():
            carichi_rows.append(["Concentrato", f"F={row.get('Carico_Totale_kN',0):.1f} kN in ({row.get('X',0):.2f}, {row.get('Y',0):.2f}), impronta {row.get('L',0):.2f}x{row.get('B',0):.2f}"])
    if len(carichi_rows) > 1:
        story.append(Table(carichi_rows, colWidths=[30*mm, 130*mm], style=tab_style))
    story.append(Spacer(1, 3*mm))

    # Risultati FEM
    if ok_status == 0 and U_z is not None:
        story.append(Paragraph("3. Risultati FEM", style_title))
        w_max = float(np.max(np.abs(U_z))) * 1000
        res_data = [
            [Paragraph("<b>Grandezza</b>", style_normal), Paragraph("<b>Valore</b>", style_normal)],
            ["Spostamento max |Uz|", f"{w_max:.4f} mm"],
        ]
        for key, unit_m in [('Mxx', 'kNm/m'), ('Myy', 'kNm/m'), ('Mxy', 'kNm/m'),
                             ('Vxz', 'kN/m'), ('Vyz', 'kN/m'),
                             ('Nxx', 'kN/m'), ('Nyy', 'kN/m'), ('Nxy', 'kN/m')]:
            if key in sforzi:
                val = float(np.max(np.abs(sforzi[key])))
                res_data.append([f"Sforzo max |{key}|", f"{val:.4f} {unit_m}"])
        story.append(Table(res_data, colWidths=[80*mm, 65*mm], style=tab_style))
        story.append(Spacer(1, 4*mm))

        # Mappa deformata
        story.append(Paragraph("4. Mappa Deformata Uz", style_title))
        buf_uz = _draw_contour_matplotlib(x_coords, y_coords, U_z * 1000, "Deformata Uz", cmap='Turbo', unit='mm')
        story.append(Image(buf_uz, width=160*mm, height=90*mm))

        # Mappe sollecitazioni (2 per pagina)
        story.append(PageBreak())
        story.append(Paragraph("5. Mappe Sollecitazioni", style_title))
        sforzi_key = ['Mxx', 'Myy', 'Vxz', 'Vyz']
        sforzi_unit = ['kNm/m', 'kNm/m', 'kN/m', 'kN/m']
        sforzi_cmap = ['RdBu_r', 'RdBu_r', 'PiYG', 'PiYG']
        i = 0
        for key, unit_m, cmap in zip(sforzi_key, sforzi_unit, sforzi_cmap):
            if key in sforzi and i < 4:
                buf_sf = _draw_contour_matplotlib(x_coords, y_coords, sforzi[key], key, cmap=cmap, unit=unit_m)
                story.append(Image(buf_sf, width=160*mm, height=80*mm))
                story.append(Spacer(1, 3*mm))
                i += 1
                if i == 2:
                    story.append(Spacer(1, 3*mm))

    else:
        story.append(Paragraph("3. Errore FEM", style_title))
        story.append(Paragraph("<b>L'analisi FEM non è andata a buon fine. Controllare i vincoli.</b>", style_normal))

    doc.build(story, onFirstPage=_header_footer_piastra, onLaterPages=_header_footer_piastra)
    buf_pdf.seek(0)
    return buf_pdf

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

        # --- DOWNLOAD PDF ---
        if 'U_z' in dir() and U_z is not None and ok == 0:
            try:
                pdf_bytes = genera_pdf_piastra(
                    int_x_str, int_y_str, dim_mesh, spessore, E_cls, nu_cls,
                    x_coords, y_coords, nodi_bloccati,
                    carichi_uni_df, carichi_conc_df,
                    U_z, sforzi, ok
                )
                st.download_button(
                    label="📄 Scarica Relazione PDF",
                    data=pdf_bytes,
                    file_name="relazione_piastra_fem.pdf",
                    mime="application/pdf"
                )
            except Exception as e_pdf:
                st.warning(f"Generazione PDF non disponibile: {e_pdf}")