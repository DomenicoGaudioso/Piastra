# 🏗️ FEM Soletta - Analisi Piastre in Calcestruzzo con OpenSeesPy e Streamlit

Un'applicazione web interattiva sviluppata in Python per la modellazione e l'analisi agli Elementi Finiti (FEM) di solette in calcestruzzo. L'applicativo unisce la potenza del solutore strutturale **OpenSeesPy** con la flessibilità dell'interfaccia utente **Streamlit**, offrendo visualizzazioni 2D e 3D avanzate grazie a **Plotly** e **Matplotlib**.

## ✨ Funzionalità Principali

* **Geometria Parametrica:** Definizione rapida della griglia strutturale inserendo semplicemente gli interassi longitudinali e trasversali.
* **Mesh Automatica:** Discretizzazione automatica della piastra in base alla dimensione massima desiderata per l'elemento finito, garantendo la coincidenza dei nodi con le linee di appoggio.
* **Vincoli Dinamici:** Piena libertà nell'assegnazione delle condizioni al contorno:
    * *Appoggi Continui:* Modellazione di setti portanti o travi infinitamente rigide vincolando intere linee della griglia (Assi X o Y).
    * *Appoggi Puntuali:* Modellazione di pilastri vincolando le singole intersezioni.
* **Gestione Carichi Multipli:** * *Aree di carico uniforme* ($q$ in kN/m²) per modellare pesi propri, permanenti e sovraccarichi accidentali.
    * *Impronte di carico concentrato* (Forza totale $F$ su un'area $L \times B$) con ripartizione automatica dei carichi nodali equivalenti.
* **Motore di Calcolo (OpenSeesPy):** Utilizzo di elementi a guscio `ShellMITC4` (spessi) formulati per evitare il *shear locking*.
* **Post-Processing Avanzato:** Calcolo rigoroso delle sollecitazioni interne derivando le curvature della superficie deformata (Teoria delle Piastre Inflesse).
* **Dashboard Interattiva:** * Modello 3D Pre-Analisi (digital twin geometrico).
    * Plot 3D interattivo della Deformata Verticale ($U_z$).
    * Contour Plot (Mappe di calore / Isolinee) per tutte le 8 componenti di sforzo: $N_{xx}$, $N_{yy}$, $N_{xy}$, $M_{xx}$, $M_{yy}$, $M_{xy}$, $V_{xz}$, $V_{yz}$.

## 🚀 Installazione e Avvio

Assicurati di avere Python 3.8+ installato. 

1. **Clona la repository:**
   ```bash
   git clone [https://github.com/tuo-username/fem-soletta-app.git](https://github.com/tuo-username/fem-soletta-app.git)
   cd fem-soletta-app
   ```

2. **Installa le dipendenze:**
   È consigliato l'uso di un ambiente virtuale (es. `venv` o `conda`).
   ```bash
   pip install -r requirements.txt
   ```
   *(Nota: se non hai un file requirements, installa manualmente: `pip install streamlit openseespy numpy matplotlib plotly pandas`)*

3. **Avvia l'applicazione:**
   ```bash
   streamlit run app.py
   ```
   L'app si aprirà automaticamente nel tuo browser predefinito all'indirizzo `http://localhost:8501`.

## 🛠️ Come utilizzare l'app

1. **Definisci la Geometria:** Usa la barra laterale per inserire gli interassi e le proprietà del materiale (Spessore, Modulo di Young, Coefficiente di Poisson).
2. **Applica i Vincoli:** Seleziona quali assi o punti fungono da appoggi. **Attenzione:** un modello senza vincoli sufficienti risulta labile e l'analisi fallirà.
3. **Applica i Carichi:** Aggiungi righe alle tabelle dei carichi uniformi e concentrati direttamente dall'interfaccia.
4. **Visualizza il Modello 3D:** Clicca su "Visualizza Modello 3D Pre-Analisi" per verificare la corretta posizione di nodi bloccati e centri di carico.
5. **Esegui l'Analisi:** Clicca su "Esegui Analisi FEM". Naviga tra le schede per visualizzare gli spostamenti e le mappe tensionali. Usa la spunta "Visualizza Isolinee" per passare da una sfumatura continua a bande di contorno nette.

## 🧮 Dettagli Tecnici (Under the Hood)

Invece di estrarre le sollecitazioni direttamente dai punti di integrazione degli elementi `ShellMITC4` (che in alcune configurazioni di OpenSeesPy possono dare problemi di output omesso), l'applicativo utilizza un approccio matematico robusto basato sugli spostamenti nodali:
* Vengono estratti gli spostamenti globali $U_x, U_y, U_z$ per ogni nodo.
* Tramite differenze finite (`numpy.gradient`), vengono calcolate le curvature della piastra.
* Si applicano le equazioni costitutive della piastra elastica isotropa per ricavare gli sforzi membranali, i momenti flettenti/torcenti e le forze di taglio.

## 🤝 Contribuire
Sentiti libero di aprire una *Issue* per segnalare bug o proporre nuove funzionalità, o inviare una *Pull Request*!

## 📄 Licenza
[Inserisci qui la licenza, es. MIT License]
