"""
PPIIMoL - PPII Helix Detector for PyMOL
---------------------------------------

This PyMOL plugin detects polyproline II (PPII) helices in protein structures using dihedral angles (phi/psi)
and helps visualize relevant interactions with non-canonical hydrogen bonds.

Main Features:
- Load PDB files or fetch proteins by ID.
- Calculate phi/psi angles for backbone atoms.
- Automatically identify candidate PPII helices.
- Visualize key atom distances and angles (CA–H–O).
- Export CSV and PDB reports of detected segments.
- User-friendly interface with Tkinter.

Usage:
- Run this script from PyMOL: `run PPIIMoL.py`
- Launch the GUI with: `lanzar_interfaz()`
- Alternatively, use core functions like: `detectar_segmentos_ppii("your_object")`

Dependencies:
- PyMOL with Python support
- Python 3.x
- Tkinter

Authors:
- Silvia Enma (2025), Instituto de Química-Física "Blas Cabrera" (CSIC)
- GitHub: https://github.com/silviaenma/ppii-detector-pymol

License:
- GNU GPL v3

Note:
- This script was originally written in Spanish for internal lab use, but has been translated for the PyMOL community.
"""

def __init_plugin__(app=None):
    """
    PyMOL plugin loader function. Called when the plugin is loaded from the Plugin Manager.
    """
    try:
        from tkinter import messagebox
        messagebox.showinfo("PPIIMoL", "PPIIMoL plugin loaded successfully.")
    except ImportError:
        print("PPIIMoL plugin loaded (no GUI message shown, Tkinter not available).")

        
from pymol import cmd, stored
import tkinter as tk
import math
from pathlib import Path
from datetime import datetime
import os
from tkinter import filedialog, messagebox
import tkinter as tk
from tkinter import ttk, messagebox
import os
import csv


# Valores por defecto (se actualizarán cuando el usuario los cambie)
tol_phi_global = 20.0
tol_psi_global = 20.0
min_ang_global = 110.0
max_ang_global = 180.0


pdb_file = None
segmentos_ppii_global = []  # Variable global para guardar segmentos PPII detectados
localizar_atomos=None
distancias=None
angulos_v=None
objetos = None
pares = None

# ==========================================================
# >>>>>>>>    FUNCIÓN: OBTENER_CARPETA_RESULTADOS    <<<<<<<<
# ==========================================================

"""
obtener_carpeta_resultados(nombre_carpeta='Resultados_PyMOL')
----------------------------------------------------------
Descripción:
    Crea una carpeta con fecha actual dentro del directorio Documentos
    (o equivalente) del usuario para guardar los resultados generados
    por el módulo.

Parámetros:
    nombre_carpeta : str, opcional
        Nombre base para la carpeta de resultados (por defecto 'Resultados_PyMOL').

Retorno:
    pathlib.Path
        Ruta completa a la carpeta creada.

Ejemplo de uso:
    >>> carpeta = obtener_carpeta_resultados("Resultados_PPIIMoL")
    >>> print(carpeta)
    /home/usuario/Documentos/Resultados_PPIIMoL_2025-07-05
----------------------------------------------------------
"""

def obtener_carpeta_resultados(nombre_carpeta="Resultados_PyMOL"):
    # Formatear fecha y hora actual (compatible con ambos sistemas)
    fecha_hora = datetime.now().strftime("%Y-%m-%d")
    
    # Usar ~/Documents o ~/Documentos según el sistema
    documentos = Path.home() / ("Documentos" if os.getenv('LANG', '').startswith('es_') else "Documents")
    
    # Si no existe, crear directorio en HOME directamente
    if not documentos.exists():
        documentos = Path.home()
    
    # Crear nombre de carpeta combinado
    nombre_completo = f"{nombre_carpeta}_{fecha_hora}"
    carpeta_resultados = documentos / nombre_completo
    
    # Crear la carpeta (con parents=True por si faltan directorios padres)
    carpeta_resultados.mkdir(parents=True, exist_ok=True)
    
    return carpeta_resultados

# ==========================================================
# >>>>>>>>    FUNCIÓN: SELECCIONAR_ARCHIVO    <<<<<<<<
# ==========================================================

"""
seleccionar_archivo()
----------------------------------------------------------
Descripción:
    Abre un cuadro de diálogo para que el usuario seleccione un archivo PDB
    de entrada. La ruta seleccionada se guarda en la variable global pdb_file.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> seleccionar_archivo()
----------------------------------------------------------
"""

def seleccionar_archivo():
    global pdb_file
    pdb_file = filedialog.askopenfilename(
        title="Selecciona un archivo PDB",
        filetypes=[("Archivos PDB", "*.pdb")]
    )
    if pdb_file:
        cmd.reinitialize()
        cmd.load(pdb_file, "proteina")
        cmd.hide("everything", "proteina")  # Ocultar todo lo visible
        cmd.show("licorice", "proteina")   # Mostrar como licorice
        messagebox.showinfo("Archivo cargado", f"Se cargó:\n{pdb_file}")

# ==========================================================
# >>>>>>>>    FUNCIÓN: DESCARGAR_MOLECULA    <<<<<<<<
# ==========================================================

"""
descargar_molecula()
----------------------------------------------------------
Descripción:
    Abre una ventana emergente para permitir al usuario introducir un ID
    de proteína del Protein Data Bank (PDB) y descargar automáticamente
    la estructura correspondiente. Una vez descargada, la molécula se 
    carga en PyMOL, se ocultan todas las representaciones gráficas 
    actuales y se muestra la proteína usando el estilo "licorice".
    La ruta o ID de la proteína descargada se guarda en la variable global
    pdb_file.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> descargar_molecula()
----------------------------------------------------------
"""

def descargar_molecula():
    global pdb_file
    def fetch_pdb():
        global pdb_file
        pdb_id = entry.get().strip()
        if pdb_id:
            try:
                cmd.reinitialize()
                cmd.fetch(pdb_id, name="proteina")
                cmd.hide("everything", "proteina")
                cmd.show("licorice", "proteina")
                # Aquí asignamos pdb_file con un valor para simular archivo cargado
                pdb_file = pdb_id  # Solo guardamos el ID para que funcione la lógica de “selección”
                messagebox.showinfo("Descarga completa", f"Molécula {pdb_id.upper()} descargada y cargada.")
                fetch_window.destroy()
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo descargar {pdb_id}: {e}")
        else:
            messagebox.showwarning("Advertencia", "Por favor ingresa un ID de PDB.")
    
    fetch_window = tk.Toplevel()
    fetch_window.title("Descargar proteína")
    tk.Label(fetch_window, text="ID de la proteína (PDB):").pack(pady=5)
    entry = tk.Entry(fetch_window, width=20)
    entry.pack(pady=5)
    tk.Button(fetch_window, text="Descargar", command=fetch_pdb).pack(pady=5)

# ==========================================================
# >>>>>>>>    FUNCIÓN: ANADIR_HIDROGENOS    <<<<<<<<
# ==========================================================

"""
anadir_hidrogenos()
----------------------------------------------------------
Descripción:
    Añade átomos de hidrógeno a todas las moléculas cargadas
    en la sesión actual de PyMOL. Una vez añadidos, la 
    estructura se reorganiza para optimizar la visualización 
    y se muestra usando el estilo "licorice". Si no hay 
    ningún archivo PDB cargado previamente, se muestra una 
    advertencia al usuario.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> anadir_hidrogenos()
----------------------------------------------------------
"""

def anadir_hidrogenos():
    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return
    cmd.h_add("all")
    cmd.sort("all extend 1")
    cmd.show("licorice", "all")
    messagebox.showinfo("Hidrógenos", "Hidrógenos añadidos y mostrados como licorice.")

# ==========================================================
# >>>>>>>>    FUNCIÓN: ELIMINAR_SOLVENTES    <<<<<<<<
# ==========================================================

"""
eliminar_solventes()
----------------------------------------------------------
Descripción:
    Elimina todas las moléculas de solvente presentes en la 
    estructura cargada en PyMOL. Tras la eliminación, se 
    muestra un mensaje informativo al usuario confirmando 
    la acción realizada.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> eliminar_solventes()
----------------------------------------------------------
"""

def eliminar_solventes():
    cmd.remove("solvent")
    messagebox.showinfo("Solventes", "Solventes eliminados.")

# ==========================================================
# >>>>>>>>    FUNCIÓN: OCULTAR_SIDE_CHAINS    <<<<<<<<
# ==========================================================

"""
ocultar_side_chains()
----------------------------------------------------------
Descripción:
    Oculta todas las cadenas laterales de los residuos de 
    la proteína cargada en PyMOL, dejando visibles únicamente 
    los átomos del esqueleto principal (backbone: N, CA, C, O). 
    Esta función permite centrar la visualización en la 
    estructura principal de la proteína sin distracciones.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> ocultar_side_chains()
----------------------------------------------------------
"""

def ocultar_side_chains():
    cmd.hide("everything", "proteina and not name N+CA+C+O")
    messagebox.showinfo("Backbone", "Cadenas laterales ocultas (solo backbone).")

# ==========================================================
# >>>>>>>>    FUNCIÓN: SEPARAR_CADENAS    <<<<<<<<
# ==========================================================

"""
separar_cadenas()
----------------------------------------------------------
Descripción:
    Separa cada cadena de la proteína cargada en PyMOL 
    en objetos individuales. Cada nueva cadena se guarda 
    como un objeto separado con un nombre identificador 
    en el formato "cadena_X", donde X es el identificador 
    de la cadena original. Esta función facilita la 
    manipulación y análisis individual de las cadenas.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> separar_cadenas()
----------------------------------------------------------
"""

def separar_cadenas():
    stored.chains = []
    cmd.iterate("proteina", "stored.chains.append(chain)")
    for cadena in set(stored.chains):
        nuevo_objeto = f"cadena_{cadena}"
        cmd.create(nuevo_objeto, f"proteina and chain {cadena}")

# ==========================================================
# >>>>>>>>    FUNCIÓN: OBTENER_ANGULOS_PHI_PSI_POR_CADENA    <<<<<<<<
# ==========================================================

"""
obtener_angulos_phi_psi_por_cadena(objeto="proteina")
----------------------------------------------------------
Descripción:
    Calcula los ángulos torsionales phi (ϕ) y psi (ψ) de cada 
    residuo en todas las cadenas del objeto especificado 
    cargado en PyMOL. Los resultados se devuelven en un 
    diccionario que asocia cada residuo con su nombre, número 
    y los valores de los ángulos phi y psi correspondientes.

Parámetros:
    objeto : str, opcional
        Nombre del objeto cargado en PyMOL sobre el cual se 
        calcularán los ángulos (por defecto "proteina").

Retorno:
    dict
        Diccionario con las claves como tuplas (cadena, número 
        de residuo) y valores como tuplas (nombre del residuo, 
        ángulo phi, ángulo psi).

Ejemplo de uso:
    >>> resultados = obtener_angulos_phi_psi_por_cadena()
    >>> print(resultados)
----------------------------------------------------------
"""

def obtener_angulos_phi_psi_por_cadena(objeto="proteina"):
    resultados = {}
    chains = cmd.get_chains(objeto)
    for chain in chains:
        sel_ca = f"{objeto} and chain {chain} and name CA"
        phipsi = cmd.get_phipsi(sel_ca)
        if not phipsi:
            continue
        for (obj, idx), (phi, psi) in sorted(phipsi.items()):
            if phi is None or psi is None:
                continue
            stored.info = []
            cmd.iterate(f"({obj}`{idx})", "stored.info.append((chain, resn, resi))", space={'stored': stored})
            if not stored.info:
                continue
            ch, resn, resi = stored.info[0]
            resultados[(ch, resi)] = (resn, phi, psi)
    return resultados

# ==========================================================
# >>>>>>>>    FUNCIÓN: GUARDAR_CSV_ANGULOS_PHI_PSI    <<<<<<<<
# ==========================================================

"""
guardar_csv_angulos_phi_psi()
----------------------------------------------------------
Descripción:
    Genera un archivo CSV con los ángulos torsionales phi (ϕ) 
    y psi (ψ) calculados para cada residuo de la proteína 
    cargada en PyMOL. El archivo se guarda en una carpeta 
    de resultados con el nombre "angulos_phi_psi.csv" e 
    incluye información de la cadena, el residuo, su número 
    y los valores de los ángulos.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> guardar_csv_angulos_phi_psi()
----------------------------------------------------------
"""

def guardar_csv_angulos_phi_psi():
    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return

    phi_map = obtener_angulos_phi_psi_por_cadena("proteina")

    stored.res_list = []
    cmd.iterate("proteina and name CA", "stored.res_list.append((chain, resn, resi))")

    datos_csv = [("Cadena", "Residuo", "Número", "Phi", "Psi")]

    for chain, resn, resi in sorted(stored.res_list, key=lambda x: (x[0], int(x[2]))):
        key = (chain, resi)
        if key in phi_map:
            resn_val, phi, psi = phi_map[key]
            datos_csv.append((chain, resn, resi, f"{phi:.2f}", f"{psi:.2f}"))
    carpeta= obtener_carpeta_resultados()
    # Guardar en CSV con separador ;
    if len(datos_csv) > 1:
        ruta_csv = os.path.join(os.getcwd(), carpeta/"angulos_phi_psi.csv")
        with open(ruta_csv, mode="w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter=";")
            writer.writerows(datos_csv)
    messagebox.showinfo("Éxito", f"CSV generado correctamente:\n{ruta_csv}")

# ==========================================================
# >>>>>>>>    FUNCIÓN: GENERAR_REPORTE_CSV    <<<<<<<<
# ==========================================================

"""
generar_reporte_csv(segmentos_ppii, max_dist=5.0, nombre_archivo="reporte_ppii.csv", min_ang=110.0, max_ang=180.0)
----------------------------------------------------------
Descripción:
    Genera un archivo CSV con un informe detallado sobre los 
    segmentos de hélices PPII detectados en la proteína cargada 
    en PyMOL. El reporte incluye pares de átomos Cα-O colindantes, 
    sus distancias, los ángulos CA-H-O calculados y los valores 
    validados según los rangos definidos. El archivo se guarda 
    en una carpeta de resultados con el nombre especificado.

Parámetros:
    segmentos_ppii : list
        Lista de segmentos PPII detectados en la proteína.
    max_dist : float, opcional
        Distancia máxima (en Ångström) para considerar los 
        pares de átomos colindantes. Por defecto es 5.0 Å.
    nombre_archivo : str, opcional
        Nombre del archivo CSV que se generará. Por defecto 
        es "reporte_ppii.csv".
    min_ang : float, opcional
        Ángulo mínimo en grados para filtrar los ángulos CA-H-O 
        válidos. Por defecto es 110.0°.
    max_ang : float, opcional
        Ángulo máximo en grados para filtrar los ángulos CA-H-O 
        válidos. Por defecto es 180.0°.

Retorno:
    None

Ejemplo de uso:
    >>> generar_reporte_csv(segmentos_ppii)
----------------------------------------------------------
"""

def generar_reporte_csv(segmentos_ppii, max_dist=5.0, nombre_archivo="reporte_ppii.csv", min_ang=110.0, max_ang=180.0):
    if not segmentos_ppii:
        messagebox.showwarning("Advertencia", "No hay segmentos PPII detectados.")
        return

    objetos_por_segmento = localizar_atomos
    pares_ca_o=distancias
    angulos_validos=angulos_v
    if objetos_por_segmento ==None:
        objetos_por_segmento = localizar_atomicos_clave_segmentos(segmentos_ppii)
    if pares_ca_o is None:
        pares_ca_o = calcular_distancias_colindantes(objetos_por_segmento, max_dist=max_dist)
    
    # PASAMOS LOS ÁNGULOS AQUÍ
    if angulos_validos is None:
        angulos_validos = calcular_angulos_ca_h_o(pares_ca_o, min_ang=min_ang, max_ang=max_ang)

    dist_dict = {}
    for ca, o in pares_ca_o:
        if not (cmd.count_atoms(ca) and cmd.count_atoms(o)):
            continue
        coord_ca = cmd.get_coords(ca)[0]
        coord_o = cmd.get_coords(o)[0]
        dist = math.sqrt(sum((a - b)**2 for a, b in zip(coord_ca, coord_o)))
        dist_dict[(ca, o)] = dist

    ang_dict = {(ca, o): ang for ca, _, o, ang in angulos_validos}

    datos_csv = [["Átomo Cα", "Átomo O", "Distancia (Å)", "Ángulo CA-H-O (°)"]]
    for (ca, o), dist in dist_dict.items():
        ang = ang_dict.get((ca, o), None)
        if ang is not None:
            datos_csv.append([ca, o, f"{dist:.2f}", f"{ang:.2f}"])

    ruta_csv = obtener_carpeta_resultados() / nombre_archivo
    with open(ruta_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerows(datos_csv)

    messagebox.showinfo("Éxito", f"CSV generado correctamente:\n{ruta_csv}")

# ==========================================================
# >>>>>>>>    FUNCIÓN: LOCALIZAR_ATOMICOS_CLAVE_SEGMENTOS    <<<<<<<<
# ==========================================================

"""
localizar_atomicos_clave_segmentos(segmentos)
----------------------------------------------------------
Descripción:
    Localiza y marca en PyMOL los átomos clave (Cα y O) de 
    cada segmento PPII detectado. Crea pseudoátomos en las 
    posiciones correspondientes para facilitar la visualización 
    y análisis. Los átomos Cα se marcan en color azul y los O 
    en color rojo, agrupándolos bajo el nombre "atomos_clave".

Parámetros:
    segmentos : list
        Lista de segmentos PPII detectados en la proteína.

Retorno:
    list
        Lista de listas, donde cada sublista contiene los 
        nombres de los pseudoátomos creados para un segmento.

Ejemplo de uso:
    >>> localizar_atomicos_clave_segmentos(segmentos_ppii)
----------------------------------------------------------
"""

def localizar_atomicos_clave_segmentos(segmentos):
    global localizar_atomos
    cmd.delete("esfera_*")
    objetos_por_segmento = []
    cmd.h_add()

    for idx, seg in enumerate(segmentos, start=1):
        objetos_segmento = []

        for (resn, resi, chain, _, _) in seg:
            sele_base = f"proteina and chain {chain} and resi {resi}"
            stored.coords = []

            # Átomo CA (carbono alfa)
            stored.ca_coords = []
            cmd.iterate_state(
                1,
                f"{sele_base} and name CA",
                "stored.ca_coords.append((x, y, z))",
                space={'stored': stored}
            )

            # Oxígeno carbonilo
            stored.o_coords = []
            cmd.iterate_state(
                1,
                f"{sele_base} and name O",
                "stored.o_coords.append((x, y, z))",
                space={'stored': stored}
            )

            # Crear pseudoatomos con nombres consistentes
            if stored.ca_coords:
                x, y, z = stored.ca_coords[0]
                esfera_name = f"esfera_s{idx}_CA_{resn}_{resi}_{chain}"
                cmd.pseudoatom(esfera_name, pos=[x, y, z])
                cmd.set("sphere_scale", 0.3, esfera_name)
                cmd.color("blue", esfera_name)
                objetos_segmento.append(esfera_name)

            if stored.o_coords:
                x, y, z = stored.o_coords[0]
                esfera_name = f"esfera_s{idx}_O_{resn}_{resi}_{chain}"
                cmd.pseudoatom(esfera_name, pos=[x, y, z])
                cmd.set("sphere_scale", 0.3, esfera_name)
                cmd.color("red", esfera_name)
                objetos_segmento.append(esfera_name)

        objetos_por_segmento.append(objetos_segmento)

    cmd.group("atomos_clave", "esfera_*")
    localizar_atomos = objetos_por_segmento
    return objetos_por_segmento

# ==========================================================
# >>>>>>>>    FUNCIÓN: CALCULAR_DISTANCIAS_COLINDANTES    <<<<<<<<
# ==========================================================

"""
calcular_distancias_colindantes(objetos_por_segmento, max_dist=5.0, archivo_salida="distancias_colindantes.txt")
----------------------------------------------------------
Descripción:
    Calcula las distancias entre los átomos Cα y O de segmentos 
    PPII consecutivos o cercanos. Identifica los pares cuya 
    distancia es menor al valor máximo especificado. Guarda los 
    resultados en un archivo de texto y devuelve la lista de 
    pares de átomos que cumplen el criterio.

Parámetros:
    objetos_por_segmento : list
        Lista de listas con los nombres de los pseudoátomos 
        (Cα y O) de cada segmento.
    max_dist : float, opcional
        Distancia máxima (en Ångström) para considerar los pares 
        Cα-O como colindantes. Por defecto es 5.0 Å.
    archivo_salida : str, opcional
        Nombre del archivo de texto donde se guardarán las 
        distancias calculadas. Por defecto es 
        "distancias_colindantes.txt".

Retorno:
    list
        Lista de tuplas con los nombres de los pares de átomos 
        (Cα, O) cuya distancia es menor a max_dist.

Ejemplo de uso:
    >>> pares = calcular_distancias_colindantes(objetos_por_segmento)
----------------------------------------------------------
"""

def calcular_distancias_colindantes(objetos_por_segmento, max_dist=5.0, archivo_salida="distancias_colindantes.txt"):
    global distancias
    pares_ca_o = []
    carpeta = obtener_carpeta_resultados()
    archivo_salida = carpeta / archivo_salida
    coordenadas = {}

    for segmento in objetos_por_segmento:
        for obj in segmento:
            coords = cmd.get_coords(obj)
            if coords is None or len(coords) == 0:
                continue
            coordenadas[obj] = coords[0]

    def es_CA(obj): return "_CA_" in obj
    def es_O(obj): return "_O_" in obj

    max_dist_sq = max_dist ** 2

    for salto in [1, 2]:
        for i in range(len(objetos_por_segmento) - salto):
            seg1 = objetos_por_segmento[i]
            seg2 = objetos_por_segmento[i + salto]

            ca1 = [obj for obj in seg1 if es_CA(obj)]
            o2 = [obj for obj in seg2 if es_O(obj)]
            ca2 = [obj for obj in seg2 if es_CA(obj)]
            o1 = [obj for obj in seg1 if es_O(obj)]

            # Comparar ca1 con o2
            for ca in ca1:
                coord_ca = coordenadas.get(ca)
                if coord_ca is None:
                    continue
                for o in o2:
                    coord_o = coordenadas.get(o)
                    if coord_o is None:
                        continue
                    dist_sq = sum((a - b) ** 2 for a, b in zip(coord_ca, coord_o))
                    if dist_sq < max_dist_sq:
                        pares_ca_o.append((ca, o))

            # Comparar ca2 con o1
            for ca in ca2:
                coord_ca = coordenadas.get(ca)
                if coord_ca is None:
                    continue
                for o in o1:
                    coord_o = coordenadas.get(o)
                    if coord_o is None:
                        continue
                    dist_sq = sum((a - b) ** 2 for a, b in zip(coord_ca, coord_o))
                    if dist_sq < max_dist_sq:
                        pares_ca_o.append((ca, o))

    with open(archivo_salida, "w") as f:
        f.write(f"Pares CA-O con distancia < {max_dist} Å:\n")
        for i, (ca, o) in enumerate(pares_ca_o):
            dist = math.sqrt(sum((a - b)**2 for a, b in zip(coordenadas[ca], coordenadas[o])))
            f.write(f"{i:03d} | {ca} - {o} : {dist:.2f} Å\n")

    print(f"[DEBUG] Total de pares CA-O guardados: {len(pares_ca_o)}")
    distancias= pares_ca_o
    return pares_ca_o

# ==========================================================
# >>>>>>>>    FUNCIÓN: CALCULAR_ANGULOS_CA_H_O    <<<<<<<<
# ==========================================================

"""
calcular_angulos_ca_h_o(pares_ca_o, min_ang=110.0, max_ang=180.0, archivo_salida="angulos_ca_h_o.txt")
----------------------------------------------------------
Descripción:
    Calcula los ángulos formados entre los átomos Cα-H-O para 
    cada par de átomos Cα-O proporcionado. Filtra y guarda 
    únicamente los ángulos que estén dentro del rango definido 
    por min_ang y max_ang. Los resultados se almacenan en un 
    archivo de texto y se devuelven como una lista.

Parámetros:
    pares_ca_o : list
        Lista de tuplas con los nombres de los pares de átomos 
        (Cα, O) sobre los que se calcularán los ángulos.
    min_ang : float, opcional
        Ángulo mínimo (en grados) para considerar válido un 
        ángulo CA-H-O. Por defecto es 110.0°.
    max_ang : float, opcional
        Ángulo máximo (en grados) para considerar válido un 
        ángulo CA-H-O. Por defecto es 180.0°.
    archivo_salida : str, opcional
        Nombre del archivo de texto donde se guardarán los 
        ángulos calculados. Por defecto es "angulos_ca_h_o.txt".

Retorno:
    list
        Lista de tuplas con la información de cada ángulo válido 
        en el formato (Cα, H, O, ángulo).

Ejemplo de uso:
    >>> angulos = calcular_angulos_ca_h_o(pares_ca_o)
----------------------------------------------------------
"""

def calcular_angulos_ca_h_o(pares_ca_o, min_ang=110.0, max_ang=180.0, archivo_salida="angulos_ca_h_o.txt"):
    global angulos_v
    angulos_validos = []
    archivo_salida = obtener_carpeta_resultados() / archivo_salida
    archivo_salida.parent.mkdir(parents=True, exist_ok=True)

    with open(archivo_salida, "w") as f:
        f.write("CA - H - O : Ángulo (grados)\n")

        for ca, o in pares_ca_o:
            if not (cmd.count_atoms(ca) and cmd.count_atoms(o)):
                continue

            # Extraer resn, resi y chain
            partes = ca.split("_")
            if len(partes) < 4:
                f.write(f"# Nombre inválido para {ca}\n")
                continue

            resn = partes[-3]
            resi = partes[-2]
            chain = partes[-1]

            # Buscar Hs sin usar "model"
            seleccion = f"resn {resn} and resi {resi} and chain {chain} and name H*"
            try:
                hidrogenos = cmd.get_model(seleccion).atom
            except Exception as e:
                f.write(f"# Error al buscar H para {ca} ({seleccion}): {e}\n")
                continue

            if not hidrogenos:
                f.write(f"# No se encontraron H para {ca} (búsqueda: {seleccion})\n")
                continue

            for h in hidrogenos:
                h_sel = f"/{h.model}//{h.chain}/{h.resi}/{h.name}"
                h_nombre = f"{h.resn}-{h.resi}{h.chain}_{h.name}"
                try:
                    ang = cmd.get_angle(ca, h_sel, o)
                    if min_ang <= abs(ang) <= max_ang:
                        angulos_validos.append((ca, h_nombre, o, ang))
                        f.write(f"{ca} - {h_nombre} - {o} : {ang:.2f}°\n")
                except Exception as e:
                    f.write(f"# Error al calcular ángulo para {ca}, {h_sel}, {o} : {e}\n")
    angulos_v= angulos_validos
    return angulos_validos


# ==========================================================
# >>>>>>>>    FUNCIÓN: CALCULAR_Y_GUARDAR_ANGULOS    <<<<<<<<
# ==========================================================

"""
calcular_y_guardar_angulos()
----------------------------------------------------------
Descripción:
    Calcula los ángulos CA-H-O para los segmentos PPII previamente 
    detectados en la proteína cargada en PyMOL. Localiza los átomos 
    clave, obtiene las distancias colindantes y calcula los ángulos 
    que cumplen con el rango de validez (110°-180°). Los resultados 
    se guardan en un archivo de texto llamado "angulos_ca_h_o.txt" 
    dentro de la carpeta de resultados.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> calcular_y_guardar_angulos()
----------------------------------------------------------
"""


def calcular_y_guardar_angulos():
    global segmentos_ppii_global
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return
    
    # Obtener los objetos atómicos clave
    objetos = localizar_atomos
    if objetos is None:
        objetos = localizar_atomicos_clave_segmentos(segmentos_ppii_global)
    
    # Calcular distancias para obtener las tripletas
    tripletas = calcular_distancias_colindantes(objetos)
    
    # Calcular y guardar ángulos
    carpeta = obtener_carpeta_resultados()
    archivo_salida = carpeta / "angulos_ca_h_o.txt"
    angulos = calcular_angulos_ca_h_o(tripletas, str(archivo_salida))  # Convertir a string para PyMOL
    
    if angulos:
        messagebox.showinfo("Éxito", f"Ángulos CA-H-O guardados en:\n{archivo_salida}")
    else:
        messagebox.showwarning("Aviso", "No se encontraron ángulos válidos (110°-180°)")


# ==========================================================
# >>>>>>>>    FUNCIÓN: VISUALIZAR_DISTANCIAS_PARES    <<<<<<<<
# ==========================================================
"""
visualizar_distancias_pares(pares_candidatos)
----------------------------------------------------------
Descripción:
    Dibuja y visualiza en PyMOL las distancias entre los pares 
    de átomos proporcionados. Cada distancia se representa como 
    una línea discontinua (dashed line) de color cian para 
    facilitar la identificación de interacciones potenciales.

Parámetros:
    pares_candidatos : list
        Lista de tuplas con los nombres de los pares de átomos 
        (atomo1, atomo2) que se desean visualizar.

Retorno:
    None

Ejemplo de uso:
    >>> visualizar_distancias_pares(pares_candidatos)
----------------------------------------------------------
"""


def visualizar_distancias_pares(pares_candidatos):
    if not pares_candidatos:
        messagebox.showinfo("Visualización", "No hay pares de átomos para visualizar.")
        return

    cmd.delete("distancia_ppii")

    for i, (at1, at2) in enumerate(pares_candidatos, start=1):
        nombre_dist = f"distancia_ppii_{i}"
        cmd.distance(nombre_dist, at1, at2)
        cmd.set("dash_width", 4, nombre_dist)
        cmd.set("dash_length", 0.5, nombre_dist)
        cmd.color("cyan", nombre_dist)

    messagebox.showinfo("Visualización", f"Visualizados {len(pares_candidatos)} pares de distancias colindantes.")

# ==========================================================
# >>>>>>>>    FUNCIÓN: DETECTAR_SEGMENTOS_PPII    <<<<<<<<
# ==========================================================

"""
detectar_segmentos_ppii(objeto="proteina", min_length=3, tol_phi=20.0, tol_psi=20.0, max_saltos=0)
----------------------------------------------------------
Descripción:
    Detecta segmentos de hélices de tipo poliprolina II (PPII) 
    en la proteína cargada en PyMOL. Analiza los ángulos 
    torsionales phi (ϕ) y psi (ψ) de cada residuo para 
    identificar regiones consecutivas compatibles con la 
    conformación PPII según los criterios de tolerancia 
    definidos. Permite configurar la longitud mínima de los 
    segmentos y el número máximo de saltos permitidos entre 
    residuos consecutivos.

Parámetros:
    objeto : str, opcional
        Nombre del objeto cargado en PyMOL sobre el que se 
        realizará la detección. Por defecto es "proteina".
    min_length : int, opcional
        Longitud mínima (en número de residuos) que debe tener 
        un segmento para ser considerado PPII. Por defecto es 3.
    tol_phi : float, opcional
        Tolerancia en grados para el ángulo phi respecto a 
        los valores característicos de la conformación PPII. 
        Por defecto es 20.0°.
    tol_psi : float, opcional
        Tolerancia en grados para el ángulo psi respecto a 
        los valores característicos de la conformación PPII. 
        Por defecto es 20.0°.
    max_saltos : int, opcional
        Número máximo de residuos consecutivos que pueden 
        incumplir los criterios de PPII dentro de un segmento 
        sin interrumpirlo. Por defecto es 0.

Retorno:
    None

Ejemplo de uso:
    >>> detectar_segmentos_ppii()
----------------------------------------------------------
"""

def detectar_segmentos_ppii(objeto="proteina", min_length=3, tol_phi=20.0, tol_psi=20.0, max_saltos=0):
    global segmentos_ppii_global

    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return

    phi_psi_map = obtener_angulos_phi_psi_por_cadena(objeto)
    lista_residuos = []
    for (chain, resi), (resn, phi, psi) in phi_psi_map.items():
        try:
            resi_num = int(resi)
        except:
            continue
        lista_residuos.append((chain, resi_num, resn, phi, psi))
    lista_residuos.sort(key=lambda x: (x[0], x[1]))

    segmentos = []
    segmento_actual = []
    saltos_restantes = max_saltos  # Contador de saltos permitidos

    def en_rango_ppii(phi, psi, tol_phi=20.0, tol_psi=20.0):
    # Valores IDEALES para PPII: φ = -75°, ψ = +145°
        return (abs(phi - (-75)) <= tol_phi) and (abs(psi - 145) <= tol_psi)

    for i, (chain, resi, resn, phi, psi) in enumerate(lista_residuos):
        if en_rango_ppii(phi, psi, tol_phi, tol_psi):
            # Si cumple, reiniciamos el contador de saltos
            saltos_restantes = max_saltos
            if not segmento_actual:
                segmento_actual.append((resn, resi, chain, phi, psi))
            else:
                _, last_resi, last_chain, _, _ = segmento_actual[-1]
                if chain == last_chain and resi == last_resi + 1:
                    segmento_actual.append((resn, resi, chain, phi, psi))
                else:
                    if len(segmento_actual) >= min_length:
                        segmentos.append(segmento_actual)
                    segmento_actual = [(resn, resi, chain, phi, psi)]
        else:
            # Si no cumple pero hay saltos disponibles
            if segmento_actual and saltos_restantes > 0:
                _, last_resi, last_chain, _, _ = segmento_actual[-1]
                if chain == last_chain and resi == last_resi + 1:
                    segmento_actual.append((resn, resi, chain, phi, psi))
                    saltos_restantes -= 1
                    continue
            
            # Si no hay saltos disponibles o no es consecutivo
            if len(segmento_actual) >= min_length:
                segmentos.append(segmento_actual)
            segmento_actual = []
            saltos_restantes = max_saltos

    # Añadir el último segmento si cumple con la longitud mínima
    if len(segmento_actual) >= min_length:
        segmentos.append(segmento_actual)

    if not segmentos:
        messagebox.showinfo("Resultado", "No se encontraron segmentos PPII.")
        return

    cmd.delete("ppii_segmento*")

    salida = f"Segmentos candidatos a hélices PPII (saltos permitidos: {max_saltos}):\n"
    for idx, seg in enumerate(segmentos, start=1):
        start_resi = seg[0][1]
        end_resi = seg[-1][1]
        chain = seg[0][2]
        salida += f"\nSegmento {idx} (Cadena {chain}, residuos {start_resi}-{end_resi}, longitud {len(seg)}):\n"
        for (resn, resi, _, phi, psi) in seg:
            salida += f"  {resn}-{resi}{chain}: (phi={phi:.1f}, psi={psi:.1f})\n"

        sel_str = f"proteina and chain {chain} and resi {start_resi}-{end_resi}"
        obj_name = f"ppii_segmento_{chain}_{start_resi}_{end_resi}"
        cmd.create(obj_name, sel_str)
        cmd.color("red", obj_name)
        cmd.show("cartoon", obj_name)

    carpeta = obtener_carpeta_resultados()
    ruta_archivo = carpeta / "segmentos_ppii.txt"
    with open(ruta_archivo, "w") as f:
        f.write(salida)

    segmentos_ppii_global = segmentos

    messagebox.showinfo("Éxito", f"{len(segmentos)} segmentos PPII detectados (con {max_saltos} saltos permitidos).\n"
                               f"Átomos clave visualizados en PyMOL.")

# ==========================================================
# >>>>>>>>    FUNCIÓN: GUARDAR_SEGMENTO_PPII_PDB    <<<<<<<<
# ==========================================================

"""
guardar_segmento_ppii_pdb(segmento, nombre_archivo="segmento_ppii.pdb")
----------------------------------------------------------
Descripción:
    Guarda en un archivo PDB un segmento específico de hélice 
    PPII previamente detectado en la proteína cargada en PyMOL. 
    El archivo se genera con el nombre especificado dentro de la 
    carpeta de resultados.

Parámetros:
    segmento : list
        Lista con los residuos que forman el segmento PPII a guardar.
    nombre_archivo : str, opcional
        Nombre del archivo PDB que se generará. Por defecto es 
        "segmento_ppii.pdb".

Retorno:
    None

Ejemplo de uso:
    >>> guardar_segmento_ppii_pdb(segmentos_ppii[0])
----------------------------------------------------------
"""

def guardar_segmentos_ppii_pdb():
    global segmentos_ppii_global
    carpeta = obtener_carpeta_resultados()
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return

    # Para cada segmento almacenado, construimos el nombre de objeto y hacemos cmd.save
    count = 0
    for seg in segmentos_ppii_global:
        start_resi = seg[0][1]
        end_resi = seg[-1][1]
        chain = seg[0][2]
        obj_name = f"ppii_segmento_{chain}_{start_resi}_{end_resi}"
        # Comprobamos que el objeto efectivamente exista en PyMOL
        if cmd.count_atoms(f"{obj_name}") > 0:
            
            filename = carpeta / f"{obj_name}.pdb"
            try:
                cmd.save(filename, obj_name)
                count += 1
            except Exception as e:
                print(f"Error guardando {obj_name}: {e}")
        else:
            print(f"Objeto {obj_name} no encontrado en la sesión de PyMOL.")

    if count > 0:
        messagebox.showinfo("Éxito", f"Se guardaron {count} archivos PDB de segmentos PPII:\n"
                                     f"{os.getcwd()}")
    else:
        messagebox.showwarning("Atención", "No se guardó ningún segmento (quizá no existan objetos en PyMOL).")

# ==========================================================
# >>>>>>>>    FUNCIÓN: CONVERTIR_A_SELECCIONES_PYMOL    <<<<<<<<
# ==========================================================

"""
convertir_a_selecciones_pymol(pares_con_distancias)
----------------------------------------------------------
Descripción:
    Convierte una lista de pares de átomos con sus distancias en 
    selecciones de PyMOL. Cada par se representa como una selección 
    que facilita la visualización y el análisis de las interacciones 
    detectadas.

Parámetros:
    pares_con_distancias : list
        Lista de tuplas en el formato (atomo1, atomo2, distancia) 
        que serán convertidas en selecciones de PyMOL.

Retorno:
    list
        Lista con los nombres de las selecciones creadas en PyMOL.

Ejemplo de uso:
    >>> selecciones = convertir_a_selecciones_pymol(pares_con_distancias)
----------------------------------------------------------
"""

def convertir_a_selecciones_pymol(pares_con_distancias):
    selecciones = []
    for at1, at2, _ in pares_con_distancias:
        sele1 = f"id {cmd.index(at1)[0][1]}"
        sele2 = f"id {cmd.index(at2)[0][1]}"
        selecciones.append((sele1, sele2))
    return selecciones

# ==========================================================
# >>>>>>>>    FUNCIÓN: CALCULAR_Y_VISUALIZAR_DISTANCIAS    <<<<<<<<
# ==========================================================

"""
calcular_y_visualizar_distancias(max_dist=5.0, min_ang=110.0, max_ang=180.0)
----------------------------------------------------------
Descripción:
    Calcula las distancias entre átomos Cα-O de segmentos PPII 
    consecutivos y los ángulos CA-H-O asociados. Visualiza en 
    PyMOL los pares de átomos que cumplen con los criterios de 
    distancia y ángulo definidos por el usuario. Los pares 
    válidos se representan con líneas discontinuas (dashed lines) 
    en color cian para facilitar el análisis.

Parámetros:
    max_dist : float, opcional
        Distancia máxima (en Ångström) para considerar los pares 
        Cα-O como colindantes. Por defecto es 5.0 Å.
    min_ang : float, opcional
        Ángulo mínimo (en grados) para considerar válido un ángulo 
        CA-H-O. Por defecto es 110.0°.
    max_ang : float, opcional
        Ángulo máximo (en grados) para considerar válido un ángulo 
        CA-H-O. Por defecto es 180.0°.

Retorno:
    None

Ejemplo de uso:
    >>> calcular_y_visualizar_distancias()
----------------------------------------------------------
"""

def calcular_y_visualizar_distancias(max_dist=5.0, min_ang=110.0, max_ang=180.0):
    global segmentos_ppii_global
    global objetos, pares
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return
    if objetos is None:
        objetos = localizar_atomicos_clave_segmentos(segmentos_ppii_global)
    if pares is None:
        pares = calcular_distancias_colindantes(objetos, max_dist=max_dist)
    
    calcular_angulos_ca_h_o(pares, min_ang=min_ang, max_ang=max_ang)

# ==========================================================
# >>>>>>>>    FUNCIÓN: DISTANCIAS_P    <<<<<<<<
# ==========================================================

"""
distancias_p()
----------------------------------------------------------
Descripción:
    Visualiza en PyMOL todas las distancias entre los pares 
    de átomos Cα-O previamente calculados y almacenados. 
    Cada distancia se dibuja como una línea discontinua 
    (dashed line) para facilitar el análisis estructural 
    de la proteína.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> distancias_p()
----------------------------------------------------------
"""

def distancias_p():
    global segmentos_ppii_global
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return

    objetos = localizar_atomicos_clave_segmentos(segmentos_ppii_global)
    pares = calcular_distancias_colindantes(objetos)
    visualizar_distancias_pares(pares)

# ==========================================================
# >>>>>>>>    FUNCIÓN: LANZAR_INTERFAZ    <<<<<<<<
# ==========================================================

"""
lanzar_interfaz()
----------------------------------------------------------
Descripción:
    Inicia la interfaz gráfica de usuario (GUI) desarrollada 
    con Tkinter para facilitar la interacción con el programa. 
    Permite acceder a las funciones principales como cargar 
    archivos PDB, detectar hélices PPII, añadir hidrógenos, 
    eliminar solventes y generar reportes, todo desde un menú 
    visual.

Parámetros:
    Ninguno

Retorno:
    None

Ejemplo de uso:
    >>> lanzar_interfaz()
----------------------------------------------------------
"""

def lanzar_interfaz():
    root = tk.Tk()
    root.title("PPIIMoL: PPII Detect")
    root.geometry("450x1000")
    root.resizable(False, False)

    style = ttk.Style(root)
    style.theme_use('classic')  # Puedes probar: 'alt', 'clam', 'default', 'classic'

    main_frame = ttk.Frame(root, padding=5)
    main_frame.pack(fill="both", expand=True)

    def wrapper_detectar_segmentos_ppii():
        try:
            tol_phi = float(entrada_tol_phi.get())
            tol_psi = float(entrada_tol_psi.get())
            max_saltos = int(entrada_saltos.get())  # Obtener el valor del campo de saltos
            if max_saltos < 0 or max_saltos > 5:    # Validar rango (0-5)
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "¡Saltos debe ser un entero entre 0 y 5!")
            return
        
        detectar_segmentos_ppii(tol_phi=tol_phi, tol_psi=tol_psi, max_saltos=max_saltos)  # Pasar el parámetro

# En la función lanzar_interfaz(), añade este frame antes del frame de ángulos CA-H-O
    # Frame para parámetros phi/psi
    phi_psi_frame = ttk.LabelFrame(main_frame, text="Parámetros de ángulos phi/psi", padding=10)
    phi_psi_frame.pack(fill="x", pady=5)
    ttk.Label(phi_psi_frame, text="Actualmente sin la tolerancia los angulos ").pack(anchor="w")
    ttk.Label(phi_psi_frame, text="por defecto que se miden son hasta: phi 75 y psi 145").pack(anchor="w")

    ttk.Label(phi_psi_frame, text="Tolerancia para phi (±°):").pack(anchor="w")
    entrada_tol_phi = ttk.Entry(phi_psi_frame)
    entrada_tol_phi.insert(0, "20.0")
    entrada_tol_phi.pack(fill="x", pady=2)

    ttk.Label(phi_psi_frame, text="Tolerancia para psi (±°): ").pack(anchor="w")
    entrada_tol_psi = ttk.Entry(phi_psi_frame)
    entrada_tol_psi.insert(0, "20.0")
    entrada_tol_psi.pack(fill="x", pady=2)
    
    
    saltos_frame = ttk.LabelFrame(main_frame, text="Parámetros de saltos", padding=5)
    saltos_frame.pack(fill="x", pady=3)

    ttk.Label(saltos_frame, text="Saltos permitidos (0-5):").pack(anchor="w")
    entrada_saltos = ttk.Entry(saltos_frame)
    entrada_saltos.insert(0, "0")  # Valor por defecto: 0 saltos
    entrada_saltos.pack(fill="x", pady=1)

    
    def wrapper_generar_reporte_csv():
        try:
            min_ang = float(entrada_min_ang.get())
            max_ang = float(entrada_max_ang.get())
        except ValueError:
            messagebox.showerror("Error", "Introduce valores numéricos válidos para los ángulos.")
            return

        if not segmentos_ppii_global:
            messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
            return

        generar_reporte_csv(segmentos_ppii_global, min_ang=min_ang, max_ang=max_ang)
    #cambiar por otra funcion envoltorio diferente   
    
    # Función envoltorio para pasar los valores
    def wrapper_calcular_y_visualizar():
        try:
            min_ang = float(entrada_min_ang.get())
            max_ang = float(entrada_max_ang.get())
        except ValueError:
            messagebox.showerror("Error", "Introduce valores numéricos válidos para los ángulos.")
            return
        calcular_y_visualizar_distancias(min_ang=min_ang, max_ang=max_ang)

    # Botones funcionales
    botones = [
        ("Seleccionar archivo PDB", seleccionar_archivo),
        ("Descargar proteína", descargar_molecula),
        ("Eliminar solventes", eliminar_solventes),
        ("Añadir hidrógenos", anadir_hidrogenos),
        ("Ocultar cadenas laterales", ocultar_side_chains),
        ("Guardar ángulos phi/psi en archivo", guardar_csv_angulos_phi_psi),
        ("Detectar segmentos PPII y resaltarlos", wrapper_detectar_segmentos_ppii),
        ("Guardar segmentos PPII en PDB", guardar_segmentos_ppii_pdb),
        ("Visualizar distancias", distancias_p),
        ("Angulos entre CA-O-H colindantes", wrapper_calcular_y_visualizar),
        ("Generar reporte completo (CSV)", wrapper_generar_reporte_csv),
    ]

    for texto, accion in botones:
        ttk.Button(main_frame, text=texto, command=accion).pack(fill="x", pady=1)

    # Entradas de ángulos
    ang_frame = ttk.LabelFrame(main_frame, text="Parámetros de ángulo CA-H-O", padding=3)
    ang_frame.pack(fill="x", pady=3)

    ttk.Label(ang_frame, text="Ángulo mínimo (°):").pack(anchor="w")
    entrada_min_ang = ttk.Entry(ang_frame)
    entrada_min_ang.insert(0, "110.0")
    entrada_min_ang.pack(fill="x", pady=1)

    ttk.Label(ang_frame, text="Ángulo máximo (°):").pack(anchor="w")
    entrada_max_ang = ttk.Entry(ang_frame)
    entrada_max_ang.insert(0, "180.0")
    entrada_max_ang.pack(fill="x", pady=1)

    root.mainloop()


lanzar_interfaz()
