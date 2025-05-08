# CIENBIO - Script Básico para Diseño de Sondas qPCR (Versión Local)

import pandas as pd
import math
import tkinter as tk
from tkinter import filedialog

# Parámetros SantaLucia 1998
NN_params = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7), 'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0), 'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
}
R = 1.987
primer_conc = 500e-9
salt_conc = 50e-3

def calcular_tm_dg(seq):
    seq = seq.upper()
    dH, dS = 0, 0
    for i in range(len(seq) - 1):
        par = seq[i:i+2]
        deltaH, deltaS = NN_params.get(par, (-7.0, -20.0))
        dH += deltaH
        dS += deltaS
    dS += -1.4
    dH += 0.2
    tm = (1000 * dH) / (dS + (R * math.log(primer_conc))) - 273.15 + 16.6 * math.log10(salt_conc)
    return round(tm, 2), round(dH, 2)

def evaluar_sondas(excel_path):
    df = pd.read_excel(excel_path)
    resultados = []
    for idx, row in df.iterrows():
        seq_id = row["ID"]
        seq = row["Secuencia"].upper().replace(" ", "").replace("\n", "")
        snp_pos = int(row["Coordenada SNP"]) - 1
        alelo_ref = row["Alelo_Ref"].upper()
        alelo_alt = row["Alelo_Alt"].upper()
        for length in range(18, 23):
            half = length // 2
            start = max(snp_pos - half, 0)
            end = min(start + length, len(seq))
            if end - start == length:
                sonda_ref = seq[start:end]
                if snp_pos in range(start, end):
                    pos_rel = snp_pos - start
                    sonda_alt = sonda_ref[:pos_rel] + alelo_alt + sonda_ref[pos_rel+1:]
                    tm_ref, dg_ref = calcular_tm_dg(sonda_ref)
                    tm_alt, dg_alt = calcular_tm_dg(sonda_alt)
                    resultados.append({
                        "ID": seq_id,
                        "Sonda_ref": sonda_ref,
                        "Sonda_alt": sonda_alt,
                        "Longitud": length,
                        "Inicio": start + 1,
                        "SNP Central": snp_pos + 1,
                        "Tm_ref (°C)": tm_ref,
                        "Tm_alt (°C)": tm_alt,
                        "ΔTm (°C)": round(tm_alt - tm_ref, 2),
                        "ΔH_ref (kcal/mol)": dg_ref,
                        "ΔH_alt (kcal/mol)": dg_alt
                    })
    return pd.DataFrame(resultados)

if __name__ == "__main__":
    print("🔍 Selecciona un archivo .xlsx con columnas: ID, Secuencia, Coordenada SNP, Alelo_Ref, Alelo_Alt")
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx")])
    if file_path:
        df_resultado = evaluar_sondas(file_path)
        output_file = file_path.replace(".xlsx", "_evaluado.xlsx")
        df_resultado.to_excel(output_file, index=False)
        print(f"✅ Evaluación completada. Archivo guardado: {output_file}")
    else:
        print("❌ No se seleccionó ningún archivo.")
