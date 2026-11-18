#!/usr/bin/env python3
"""
Diseñado para explorar el espacio de parámetros variando mZL11 como candidato a DM (~ZL11)
y manteniendo mZL22 controlada por encima para asegurar la jerarquía de masas.

Uso:
    python3 version_sabado_busqueda_flavons.py --regime high_mass --samples "n"
"""

import os
import subprocess
import re
import csv
import time
import random
import argparse
import multiprocessing
from datetime import datetime
import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# --- CONFIGURACIÓN GLOBAL ---
DEFAULT_MAIN_SRC = "main.cpp"
DEFAULT_PAR_FILE = "data.par"
DEFAULT_OUTPUT_VALID = "valid_dm_dataset.csv"
DEFAULT_OUTPUT_ALL = "all_scan_results.csv"
DEFAULT_OUTPUT_EXCLUDED = "excluded_regions.csv"
TARGET_DM = "~ZL11"
OMEGA_BOUNDS = (0.01, 0.15)  # Rango compatible con subabundancia y Planck
SCAN_TIMEOUT = 180

# --- FUNCIONES AUXILIARES PARA data.par PLANO ---

def read_baseline_parameters(par_file):
    """Lee los parámetros base desde data.par (formato clave-valor plano)."""
    params = {}
    if not os.path.exists(par_file):
        print(f"[!] Advertencia: {par_file} no encontrado. Usando diccionario vacío.")
        return params

    with open(par_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                params[parts[0]] = parts[1]
    return params

def write_par_file(filename, params):
    """Escribe el archivo temporal .par plano manteniendo las claves exactas."""
    with open(filename, 'w') as f:
        for key, value in params.items():
            f.write(f"{key} {value}\n")

def compile_micromegas(src=DEFAULT_MAIN_SRC):
    """Compila el ejecutable de micrOMEGAs."""
    print(f"[*] Compilando el ejecutable de micrOMEGAs desde {src}...")
    try:
        result = subprocess.run(
            ["make", f"main={src}"],
            capture_output=True,
            text=True,
            check=True
        )
        print("[+] Compilación exitosa.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[!] Error en la compilación:\n{e.stderr}")
        return False

def parse_micromegas_output(stdout_text):
    """
    Analiza la salida estándar de micrOMEGAs extrayendo candidato, masa, 
    densidad relicua, sección eficaz SI y estado de exclusión por DD.
    """
    results = {
        'dm_candidate': None,
        'dm_mass': None,
        'omega': None,
        'sigma_si': None,
        'excluded_by_dd': False,
        'dd_experiment': "None"
    }

    # 1. Candidato y Masa
    dm_match = re.search(r"Dark\s+matter\s+candidate\s+is\s+['\"]([^'\"]+)['\"]\s+with\s+spin=[^\s]+\s+mass=([0-9\.eE\+\-]+)", stdout_text)
    if dm_match:
        results['dm_candidate'] = dm_match.group(1)
        try:
            results['dm_mass'] = float(dm_match.group(2))
        except ValueError:
            pass

    # 2. Densidad Relicua
    omega_match = re.search(r"Xf=[^\s]+\s+Omega=([0-9\.eE\+\-]+)", stdout_text)
    if omega_match:
        try:
            results['omega'] = float(omega_match.group(1))
        except ValueError:
            pass

    # 3. Sección eficaz Spin-Independent (Protón)
    si_match = re.search(r"proton\s+SI\s+([0-9\.eE\+\-]+)", stdout_text)
    if si_match:
        try:
            results['sigma_si'] = float(si_match.group(1))
        except ValueError:
            pass

    # 4. Exclusión por Detección Directa
    excl_match = re.search(r"Excluded\s+by\s+([A-Za-z0-9_\-]+)", stdout_text)
    if excl_match:
        results['excluded_by_dd'] = True
        results['dd_experiment'] = excl_match.group(1)

    not_excl_match = re.search(r"Not\s+excluded\s+by\s+DD\s+experiments", stdout_text, re.IGNORECASE)
    if not_excl_match:
        results['excluded_by_dd'] = False
        results['dd_experiment'] = "None"

    return results

# --- PROCESAMIENTO PARALELO ---

def run_point_worker(args):
    """Evalúa un punto individual del espacio de parámetros de forma aislada."""
    worker_id, params, timeout = args
    temp_par = f"data_worker_{worker_id}.par"
    write_par_file(temp_par, params)

    try:
        proc = subprocess.run(
            ["./main", temp_par],
            capture_output=True,
            text=True,
            timeout=timeout
        )

        if proc.returncode != 0:
            return {
                'status': 'error',
                'params': params,
                'error': f"Exit code {proc.returncode}"
            }

        parsed = parse_micromegas_output(proc.stdout)

        # Criterios de validación física
        valid = True
        reasons = []

        if parsed['dm_candidate'] != TARGET_DM:
            valid = False
            reasons.append(f"DM_is_{parsed['dm_candidate']}")

        if parsed['omega'] is None:
            valid = False
            reasons.append("Omega_NaN")
        elif not (OMEGA_BOUNDS[0] <= parsed['omega'] <= OMEGA_BOUNDS[1]):
            valid = False
            reasons.append(f"Omega_OOB({parsed['omega']:.2e})")

        if parsed['excluded_by_dd']:
            valid = False
            reasons.append(f"DD_Excl_{parsed['dd_experiment']}")

        parsed['is_valid'] = valid
        parsed['validity_reason'] = ";".join(reasons) if reasons else "Valid"
        parsed['params'] = params

        return {'status': 'success', 'data': parsed}

    except subprocess.TimeoutExpired:
        return {'status': 'error', 'params': params, 'error': 'Timeout'}
    except Exception as e:
        return {'status': 'error', 'params': params, 'error': str(e)}
    finally:
        if os.path.exists(temp_par):
            try:
                os.remove(temp_par)
            except:
                pass

# --- CLASE ESCÁNER ---

class MicromegasScanner:
    def __init__(self, n_workers=None):
        self.n_workers = n_workers or max(1, multiprocessing.cpu_count() - 1)
        self.baseline = read_baseline_parameters(DEFAULT_PAR_FILE)
        self.results_valid = []
        self.results_all = []
        self.results_excluded = []

    def sample_parameters(self, regime, n_samples):
        """
        Genera puntos de muestreo variando mZL11 como candidato y ajustando mZL22.
        """
        points = []
        for _ in range(n_samples):
            p = self.baseline.copy()
            
            if regime == "funnel":
                m_val = random.uniform(60.0, 65.0)  # Cerca de la resonancia
                g_val = 10**random.uniform(-4.0, -2.0)
            elif regime == "high_mass":
                m_val = random.uniform(500.0, 1000.0)  # Régimen de alta masa
                g_val = 10**random.uniform(-2.0, 0.5)
            else:  # mixed
                m_val = random.uniform(60.0, 80.0) if random.random() < 0.5 else random.uniform(500.0, 1000.0)
                g_val = 10**random.uniform(-3.0, 0.0)

            # 1. Variar la masa del candidato principal mZL11
            p["mZL11"] = f"{m_val:.6g}"

            # 2. Asegurar que mZL22 sea estrictamente mayor (para que mZL11 sea el candidato LOP)
            if "mZL22" in p:
                delta = random.uniform(1.0, 10.0) if regime == "high_mass" else random.uniform(20.0, 50.0)
                p["mZL22"] = f"{m_val + delta:.6g}"

            # 3. Variar el acoplamiento gLL
            if "gLL" in p:
                p["gLL"] = f"{g_val:.6g}"

            points.append(p)
        return points

    def run_scan(self, points, timeout=SCAN_TIMEOUT):
        print(f"[*] Iniciando escaneo de {len(points)} puntos con {self.n_workers} procesos...")
        task_args = [(i, p, timeout) for i, p in enumerate(points)]

        if tqdm:
            pbar = tqdm(total=len(points), desc="Escaneando")

        with multiprocessing.Pool(processes=self.n_workers) as pool:
            for res in pool.imap_unordered(run_point_worker, task_args):
                if res['status'] == 'success':
                    data = res['data']
                    entry = {**data['params'], **{k: v for k, v in data.items() if k != 'params'}}
                    self.results_all.append(entry)
                    if data['is_valid']:
                        self.results_valid.append(entry)
                    else:
                        self.results_excluded.append(entry)

                if tqdm:
                    pbar.update(1)
                    pbar.set_description(f"Escaneando (Válidos: {len(self.results_valid)})")

        if tqdm:
            pbar.close()

        print(f"\n[+] Escaneo finalizado.")
        print(f"    Evaluados: {len(self.results_all)} | Válidos: {len(self.results_valid)} | Excluidos/No viables: {len(self.results_excluded)}")

    def save_results(self):
        if not self.results_all:
            return
        headers = list(self.results_all[0].keys())

        # Guardar todo
        with open(DEFAULT_OUTPUT_ALL, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(self.results_all)

        # Guardar válidos
        if self.results_valid:
            with open(DEFAULT_OUTPUT_VALID, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                writer.writerows(self.results_valid)
            print(f"[+] Puntos válidos guardados en {DEFAULT_OUTPUT_VALID}")

        # Guardar regiones excluidas
        if self.results_excluded:
            with open(DEFAULT_OUTPUT_EXCLUDED, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                writer.writerows(self.results_excluded)
            print(f"[+] Dataset de exclusión guardado en {DEFAULT_OUTPUT_EXCLUDED}")

# --- EJECUCIÓN PRINCIPAL ---

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Escáner acotado por literatura IDM/Singlete.")
    parser.add_argument("--regime", choices=["funnel", "high_mass", "mixed"], default="funnel", help="Régimen físico a explorar.")
    parser.add_argument("--samples", type=int, default=200, help="Número de muestras.")
    parser.add_argument("--workers", type=int, default=None, help="Procesos paralelos.")
    parser.add_argument("--no-compile", action="store_true", help="Omitir compilación.")
    args = parser.parse_args()

    if not args.no_compile or not os.path.exists("./main"):
        if not compile_micromegas():
            exit(1)

    scanner = MicromegasScanner(n_workers=args.workers)
    points = scanner.sample_parameters(args.regime, args.samples)

    try:
        scanner.run_scan(points)
        scanner.save_results()
    except KeyboardInterrupt:
        print("\n[!] Interrumpido. Guardando resultados parciales...")
        scanner.save_results()
