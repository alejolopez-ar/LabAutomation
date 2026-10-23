"""Distorsión armónica (THD y THD+N) con el generador UTG932E y el osciloscopio DHO804.

El generador inyecta una senoidal en la entrada del DUT. El osciloscopio adquiere
la salida (CH1, sonda 10X) y la entrada (CH2, BNC 1X) en una sola adquisición
larga (1M puntos por defecto) y el análisis se hace en numpy:

    THD   = sqrt(sum P_h, h = 2..N) / sqrt(P_1)   (FFT con ventana Blackman-Harris)
    THD+N = Vrms(residuo en 20 Hz..--bw) / Vrms(fundamental)
    residuo = señal - senoidal ajustada (ajuste de 4 parámetros, IEEE 1057)

La entrada se analiza igual que la salida: su THD es el piso del generador y del
osciloscopio en esa misma adquisición. Una THD de salida parecida a la de
entrada no es del amplificador.

Conexionado:
    Generador CH1  -->  entrada del DUT  +  osciloscopio CH2 (BNC directo, 1X)
    salida del DUT -->  osciloscopio CH1 (punta de sonda, 10X)

Uso (con el venv del proyecto, sin tocar PYTHONPATH):
    venv/bin/python examples/thd.py --freq 1000 --vpp 0.1,0.2,0.5,1 --rl 8.2

Referencia completa de argumentos y recetas: ver examples/thd.md.
"""

import argparse
import csv
import math
import os
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Libraries"))

import numpy as np
import matplotlib.pyplot as plt

from rigol_dho804.dho804 import DHO804
from unit_utg932e.utg962 import Utg962
from bode import autorange, clamp, valid, H_DIVS, TB_MIN, TB_MAX, V_DIVS, V_SCALE_MIN, V_SCALE_MAX

HALF_LOBE = 5        # bins a cada lado de un armónico que se suman (lóbulo de Blackman-Harris: ±4)
F_LOW = 20.0         # límite inferior de la banda de THD+N (Hz)
VPP_TOL = 0.05       # discrepancia máxima aceptada contra :MEASure:VPP?


def blackman_harris(n):
    k = np.arange(n) * (2 * np.pi / (n - 1))
    return 0.35875 - 0.48829 * np.cos(k) + 0.14128 * np.cos(2 * k) - 0.01168 * np.cos(3 * k)


def sine_fit(v, t, f0, iters=10):
    """Ajuste de senoidal de 4 parámetros (amplitud, fase, continua y frecuencia).

    Arranca con la frecuencia estimada por FFT y la refina por Gauss-Newton.
    Hace falta ajustar también la frecuencia: generador y osciloscopio tienen
    relojes distintos, y un error de 0,01 bin ya deja un residuo de la fundamental
    del orden del 1 %, que taparía cualquier distorsión.
    Devuelve (f, amplitud de pico, continua, senoidal ajustada).
    """
    t = t - t.mean()                  # mejor condicionado
    w = 2 * np.pi * f0
    c_, s_ = np.cos(w * t), np.sin(w * t)
    a, b, c = np.linalg.lstsq(np.column_stack([c_, s_, np.ones_like(t)]), v, rcond=None)[0]
    for _ in range(iters):
        c_, s_ = np.cos(w * t), np.sin(w * t)
        m = np.column_stack([c_, s_, np.ones_like(t), t * (-a * s_ + b * c_)])
        a, b, c, dw = np.linalg.lstsq(m, v, rcond=None)[0]
        w += dw
        if abs(dw) < 1e-12 * w:
            break
    c_, s_ = np.cos(w * t), np.sin(w * t)
    return w / (2 * np.pi), math.hypot(a, b), c, a * c_ + b * s_ + c


def analyze(v, dt, f_hint, nharm, bw):
    """THD, THD+N y residuo de una senoidal muestreada.

    v: muestras en V, dt: período de muestreo, f_hint: frecuencia aproximada.
    """
    n = len(v)
    df = 1.0 / (n * dt)
    spec = np.abs(np.fft.rfft((v - v.mean()) * blackman_harris(n))) ** 2

    # Fundamental: pico cerca de f_hint, interpolado (parábola sobre el log).
    lo, hi = int(0.8 * f_hint / df), int(1.2 * f_hint / df) + 1
    k = lo + int(np.argmax(spec[lo:hi]))
    y0, y1, y2 = np.log(spec[k - 1:k + 2])
    k_frac = k + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)

    t = np.arange(n) * dt
    f1, amp, dc, fit = sine_fit(v, t, k_frac * df)

    def band_power(f):
        kc = int(round(f / df))
        if kc + HALF_LOBE >= len(spec):
            return math.nan
        return spec[max(kc - HALF_LOBE, 1):kc + HALF_LOBE + 1].sum()

    p1 = band_power(f1)
    harm = [band_power(h * f1) / p1 for h in range(2, nharm + 1)]
    thd = math.sqrt(np.nansum(harm)) * 100.0

    # Residuo limitado a F_LOW..bw: sin la fundamental (la quita el ajuste) y sin
    # el ruido de banda ancha del osciloscopio. Es la forma de onda que muestra Self.
    r = np.fft.rfft(v - fit)
    f = np.fft.rfftfreq(n, dt)
    r[(f < F_LOW) | (f > bw)] = 0
    residual = np.fft.irfft(r, n)
    thdn = residual.std() / (amp / math.sqrt(2)) * 100.0

    # Picos de la señal limitada a la misma banda (ajuste + residuo): sin el ruido
    # de banda ancha, que con 1M de muestras inflaría el máximo y el mínimo.
    smooth = fit + residual

    return {
        "f1": f1, "amp": amp, "dc": dc, "vrms": math.sqrt(np.mean(v ** 2)),
        "vmax": smooth.max(), "vmin": smooth.min(),
        "thd": thd, "thdn": thdn,
        "harm_dbc": [10 * math.log10(h) if h > 0 else math.nan for h in harm],
        "fit": fit, "residual": residual,
        # Para graficar: relativo al bin pico de la fundamental, así se lee en dBc.
        "spec_f": f, "spec_dbc": 10 * np.log10(spec / spec[k] + 1e-30),
    }


def plot_point(path, out, inp, v_out, dt, title, nharm):
    """Figura de un punto: salida y residuo (dos períodos) y espectro de la salida."""
    n_show = min(len(v_out), int(round(2 / (out["f1"] * dt))))
    t_ms = np.arange(n_show) * dt * 1e3
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 9))
    ax1.plot(t_ms, v_out[:n_show], lw=0.8)
    ax1.set_ylabel("Salida [V]")
    ax1.set_title(title)
    ax1.grid(True, ls=":")

    ax2.plot(t_ms, out["residual"][:n_show] * 1e3, lw=0.8, color="tab:red", label="salida")
    ax2.plot(t_ms, inp["residual"][:n_show] * 1e3 * out["amp"] / inp["amp"], lw=0.8,
             color="tab:gray", alpha=0.7, label="entrada (escalada a la ganancia)")
    ax2.set_ylabel("Residuo [mV]")
    ax2.set_xlabel("Tiempo [ms]")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, ls=":")

    fmax = (nharm + 1) * out["f1"]
    m = out["spec_f"] <= fmax
    ax3.plot(out["spec_f"][m], out["spec_dbc"][m], lw=0.7)
    ax3.set_xlim(0, fmax)
    ax3.set_ylim(max(out["spec_dbc"][m].min(), -160), 5)
    ax3.set_ylabel("Salida [dBc]")
    ax3.set_xlabel("Frecuencia [Hz]")
    ax3.grid(True, ls=":")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def excursion_report(found, target, args):
    """Texto con Pmax y el límite de excursión en el punto encontrado por --find-thd.

    Para saber qué semiciclo recorta se compara cada pico medido con el de la
    senoidal ajustada (continua ± amplitud): el lado que recorta queda achatado,
    por debajo del ajuste.
    """
    vpp, out, inp = found
    vout_rms = out["amp"] / math.sqrt(2)
    pmax = vout_rms ** 2 / args.rl
    gain = out["amp"] / inp["amp"]
    short_pos = out["dc"] + out["amp"] - out["vmax"]   # cuánto le falta al pico positivo
    short_neg = out["vmin"] - (out["dc"] - out["amp"])  # ídem al negativo
    worst = max(short_pos, short_neg)
    if worst < 0.005 * out["amp"]:
        side = "No se ve recorte en los picos: la distorsión no viene del límite de excursión."
    elif min(short_pos, short_neg) > 0.75 * worst:
        side = "Recortan los dos semiciclos a la vez (excursión simétrica)."
    else:
        which = "negativo" if short_neg > short_pos else "positivo"
        side = (f"Recorta primero el semiciclo {which}: su pico queda {worst:.2f} V por debajo "
                f"de la senoidal ajustada (el otro, {min(short_pos, short_neg):.2f} V).")
    return "\n".join([
        f"Pmax al {target} % de THD ({out['f1']:.0f} Hz, R_L = {args.rl:g} Ω): {pmax:.2f} W",
        f"  THD medida {out['thd']:.3f} %, THD+N {out['thdn']:.3f} %; fundamental {vout_rms:.3f} Vrms",
        f"  Generador {vpp:.4f} Vpp; ganancia entrada-salida {gain:.2f} ({20 * math.log10(gain):.1f} dB)",
        f"Límite de excursión:",
        f"  Vp+ = {out['vmax']:+.2f} V  (margen al riel +{args.rail:g} V: {args.rail - out['vmax']:.2f} V)",
        f"  Vp- = {out['vmin']:+.2f} V  (margen al riel -{args.rail:g} V: {args.rail + out['vmin']:.2f} V)",
        f"  {side}",
    ])


def main():
    p = argparse.ArgumentParser(description="THD y THD+N con UTG932E + DHO804")
    p.add_argument("--freq", type=float, default=1000.0, help="Frecuencia de prueba (Hz, def. 1000)")
    p.add_argument("--vpp", required=True, help="Amplitud(es) del generador en Vpp, separadas por coma (barrido de potencia)")
    p.add_argument("--rl", type=float, default=8.0, help="Resistencia de carga medida (Ω) para calcular la potencia (def. 8)")
    p.add_argument("--cycles", type=float, default=200, help="Períodos en la adquisición (def. 200)")
    p.add_argument("--mdepth", default="1M", help="Profundidad de memoria (def. 1M; con 2 canales hasta 10M)")
    p.add_argument("--nharm", type=int, default=10, help="Último armónico que entra en la THD (def. 10)")
    p.add_argument("--bw", type=float, default=80e3, help="Ancho de banda de THD+N y del residuo (Hz, def. 80k, como Self)")
    p.add_argument("--stop-thd", type=float, default=3.0, help="Corta el barrido cuando la THD de salida supera este %% (def. 3)")
    p.add_argument("--find-thd", type=float, default=None, help="Busca la amplitud que da esta THD (%%, p. ej. 1) por bisección entre los puntos de --vpp que la encierran, y reporta Pmax y el límite de excursión")
    p.add_argument("--thd-tol", type=float, default=0.05, help="Tolerancia de --find-thd en puntos de %% (def. 0.05: 0,95-1,05 %%)")
    p.add_argument("--find-iter", type=int, default=10, help="Máximo de mediciones de la bisección (def. 10)")
    p.add_argument("--rail", type=float, default=15.0, help="Tensión de los rieles (V, ±) para informar el margen de excursión (def. 15)")
    p.add_argument("--hold", action="store_true", help="Con --find-thd: deja el generador en el punto encontrado hasta apretar Enter (para leer las corrientes de la fuente)")
    p.add_argument("--tag", default="", help="Etiqueta de la corrida (p. ej. iq0, iq_nominal); va al CSV y a los nombres de archivo")
    p.add_argument("--gen-channel", type=int, default=1, help="Canal del generador (def. 1)")
    p.add_argument("--ch-in", type=int, default=2, help="Canal osc. de entrada/referencia (def. 2, BNC directo)")
    p.add_argument("--ch-out", type=int, default=1, help="Canal osc. de salida del DUT (def. 1, sonda)")
    p.add_argument("--probe-in", type=float, default=1, help="Relación de sonda del canal de entrada (def. 1X)")
    p.add_argument("--probe-out", type=float, default=10, help="Relación de sonda del canal de salida (def. 10X)")
    p.add_argument("--scope-ip", default="192.168.68.68", help="IP del osciloscopio")
    p.add_argument("--scope-adb-port", default="55555", help="Puerto ADB del osciloscopio")
    p.add_argument("--scope-port", type=int, default=None, help="SCPI por socket TCP en este puerto (p. ej. 5555) en vez de VXI-11; usar si VXI-11 queda colgado")
    p.add_argument("--gen-host", default=None, help="Host del puente TCP del generador")
    p.add_argument("--gen-port", default=None, help="Puerto del puente TCP del generador")
    p.add_argument("--gen-local", action="store_true", help="Generador por USB directo (pyvisa/USBTMC), sin Pi. Ignora --gen-host/--gen-port")
    p.add_argument("--settle", type=float, default=1.0, help="Espera tras cambiar la amplitud (s, def. 1)")
    p.add_argument("--outdir", default=None, help="Carpeta de salida (def. thd_<tag>_<timestamp>)")
    p.add_argument("--no-show", action="store_true", help="No abrir la ventana del gráfico")
    args = p.parse_args()

    vpps = [float(x) for x in args.vpp.split(",")]
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or f"thd_{args.tag + '_' if args.tag else ''}{stamp}"
    os.makedirs(outdir, exist_ok=True)

    if args.gen_local:
        from unit_utg932e.utg962_local import Utg962Local
        gen = Utg962Local()
    else:
        gen = Utg962(host=args.gen_host, port=args.gen_port)
    scope = DHO804(args.scope_ip, args.scope_adb_port, args.scope_port)
    scope.connect()

    rows = []
    found = None
    try:
        # Sin promediado: promediar en el tiempo una señal que no está sincronizada
        # con el disparo deforma el espectro (ver TP07 de ME1).
        scope.set_acquire_average(1)
        for ch, probe in ((args.ch_in, args.probe_in), (args.ch_out, args.probe_out)):
            scope.set_channel_probe(ch, probe)
            scope.set_channel_display(ch, True)
            scope.set_channel_offset(ch, 0.0)
            scope.set_channel_coupling(ch, "DC")
        scope.set_channel_scale(args.ch_in, clamp(max(vpps) / (0.75 * V_DIVS),
                                                  V_SCALE_MIN * args.probe_in, V_SCALE_MAX * args.probe_in))
        scope.set_trigger_edge(args.ch_in, level=0.0)
        scope.set_timebase(clamp(args.cycles / args.freq / H_DIVS, TB_MIN, TB_MAX))
        # La memoria solo se puede cambiar en RUN; en STOP el DHO804 lo ignora sin avisar
        # (pasó: quedaron 100k puntos, bins de 25 Hz y un piso de THD de ~0,1 % de ruido).
        scope.run()
        time.sleep(0.5)
        scope.set_memory_depth(args.mdepth)
        time.sleep(0.5)
        mdep = float(scope.scpi_bridge.query(':ACQuire:MDEPth?'))
        if abs(mdep - float(args.mdepth.replace("M", "e6").replace("k", "e3"))) > 1:
            raise RuntimeError(f"La memoria quedó en {mdep:g} puntos, no en {args.mdepth}")
        # Límite de 20 MHz: saca ruido de banda ancha que se pliega en el muestreo y no
        # afecta los armónicos de audio.
        for ch in (args.ch_in, args.ch_out):
            scope.scpi_bridge.write(f':CHANnel{ch}:BWLimit 20M')

        csv_path = os.path.join(outdir, "thd.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Tag", "Etapa", "Freq_Hz", "Gen_Vpp", "Vin_rms", "Vout_rms", "Pout_W",
                             "Vpk_pos_V", "Vpk_neg_V",
                             "THD_out_pct", "THDN_out_pct", "THD_in_pct", "THDN_in_pct",
                             "H2_dBc", "H3_dBc", "H4_dBc", "H5_dBc"])
            print(f"{'Gen[Vpp]':>9} {'Vout[Vrms]':>10} {'Pout[W]':>9} {'Vp+[V]':>7} {'Vp-[V]':>7} "
                  f"{'THD out[%]':>11} {'THD+N out[%]':>13} {'THD in[%]':>10} {'H2[dBc]':>8} {'H3[dBc]':>8}")

            gain_est = []      # ganancia del último punto medido (lista para poder mutarla)

            def measure(vpp, stage):
                """Mide un punto, lo agrega al CSV y a las figuras y devuelve (out, inp)."""
                i = len(rows)
                gen.set_sine(args.gen_channel, args.freq, -vpp / 2.0, vpp / 2.0)
                # Escala anticipada: con la ganancia del punto anterior, la salida nueva
                # entra en pantalla con 30 % de margen. Si la traza se sale de pantalla,
                # el ADC la recorta y el DHO igual devuelve una Vpp "válida", así que el
                # autorango solo no alcanza para detectarlo.
                if gain_est:
                    scope.set_channel_scale(args.ch_out, clamp(
                        1.3 * gain_est[0] * vpp / (0.75 * V_DIVS),
                        V_SCALE_MIN * args.probe_out, V_SCALE_MAX * args.probe_out))
                scope.run()
                time.sleep(args.settle)
                autorange(scope, args.ch_in, args.settle, args.probe_in)
                autorange(scope, args.ch_out, args.settle, args.probe_out)

                for intento in range(4):
                    vpp_ref = scope.measure_vpp(args.ch_out)
                    try:
                        scope.single(force_after=5.0)
                        v_out, dt = scope.read_waveform(args.ch_out)
                        v_in, _ = scope.read_waveform(args.ch_in)
                    except (IOError, TimeoutError) as e:
                        print(f"  reintento: {e}")
                        scope.run()
                        time.sleep(args.settle)
                        continue
                    scope.run()
                    # Sobre-rango del ADC: la pantalla cubre ±4 divisiones. Si la señal
                    # llega al borde, alejar la escala y adquirir de nuevo.
                    over = [ch for ch, v in ((args.ch_out, v_out), (args.ch_in, v_in))
                            if np.abs(v).max() >= 0.97 * (V_DIVS / 2) * scope.get_channel_scale(ch)]
                    if not over:
                        break
                    for ch in over:
                        scope.set_channel_scale(ch, 2.0 * scope.get_channel_scale(ch))
                    print(f"  señal fuera de pantalla en CH{over}: escala ×2 y repito")
                    time.sleep(args.settle)
                else:
                    raise RuntimeError("No se pudo adquirir el punto sin sobre-rango")

                out = analyze(v_out, dt, args.freq, args.nharm, args.bw)
                inp = analyze(v_in, dt, args.freq, args.nharm, args.bw)
                # Control de la decodificación (orden de bytes, escala, sonda): la Vpp
                # de los datos leídos tiene que coincidir con la que mide el osciloscopio.
                # Se usa VPP y no VRMS porque el DHO804 a veces devuelve VRMS inválido
                # (9.9E37) con esta base de tiempo; si la referencia es inválida, se saltea.
                vpp_raw = float(v_out.max() - v_out.min())
                if valid(vpp_ref) and abs(vpp_raw - vpp_ref) > VPP_TOL * vpp_ref:
                    print(f"  ¡OJO! Vpp decodificada {vpp_raw:.4f} V vs. osciloscopio "
                          f"{vpp_ref:.4f} V: revisar la lectura de la forma de onda")

                vout_rms = out["amp"] / math.sqrt(2)
                pout = vout_rms ** 2 / args.rl
                h = out["harm_dbc"]
                writer.writerow([args.tag, stage, f"{out['f1']:.3f}", f"{vpp:.5f}",
                                 f"{inp['amp'] / math.sqrt(2):.6f}", f"{vout_rms:.6f}", f"{pout:.5f}",
                                 f"{out['vmax']:.4f}", f"{out['vmin']:.4f}",
                                 f"{out['thd']:.6f}", f"{out['thdn']:.6f}",
                                 f"{inp['thd']:.6f}", f"{inp['thdn']:.6f}"]
                                + [f"{x:.2f}" for x in h[:4]])
                f.flush()
                rows.append((pout, out["thd"], out["thdn"], inp["thd"]))
                gain_est[:] = [2 * out["amp"] / vpp]       # Vpp de salida por Vpp de generador
                print(f"{vpp:9.4f} {vout_rms:10.4f} {pout:9.4f} {out['vmax']:7.2f} {out['vmin']:7.2f} "
                      f"{out['thd']:11.5f} {out['thdn']:13.5f} {inp['thd']:10.5f} {h[0]:8.1f} {h[1]:8.1f}")

                name = f"{args.tag + '_' if args.tag else ''}{i:02d}_{vpp:g}Vpp"
                np.savez_compressed(os.path.join(outdir, f"raw_{name}.npz"),
                                    v_out=v_out.astype(np.float32), v_in=v_in.astype(np.float32), dt=dt)
                plot_point(os.path.join(outdir, f"residuo_{name}.png"), out, inp, v_out, dt,
                           f"{args.tag} {out['f1']:.0f} Hz, {pout:.3g} W: THD {out['thd']:.4f} %, "
                           f"THD+N {out['thdn']:.4f} % ({args.bw / 1e3:g} kHz)", args.nharm)
                return out, inp

            # Barrido: en modo --find-thd corta en el primer punto que pasa el objetivo.
            limit = args.find_thd if args.find_thd is not None else args.stop_thd
            below = above = None
            for vpp in vpps:
                out, inp = measure(vpp, "barrido")
                if out["thd"] > limit:
                    above = (vpp, out, inp)
                    print(f"THD > {limit} %: fin del barrido")
                    break
                below = (vpp, out, inp)

            if args.find_thd is not None:
                target = args.find_thd
                if above is None:
                    print(f"\nNinguna amplitud de --vpp llegó al {target} % de THD: agregar valores más altos.")
                elif below is None:
                    print(f"\nLa primera amplitud de --vpp ya supera el {target} %: agregar valores más bajos.")
                else:
                    # Bisección sobre la amplitud del generador. Cerca del recorte la THD
                    # crece muy rápido (de 0,01 % a varios % en menos de 1 dB), por eso
                    # interpolar entre puntos de un barrido fijo erra mucho.
                    print(f"\nBuscando {target} ± {args.thd_tol} % entre {below[0]:g} y {above[0]:g} Vpp")
                    lo, hi = below[0], above[0]
                    found = min((below, above), key=lambda x: abs(x[1]["thd"] - target))
                    for _ in range(args.find_iter):
                        mid = (lo + hi) / 2
                        out, inp = measure(mid, "busqueda")
                        if abs(out["thd"] - target) < abs(found[1]["thd"] - target):
                            found = (mid, out, inp)
                        if abs(out["thd"] - target) <= args.thd_tol:
                            break
                        lo, hi = (mid, hi) if out["thd"] < target else (lo, mid)
                    report = excursion_report(found, target, args)
                    print("\n" + report)
                    with open(os.path.join(outdir, "pmax.txt"), "w") as fr:
                        fr.write(report + "\n")
                    if args.hold:
                        gen.set_sine(args.gen_channel, args.freq, -found[0] / 2.0, found[0] / 2.0)
                        input("\nGenerador sostenido en el punto de Pmax: leer las corrientes de "
                              "los dos rieles en la fuente y apretar Enter para apagar.")
    finally:
        gen.set_output(args.gen_channel, False)
        gen.close()
        scope.run()
        scope.close()

    print(f"\nResultados en {outdir}/")
    if len(rows) > 1:
        # La bisección agrega puntos fuera de orden: ordenar por potencia.
        pw, thd, thdn, thd_in = map(np.array, zip(*sorted(rows)))
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.loglog(pw, thd, "-o", markersize=4, label="THD salida")
        ax.loglog(pw, thdn, "-s", markersize=4, label=f"THD+N salida ({args.bw / 1e3:g} kHz)")
        ax.loglog(pw, thd_in, "--", color="tab:gray", label="THD entrada (piso generador + osc.)")
        if found is not None:
            pmax = (found[1]["amp"] / math.sqrt(2)) ** 2 / args.rl
            ax.loglog([pmax], [found[1]["thd"]], "*", color="tab:red", markersize=14,
                      label=f"Pmax = {pmax:.2f} W al {found[1]['thd']:.2f} %")
        ax.set_xlabel(f"Potencia de salida [W] (R_L = {args.rl:g} Ω)")
        ax.set_ylabel("Distorsión [%]")
        ax.set_title(f"THD vs. potencia, {args.freq:g} Hz {args.tag}")
        ax.grid(True, which="both", ls=":")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "thd_vs_potencia.png"), dpi=120)
        if not args.no_show:
            plt.show()


if __name__ == "__main__":
    main()
