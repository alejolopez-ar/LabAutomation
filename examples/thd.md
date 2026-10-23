# thd.py — Distorsión armónica (THD y THD+N) automatizada

Mide la distorsión de un DUT (pensado para el amplificador clase AB del TP04 de EA2)
inyectando una senoidal con el **UTG932E** y adquiriendo entrada y salida con el
**DHO804** en una sola adquisición larga. El análisis lo hace numpy, no la FFT del
osciloscopio:

- **THD** = √(Σ P_h, h = 2..N) / √P₁ con FFT de 1M puntos y ventana Blackman-Harris.
  Cada armónico es la suma de ±5 bins alrededor de h·f₁.
- **THD+N** = Vrms del residuo en 20 Hz–`--bw` / Vrms de la fundamental.
- **Residuo** = señal − senoidal ajustada. El ajuste es de 4 parámetros (IEEE 1057):
  también ajusta la frecuencia, porque el generador y el osciloscopio tienen relojes
  distintos y un error de 0,01 bin ya deja ~1 % de fundamental en el residuo.
  Es la forma de onda que muestra Self (Figs. 10.19–10.22): con Iq baja aparecen
  **picos en cada cruce por cero**.
- **La entrada se analiza igual que la salida.** Su THD es el **piso** del generador más
  el osciloscopio en esa misma adquisición. Si la salida da parecido a la entrada, esa
  distorsión no es del amplificador.

## Conexionado

Igual que `bode.py`:

```
Generador CH1 ──► entrada del DUT ──┬──► Osc CH2  (BNC directo, 1X)   → entrada (piso)
                                    │
       salida DUT (sobre la carga) ─┴ ... ──► Osc CH1 (sonda 10X) → salida
```

En el TP04: **carga de 8,2 Ω 10 W** en J3 (no el parlante). La masa de la sonda va a
GNDPWR, en el borne de J3.

## Uso

```bash
venv/bin/python examples/thd.py --freq 1000 --vpp 0.1,0.2,0.5,1 --rl 8.2 --tag prueba

# Como se usó en el banco del TP04: DHO804 por socket, UTG por USB directo
venv/bin/python examples/thd.py --gen-local --scope-ip 10.229.209.41 --scope-port 5555 \
  --freq 1000 --vpp 0.1,0.2,0.5,1 --rl 8.2 --tag prueba
```

Cada corrida crea `thd_<tag>_<timestamp>/` con:
| archivo | qué tiene |
|---|---|
| `thd.csv` | por punto: etapa (`barrido` o `busqueda`), Vin, Vout, Pout, **Vp+ y Vp−**, THD y THD+N de salida y de entrada, H2–H5 en dBc |
| `pmax.txt` | con `--find-thd`: Pmax, ganancia y límite de excursión (ver abajo) |
| `residuo_<tag>_NN_<vpp>Vpp.png` | dos períodos de la salida, el residuo (salida en rojo; entrada en gris, escalada a la ganancia) y el espectro de la salida en dBc |
| `thd_vs_potencia.png` | THD y THD+N vs Pout en log-log, con el piso de entrada (si hay más de un punto) |
| `raw_*.npz` | formas de onda crudas (`v_out`, `v_in`, `dt`) para volver a analizarlas sin el banco |

Para volver a analizar un `.npz`:
```python
from thd import analyze
d = np.load("raw_....npz"); r = analyze(d["v_out"].astype(float), float(d["dt"]), 1000, 10, 80e3)
```

## Argumentos

| flag | def | qué hace |
|---|---|---|
| `--freq` | 1000 | frecuencia de prueba (Hz) |
| `--vpp` | **obligatorio** | amplitud(es) del generador, separadas por coma. Varias = barrido de potencia |
| `--rl` | 8 | resistencia de carga **medida** (Ω), solo para calcular Pout |
| `--cycles` | 200 | períodos en la adquisición. Con 1M puntos a 1 kHz son 200 ms, 5 MSa/s y bins de 5 Hz |
| `--mdepth` | 1M | memoria. Más puntos bajan el piso de ruido por bin, pero la lectura es más lenta (~2 MB por canal y por millón) |
| `--nharm` | 10 | último armónico de la THD |
| `--bw` | 80000 | banda de THD+N y del residuo (Self usa 80 kHz). Los picos de cruce son de alto orden: con más banda se ven más agudos |
| `--stop-thd` | 3 | corta el barrido cuando la THD de salida supera este %: no se sigue recortando sobre la carga |
| `--find-thd` | — | busca la amplitud que da esta THD (p. ej. `1`): barre `--vpp` hasta pasarla y después hace bisección entre los dos puntos que la encierran |
| `--thd-tol` | 0.05 | tolerancia de la búsqueda, en puntos de % (1 ± 0,05 %) |
| `--find-iter` | 10 | máximo de mediciones de la bisección (suelen alcanzar 5–7) |
| `--rail` | 15 | rieles ± (V), solo para informar el margen de excursión |
| `--hold` | — | con `--find-thd`: deja el generador en el punto encontrado hasta apretar Enter, para leer las corrientes de los rieles |
| `--tag` | — | etiqueta de la corrida (al CSV y a los nombres de archivo) |
| `--settle` | 1 | espera tras cambiar la amplitud (s) |
| `--outdir` | auto | carpeta de salida |
| `--no-show` | — | no abrir el gráfico |

Canales, sondas y red: los mismos flags que `bode.py` (`--ch-in/--ch-out`,
`--probe-in/--probe-out`, `--gen-channel`, `--scope-ip`, `--scope-port`, `--gen-host/--gen-port`,
`--gen-local`). **Usar `--scope-port 5555`** (SCPI por socket TCP): si una transferencia se corta a
mitad de camino, el servicio VXI-11 del DHO804 queda colgado hasta reiniciar su LAN, y el socket no.

## Cómo funciona (resumen)

1. Promediado de adquisición **apagado**: promediar en el tiempo una señal que no está
   sincronizada con el disparo deforma el espectro (hallazgo del TP07 de ME1).
2. **Memoria fijada en RUN y verificada.** En STOP el DHO804 ignora `:ACQuire:MDEPth` sin
   avisar: en el banco quedaron 100k puntos (bins de 25 Hz) y el piso de THD subió a un
   falso ~0,1 %. Si la memoria no queda en el valor pedido, el script se detiene.
3. **Límite de ancho de banda de 20 MHz** en los dos canales: saca ruido de banda ancha
   que se pliega en el muestreo y no afecta los armónicos de audio.
4. Por punto: fija la amplitud, **escala anticipada** de la salida (ganancia del punto
   anterior, con 30 % de margen), autorango vertical (el de `bode.py`, que llena ~6 de 8
   divisiones para aprovechar los 12 bits), `:SINGle`, lectura de la memoria en modo RAW
   y formato WORD por bloques de 100k puntos, y vuelve a `:RUN`.
5. **Sobre-rango del ADC:** si los datos llegan al borde de la pantalla (±4 divisiones),
   duplica la escala y repite. Hace falta porque con la traza recortada por el ADC el DHO804
   igual devuelve una Vpp "válida" y el autorango no se entera (en el banco dio un punto
   falso de 26 % de THD en medio de un barrido limpio).
6. **Reintentos:** si una lectura vuelve vacía o la adquisición no termina, reintenta hasta
   4 veces.
7. **Control de la lectura:** compara la Vpp de los datos leídos con `:MEASure:VPP?` del
   osciloscopio (VRMS no sirve: con esta base de tiempo el DHO804 a veces la da inválida). Si difieren más de un 5 % imprime `¡OJO!`. La guía de programación no
   dice el orden de bytes del WORD (se asume little-endian): este control lo detecta,
   igual que un error de escala o de sonda.

## Pmax y límite de excursión (`--find-thd`)

Self define la potencia máxima en el punto donde el recorte lleva la THD al 0,1 % o al
1 % (cap. 1, p. 7). Cerca del recorte la THD pasa de ~0,01 % a varios % en menos de 1 dB,
así que interpolar entre puntos de un barrido fijo erra mucho. La bisección converge en
~5–7 mediciones, lo que deja poco tiempo a plena potencia sobre la carga.

- **Picos:** Vp+ y Vp− son el máximo y el mínimo de la señal limitada a `--bw` (ajuste +
  residuo), así el ruido de banda ancha no infla los picos.
- **Qué lado recorta:** se compara cada pico con el de la senoidal ajustada. El semiciclo
  que recorta queda achatado, por debajo del ajuste. El reporte dice si recorta el positivo,
  el negativo o los dos a la vez.
- **Potencia:** la de la fundamental (Vrms² / R_L), igual que el resto del script.

`pmax.txt` medido en el TP04 (2026-10-22, amplificador CFP con NE5532, ±15 V):
```
Pmax al 1.0 % de THD (1000 Hz, R_L = 8.2 Ω): 8.60 W
  THD medida 1.007 %, THD+N 1.064 %; fundamental 8.396 Vrms
  Generador 1.4937 Vpp; ganancia entrada-salida 16.31 (24.2 dB)
Límite de excursión:
  Vp+ = +11.80 V  (margen al riel +15 V: 3.20 V)
  Vp- = -11.74 V  (margen al riel -15 V: 3.26 V)
  Recorta primero el semiciclo negativo: su pico queda 0.29 V por debajo de la senoidal ajustada (el otro, -0.08 V).
```
`--hold` espera un Enter en la terminal: si el script lo corre otra herramienta sin terminal,
no usarlo y sostener la señal aparte con `Utg962Local().set_sine(...)`.

## Recetas para el TP04

La ganancia depende de la frecuencia (la "sonrisa") y del pote de volumen: **dejar el
volumen al máximo** y anotarlo. Primero una corrida corta para conocer la ganancia:
`Vpp_gen = 2·√(2·P·R_L) / G`.

```bash
# THD vs potencia a 1 kHz, hasta el recorte (corta solo al pasar el 3 %)
venv/bin/python examples/thd.py --freq 1000 --rl 8.2 --tag p1k \
  --vpp 0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.1,1.2,1.3,1.4,1.5

# Pmax al 1 % de THD + límite de excursión. --hold deja la salida en ese punto para leer
# las corrientes de los dos rieles en la fuente (rendimiento = P / (15·(I+ + I−))).
venv/bin/python examples/thd.py --freq 1000 --rl 8.2 --tag pmax1k --find-thd 1 --hold \
  --vpp 0.2,0.5,0.8,1.0,1.2,1.4,1.6,1.8,2.0

# Lo mismo a 10 kHz (la distorsión de cruce sube 6 dB/octava, según Self)
venv/bin/python examples/thd.py --freq 10000 --rl 8.2 --tag p10k \
  --vpp 0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.1,1.2,1.3,1.4,1.5

# Cruce vs Iq: una amplitud fija (~1–2 Vp a la salida); una corrida por cada ajuste de RV2.
# ANOTAR Vq antes de mover RV2 para volver al nominal.
venv/bin/python examples/thd.py --freq 1000  --rl 8.2 --vpp 0.2 --tag iq0_1k
venv/bin/python examples/thd.py --freq 10000 --rl 8.2 --vpp 0.2 --tag iq0_10k
#   ... repetir con --tag iq_mitad_* e iq_nominal_*
```

## Limitaciones

- **Piso real medido (TP04, 1 kHz, 1M puntos, sonda ×10 en la salida):** THD de
  0,03–0,06 % entre 0,01 y 6,5 W, sin ningún armónico sobre el ruido del canal (≈ −85 a
  −90 dBc por bin de 2,5 Hz). Ese número es ruido sumado en ±5 bins alrededor de cada
  armónico, no distorsión: hay que mirar el espectro y el residuo antes de informarlo. Con
  cuantización ideal el análisis llegaría a ~0,0006 % (señales sintéticas), así que el
  límite es el ruido del front-end del DHO804. Por eso se informa la THD de entrada como
  piso, punto por punto. La THD+N (0,2–0,5 %) también es ruido del osciloscopio.
- **El piso de entrada y el de salida no son idénticos:** se miden en canales distintos,
  con escalas distintas y con sonda en la salida. Es una cota, no una resta.
- **Validado** offline con señales sintéticas (THD conocida de 0,0105 % y de 1 %: error
  relativo < 0,2 %) y con instrumentos simulados, y en el banco con el DHO804 y el UTG932E
  (TP04, 2026-10-22): la Vpp leída coincide con la del osciloscopio dentro del ~1 % (ruido
  sobre los picos).
- **WiFi:** transferir 1M de puntos por canal tarda varios segundos y una transferencia
  cortada deja colgado el SCPI del osciloscopio (VXI-11 y, una vez, también el socket).
  Remedio: reiniciar la LAN del DHO (Utility → IO). Los bloques de 100k bajan la probabilidad.
- **Zumbido de red:** con el generador y el osciloscopio a tierra aparece 50/150/250 Hz a
  ≈ −58 dBc en señales chicas (lazo de masa). No entra en la THD (no son armónicos de la
  prueba) pero sí en la THD+N.

## Archivos relacionados

- `Libraries/rigol_dho804/dho804/dho804.py`: `set_memory_depth`, `single` (con
  `force_after` para forzar el disparo), `run` y `read_waveform`, agregados para este script;
  `DHO804(ip, adbport, scpi_port=5555)` para SCPI por socket.
- `examples/bode.py`: de ahí salen `autorange` y las constantes de pantalla.
